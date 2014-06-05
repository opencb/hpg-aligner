/*
 * wanderer.c
 *
 *  Created on: 24/04/2014
 *      Author: rmoreno
 */

#include "wanderer.h"

/**
 * STATIC VARS
 */
//bwander_obtain_region
static bam1_t *last_read = NULL;
static int last_read_bytes = 0;

/**
 * STATIC FUNCTIONS
 */
static int bwander_region_insert(bam_wanderer_t *wanderer, bam_region_t *region);
static int bwander_obtain_region(bam_wanderer_t *wanderer, bam_region_t *current_region);

void
bwander_init(bam_wanderer_t *wanderer)
{
	int i;
	bam_region_t *region;

	assert(wanderer);

	//Set all to zero
	memset(wanderer, 0, sizeof(bam_wanderer_t));

	//Create regions
	wanderer->regions_list = linked_list_new(COLLECTION_MODE_SYNCHRONIZED);

	//Init locks
	omp_init_lock(&wanderer->regions_lock);
	omp_init_lock(&wanderer->free_slots);
	omp_init_lock(&wanderer->output_file_lock);
	omp_init_lock(&wanderer->reference_lock);
}

void
bwander_destroy(bam_wanderer_t *wanderer)
{
	int i;
	bam_region_t *region;
	linked_list_t *list;
	size_t list_l;

	assert(wanderer);
	assert(wanderer->regions_list);

	//Handle to list
	list = wanderer->regions_list;

	//Regions exists?
	if(wanderer->regions_list)
	{
		//for(i = 0; i < wanderer->regions_l; i++)
		list_l = linked_list_size(list);
		for(i = 0; i < list_l; i++)
		{
			//Get region
			region = linked_list_get(i, list);
			breg_destroy(region, 1);
			free(region);
		}
		linked_list_free(list, NULL);
	}

	//Destroy lock
	omp_destroy_lock(&wanderer->regions_lock);
	omp_destroy_lock(&wanderer->output_file_lock);
	omp_destroy_lock(&wanderer->reference_lock);
}

void 
bwander_configure(bam_wanderer_t *wanderer, bam_file_t *in_file, bam_file_t *out_file, genome_t *reference, wanderer_function wf, processor_function pf)
{
	assert(wanderer);
	assert(in_file);
	assert(wf);
	assert(pf);

	//Set I/O
	wanderer->input_file = in_file;
	omp_set_lock(&wanderer->output_file_lock);
	wanderer->output_file = out_file;		
	omp_unset_lock(&wanderer->output_file_lock);
	omp_set_lock(&wanderer->reference_lock);
	wanderer->reference = reference;
	omp_unset_lock(&wanderer->reference_lock);

	//Set functions
	wanderer->wander_f = wf;
	wanderer->processing_f = pf;

	//Logging
	LOG_INFO("Wanderer configured\n");
}

static int
bwander_run_sequential(bam_wanderer_t *wanderer)
{
	int i, err;
	size_t reads, reads_to_write;
	double times;
	bam_region_t *region;

	err = WANDER_REGION_CHANGED;
	reads = 0;
	while(err)
	{
		//Create new current region
		region = (bam_region_t *)malloc(sizeof(bam_region_t));
		breg_init(region);

		//Fill region
#ifdef D_TIME_DEBUG
		times = omp_get_wtime();
#endif
		err = bwander_obtain_region(wanderer, region);
#ifdef D_TIME_DEBUG
		times = omp_get_wtime() - times;
		if(region->size != 0)
			time_add_time_slot(D_FWORK_READ, TIME_GLOBAL_STATS, times / (double)region->size);
#endif
		if(err)
		{
			if(err == WANDER_REGION_CHANGED || err == WANDER_READ_EOF)
			{
				//Add region to wanderer regions
				bwander_region_insert(wanderer, region);

#ifdef D_TIME_DEBUG
				times = omp_get_wtime();
#endif
				//Process region
				wanderer->processing_f(wanderer, region);
#ifdef D_TIME_DEBUG
				times = omp_get_wtime() - times;
				if(region->size != 0)
				{
					time_add_time_slot(D_FWORK_PROC,  TIME_GLOBAL_STATS, times / (double)region->size);
					time_add_time_slot(D_FWORK_PROC_FUNC, TIME_GLOBAL_STATS, times / (double)region->size);
				}
				times = omp_get_wtime();
#endif

				reads_to_write = region->size;
				reads += reads_to_write;
				printf("Reads processed: %d\r", reads);

				//Write region
				breg_write_n(region, reads_to_write, wanderer->output_file);

				//Remove region from list
				linked_list_remove(region, wanderer->regions_list);

				//Free region
				breg_destroy(region, 1);
				free(region);

#ifdef D_TIME_DEBUG
				times = omp_get_wtime() - times;
				if(reads_to_write != 0)
					time_add_time_slot(D_FWORK_WRITE, TIME_GLOBAL_STATS, times / (double)reads_to_write);
#endif

				//End readings
				if(err == WANDER_READ_EOF)
					 break;
			}
			else
			{
				LOG_FATAL_F("Failed to read next region, error code: %d\n", err);
				break;
			}
		}
		else
		{
			//No more regions, end loop
			LOG_INFO("No more regions to read");
		}
	}

	printf("\n");
	return err;
}

static  int
bwander_run_threaded(bam_wanderer_t *wanderer)
{
	int err;
	bam_region_t *region;
	linked_list_t *regions;
	double times;

	omp_lock_t end_condition_lock;
	int end_condition;

	omp_lock_t reads_lock;
	size_t reads;
	size_t reads_to_write;

	//Init lock
	omp_init_lock(&end_condition_lock);
	omp_init_lock(&reads_lock);
	//#pragma omp parallel private(err, region, regions, times, reads_to_write)
	{
		//#pragma omp single
		{
			printf("Running in multithreading mode with %d threads\n", omp_get_max_threads());
			end_condition = 1;
			reads = 0;
		}

		#pragma omp parallel sections private(err, region, regions, times, reads_to_write)
		{
			//Region read
			#pragma omp section
			{
				regions = wanderer->regions_list;
				while(1)
				{
					//Create new current region
					region = (bam_region_t *)malloc(sizeof(bam_region_t));
					breg_init(region);

					//Fill region
#ifdef D_TIME_DEBUG
					times = omp_get_wtime();
#endif
					err = bwander_obtain_region(wanderer, region);
#ifdef D_TIME_DEBUG
					times = omp_get_wtime() - times;
					omp_set_lock(&region->lock);
					if(region->size != 0)
						time_add_time_slot(D_FWORK_READ, TIME_GLOBAL_STATS, times / (double)region->size);
					omp_unset_lock(&region->lock);
#endif
					if(err)
					{
						if(err == WANDER_REGION_CHANGED || err == WANDER_READ_EOF)
						{
							//Until process, this region cant be writed
							omp_test_lock(&region->write_lock);

							//Add region to wanderer regions
							bwander_region_insert(wanderer, region);

							#pragma omp task untied firstprivate(region) private(err)
							{
								double aux_time;

								//Process region
								omp_set_lock(&region->lock);
#ifdef D_TIME_DEBUG
								times = omp_get_wtime();
#endif
								wanderer->processing_f(wanderer, region);
#ifdef D_TIME_DEBUG
								times = omp_get_wtime() - times;
								if(region->size != 0)
									time_add_time_slot(D_FWORK_PROC_FUNC, TIME_GLOBAL_STATS, times / (double)region->size);
								aux_time = omp_get_wtime();
#endif
								omp_unset_lock(&region->lock);

								omp_set_lock(&reads_lock);
								reads += region->size;
								printf("Reads processed: %d\r", reads);
								omp_unset_lock(&reads_lock);

#ifdef D_TIME_DEBUG
								aux_time = omp_get_wtime() - aux_time;
								omp_set_lock(&region->lock);
								if(region->size != 0)
									time_add_time_slot(D_FWORK_PROC, TIME_GLOBAL_STATS, (times + aux_time) / (double)region->size);
								omp_unset_lock(&region->lock);
#endif

								//Set this region as writable
								omp_unset_lock(&region->write_lock);
							}

							//End readings
							if(err == WANDER_READ_EOF)
								 break;
						}
						else
						{
							LOG_FATAL_F("Failed to read next region, error code: %d\n", err);
							break;
						}
					}
					else
					{
						//No more regions, end loop
						LOG_INFO("No more regions to read");
						break;
					}
				}

				omp_set_lock(&end_condition_lock);
				end_condition = 0;
				omp_unset_lock(&end_condition_lock);
				//LOG_WARN("Read thread exit\n");
			}//End read section

			//Write section
			#pragma omp section
			{
				regions = wanderer->regions_list;
				omp_set_lock(&end_condition_lock);
				while(end_condition || linked_list_size(regions) > 0)
				{
					omp_unset_lock(&end_condition_lock);
#ifdef D_TIME_DEBUG
					times = omp_get_wtime();
#endif

					//Get next region
					omp_set_lock(&wanderer->regions_lock);
					region = linked_list_get_first(regions);
					omp_unset_lock(&wanderer->regions_lock);
					if(region == NULL)
					{
						omp_set_lock(&end_condition_lock);
						continue;
					}

					//Wait region to be writable
					omp_set_lock(&region->write_lock);

					//Write region
					omp_set_lock(&wanderer->output_file_lock);
					reads_to_write = region->size;
					breg_write_n(region, reads_to_write, wanderer->output_file);
					omp_unset_lock(&wanderer->output_file_lock);

					//Remove from list
					omp_set_lock(&wanderer->regions_lock);
					if(linked_list_size(regions) == 1)	//Possible bug?
						linked_list_clear(regions, NULL);
					else
						linked_list_remove_first(regions);

					//Signal read section if regions list is full
					if(linked_list_size(regions) < (WANDERER_REGIONS_MAX / 2) )
						omp_unset_lock(&wanderer->free_slots);

					omp_unset_lock(&wanderer->regions_lock);

#ifdef D_TIME_DEBUG
					times = omp_get_wtime() - times;
					omp_set_lock(&region->lock);
					if(reads_to_write != 0)
						time_add_time_slot(D_FWORK_WRITE, TIME_GLOBAL_STATS, times / (double)reads_to_write);
					omp_unset_lock(&region->lock);
#endif

					//Free region
					breg_destroy(region, 1);
					free(region);

					omp_set_lock(&end_condition_lock);
				}
				omp_unset_lock(&end_condition_lock);

				//LOG_WARN("Write thread exit\n");
			}//End write section

		}//End sections

	}//End parallel

	//Free
	omp_destroy_lock(&end_condition_lock);

	return NO_ERROR;
}

int
bwander_run(bam_wanderer_t *wanderer)
{
	int err;
	double times;
	bam1_t *read;

	assert(wanderer);
	assert(wanderer->input_file);
	assert(wanderer->regions_list);

	//Logging
	LOG_INFO("Wanderer is now running\n");

	if(omp_get_max_threads() > 1)
	{
		//Run in multithreaded mode
		err = bwander_run_threaded(wanderer);
	}
	else
	{
		//Run in sequential mode
		err = bwander_run_sequential(wanderer);
	}

	//Logging
	LOG_INFO("Wanderer SUCCESS!\n");

	return err;
}

/**
 * STATIC FUNCTIONS
 */

static inline int
bwander_region_insert(bam_wanderer_t *wanderer, bam_region_t *region)
{
	linked_list_t *list;
	size_t list_l;

	assert(wanderer);
	assert(region);

	omp_set_lock(&wanderer->regions_lock);

	//List handle
	list = wanderer->regions_list;
	list_l = linked_list_size(list);

	if(list_l >= WANDERER_REGIONS_MAX)
	{
		omp_unset_lock(&wanderer->regions_lock);

		//Wait for free slots
		if(!omp_test_lock(&wanderer->free_slots))
		{
			if(omp_get_num_threads() == 2)
			{
				#pragma omp taskwait	//Force processing
			}
			omp_set_lock(&wanderer->free_slots);
		}

		//LOG_FATAL_F("Not enough region slots, current: %d\n", list_l);
	}

	//This lock must be always locked until regions buffer have regions free
	omp_test_lock(&wanderer->free_slots);

	//Add region to list
	linked_list_insert_last(region, list);

	omp_set_lock(&region->lock);
	LOG_INFO_F("Inserting region %d:%d-%d with %d reads\n",
				region->chrom + 1, region->init_pos + 1,
				region->end_pos + 1, region->size);
	LOG_INFO_F("Regions to process %d\n", linked_list_size(list));
	omp_unset_lock(&region->lock);

	omp_unset_lock(&wanderer->regions_lock);

	return NO_ERROR;
}

static inline int
bwander_obtain_region(bam_wanderer_t *wanderer, bam_region_t *region)
{
	int i, err, bytes;
	bam1_t *read;
	double times;

	//Get first read
	if(last_read != NULL)
	{
		read = last_read;
		bytes = last_read_bytes;
		last_read = NULL;
	}
	else
	{
		//Get first read from file
		read = bam_init1();
		assert(read);
		bytes = bam_read1(wanderer->input_file->bam_fd, read);
	}

	//Iterate reads
	while(bytes > 0)
	{
		//Wander this read
		omp_set_lock(&region->lock);
#ifdef D_TIME_DEBUG
			times = omp_get_wtime();
#endif
		err = wanderer->wander_f(wanderer, region, read);
#ifdef D_TIME_DEBUG
			times = omp_get_wtime() - times;
			time_add_time_slot(D_FWORK_WANDER_FUNC, TIME_GLOBAL_STATS, times);
#endif
		omp_unset_lock(&region->lock);
		switch(err)
		{
		case WANDER_READ_FILTERED:
			//This read dont pass the filters
		case NO_ERROR:
			//Add read to region
			omp_set_lock(&region->lock);
			region->reads[region->size]  = read;
			region->size++;

			//Region is full?
			if(region->size >= region->max_size)
			{
				omp_unset_lock(&region->lock);
				return WANDER_REGION_CHANGED;
			}
			omp_unset_lock(&region->lock);

			//Get next read from file
			read = bam_init1();
			assert(read);
			bytes = bam_read1(wanderer->input_file->bam_fd, read);
			break;

		case WANDER_REGION_CHANGED:
			//The region have changed
			last_read = read;
			last_read_bytes = bytes;
			return err;

		default:
			//Unknown error
			LOG_ERROR_F("Wanderer fails with error code: %d\n", err);
			return err;
		}
	}

	//Check read error
	if(bytes <= 0)
	{
		//Destroy bam
		bam_destroy1(read);

		//End of file
		if(bytes == -1)
		{
			LOG_INFO("EOF\n");
			return WANDER_READ_EOF;
		}
		else
		{
			LOG_INFO("TRUNCATED\n");
			return WANDER_READ_TRUNCATED;
		}
	}

	return NO_ERROR;
}
