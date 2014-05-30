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

static INLINE int bwander_obtain_region(bam_wanderer_t *wanderer, bam_region_t *current_region);

void
bwander_init(bam_wanderer_t *wanderer)
{
	int i;
	bam_region_t *region;

	assert(wanderer);

	//Set all to zero
	memset(wanderer, 0, sizeof(bam_wanderer_t));

	//Create regions
	wanderer->regions = (bam_region_t **)malloc(WANDERER_REGIONS_MAX * sizeof(bam_region_t *));
	/*for(i = 0; i < WANDERER_REGIONS_MAX; i++)
	{
		region = wanderer->regions[i];
		breg_init(region);
	}*/

	//Init locks
	omp_init_lock(&wanderer->regions_lock);
	omp_init_lock(&wanderer->output_file_lock);
}

void
bwander_destroy(bam_wanderer_t *wanderer)
{
	int i;
	bam_region_t *region;

	assert(wanderer);

	//Regions exists?
	if(wanderer->regions)
	{
		for(i = 0; i < wanderer->regions_l; i++)
		{
			//Get region
			region = wanderer->regions[i];
			breg_destroy(region, 0);
			free(region);
		}
		free(wanderer->regions);
	}

	//Destroy lock
	omp_destroy_lock(&wanderer->regions_lock);
	omp_destroy_lock(&wanderer->output_file_lock);
}

void 
bwander_configure(bam_wanderer_t *wanderer, bam_file_t *in_file, bam_file_t *out_file, wanderer_function wf, processor_function pf)
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

	//Set functions
	wanderer->wander_f = wf;
	wanderer->processing_f = pf;

	//Logging
	LOG_INFO("Wanderer configured\n");
}

int
bwander_run(bam_wanderer_t *wanderer)
{
	int i, err;
	bam_region_t *read_region;
	bam_region_t *process_region;
	bam_region_t *write_region;
	size_t total;
	int bytes;
	int while_cond = 1;

	bam1_t *read;

	assert(wanderer);
	assert(wanderer->input_file);
	assert(wanderer->regions);

	//Logging
	LOG_INFO("Wanderer is now running\n");

	#pragma omp parallel //sections
	{
		//Region read
		#pragma omp single
		{
			while(1)
			{
				//Create new current region
				read_region = (bam_region_t *)malloc(sizeof(bam_region_t));
				breg_init(read_region);

				//Fill region
				err = bwander_obtain_region(wanderer, read_region);
				if(err)
				{
					if(err == WANDER_REGION_CHANGED)
					{
						//Add region to wanderer regions
						bwander_region_insert(wanderer, read_region);
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
		}//Read section

		//Region process section
		/*#pragma omp section
		{
			while(wanderer->regions_l < 6)
			{
				//Iterate regions
				for(i = 0; i < wanderer->regions_l; i++)
				{
					//Get next region
					region = wanderer->regions[i];

					//Process region
					wanderer->processing_f(wanderer, region);
				}
			}
		} //End process section*/

		//Region writing (in region order)
		/*#pragma omp section
		{
			//Iterate regions
			for(i = 0; i < wanderer->regions_l; i++)
			{
				//Get next region
				region = wanderer->regions[i];

				//This region is processed?
				//TODO if()
				{
					//Write region to disk
					breg_write_n(region, region->size, wanderer->output_file);
					wanderer->regions_l--;
				}
				//else break;
			}
		} //End write section*/

	}//End parallel

	//Logging
	LOG_INFO("Wanderer SUCCESS!\n");

	return NO_ERROR;
}

EXTERNC INLINE int
filter_read(bam1_t *read, uint8_t filters)
{
	assert(read);

	//Filter read
	if(filters != 0)
	{
		if(filters & FILTER_ZERO_QUAL)
		{
			if(read->core.qual == 0)
				return 1;
		}

		if(filters & FILTER_DIFF_MATE_CHROM)
		{
			if(read->core.tid != read->core.mtid)
				return 1;
		}

		if(filters & FILTER_NO_CIGAR)
		{
			if(read->core.n_cigar == 0)
				return 1;
		}

		if(filters & FILTER_DEF_MASK)
		{
			if(read->core.flag & BAM_DEF_MASK)
				return 1;
		}
	}

	return NO_ERROR;
}

EXTERNC INLINE int
bwander_region_insert(bam_wanderer_t *wanderer, bam_region_t *region)
{
	assert(wanderer);
	assert(region);

	omp_set_lock(&wanderer->regions_lock);

	if(wanderer->regions_l >= WANDERER_REGIONS_MAX)
	{
		omp_unset_lock(&wanderer->regions_lock);
		LOG_FATAL("Not enough region slots\n");
	}

	wanderer->regions[wanderer->regions_l] = region;
	wanderer->regions_l++;

	omp_set_lock(&region->lock);
	LOG_INFO_F("Inserting region %d:%d-%d with %d reads\n",
				region->chrom + 1, region->init_pos + 1,
				region->end_pos + 1, region->size);
	LOG_INFO_F("Regions to process %d\n", wanderer->regions_l);
	omp_unset_lock(&region->lock);

	omp_unset_lock(&wanderer->regions_lock);

	#pragma omp task
	{
		omp_set_lock(&region->lock);
		wanderer->processing_f(wanderer, region);
		omp_unset_lock(&region->lock);

		//Write region
		omp_set_lock(&wanderer->output_file_lock);
		breg_write_n(region, region->size, wanderer->output_file);
		omp_unset_lock(&wanderer->output_file_lock);
	}

	return NO_ERROR;
}

static INLINE int
bwander_obtain_region(bam_wanderer_t *wanderer, bam_region_t *region)
{
	int i, err, bytes;
	bam1_t *read;

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
		err = wanderer->wander_f(wanderer, region, read);
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

/**
 * REGISTER WINDOW
 */
/*int
bwander_window_register(bam_wanderer_t *wanderer, bam_region_window_t *window)
{
	bam_region_window_t *dest_window;

	assert(wanderer);
	assert(wanderer->current_region);
	assert(window);

	//Logging
	LOG_INFO_F("Registering window %d:%d-%d with %d reads\n",
			window->region->chrom + 1, window->init_pos + 1,
			window->end_pos + 1, window->size);

	//Slots free?
	if(wanderer->windows_l >= WANDERER_WINDOWS_MAX)
	{
		LOG_WARN("Wanderer bam_region_window_t buffer is full\n");
		return WANDER_WINDOW_BUFFER_FULL;
	}

	//Duplicate and register
	dest_window = &wanderer->windows[wanderer->windows_l];
	dest_window->init_pos = window->init_pos;
	dest_window->end_pos = window->end_pos;
	dest_window->filter_flags = window->filter_flags;
	dest_window->size = window->size;
	dest_window->region = window->region;
	memcpy(dest_window->filter_reads, window->filter_reads, window->size * sizeof(bam1_t *));

	//Increase windows counter
	wanderer->windows_l++;

	return NO_ERROR;
}*/

/*void
bwander_window_clear(bam_wanderer_t *wanderer)
{
	int i;
	bam_region_window_t *window;

	assert(wanderer);
	assert(wanderer->current_region);

	//Logging
	LOG_INFO_F("Clearing %d windows\n", wanderer->windows_l);

	//Iterate windows
	for(i = 0; i < wanderer->windows_l; i++)
	{
		//Get next window
		window = &wanderer->windows[i];

		//Clear window
		breg_window_clear(window);
	}

	//Set size to zero
	wanderer->windows_l = 0;
}*/

