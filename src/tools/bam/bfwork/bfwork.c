/**
* Copyright (C) 2013 Raúl Moreno Galdón
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <errno.h>
#include <string.h>

#include "bfwork.h"

/**
 * STATIC VARS
 */
//bfwork_obtain_region
static bam1_t *last_read = NULL;
static int last_read_bytes = 0;

/**
 * STATIC FUNCTIONS
 */

/**
 * \brief PRIVATE. Insert a region in a context.
 *
 * \param[in] fwork Target framework.
 * \param[in] region Region to insert.
 */
static int bfwork_region_insert(bam_fwork_t *fwork, bam_region_t *region);

/**
 * \brief PRIVATE. Wander for a region.
 *
 * \param[in] fwork Target framework.
 * \param[out] region Region to fill using current context wandering function.
 */
static int bfwork_obtain_region(bam_fwork_t *fwork, bam_region_t *current_region);

/**
 * Initialize empty BAM framework data structure.
 */
void
bfwork_init(bam_fwork_t *fwork)
{
	int i, threads;
	bam_region_t *region;

	assert(fwork);

	//Set all to zero
	memset(fwork, 0, sizeof(bam_fwork_t));

	//Create regions
	fwork->regions_list = linked_list_new(COLLECTION_MODE_SYNCHRONIZED);

	//Init locks
	omp_init_lock(&fwork->regions_lock);
	omp_init_lock(&fwork->free_slots);
	omp_init_lock(&fwork->output_file_lock);
	omp_init_lock(&fwork->reference_lock);
}

/**
 * Destroy BAM framework data structure.
 */
void
bfwork_destroy(bam_fwork_t *fwork)
{
	int i;
	bam_region_t *region;
	linked_list_t *list;
	size_t list_l;

	assert(fwork);
	assert(fwork->regions_list);

	//Handle to list
	list = fwork->regions_list;

	//Regions exists?
	if(fwork->regions_list)
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
	omp_destroy_lock(&fwork->regions_lock);
	omp_destroy_lock(&fwork->output_file_lock);
	omp_destroy_lock(&fwork->reference_lock);
}

/**
 * Configure BAM framework data structure for wandering.
 */
int
bfwork_configure(bam_fwork_t *fwork, const char *in_file, const char *out_file, const char *reference, bfwork_context_t *context)
{
	assert(fwork);
	assert(in_file);

	//Set I/O
	fwork->input_file_str = (char *)in_file;
	omp_set_lock(&fwork->output_file_lock);
	fwork->output_file_str = (char *)out_file;
	omp_unset_lock(&fwork->output_file_lock);
	omp_set_lock(&fwork->reference_lock);
	fwork->reference_str = (char *)reference;
	omp_unset_lock(&fwork->reference_lock);

	//Initial context have been specified?
	if(context)
	{
		//Set context
		fwork->context = context;

		//Add to context list
		fwork->v_context[0] = context;
		fwork->v_context_l = 1;
	}

	//Logging
	LOG_INFO("Framework configured\n");

	return NO_ERROR;
}

/**
 * Add additional context to execute in framework.
 */
int
bfwork_add_context(bam_fwork_t *fwork, bfwork_context_t *context, uint8_t flags)
{
	assert(fwork);
	assert(context);
	assert(flags != 0);

	//What kind of execution queue have this context?
	if(flags & FWORK_CONTEXT_PARALLEL)
	{
		//Not supported yet
		LOG_WARN("FWORK_CONTEXT_QUEUE_PARALLEL is not supported yet, changing to FWORK_CONTEXT_QUEUE_SEQUENTIAL\n");
		flags = FWORK_CONTEXT_SEQUENTIAL;
	}
	if(flags & FWORK_CONTEXT_SEQUENTIAL)
	{
		//Add to sequential execution queue
		fwork->v_context[fwork->v_context_l] = context;
		fwork->v_context_l++;
	}
	else
	{
		LOG_FATAL_F("Trying to add a context with not known flags: %x\n", flags);
	}

	return NO_ERROR;
}

static int
bfwork_run_sequential(bam_fwork_t *fwork)
{
	int i, err;
	size_t reads, reads_to_write;
	double times;
	bam_region_t *region;

	//Context
	bfwork_context_t *context;
	size_t pf_l;

	err = WANDER_REGION_CHANGED;
	reads = 0;
	context = fwork->context;
	pf_l = context->processing_f_l;
	while(err)
	{
		//Create new current region
		region = (bam_region_t *)malloc(sizeof(bam_region_t));
		breg_init(region);

		//Fill region
#ifdef D_TIME_DEBUG
		times = omp_get_wtime();
#endif
		err = bfwork_obtain_region(fwork, region);
#ifdef D_TIME_DEBUG
		times = omp_get_wtime() - times;
		if(region->size != 0)
			if(context->time_stats)
			time_add_time_slot(D_FWORK_READ, context->time_stats, times / (double)region->size);
#endif
		if(err)
		{
			if(err == WANDER_REGION_CHANGED || err == WANDER_READ_EOF)
			{
				//Add region to framework regions
				bfwork_region_insert(fwork, region);

#ifdef D_TIME_DEBUG
				times = omp_get_wtime();
#endif
				//Process region
				for(i = 0; i < pf_l; i++)
				{
					context->processing_f[i](fwork, region);
				}
#ifdef D_TIME_DEBUG
				times = omp_get_wtime() - times;
				if(context->time_stats)
				if(region->size != 0)
				{
					time_add_time_slot(D_FWORK_PROC,  context->time_stats, times / (double)region->size);
					time_add_time_slot(D_FWORK_PROC_FUNC, context->time_stats, times / (double)region->size);
				}
				times = omp_get_wtime();
#endif

				reads_to_write = region->size;
				reads += reads_to_write;
				printf("Reads processed: %d\r", reads);

				//Write region
				breg_write_n(region, reads_to_write, fwork->output_file);

				//Remove region from list
				linked_list_remove(region, fwork->regions_list);

				//Free region
				breg_destroy(region, 1);
				free(region);

#ifdef D_TIME_DEBUG
				times = omp_get_wtime() - times;
				if(context->time_stats)
				if(reads_to_write != 0)
					time_add_time_slot(D_FWORK_WRITE, context->time_stats, times / (double)reads_to_write);
#endif

				//End readings
				if(err == WANDER_READ_EOF)
					 break;
			}
			else
			{
				if(err == WANDER_READ_TRUNCATED)
				{
					LOG_WARN("Readed truncated read\n");
				}
				else
				{
					LOG_FATAL_F("Failed to read next region, error code: %d\n", err);
				}
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
bfwork_run_threaded(bam_fwork_t *fwork)
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
				regions = fwork->regions_list;
				while(1)
				{
					//Create new current region
					region = (bam_region_t *)malloc(sizeof(bam_region_t));
					breg_init(region);

					//Fill region
#ifdef D_TIME_DEBUG
					times = omp_get_wtime();
#endif
					err = bfwork_obtain_region(fwork, region);
#ifdef D_TIME_DEBUG
					times = omp_get_wtime() - times;
					omp_set_lock(&region->lock);
					if(fwork->context->time_stats)
					if(region->size != 0)
						time_add_time_slot(D_FWORK_READ, fwork->context->time_stats, times / (double)region->size);
					omp_unset_lock(&region->lock);
#endif
					if(err)
					{
						if(err == WANDER_REGION_CHANGED || err == WANDER_READ_EOF)
						{
							//Until process, this region cant be writed
							omp_test_lock(&region->write_lock);

							//Add region to framework regions
							bfwork_region_insert(fwork, region);

							#pragma omp task untied firstprivate(region) private(err)
							{
								int i;
								size_t pf_l;
								double aux_time;

								//Process region
								omp_set_lock(&region->lock);
#ifdef D_TIME_DEBUG
								times = omp_get_wtime();
#endif
								//Process region
								pf_l = fwork->context->processing_f_l;
								for(i = 0; i < pf_l; i++)
								{
									fwork->context->processing_f[i](fwork, region);
								}
#ifdef D_TIME_DEBUG
								times = omp_get_wtime() - times;
								if(fwork->context->time_stats)
								if(region->size != 0)
									time_add_time_slot(D_FWORK_PROC_FUNC, fwork->context->time_stats, times / (double)region->size);
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
								if(fwork->context->time_stats)
								if(region->size != 0)
									time_add_time_slot(D_FWORK_PROC, fwork->context->time_stats, (times + aux_time) / (double)region->size);
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
							if(err == WANDER_READ_TRUNCATED)
							{
								LOG_WARN("Readed truncated read\n");
							}
							else
							{
								LOG_FATAL_F("Failed to read next region, error code: %d\n", err);
							}
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
				regions = fwork->regions_list;
				omp_set_lock(&end_condition_lock);
				while(end_condition || linked_list_size(regions) > 0)
				{
					omp_unset_lock(&end_condition_lock);
#ifdef D_TIME_DEBUG
					times = omp_get_wtime();
#endif

					//Get next region
					omp_set_lock(&fwork->regions_lock);
					region = linked_list_get_first(regions);
					omp_unset_lock(&fwork->regions_lock);
					if(region == NULL)
					{
						omp_set_lock(&end_condition_lock);
						continue;
					}

					//Wait region to be writable
					omp_set_lock(&region->write_lock);

					//Write region
					omp_set_lock(&fwork->output_file_lock);
					reads_to_write = region->size;
					breg_write_n(region, reads_to_write, fwork->output_file);
					omp_unset_lock(&fwork->output_file_lock);

					//Remove from list
					omp_set_lock(&fwork->regions_lock);
					if(linked_list_size(regions) == 1)	//Possible bug?
						linked_list_clear(regions, NULL);
					else
						linked_list_remove_first(regions);

					//Signal read section if regions list is full
					if(linked_list_size(regions) < (FWORK_REGIONS_MAX / 2) )
						omp_unset_lock(&fwork->free_slots);

					omp_unset_lock(&fwork->regions_lock);

#ifdef D_TIME_DEBUG
					times = omp_get_wtime() - times;
					omp_set_lock(&region->lock);
					if(fwork->context->time_stats)
					if(reads_to_write != 0)
						time_add_time_slot(D_FWORK_WRITE, fwork->context->time_stats, times / (double)reads_to_write);
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

	//Lineskip
	printf("\n");

	//Free
	omp_destroy_lock(&end_condition_lock);

	return NO_ERROR;
}

/**
 * Run framework contexts.
 */
int
bfwork_run(bam_fwork_t *fwork)
{
	int err, c;
	double times;
	bam1_t *read;

	//Reference
	char *ref_path;
	char *ref_name;

	assert(fwork);
	assert(fwork->input_file_str);
	assert(fwork->regions_list);

	printf("============== BEGIN RUN  ==============\n");

	//Check if contexts present
	if(fwork->v_context_l == 0)
	{
		LOG_WARN("No contexts have been specified to run!\n");
		printf("============== END RUN  ==============\n\n");

		return NO_ERROR;
	}

	//Open reference
	if(fwork->reference_str)
	{
		//Obtain reference filename and dirpath from full path
		ref_path = strdup(fwork->reference_str);
		ref_path = dirname(ref_path);
		ref_name = strrchr(fwork->reference_str, '/');
		printf("Reference path: %s\n", ref_path);
		printf("Reference name: %s\n", ref_name);
		printf("Opening reference genome from \"%s%s\" ...\n", ref_path, ref_name);
		fwork->reference = genome_new(ref_name, ref_path, BWT_MODE);
		assert(fwork->reference);
		printf("Reference opened!...\n");
	}

	printf("--------------------------------------\n");

	for(c = 0; c < fwork->v_context_l; c++)
	{
		//Select next context
		fwork->context = fwork->v_context[c];
		assert(fwork->context);

#ifdef D_TIME_DEBUG
		times = omp_get_wtime();
#endif

		//Open input bam
		{
			//If last context had no output
			if(!fwork->last_temp_file_str)
			{
				//Open initial input file
				printf("Opening BAM from \"%s\" ...\n", fwork->input_file_str);
				fwork->input_file = bam_fopen(fwork->input_file_str);
				assert(fwork->input_file);
				printf("BAM opened!...\n");
			}
			else
			{
				//Open last context output
				printf("Opening intermediate BAM from \"%s\" ...\n", fwork->last_temp_file_str);
				fwork->input_file = bam_fopen(fwork->last_temp_file_str);
				assert(fwork->input_file);
				printf("Intermediate BAM opened!...\n");
			}
		}

		//Create new output bam if last context
		fwork->erase_tmp = 0;
		if(c == fwork->v_context_l - 1)
		{
			if(fwork->output_file_str)
			{
				if(fwork->context->output_file_str == NULL)
				{
					//Allocate
					fwork->context->output_file_str = (char *)malloc(256 * sizeof(char));
				}

				//Set final output
				strncpy(fwork->context->output_file_str, fwork->output_file_str, 256);
			}
		}
		else
		{
			//Temporary file?
			if(fwork->context->output_temp)
			{
				//Allocate
				fwork->context->output_file_str = malloc(256 * sizeof(char));

				//Set output temp path
				sprintf(fwork->context->output_file_str, "/tmp/bfwork_tmp_%d.tmp", c);
				fwork->erase_tmp = 1;
			}
		}

		//Create new temporary bam if context have output
		if(fwork->context->output_file_str != NULL)
		{
			printf("Creating new intermediate bam file in \"%s\"...\n", fwork->context->output_file_str);
			fwork->output_file = bam_fopen_mode(fwork->context->output_file_str, fwork->input_file->bam_header_p, "w");
			assert(fwork->output_file);
			bam_fwrite_header(fwork->output_file->bam_header_p, fwork->output_file);
			fwork->output_file->bam_header_p = NULL;
			printf("New intermediate BAM initialized!...\n");
		}

#ifdef D_TIME_DEBUG
		times = omp_get_wtime() - times;
		if(fwork->context->time_stats)
			time_add_time_slot(D_FWORK_INIT, fwork->context->time_stats, times);
#endif

		//Logging
		if(fwork->context->tag != NULL)
		{
			printf("Context %s is now running\n", fwork->context->tag);
		}
		else
		{
			printf("Context %d is now running\n", c);
		}

#ifdef D_TIME_DEBUG
		times = omp_get_wtime();
#endif

		//Run this context
		if(omp_get_max_threads() > 1)
		{
			//Run in multithreaded mode
			err = bfwork_run_threaded(fwork);
		}
		else
		{
			//Run in sequential mode
			err = bfwork_run_sequential(fwork);
		}

		//Reduce needed?
		if(fwork->context->reduce != NULL && fwork->context->reduce_dest != NULL)
		{
			//Reduce into context reduce data
			bfwork_context_local_user_data_reduce(fwork->context, fwork->context->reduce_dest, fwork->context->reduce);
		}

#ifdef D_TIME_DEBUG
		times = omp_get_wtime() - times;
		if(fwork->context->time_stats)
			time_add_time_slot(D_FWORK_TOTAL, fwork->context->time_stats, times);
#endif

		//Close input BAM
		printf("\nClosing BAM file...\n");
		bam_fclose(fwork->input_file);
		fwork->input_file = NULL;
		printf("BAM closed.\n");

		//Close output file
		if(fwork->output_file != NULL)
		{
			printf("Closing \"%s\" BAM file...\n", fwork->output_file->filename);
			bam_fclose(fwork->output_file);
			fwork->output_file = NULL;
			printf("BAM closed.\n");
		}

		//Remove last temporary file
		if(fwork->last_temp_file_str != NULL && fwork->erase_tmp)
		{
			//Delete file
			printf("Deleting %s...\n", fwork->last_temp_file_str);
			remove(fwork->last_temp_file_str);
			fwork->last_temp_file_str = NULL;
			fwork->erase_tmp = 0;
		}

		//Set last file
		if(fwork->context->output_file_str)
		{
			//Set last temporary file if not the last context
			if(c < fwork->v_context_l - 1)
				fwork->last_temp_file_str = fwork->context->output_file_str;
		}

		//Logging
		printf("--------------------------------------\n");
		LOG_INFO("Context SUCCESS!\n");
	}

	//Remove last temporary file
	if(fwork->last_temp_file_str != NULL && fwork->erase_tmp)
	{
		//Delete file
		printf("Deleting %s...\n", fwork->last_temp_file_str);
		remove(fwork->last_temp_file_str);
	}

	//Close reference
	if(fwork->reference != NULL)
	{
		printf("\nClosing reference file...\n");
		genome_free(fwork->reference);
		printf("Reference closed.\n");
	}

	//Logging
	LOG_INFO("Framework SUCCESS!\n");

	printf("============== END RUN  ==============\n\n");

	return err;
}

/**
 * WANDERING CONTEXT OPERATIONS
 */

/**
 * Initialize empty BAM context data structure.
 */
void
bfwork_context_init(bfwork_context_t *context, wanderer_function wf, processor_function pf, reducer_function rf, void *reduce_dest)
{
	int threads;

	assert(context);
	assert(wf);
	assert(pf);

	//Set context to zeros
	memset(context, 0, sizeof(bfwork_context_t));

	//Create local data
	threads = omp_get_max_threads();
	context->local_user_data = (void **)malloc(threads * sizeof(void*));
	memset(context->local_user_data, 0, threads * sizeof(void*));

	//Assign functions
	context->wander_f = wf;
	context->processing_f[0] = pf;
	context->processing_f_l = 1;
	context->reduce = rf;
	context->reduce_dest = reduce_dest;

	//Init locks
	omp_init_lock(&context->user_data_lock);
}

/**
 * Destroy BAM context data structure.
 */
void
bfwork_context_destroy(bfwork_context_t *context)
{
	assert(context);

	//Free
	if(context->output_file_str)
		free(context->output_file_str);

	//Destroy locks
	omp_destroy_lock(&context->user_data_lock);
}

/**
 * Add additional processing function to this context.
 */
int
bfwork_context_add_proc(bfwork_context_t *context, processor_function pf)
{
	assert(context);
	assert(pf);

	//Check max functions
	if(context->processing_f_l >= FWORK_PROC_FUNC_MAX)
	{
		LOG_ERROR("Trying to add processor function, maximun number of processor functions reached\n");
		return WANDER_PROC_FUNC_FULL;
	}

	//Add function to list
	context->processing_f[context->processing_f_l] = pf;
	context->processing_f_l++;

	//Logging
	LOG_INFO("Added processor function\n");

	return NO_ERROR;
}

/**
 * Define if this context modify original input file for next contexts.
 */
int
bfwork_context_set_output(bfwork_context_t *context, const char *file_str)
{
	assert(context);

	if(context->output_file_str)
		free(context->output_file_str);

	//Set intermediate file name
	if(file_str)
	{
		//Allocate
		context->output_file_str = malloc(256 * sizeof(char));

		//Set to string file
		strncpy(context->output_file_str, file_str, 256);

		//Set temp
		context->output_temp = 0;
	}
	else
	{
		//Set to temp
		context->output_temp = 1;
	}

	return NO_ERROR;
}


/**
 * USER DATA
 */

/**
 * Set a pointer to be shared as a global user data among all threads in a context.
 */
int
bfwork_context_set_user_data(bfwork_context_t *context, void *user_data)
{
	assert(context);
	assert(user_data);

	//Set user data
	context->user_data = user_data;

	return NO_ERROR;
}

/**
 * TIMING
 */

/**
 * Initialize and activate timing of a context.
 */
int
bfwork_context_init_timing(bfwork_context_t *context, const char *tag, const char *path_folder)
{
	int err;
	char *sched;
	char filename[100];
	char intaux[20];
	char cwd[1024];

	assert(context);

	//Set tag
	context->tag = (char *) tag;

	//Create timing
	if(time_new_stats(20, &context->time_stats))
	{
		LOG_ERROR("Failed to initialize time stats\n");
	}

	//Get OMP schedule
	sched = getenv("OMP_SCHEDULE");

	//Target folder for stats?
	if(path_folder != NULL)
	{
		strcpy(filename, path_folder);

		//Create stats directory
		err = mkdir(filename, S_IRWXU);
		if(err)
		{
			err = errno;
			LOG_WARN_F("Failed to create stats directory \"%s\", error code: %s\n", filename, strerror(err));
		}

		strcat(filename, "/");
		if(tag)
		{
			strcat(filename,tag);
		}
		if(sched)
		{
			strcat(filename,"_");
			strcat(filename,sched);
		}
		else
		{
			printf("ERROR: Obtaining OMP_SCHEDULE environment value\n");
		}

		//strcat(filename,"_");
		//sprintf(intaux, "%d", MAX_BATCH_SIZE);
		//strcat(filename, intaux);
		strcat(filename, "_");
		sprintf(intaux, "%d", omp_get_max_threads());
		strcat(filename, intaux);
		strcat(filename, ".stats");

		//Set output file
		if(time_set_output_file(filename, context->time_stats))
		{
			LOG_ERROR_F("Failed to set timing file output to \"%s\"\n", filename);
		}
		else
		{
			printf("STATISTICS ACTIVATED, output file: %s\n\n", filename);
		}

	}

	return NO_ERROR;
}

/**
 * Destroy timing of a context.
 */
void
bfwork_context_destroy_timing(bfwork_context_t *context)
{
	assert(context);

	//Destroy time stats
	if(context->time_stats)
	{
		time_destroy_stats(&context->time_stats);
	}
}

/**
 * Output timings in standard output.
 */
int
bfwork_context_print_times(bfwork_context_t *context)
{
	//Print times
	double min, max, mean;

#ifdef D_TIME_DEBUG
	if(context->time_stats)
	{
		//Print time stats
		printf("----------------------------\nTime stats for %s: \n", context->tag != NULL ? context->tag : "context");

		printf("\n====== General times ======\n");
		time_get_mean_slot(D_FWORK_TOTAL, context->time_stats, &mean);
		time_get_min_slot(D_FWORK_TOTAL, context->time_stats, &min);
		time_get_max_slot(D_FWORK_TOTAL, context->time_stats, &max);
		printf("Total time to process -> %.2f s - min/max = %.2f/%.2f\n",
				mean, min, max);

		time_get_mean_slot(D_FWORK_INIT, context->time_stats, &mean);
		time_get_min_slot(D_FWORK_INIT, context->time_stats, &min);
		time_get_max_slot(D_FWORK_INIT, context->time_stats, &max);
		printf("Time used to initialize framework -> %.2f s - min/max = %.2f/%.2f\n",
				mean, min, max);

		printf("\n====== Wandering function ======\n");
		time_get_mean_slot(D_FWORK_WANDER_FUNC, context->time_stats, &mean);
		time_get_min_slot(D_FWORK_WANDER_FUNC, context->time_stats, &min);
		time_get_max_slot(D_FWORK_WANDER_FUNC, context->time_stats, &max);
		printf("Time of wandering function per alignment (inside framework read time) -> %.2f us - min/max = %.2f/%.2f\n",
					mean*1000000.0, min*1000000.0, max*1000000.0);

		printf("\n====== Processing function ======\n");
		time_get_mean_slot(D_FWORK_PROC_FUNC, context->time_stats, &mean);
		time_get_min_slot(D_FWORK_PROC_FUNC, context->time_stats, &min);
		time_get_max_slot(D_FWORK_PROC_FUNC, context->time_stats, &max);
		printf("Time of processing function per alignment (inside framework process time) -> %.2f us - min/max = %.2f/%.2f\n",
				mean*1000000.0, min*1000000.0, max*1000000.0);

		printf("\n====== Framework ======\n");

		time_get_mean_slot(D_FWORK_PROC, context->time_stats, &mean);
		time_get_min_slot(D_FWORK_PROC, context->time_stats, &min);
		time_get_max_slot(D_FWORK_PROC, context->time_stats, &max);
		printf("Time used for process per alignment -> %.2f us - min/max = %.2f/%.2f\n",
				mean*1000000.0, min*1000000.0, max*1000000.0);

		time_get_mean_slot(D_FWORK_READ, context->time_stats, &mean);
		time_get_min_slot(D_FWORK_READ, context->time_stats, &min);
		time_get_max_slot(D_FWORK_READ, context->time_stats, &max);
		printf("Time used for read per alignment -> %.2f us - min/max = %.2f/%.2f\n",
				mean*1000000.0, min*1000000.0, max*1000000.0);

		time_get_mean_slot(D_FWORK_WRITE, context->time_stats, &mean);
		time_get_min_slot(D_FWORK_WRITE, context->time_stats, &min);
		time_get_max_slot(D_FWORK_WRITE, context->time_stats, &max);
		printf("Time used for write per alignment -> %.2f us - min/max = %.2f/%.2f\n",
				mean*1000000.0, min*1000000.0, max*1000000.0);

		printf("-------------\n\n");
	}
#endif

	return NO_ERROR;
}

/**
 * STATIC FUNCTIONS
 */

/**
 * PRIVATE. Insert a region in a context.
 */
static inline int
bfwork_region_insert(bam_fwork_t *fwork, bam_region_t *region)
{
	linked_list_t *list;
	size_t list_l;

	assert(fwork);
	assert(region);

	omp_set_lock(&fwork->regions_lock);

	//List handle
	list = fwork->regions_list;
	list_l = linked_list_size(list);

	if(list_l >= FWORK_REGIONS_MAX)
	{
		omp_unset_lock(&fwork->regions_lock);

		//Wait for free slots
		if(!omp_test_lock(&fwork->free_slots))
		{
			if(omp_get_num_threads() == 2)
			{
				#pragma omp taskwait	//Force processing
			}
			omp_set_lock(&fwork->free_slots);
		}

		//LOG_FATAL_F("Not enough region slots, current: %d\n", list_l);
	}

	//This lock must be always locked until regions buffer have regions free
	omp_test_lock(&fwork->free_slots);

	//Add region to list
	linked_list_insert_last(region, list);

	omp_set_lock(&region->lock);
	LOG_INFO_F("Inserting region %d:%lu-%lu with %d reads\n",
				region->chrom + 1, region->init_pos + 1,
				region->end_pos + 1, region->size);
	LOG_INFO_F("Regions to process %lu\n", linked_list_size(list));
	omp_unset_lock(&region->lock);

	omp_unset_lock(&fwork->regions_lock);

	return NO_ERROR;
}

/**
 * PRIVATE. Wander for a region.
 */
static inline int
bfwork_obtain_region(bam_fwork_t *fwork, bam_region_t *region)
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
		bytes = bam_read1(fwork->input_file->bam_fd, read);
	}

	//Iterate reads
	while(bytes > 0)
	{
		//Wander this read
		omp_set_lock(&region->lock);
#ifdef D_TIME_DEBUG
			times = omp_get_wtime();
#endif
		err = fwork->context->wander_f(fwork, region, read);
#ifdef D_TIME_DEBUG
			times = omp_get_wtime() - times;
			if(fwork->context->time_stats)
				time_add_time_slot(D_FWORK_WANDER_FUNC, fwork->context->time_stats, times);
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
			bytes = bam_read1(fwork->input_file->bam_fd, read);
			break;

		case WANDER_REGION_CHANGED:
			//The region have changed
			last_read = read;
			last_read_bytes = bytes;
			return err;

		default:
			//Unknown error
			LOG_ERROR_F("Framework fails with error code: %d\n", err);
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
			return WANDER_READ_EOF;
		}
		else
		{
			return WANDER_READ_TRUNCATED;
		}
	}

	return NO_ERROR;
}
