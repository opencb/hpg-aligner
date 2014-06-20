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

#include "bam_recal.h"

/***********************************************
 * FRAMEWORK RECALIBRATOR
 **********************************************/

int
recalibrate_wanderer(bam_fwork_t *fwork, bam_region_t *region, bam1_t *read)
{
	assert(fwork);
	assert(region);
	assert(read);

	//Filter read
	if(filter_read(read, FILTER_ZERO_QUAL | FILTER_DIFF_MATE_CHROM | FILTER_NO_CIGAR | FILTER_DEF_MASK))
	{
		//Read is not valid for process
		return WANDER_READ_FILTERED;
	}

	//Update region bounds
	if(region->init_pos > read->core.pos)
	{
		region->init_pos = read->core.pos;
		region->chrom = read->core.tid;
	}
	if(region->end_pos < read->core.pos)
	{
		region->end_pos = read->core.pos;
	}

	return NO_ERROR;
}

int
recalibrate_collect_processor(bam_fwork_t *fwork, bam_region_t *region)
{
	int err, i;
	recal_info_t *data;
	recal_data_collect_env_t *collect_env;
	bam1_t *read;
	size_t *cycles;

	//Get data
	bfwork_local_user_data(fwork, (void **)&data);
	if(data == NULL)
	{
		//Local data is not initialized
		data = (recal_info_t *)malloc(sizeof(recal_info_t));

		//Lock cycles
		bfwork_lock_user_data(fwork, (void **)&cycles);
		recal_init_info(*cycles, data);
		bfwork_unlock_user_data(fwork);

		//Set local data
		bfwork_local_user_data_set(fwork, data);
	}

	//Initialize get data environment
	collect_env = (recal_data_collect_env_t *) malloc(sizeof(recal_data_collect_env_t));
	recal_get_data_init_env(data->num_cycles, collect_env);

	//Obtain data from all reads in region
	for(i = 0; i < region->size; i++)
	{
		//Get next read
		read = region->reads[i];
		assert(read);

		//Get data
		omp_set_lock(&fwork->reference_lock);
		recal_get_data_from_bam_alignment(read, fwork->reference, data, collect_env);
		omp_unset_lock(&fwork->reference_lock);
	}

	//Destroy environment
	recal_get_data_destroy_env(collect_env);

	return NO_ERROR;
}

int
recalibrate_recalibrate_processor(bam_fwork_t *fwork, bam_region_t *region)
{
	int err, i;
	recal_info_t *gdata;
	recal_info_t *data;
	recal_recalibration_env_t *recal_env;
	bam1_t *read;

	//Get data
	bfwork_local_user_data(fwork, (void **)&data);
	if(data == NULL)
	{
		//Local data is not initialized
		data = (recal_info_t *)malloc(sizeof(recal_info_t));

		//Lock data
		bfwork_lock_user_data(fwork, (void **)&gdata);

		//Init struct
		recal_init_info(gdata->num_cycles, data);

		//Clone global data
		if(recal_reduce_info(data, gdata))	//Copy data into local
		{
			LOG_FATAL_F("In local recalibration data copy\nLOCAL: Cycles: %d, Min Q: %d, Num Q: %d, Dinuc: %d\nGLOBAL: Cycles: %d, Min Q: %d, Num Q: %d, Dinuc: %d\n",
					data->num_cycles, data->min_qual, data->num_quals, data->num_dinuc, gdata->num_cycles, gdata->min_qual, gdata->num_quals, gdata->num_dinuc);
		}

		//Recalculate deltas (reduce only merge bases and misses)
		recal_calc_deltas(data);

		bfwork_unlock_user_data(fwork);

		//Set local data
		bfwork_local_user_data_set(fwork, data);
	}

	//Initialize get data environment
	recal_env = (recal_recalibration_env_t *) malloc(sizeof(recal_recalibration_env_t));
	recal_recalibration_init_env(data->num_cycles, recal_env);

	//Recalibrate region
	for(i = 0; i < region->size; i++)
	{
		//Get next read
		read = region->reads[i];
		assert(read);

		//Recalibrate read
		recal_recalibrate_alignment(read, data, recal_env);
	}

	//Destroy environment
	recal_recalibration_destroy_env(recal_env);

	return NO_ERROR;
}

void
reduce_data(void *dest, void *data)
{
	recal_info_t *data_ptr, *dest_ptr;
	assert(data);
	assert(dest);

	//Cast pointers
	data_ptr = (recal_info_t *)data;
	dest_ptr = (recal_info_t *)dest;

	//Combine
	recal_reduce_info(dest_ptr, data_ptr);

	//Delta processing
	recal_calc_deltas(dest_ptr);
	//printf("Estimated %.2f \tEmpirical %.2f \t TotalDelta %.2f\n", dest_ptr->total_estimated_Q, dest_ptr->total_delta + dest_ptr->total_estimated_Q, dest_ptr->total_delta);
}

void
destroy_data(void *data)
{
	recal_info_t *aux;

	assert(data);

	aux = (recal_info_t *)data;
	recal_destroy_info(aux);
}

/**
 * Recalibrate BAM file
 */
ERROR_CODE
recal_bam_file(uint8_t flags, const char *bam_path, const char *ref, const char *data_file, const char *info_file, const char *outbam, int cycles, const char *stats_path)
{
	int bytes;

	//Data
	recal_info_t info;
	U_CYCLES aux_cycles;
	int cycles_param;

	//Times
	double times;

	//Wanderer
	bam_fwork_t fwork;
	bfwork_context_t collect_context;
	bfwork_context_t recal_context;

	assert(bam_path);

	//Get phases
	if((flags & RECALIBRATE_RECALIBRATE) && !(flags & RECALIBRATE_COLLECT))
	{
		//Second phase only, needs data file
		assert(data_file);
	}

	//Create new data
	aux_cycles = cycles;
	recal_init_info(aux_cycles, &info);

	//Init wandering
	bfwork_init(&fwork);

	//Configure framework
	bfwork_configure(&fwork, bam_path, outbam, ref, NULL);

	//Collection is needed?
	if(flags & RECALIBRATE_COLLECT)
	{
		assert(ref);

		//Create data collection context
		bfwork_context_init(&collect_context,
				(int (*)(void *, bam_region_t *, bam1_t *))recalibrate_wanderer,
				(int (*)(void *, bam_region_t *))recalibrate_collect_processor,
				(int (*)(void *, void *))reduce_data,
				&info
		);

#ifdef D_TIME_DEBUG
		//Init timing
		bfwork_context_init_timing(&collect_context, "collect", stats_path);
#endif

		//Set user data
		cycles_param = cycles;
		bfwork_context_set_user_data(&collect_context, &cycles_param);

		//Add context for data collection
		bfwork_add_context(&fwork, &collect_context, FWORK_CONTEXT_SEQUENTIAL);

		printf("Cycles: %d\n",cycles);
	}

	//Recalibration is needed?
	if(flags & RECALIBRATE_RECALIBRATE)
	{
		//Create recalibration context
		bfwork_context_init(&recal_context,
						(int (*)(void *, bam_region_t *, bam1_t *))recalibrate_wanderer,
						(int (*)(void *, bam_region_t *))recalibrate_recalibrate_processor,
						NULL,	//No reduction needed
						NULL
		);

#ifdef D_TIME_DEBUG
		//Init timing
		bfwork_context_init_timing(&recal_context, "recalibrate", stats_path);
#endif

		//Set context user data
		bfwork_context_set_user_data(&recal_context, &info);

		//Previous collect?
		if(!(flags & RECALIBRATE_COLLECT))
		{
			//No previous collect, load data file
			if(data_file)
			{
				//Load data from disk
				recal_load_recal_info(data_file, &info);
			}
			else
			{
				//If only recalibrate input data must be present
				LOG_FATAL("Only recalibration phase 2 requested but no input data specified\n");
			}
		}

		//Add context for recalibration
		bfwork_add_context(&fwork, &recal_context, FWORK_CONTEXT_SEQUENTIAL);
	}

	//Run wander
	bfwork_run(&fwork);

	//Collection last operations
	if(flags & RECALIBRATE_COLLECT)
	{
		//Save data file
		if(data_file)
		{
			recal_save_recal_info(&info, data_file);
		}

		//Save info file
		if(info_file)
		{
			recal_fprint_info(&info, info_file);
		}

		bfwork_context_local_user_data_free(&collect_context, destroy_data);

#ifdef D_TIME_DEBUG
		//Print times
		bfwork_context_print_times(&collect_context);

		//Destroy wanderer time
		bfwork_context_destroy_timing(&collect_context);
#endif

		//Destroy context
		bfwork_context_destroy(&collect_context);
	}

	//Recalibration last operations
	if(flags & RECALIBRATE_RECALIBRATE)
	{
		//Free local data
		bfwork_context_local_user_data_free(&recal_context, destroy_data);

#ifdef D_TIME_DEBUG
		//Print times
		bfwork_context_print_times(&recal_context);

		//Destroy wanderer time
		bfwork_context_destroy_timing(&recal_context);
#endif

		//Destroy context
		bfwork_context_destroy(&recal_context);
	}

	//Destroy wanderer
	bfwork_destroy(&fwork);

	//Free data memory
	recal_destroy_info(&info);

	return NO_ERROR;
}

/**
 * Recalibrate BAM file from path and store in file.
 */
ERROR_CODE
recal_recalibrate_bam_file(const char *orig_bam_path, const recal_info_t *bam_info, const char *recal_bam_path)
{
	bam_file_t *orig_bam_f, *recal_bam_f;
	bam_header_t *recal_bam_header;

	//Open bam
	printf("Opening BAM from \"%s\" to being recalibrated ...\n", orig_bam_path);
	orig_bam_f = bam_fopen((char *)orig_bam_path);
	printf("BAM opened!...\n");

	//Allocate
	//recal_bam_header = (bam_header_t *) calloc(1, sizeof(bam_header_t));

	//Create new bam
	printf("Creating new bam file in \"%s\"...\n", recal_bam_path);
	//init_empty_bam_header(orig_bam_f->bam_header_p->n_targets, recal_bam_header);
	recal_bam_f = bam_fopen_mode((char *)recal_bam_path, orig_bam_f->bam_header_p, "w");
	bam_fwrite_header(recal_bam_f->bam_header_p, recal_bam_f);
	recal_bam_f->bam_header_p = NULL;
	printf("New BAM initialized!...\n");

	//Recalibrate bams
	recal_recalibrate_bam(orig_bam_f, bam_info, recal_bam_f);

	//Memory free
	printf("Closing \"%s\" BAM file...\n", recal_bam_path);
	bam_fclose(recal_bam_f);
	printf("Closing \"%s\" BAM file...\n", orig_bam_path);
	bam_fclose(orig_bam_f);

	printf("BAMs closed.\n");
	printf("Recalibration DONE.\n");

	return NO_ERROR;
}

/**
 * Recalibrate BAM file and store in file.
 */
ERROR_CODE
recal_recalibrate_bam(const bam_file_t *orig_bam_f, const recal_info_t *bam_info, bam_file_t *recal_bam_f)
{
	bam_batch_t *batch;
	bam_batch_t *rdy_batch;
	bam_batch_t *read_batch;
	int count = 0, countb = 0;
	ERROR_CODE err;

	//Thread output
	//pthread_t out_thread;
	//pthread_attr_t out_thread_attr;
	//batch_out_t out_args;
	//void *status;

	// Set thread attributes
	//pthread_attr_init(&out_thread_attr);
	//pthread_attr_setdetachstate(&out_thread_attr, PTHREAD_CREATE_JOINABLE);

	//Allocate memory for batchs
	batch = bam_batch_new(MAX_BATCH_SIZE, MULTIPLE_CHROM_BATCH);
	batch->num_alignments = 0;
	rdy_batch = bam_batch_new(1, MULTIPLE_CHROM_BATCH);
	rdy_batch->num_alignments = 0;

	printf("---------------------\n", count);

	//Set thread num
	//omp_set_num_threads(NUM_THREADS);
	omp_set_nested(1);

	//OMP PARALLEL
	//#pragma omp parallel num_threads(omp_get_num_procs())
	#pragma omp parallel
	{

		#pragma omp single
		printf("Using %d threads\n", omp_get_num_threads());

		do
		{
			#pragma omp sections
			{

				//Read batch
				#pragma omp section
				{
					//Free memory and take a new batch
					//bam_batch_free(batch, 1);
					read_batch = bam_batch_new(MAX_BATCH_SIZE, MULTIPLE_CHROM_BATCH);

					//Read batch
					bam_fread_max_size(read_batch, MAX_BATCH_SIZE, 0, (bam_file_t *)orig_bam_f);
				}

				//Recalibrate batch
				#pragma omp section
				{
					//Recalibrate batch
					if(batch->num_alignments != 0)
					{
						err = recal_recalibrate_batch(batch, bam_info);
						if(err)
							printf("ERROR (recal_recalibrate_batch): %d\n", err);
					}
				}


				//Write batch task
				#pragma omp section
				{
					if(rdy_batch)
					{
						if(rdy_batch->num_alignments != 0)
						{

							bam_fwrite_batch(rdy_batch, recal_bam_f);

							//Update read counter
							count += rdy_batch->num_alignments;

							//Update batch counter
							countb++;

							//Show total progress
							printf("Total alignments recalibrated: %d\r", count);
							fflush(stdout);
						}

						//Free batch
						bam_batch_free(rdy_batch, 1);
						rdy_batch = NULL;
					}
				}
			}

			//#pragma omp barrier

			#pragma omp single
			{
				//Setup next iteration
				rdy_batch = batch;
				batch = read_batch;
				read_batch = NULL;
			}

			//#pragma omp barrier

			//printf("READY BATCH %d\n", rdy_batch->num_alignments);
			//printf("BATCH %d\n", batch->num_alignments);

		}while( (rdy_batch && rdy_batch->num_alignments != 0)
				|| (batch && batch->num_alignments != 0)
			#ifdef D_MAX_READS_W
				&& count < D_MAX_READS_W
			#endif
			);
	}

	printf("\nBatchs writed: %d\n", countb);

	bam_batch_free(batch, 1);
	bam_batch_free(read_batch, 1);

	printf("\n---------------------\n", count);

	return NO_ERROR;
}

/**
 * Recalibrate BAM batch of alignments and store in file.
 */
ERROR_CODE
recal_recalibrate_batch(const bam_batch_t* batch, const recal_info_t *bam_info)
{
	int i;
	ERROR_CODE err;
	//int num_thr;

	//Measures
	double init_time, end_time;

	//Get data environment
	recal_recalibration_env_t *recalibration_env;

	//CHECK ARGUMENTS
	{
		//Check nulls
		if(!batch || !bam_info)
		{
			return INVALID_INPUT_PARAMS_NULL;
		}
	}

	//num_thr = omp_get_num_procs();
	#pragma omp parallel /*num_threads(num_thr)*/ private(recalibration_env, err, init_time, end_time)
	{

		//Initialize get data environment
		recalibration_env = (recal_recalibration_env_t *) malloc(sizeof(recal_recalibration_env_t));
		recal_recalibration_init_env(bam_info->num_cycles, recalibration_env);

		//Process all alignments of the batchs
		#pragma omp for schedule(runtime)
		for(i = 0; i < batch->num_alignments; i++)
		{
			/*#pragma omp critical
			{
				printf("INSIDE ALIGMENTS SEC: %d - %d - %d\n", omp_get_num_threads(), omp_get_thread_num(), i);
				fflush(stdout);
			}*/

			//Process every alignment
			/*err = */recal_recalibrate_alignment_priv(batch->alignments_p[i], bam_info, recalibration_env);
			//if(err)
			//	printf("ERROR (recal_recalibrate_alignment): %d\n", err);
		}

		//Destroy environment
		recal_recalibration_destroy_env(recalibration_env);
	}

	return NO_ERROR;
}

/**
 * Recalibrate alignment and store in file.
 */
ERROR_CODE
recal_recalibrate_alignment(const bam1_t* alig, const recal_info_t *bam_info, recal_recalibration_env_t *recalibration_env)
{
	//Lengths
	uint32_t bam_seq_l;
	uint32_t bam_seq_max_l;

	//CHECK ARGUMENTS (Assuming this function is called always from recal_recalibrate_batch)
	{
		//Check nulls
		if(!alig || !bam_info || !recalibration_env)
			return INVALID_INPUT_PARAMS_NULL;
	}

	//SET VARS
	{
		bam_seq_max_l = recalibration_env->bam_seq_max_l;
	}

	//Sequence length
	bam_seq_l = alig->core.l_qseq;
	if(bam_seq_l > bam_seq_max_l || bam_seq_l == 0)
	{
		return INVALID_SEQ_LENGTH;
	}

	//Process
	recal_recalibrate_alignment_priv(alig, bam_info, recalibration_env);

	return NO_ERROR;
}

//Function for library internal use
static INLINE void
recal_recalibrate_alignment_priv(const bam1_t* alig, const recal_info_t *bam_info, recal_recalibration_env_t *recalibration_env)
{
	U_QUALS qual_index;
	unsigned int matrix_index;
	unsigned int i;
	U_DINUC dinuc;

	//Sequence
	char *bam_seq;
	char *bam_quals;
	char *res_quals;
	U_CYCLES bam_seq_l;
	U_CYCLES bam_seq_max_l;

	//Recalibration
	double delta_r, delta_rc, delta_rd;

	alignment_t* aux_alig;
	bam1_t *aux_alig1;

	//CHECK ARGUMENTS (Assuming this function is called always from recal_recalibrate_batch)
	{
		//Check nulls
		ASSERT(alig);
		ASSERT(bam_info);
		ASSERT(recalibration_env);
	}

	//FILTERS
	{

	}

	//SET VARS
	{
		bam_quals = recalibration_env->bam_quals;
		bam_seq_l = 0;
		bam_seq_max_l = recalibration_env->bam_seq_max_l;
	}

	//Get sequence length
	bam_seq_l = alig->core.l_qseq;

	//Sequence length check
	ASSERT(alig->core.l_qseq <= bam_seq_max_l);
	ASSERT(bam_seq_l != 0);

	//Get sequence
	bam_seq = new_sequence_from_bam((bam1_t *)alig);

	//Get quals
	new_quality_from_bam_ref((bam1_t *)alig, 0, bam_quals, bam_seq_max_l);

	//Allocate for result
	res_quals = (char *)malloc(bam_seq_l * sizeof(char));
	memcpy(res_quals, bam_quals, bam_seq_l * sizeof(char));

	//Iterates nucleotides in this read
	dinuc = 0;
	for(i = 0; i < bam_seq_l; i++)
	{
		//Compare only if the nucleotide is not "N"
		#ifdef NOT_COUNT_NUCLEOTIDE_N
		if(bam_seq[i] != 'N')
		#endif
		{
			//Recalibrate quality
			qual_index = bam_quals[i] - bam_info->min_qual;
			delta_r = bam_info->qual_delta[qual_index];

			matrix_index = qual_index * bam_info->num_cycles + i;
			delta_rc = bam_info->qual_cycle_delta[matrix_index];

			//dont take prev dinuc in first cycle (delta = 0)
			if(i > 0)
			{
				recal_get_dinuc(bam_seq[i-1], bam_seq[i], &dinuc);
				matrix_index = qual_index * bam_info->num_dinuc + i;
				delta_rd = bam_info->qual_dinuc_delta[matrix_index];
			}
			else
			{
				delta_rd = 0.0;
			}

			//Recalibration formula
			double global_delta = bam_info->total_delta;
			double calidad = (double)bam_quals[i];
			if(calidad > MIN_QUALITY_TO_STAT)
			{
				double res = bam_info->total_estimated_Q + bam_info->total_delta + delta_r + delta_rc + delta_rd;
				res_quals[i] = (char)res;
			}
			/*else
			{
				res_quals[i] = (char)calidad;
			}*/
		}
	}

	//Convert bam1_t to alignment_t
	//aux_alig = alignment_new_by_bam(alig, 0);

	//Set qualities in alignment
	//memcpy(aux_alig->quality, res_quals, bam_seq_l);

	//Fix reads in alignment (sequence conversion to string is bug)
	//memcpy(aux_alig->sequence, bam_seq, bam_seq_l);

	//Convert alig to bam1 format
	//aux_alig1 = convert_to_bam(aux_alig, 0);

	//Overwrite batch readings
	//free(alig->data);
	//memcpy(alig, aux_alig1, sizeof(bam1_t));

	memcpy(bam1_qual(alig), res_quals, bam_seq_l);

	//Memory free
#ifdef __SSE2__
	_mm_free(bam_seq);
#else
	free(bam_seq);
#endif
	free(res_quals);
}

