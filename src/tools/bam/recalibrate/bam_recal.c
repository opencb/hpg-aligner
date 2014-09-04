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
recal_recalibrate_alignment_priv(bam1_t* alig, const recal_info_t *bam_info, recal_recalibration_env_t *recalibration_env)
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

