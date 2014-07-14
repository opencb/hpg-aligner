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

#include "data_collector.h"

#ifdef CHECK_DUPLICATES
	static char *ult_seq = NULL;
	static int l_ult_seq = 0;
	static int pos_ult_seq;
#endif

/**
 * Get recalibration data from BAM path.
 */
ERROR_CODE
recal_get_data_from_file(const char *bam_path, const char *ref_name, const char *ref_path, recal_info_t *out_info)
{
	genome_t* ref;
	bam_file_t *bam_f;

	//Open bam
	printf("Opening BAM from \"%s\" ...\n", bam_path);
	bam_f = bam_fopen(bam_path);
	printf("BAM opened!...\n");

	//Open reference genome
	printf("Opening reference genome from \"%s%s\" ...\n", ref_path, ref_name);
	ref = genome_new(ref_name, ref_path, BWT_MODE);
	printf("Reference opened!...\n");

	//Fill data
	#ifdef D_TIME_DEBUG
		time_init_slot(D_SLOT_PH1_COLLECT_BAM, TIME_GLOBAL_STATS);
	#endif
	recal_get_data_from_bam(bam_f, ref, out_info);
	#ifdef D_TIME_DEBUG
		time_set_slot(D_SLOT_PH1_COLLECT_BAM, TIME_GLOBAL_STATS);
	#endif

	//Memory free
	printf("\nClosing BAM file...\n");
	bam_fclose(bam_f);
	genome_free(ref);
	printf("BAM closed.\n");

	return NO_ERROR;
}

/**
 * Get recalibration data from BAM file.
 */
ERROR_CODE
recal_get_data_from_bam(const bam_file_t *bam, const genome_t* ref, recal_info_t* output_data)
{
	bam_batch_t *read_batch;
	bam_batch_t *collect_batch;
	ERROR_CODE err;

	//Duplicate check
	char *last_seq;
	U_CYCLES l_last_seq;
	U_CYCLES pos_last_seq;
	bam1_t *last_alig;

	//Time measure
	double init_read = 0.0, init_collect = 0.0;
	double end_read = 0.0, end_collect = 0.0;
	double init_it = 0.0, end_it = 0.0;

	//Number alignment readed
	int count = 0;

	unmapped = 0;
	duplicated = 0;
	#ifdef NOT_PRIMARY_ALIGNMENT
		notprimary = 0;
	#endif
	#ifdef NOT_MAPPING_QUAL_ZERO
		mapzero = 0;
	#endif

	//DEBUG
	#ifdef D_SAMPLES_OUTPUT
	FILE *fp;
	out_path = "cadenas.out";
	fp = fopen(out_path, "w");
	fclose(fp);
	#endif

	printf("\n----------------\nProcessing \"%s\" file...\n----------------", bam->filename);

	//Allocate first iteration empty batch
	collect_batch = bam_batch_new(1, MULTIPLE_CHROM_BATCH);
	collect_batch->num_alignments = 0;

	//Duplicate checking
	{
		#ifdef CHECK_DUPLICATES
			ult_seq = (char *)malloc(sizeof(char) * output_data->num_cycles);
		#endif
		last_seq = (char *)malloc(sizeof(char));	// Avoid comprobation in bucle and always use free
	}

	omp_set_nested(1);

	//OpenMP parallel, principal proccess
	#pragma omp parallel
	{

		#pragma omp single
		printf("\nUsing %d threads\n", omp_get_num_threads());

		do
		{
			#ifdef D_TIME_DEBUG
				#pragma omp single
				init_it = omp_get_wtime();
			#endif

			#pragma omp sections
			{
				//Read next batch
				#pragma omp section
				{
					read_batch = bam_batch_new(MAX_BATCH_SIZE, SINGLE_CHROM_BATCH);

					//Read batch
					#ifdef D_TIME_DEBUG
						init_read = omp_get_wtime();
					#endif
					//bam_fread_max_size(batch, MAX_BATCH_SIZE, 1, bam);
					bam_fread_max_size_no_duplicates(read_batch, MAX_BATCH_SIZE, 0, bam, last_seq, &l_last_seq, &pos_last_seq);
					#ifdef D_TIME_DEBUG
						end_read = omp_get_wtime();
					#endif
				}

				//Collect batch
				#pragma omp section
				{
					if(collect_batch->num_alignments != 0)
					{
						//Process batch
						#ifdef D_TIME_DEBUG
							init_collect = omp_get_wtime();
						#endif
						err = recal_get_data_from_bam_batch(collect_batch, ref, output_data);
						#ifdef D_TIME_DEBUG
							end_collect = omp_get_wtime();
						#endif

						if(err)
							printf("ERROR (recal_get_data_from_bam_batch): %d\n", err);

						//Update read counter
						count += collect_batch->num_alignments;

						//Show total progress
						printf("Total alignments readed: %d\r", count);
						fflush(stdout);

						//Get last alignment
#ifdef __SSE2__
						_mm_free(last_seq);
#else
						free(last_seq);
#endif
						last_alig = collect_batch->alignments_p[collect_batch->num_alignments-1];
						last_seq = new_sequence_from_bam(last_alig);
						l_last_seq = last_alig->core.l_qseq;
						pos_last_seq = last_alig->core.pos;

						//Free memory and take a new batch
						bam_batch_free(collect_batch, 1);
					}
				}

			}	/* END SECTIONS */

			//#pragma omp barrier

			#pragma omp single
			{
				#ifdef D_TIME_DEBUG
					#ifdef D_TIME_OPENMP_VERBOSE
						fflush(stdout);
						printf("Times:\n");
						printf("Read %.2f ms\n", (end_read - init_read) * 1000.0);
						printf("Collect %.2f ms\n", (end_collect - init_collect) * 1000.0);
						printf("NEW ITERATION\n");
						fflush(stdout);
					#endif

					end_it = omp_get_wtime();

					//Add times
					if(end_read != 0.0) time_add_time_slot(D_SLOT_PH1_READ_BATCH, TIME_GLOBAL_STATS, end_read - init_read);
					if(end_collect != 0.0) time_add_time_slot(D_SLOT_PH1_COLLECT_BATCH, TIME_GLOBAL_STATS, end_collect - init_collect);
					time_add_time_slot(D_SLOT_PH1_ITERATION, TIME_GLOBAL_STATS, end_it - init_it);

					//Reset counters to avoid resample
					init_read = 0.0;
					end_read = 0.0;
					init_collect = 0.0;
					end_collect = 0.0;
				#endif

				//Setup next iteration
				collect_batch = read_batch;
				read_batch = NULL;
			}

			//#pragma omp barrier

		} while( collect_batch->num_alignments != 0
			#ifdef D_MAX_READS
				&& count < D_MAX_READS
			#endif
			);

	}	/* END PARALLEL */

	printf("\n----------------\n%d alignments readed.", count);
	printf("\n%d alignments processed.", count - unmapped - mapzero - duplicated - notprimary);
	printf("\n%d alignments duplicated.", duplicated);
	#ifdef NOT_PRIMARY_ALIGNMENT
		printf("\n%d not primary alignments.", notprimary);
	#endif
	#ifdef NOT_MAPPING_QUAL_ZERO
	printf("\n%d alignments with map quality zero.", mapzero);
	#endif
	printf("\n%d alignments unmapped.", unmapped);

	//Last free
	free(last_seq);
	#ifdef CHECK_DUPLICATES
		free(ult_seq);
	#endif

	return NO_ERROR;
}

/**
 * Get recalibration data from BAM batch of alignments.
 */
ERROR_CODE
recal_get_data_from_bam_batch(const bam_batch_t* batch, const genome_t* ref, recal_info_t* output_data)
{
	int i, j;
	ERROR_CODE err;

	//Batch splitting
	bam_batch_t *current_batch;
	bam_batch_t *v_batchs;
	size_t batchs_l;
	size_t num_chroms;

	//Time measures
	double init_time, end_time;

	//Get data environment
	recal_data_collect_env_t *collect_env;
	recal_info_t *data;

	//CHECK ARGUMENTS
	{
		//Check nulls
		if(!batch || !ref || !output_data)
		{
			return INVALID_INPUT_PARAMS_NULL;
		}
	}

	//Parallel zone
	#pragma omp parallel private(collect_env, data, init_time, end_time, err)
	{

		//Initialize get data environment
		collect_env = (recal_data_collect_env_t *) malloc(sizeof(recal_data_collect_env_t));
		recal_get_data_init_env(output_data->num_cycles, collect_env);

		//Initialize output data
		recal_init_info(output_data->num_cycles, &data);

		//Current is general batch
		current_batch = batch;

		#ifdef SPLIT_BATCHS_BY_CHROM
			//printf("Number alignments in original batch: %d\n", batch->num_alignments);

			//Get number of chroms in this batch
			batch_count_chroms(batch, &num_chroms);

			//Split batchs
			v_batchs = (bam_batch_t *) malloc(sizeof(bam_batch_t) * num_chroms);
			batch_split_by_chrom(batch, v_batchs, &batchs_l, num_chroms);

			//printf("Number of splitted batchs: %d\n", batchs_l);

			//Process every batch
			for(j = 0; j < batchs_l; j++)
			{
				//Get next branch
				current_batch = &v_batchs[j];
				//printf("BATCH %d: Chrom = %d, Aligs = %d, Startpos: %d\n", j, current_batch->alignments_p[0]->core.tid, current_batch->num_alignments, current_batch->alignments_p[0]->core.pos);
		#endif

				//Process all alignments of the batch
				#pragma omp for schedule(runtime)
				for(i = 0; i < current_batch->num_alignments; i++)
				{
					#ifdef D_TIME_DEBUG
					init_time = omp_get_wtime();
					#endif
					//Recollection
					recal_get_data_from_bam_alignment(current_batch->alignments_p[i], ref, data, collect_env);
					#ifdef D_TIME_DEBUG
						end_time = omp_get_wtime();
						#pragma omp critical
							time_add_time_slot(D_SLOT_PH1_COLLECT_ALIG, TIME_GLOBAL_STATS, end_time - init_time);
					#endif
				}

		#ifdef SPLIT_BATCHS_BY_CHROM
				free(current_batch->alignments_p);
			}

			free(v_batchs);
		#endif

		//Reduce data
		#pragma omp critical
		{
			#ifdef D_TIME_DEBUG
				init_time = omp_get_wtime();
			#endif

			err = recal_reduce_info(output_data, data);

			if(err)
				printf("ERROR: Failed to reduce collection data!, error code: %d\n", err);

			#ifdef D_TIME_DEBUG
				end_time = omp_get_wtime();
				time_add_time_slot(D_SLOT_PH1_COLLECT_REDUCE_DATA, TIME_GLOBAL_STATS, end_time - init_time);
			#endif
		}

		//Free data memory
		recal_destroy_info(&data);

		//Destroy environment
		recal_get_data_destroy_env(collect_env);

	} /* END PARALLEL */

	return NO_ERROR;
}

/**
 * Get recalibration data from alignment.
 */
ERROR_CODE
recal_get_data_from_bam_alignment(const bam1_t* alig, const genome_t* ref, recal_info_t* output_data, recal_data_collect_env_t *collect_env)
{
	char *ref_seq;
	char aux_comp[16];
	size_t init_pos, end_pos;
	char *comp_res;
	char *dinucs;
	uint32_t flag;

	//Enviroment
	char *bam_seq;
	char *bam_quals;
	U_CYCLES bam_seq_l;
	char *aux_res_seq;
	char *aux_res_qual;
	U_CYCLES aux_res_seq_l;
	U_CYCLES bam_seq_max_l;

	//SSE
	#ifdef __SSE2__
	__m128i v_ref, v_seq, v_comp;
	#endif

	unsigned int i, j;

	//CHECK ARGUMENTS (Assuming this function is called always from recal_get_data_from_bam_batch)
	{
		//Check nulls
		if(!alig)
			return INVALID_INPUT_PARAMS_NULL;
	}

	//FILTERS
	{
		//Map quality is zero
		#ifdef NOT_MAPPING_QUAL_ZERO
		if(alig->core.qual == 0)
		{
			mapzero++;
			return NO_ERROR;
		}
		#endif

		//Not primary
		#ifdef NOT_PRIMARY_ALIGNMENT
		if(alig->core.flag & 256)	//Not primary alignment flag
		{
			notprimary++;
			return NO_ERROR;
		}
		#endif

		//Unmapped readings
		if(alig->core.tid < 0)
		{
			unmapped++;
			return NO_ERROR;
		}
	}

	//SET VARS
	{
		bam_seq = collect_env->bam_seq;
		bam_quals = collect_env->bam_quals;
		bam_seq_l = 0;
		aux_res_seq = collect_env->aux_res_seq;
		aux_res_qual = collect_env->aux_res_qual;
		aux_res_seq_l = 0;
		bam_seq_max_l = collect_env->bam_seq_max_l;
	}

	//Get sequence
	new_sequence_from_bam_ref(alig, bam_seq, bam_seq_max_l);

	//Get quals
	new_quality_from_bam_ref(alig, 0, bam_quals, bam_seq_max_l);

	//Indel suppression
 	supress_indels_from_32_cigar(bam_seq, bam_quals, alig->core.l_qseq, bam1_cigar(alig), alig->core.n_cigar,
 			aux_res_seq, aux_res_qual, &aux_res_seq_l, bam_seq_max_l);

	//Check if sequence length is valid
	if(aux_res_seq_l == 0)
	{
		return INVALID_SEQ_LENGTH;
	}

	//Save sequence to primary array
	memcpy(bam_seq, aux_res_seq, aux_res_seq_l * sizeof(char));
	memcpy(bam_quals, aux_res_qual, aux_res_seq_l * sizeof(char));

	//Get cycles and positions
	//cycles = alig->core.l_qseq;
	bam_seq_l = aux_res_seq_l;
	init_pos = alig->core.pos + 1;
	end_pos = alig->core.pos + bam_seq_l;

	//Duplicates check
	#ifdef CHECK_DUPLICATES
	{
		if(l_ult_seq)
		{
			if(pos_ult_seq == init_pos && l_ult_seq == bam_seq_l && strcmp(bam_seq, ult_seq) == 0)
			{
				//printf("\nDUPLICATE POS: %d CYCLES: %d\n\tSEQ:  %s\n\tLAST: %s", init_pos, cycles, bam_seq, ult_seq);
				duplicated++;
				return NO_ERROR;
			}
		}
	}
	#endif

	//Allocations
#ifdef __SSE2__
	ref_seq = (char *)_mm_malloc((bam_seq_l + 1) * sizeof(char), MEM_ALIG_SSE_SIZE);
	comp_res = (char *)_mm_malloc(bam_seq_l * sizeof(char), MEM_ALIG_SSE_SIZE);
	dinucs = (char *)_mm_malloc(bam_seq_l * sizeof(char), MEM_ALIG_SSE_SIZE);
#else
	ref_seq = (char *)malloc((bam_seq_l + 1) * sizeof(char));
	comp_res = (char *)malloc(bam_seq_l * sizeof(char));
	dinucs = (char *)malloc(bam_seq_l * sizeof(char));
#endif

	//Obtain reference for this 100 nucleotides
	flag = (uint32_t) alig->core.flag;

	genome_read_sequence_by_chr_index(ref_seq, (flag & BAM_FREVERSE) ? 1 : 0, (unsigned int)alig->core.tid, &init_pos, &end_pos, ref);

	//Iterates nucleotides in this read
	for(i = 0; i < bam_seq_l; i++)
	{
#ifdef __SSE2__
		//#ifdef __SSE2__ //SSE Block
		if( (i + 16) < bam_seq_l)
		{
			//Use SSE
			//_mm_prefetch(&ref_seq[i + 16], _MM_HINT_T0);
			//_mm_prefetch(&bam_seq[i + 16], _MM_HINT_T0);

			//Pack sequences
			v_ref = _mm_load_si128(&ref_seq[i]);
			v_seq = _mm_load_si128(&bam_seq[i]);

			//Compare sequences
			v_comp = _mm_cmpeq_epi8(v_ref, v_seq);

			//Store comparation values
			_mm_store_si128(&comp_res[i], v_comp);

			i += 15;
		}
		else
#endif //SSE Block
		{
			if(ref_seq[i] != bam_seq[i])
			{
				comp_res[i] = 0x00;	//Diff
			}
			else
			{
				comp_res[i] = 0xFF;	//Equals
			}
		}
	}

	//Dinucs
	for(i = 0; i < bam_seq_l; i++)
	{
		if(i > 0)
		{
			recal_get_dinuc(bam_seq[i-1], bam_seq[i], &dinucs[i]);
		}
		else
		{
			dinucs[i] = d_X;
		}
	}

	//Add data
	recal_add_base_v(output_data, bam_seq, bam_quals, 0, bam_seq_l, dinucs, comp_res);

	//Set last sequence for duplicates
	#ifdef CHECK_DUPLICATES
	{
		strcpy(ult_seq, bam_seq);
		l_ult_seq = alig->core.l_qseq;
		pos_ult_seq = init_pos;
	}
	#endif

	//Free resources
	{
#ifdef __SSE2__
		_mm_free(ref_seq);
		_mm_free(comp_res);
		_mm_free(dinucs);
#else
		free(ref_seq);
		free(comp_res);
		free(dinucs);
#endif
	}

	return NO_ERROR;
}
