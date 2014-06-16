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
	bam_f = bam_fopen((char *)bam_path);
	printf("BAM opened!...\n");

	//Open reference genome
	printf("Opening reference genome from \"%s%s\" ...\n", ref_path, ref_name);
	ref = genome_new((char *)ref_name, (char *)ref_path);
	printf("Reference opened!...\n");

	//Fill data
	recal_get_data_from_bam(bam_f, ref, out_info);

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

			#pragma omp sections
			{
				//Read next batch
				#pragma omp section
				{
					read_batch = bam_batch_new(MAX_BATCH_SIZE, SINGLE_CHROM_BATCH);

					//Read batch
					//bam_fread_max_size(batch, MAX_BATCH_SIZE, 1, bam);
					bam_fread_max_size_no_duplicates(read_batch, MAX_BATCH_SIZE, 0, (bam_file_t *)bam, last_seq, &l_last_seq, &pos_last_seq);
				}

				//Collect batch
				#pragma omp section
				{
					if(collect_batch->num_alignments != 0)
					{
						//Process batch
						err = recal_get_data_from_bam_batch(collect_batch, ref, output_data);

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
		data = (recal_info_t *)malloc(sizeof(recal_info_t));
		recal_init_info(output_data->num_cycles, data);

		//Current is general batch
		current_batch = (bam_batch_t *)batch;

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
					//Recollection
					recal_get_data_from_bam_alignment(current_batch->alignments_p[i], ref, data, collect_env);
				}

		#ifdef SPLIT_BATCHS_BY_CHROM
				free(current_batch->alignments_p);
			}

			free(v_batchs);
		#endif

		//Reduce data
		#pragma omp critical
		{
			err = recal_reduce_info(output_data, data);

			if(err)
				printf("ERROR: Failed to reduce collection data!, error code: %d\n", err);
		}

		//Free data memory
		recal_destroy_info(data);
		free(data);

		//Destroy environment
		recal_get_data_destroy_env(collect_env);

	} /* END PARALLEL */

	return NO_ERROR;
}

/**
 * Get recalibration data from alignment.
 */
ERROR_CODE
recal_get_data_from_bam_alignment(const bam1_t* read, const genome_t* ref, recal_info_t* output_data, recal_data_collect_env_t *collect_env)
{
	char *ref_seq;
	char aux_comp[16];
	size_t init_pos, end_pos;
	size_t init_pos_ref, end_pos_ref;
	char *comp_res;
	char *comp_mask;
	char *dinucs;
	uint32_t flag;

	//Enviroment
	char *bam_seq;
	char *bam_quals;
	U_CYCLES bam_seq_l;
	char *aux_res_seq;
	char *aux_res_qual;
	//U_CYCLES aux_res_seq_l;
	//U_CYCLES bam_seq_max_l;
	uint32_t *read_left_cigar;
	uint32_t *read_cigar;
	size_t read_disp_clip;
	char *read_seq_ref;
	size_t read_seq_ref_l;
	uint32_t misses;
	size_t read_l;

	//SSE
	#ifdef __SSE2__
	__m128i v_ref, v_seq, v_comp;
	#endif

	unsigned int i, j;

	//CHECK ARGUMENTS (Assuming this function is called always from recal_get_data_from_bam_batch)
	{
		//Check nulls
		if(!read)
			return INVALID_INPUT_PARAMS_NULL;
	}

	//FILTERS
	{
		//Map quality is zero
		#ifdef NOT_MAPPING_QUAL_ZERO
		if(read->core.qual == 0)
		{
			mapzero++;
			return NO_ERROR;
		}
		#endif

		//Not primary
		#ifdef NOT_PRIMARY_ALIGNMENT
		if(read->core.flag & 256)	//Not primary alignment flag
		{
			notprimary++;
			return NO_ERROR;
		}
		#endif

		//Unmapped readings
		if(read->core.tid < 0)
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
		//aux_res_seq_l = 0;
		//bam_seq_max_l = collect_env->bam_seq_max_l;
	}

	//Get sequence
	new_sequence_from_bam_ref((bam1_t *)read, bam_seq, read->core.l_qseq + 1);

	//Get quals
	new_quality_from_bam_ref((bam1_t *)read, 0, bam_quals, read->core.l_qseq + 1);

	//Indel suppression
 	//supress_indels_from_32_cigar(bam_seq, bam_quals, read->core.l_qseq, bam1_cigar(read), read->core.n_cigar,
 	//		aux_res_seq, aux_res_qual, &aux_res_seq_l, bam_seq_max_l);

	//Check if sequence length is valid
	/*if(aux_res_seq_l == 0)
	{
		return INVALID_SEQ_LENGTH;
	}*/

	//Save sequence to primary array
	//memcpy(bam_seq, aux_res_seq, aux_res_seq_l * sizeof(char));
	//memcpy(bam_quals, aux_res_qual, aux_res_seq_l * sizeof(char));

	//Get cycles and positions
	//cycles = alig->core.l_qseq;
	//bam_seq_l = aux_res_seq_l;
	bam_seq_l = read->core.l_qseq;
	init_pos = read->core.pos;
	end_pos = read->core.pos + (bam_seq_l  * 2);
	init_pos_ref = init_pos + RECAL_REFERENCE_CORRECTION_OFFSET;
	end_pos_ref = end_pos + RECAL_REFERENCE_CORRECTION_OFFSET;

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
/*#ifdef __SSE2__
	ref_seq = (char *)_mm_malloc((read->core.l_qseq + 1) * sizeof(char), MEM_ALIG_SSE_SIZE);
	comp_res = (char *)_mm_malloc((read->core.l_qseq + 1) * sizeof(char), MEM_ALIG_SSE_SIZE);
	dinucs = (char *)_mm_malloc((read->core.l_qseq + 1) * sizeof(char), MEM_ALIG_SSE_SIZE);
#else*/
	ref_seq = (char *)malloc(((end_pos - init_pos) + 2) * sizeof(char));
	comp_res = (char *)malloc((read->core.l_qseq + 1) * sizeof(char));
	comp_mask = (char *)malloc((read->core.l_qseq + 1) * sizeof(char));
	dinucs = (char *)malloc((read->core.l_qseq + 1) * sizeof(char));
	memset(comp_res, 0, (read->core.l_qseq + 1) * sizeof(char));
	memset(comp_mask, 0, (read->core.l_qseq + 1) * sizeof(char));
//#endif

	//Obtain reference for this 100 nucleotides
	flag = (uint32_t) read->core.flag;

	genome_read_sequence_by_chr_index(ref_seq, (flag & BAM_FREVERSE) ? 1 : 0, (unsigned int)read->core.tid, &init_pos_ref, &end_pos_ref, (genome_t *)ref);

	//Get initial clip displacement
	//cigar32_count_clip_displacement(read_cigar, read->core.n_cigar, &read_disp_clip);
	cigar32_count_nucleotides_not_clip(bam1_cigar(read), read->core.n_cigar, &read_l);

	//Create sequence to compare with reference
	read_seq_ref = (char *) malloc((read->core.l_qseq + 1) * sizeof(char));
	cigar32_create_ref(bam1_cigar(read), read->core.n_cigar,
			ref_seq, (end_pos_ref - init_pos_ref),
			bam_seq, read->core.l_qseq,
			read_seq_ref, &read_seq_ref_l, comp_mask);

	//Correct comparation?
	if(read_seq_ref_l != read->core.l_qseq)
	{
		LOG_WARN_F("Read-Ref: %d, Read: %d, Not enough reference sequence length\n", read_seq_ref_l, read->core.l_qseq);
	}

	//Get raw score with reference
	nucleotide_compare(read_seq_ref, bam_seq, read_seq_ref_l, comp_res, &misses);

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
	recal_add_base_v(output_data, bam_seq, bam_quals, 0, bam_seq_l, dinucs, comp_res, comp_mask);

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
/*#ifdef __SSE2__
		_mm_free(ref_seq);
		_mm_free(comp_res);
		_mm_free(dinucs);
#else*/
		free(ref_seq);
		free(comp_res);
		free(comp_mask);
		free(dinucs);
		free(read_seq_ref);
//#endif
	}

	return NO_ERROR;
}
