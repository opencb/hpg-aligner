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

#include "alig.h"

uint64_t cigar_changed = 0;

char log_msg[1024];
char aux_msg[512];

/**
 * CONTEXT
 */

/**
 * Initialize empty realignment data structure.
 */
ERROR_CODE
alig_init(alig_context_t *context, linked_list_t *in_list, genome_t *genome, uint8_t flags)
{
	ERROR_CODE err;

	assert(context);
	assert(in_list);
	assert(genome);

	//Set all to 0 / NULL
	memset(context, 0, sizeof(alig_context_t));

	//Set fields
	context->in_list = in_list;
	context->genome = genome;

	//Create filtered list
	context->filtered_list = array_list_new(ALIG_LIST_NEXT_SIZE, 1.5, COLLECTION_MODE_SYNCHRONIZED);
	if(context->filtered_list == NULL)
	{
		return ALIG_INIT_FAIL;
	}

	//Create realign list
	context->realign_list = array_list_new(ALIG_LIST_NEXT_SIZE, 1.5, COLLECTION_MODE_SYNCHRONIZED);
	if(context->realign_list == NULL)
	{
		return ALIG_INIT_FAIL;
	}

	//Create haplotype list
	context->haplo_list = array_list_new(100, 1.2f, COLLECTION_MODE_SYNCHRONIZED);
	if(context->haplo_list == NULL)
	{
		return ALIG_INIT_FAIL;
	}

	//Init reference
	memset(&context->reference, 0, sizeof(alig_reference_t));

	//Init scores
	memset(&context->scores, 0, sizeof(alig_scores_t));

	//Set flags
	context->flags = flags;

	return NO_ERROR;
}

/**
 * Free resources from realigner.
 */
ERROR_CODE
alig_destroy(alig_context_t *context)
{
	int i;
	bam1_t *read;
	aux_indel_t *haplo;
	size_t list_l;

	assert(context);

	//Free filtered list
	if(context->filtered_list)
	{
		//Free list
		array_list_free(context->filtered_list, NULL);
	}

	//Free filtered list
	if(context->realign_list)
	{
		//Free list
		array_list_free(context->realign_list, NULL);
	}

	//Free haplotype list
	if(context->haplo_list)
	{
		//Free haplotypes
		for(i = 0; i < array_list_size(context->haplo_list); i++)
		{
			//Get haplo
			haplo = array_list_get(i, context->haplo_list);

			//Free haplo
			if(haplo)
				free(haplo);
		}

		//Free list
		array_list_free(context->haplo_list, free);
	}

	//Free reference
	if(context->reference.reference)
	{
		free(context->reference.reference);
	}

	//Free scores
	if(context->scores.m_positions)
		free(context->scores.m_positions);
	if(context->scores.m_positions)
		free(context->scores.m_positions);

	//Set all to NULL
	memset(context, 0, sizeof(alig_context_t));

	return NO_ERROR;
}

/**
 * Check if a realigner context is valid.
 */
ERROR_CODE
alig_validate(alig_context_t *context)
{
	if(context == NULL)
		return ALIG_INVALID_CONTEXT;

	if(context->in_list == NULL)
		return INVALID_INPUT_BAM;
	if(context->genome == NULL)
		return INVALID_GENOME;

	//Have filtered list?
	if(context->filtered_list == NULL)
		return ALIG_INVALID_PROCESS_LIST;

	//Have realign list?
	if(context->realign_list == NULL)
		return ALIG_INVALID_PROCESS_LIST;

	//Have process list?
	if(context->haplo_list == NULL)
		return ALIG_INVALID_HAPLO_LIST;

	//Valid reference?
	if(context->reference.length != 0 && !context->reference.reference)
		return ALIG_INVALID_REFERENCE;

	//Valid scores?
	if(context->scores.m_total != 0 && (!context->scores.m_positions || !context->scores.m_scores))
		return ALIG_INVALID_SCORES;

	return NO_ERROR;
}

/**
 * REGION OPERATIONS
 */

/**
 * Get next region of reads to process.
 */
ERROR_CODE
alig_region_next(alig_context_t *context)
{
	ERROR_CODE err;
	int i;

	//List
	linked_list_iterator_t *it;
	array_list_t *aux_list;
	size_t aux_list_l;

	//Read
	bam1_t* read = NULL;
	int32_t read_pos;
	int32_t last_read_chrom = -1;
	size_t indels;
	size_t interval_read_begin;
	size_t interval_read_end;

	assert(context);

	//Validate context
	err = alig_validate(context);
	if(err)
	{
		return err;
	}

	//Create list
	aux_list = array_list_new(ALIG_LIST_NEXT_SIZE, 1.5, COLLECTION_MODE_SYNCHRONIZED);

	//Set last readed to 0
	context->last_readed_count = 0;

	//Get first read
	it = linked_list_iterator_new(context->in_list);
	assert(it);
	read = (bam1_t *)linked_list_iterator_curr(it);
	if(read == NULL)
	{
		//No more input reads
		linked_list_iterator_free(it);
		return NO_ERROR;
	}

	last_read_chrom = read->core.tid;

	//Reset interval
	memset(&context->region, 0, sizeof(alig_region_t));

	//Read alignments until interval reached
	while(read != NULL)
	{
		//FILTERS: MAP QUALITY = 0, BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP, DIFFERENT MATE CHROM, 1 INDEL, CIGAR PRESENT
		if( read->core.qual != 0
			&& !(read->core.flag & BAM_DEF_MASK)
			&& read->core.mtid == read->core.tid
			&& read->core.n_cigar != 0
			)
		{
			//Count indels
			cigar32_count_indels(bam1_cigar(read), read->core.n_cigar, &indels);

			//Check if chrom is different
			if(last_read_chrom != read->core.tid)
			{
				//Interval end!

				//Skip while
				break;
			}

			//Get read position
			read_pos = read->core.pos;

			//Cases
			if(context->region.valid == 0)
			{
				//NO INTERVAL CASE

				//Get interval for this alignment
				err = region_get_from_bam1(read, &interval_read_begin, &interval_read_end);
				if(err)
				{
					sprintf(log_msg, "Trying to get region from invalid read: %s\n", bam1_qname(read));
					LOG_ERROR(log_msg);
					return INVALID_INPUT_BAM;
				}

				//This alignment have an interval?
				if(interval_read_begin != SIZE_MAX)
				{
					//Interval found
					//Set interval values
					context->region.start_pos = interval_read_begin;
					context->region.end_pos = interval_read_end;
					context->region.chrom = read->core.tid;
					context->region.valid = 1;

					assert(context->region.start_pos <= context->region.end_pos);
				}

			}
			else if (context->region.end_pos > read_pos)
			{
				//ALIGNMENT INSIDE INTERVAL CASE

				//Get interval for this alignment
				region_get_from_bam1(read, &interval_read_begin, &interval_read_end);

				//This alignment have an interval?
				if(interval_read_begin != SIZE_MAX)
				{
					//Interval found
					//Set interval begin
					if(interval_read_begin < context->region.start_pos)
						context->region.start_pos = interval_read_begin;

					//Set interval end
					if(interval_read_end > context->region.end_pos)
						context->region.end_pos = interval_read_end;

					//Check
					assert(context->region.start_pos <= context->region.end_pos);
				}
			}
			else //if (context->region.end_pos < read_pos)
			{
				//ALIGNMENT PAST INTERVAL
				//Interval end!

				//Filter out of region reads
				break;

			}	//Cases if

			//Save last read chrom
			last_read_chrom = read->core.tid;

			//Add filtered read to list
			array_list_insert(read, aux_list);

		}	// Filters if

		//Increment progress
		context->read_count++;
		context->last_readed_count++;

		//Read next alignment
		read = (bam1_t *)linked_list_iterator_next(it);
	}

	//There is an valid interval?
	if(context->region.valid)
	{
		//Logging
		LOG_INFO("************************************\n");
		sprintf(log_msg, "INTERVAL %d - %d\n", context->region.start_pos + 1, context->region.end_pos + 1);
		LOG_INFO(log_msg);

		//Save interval reads in filtered list
		aux_list_l = array_list_size(aux_list);
		for(i = 0; i < aux_list_l; i++)
		{
			//Get read
			read = array_list_get(i, aux_list);

			if(region_bam_overlap(read, &context->region))
			{
				//Logging
				sprintf(log_msg, "%s \t%d:%d\n", bam1_qname(read), read->core.tid + 1, read->core.pos + 1);
				LOG_INFO(log_msg);

				//Add to filtered read list
				array_list_insert(read, context->filtered_list);
			}
		}

	}

	//Free
	linked_list_iterator_free(it);
	array_list_free(aux_list, NULL);

	return NO_ERROR;
}

/**
 * Load reference sequence for present region in context.
 */
ERROR_CODE
alig_region_load_reference(alig_context_t *context)
{
	ERROR_CODE err;

	//Read
	bam1_t *read;

	//Reference
	char *ref_seq = NULL;
	uint32_t flag;
	int64_t ref_pos_begin = SIZE_MAX;
	int64_t ref_pos_end = SIZE_MAX;
	size_t aux_pos_begin = SIZE_MAX;
	size_t aux_pos_end = SIZE_MAX;
	size_t ref_length = SIZE_MAX;

	assert(context);

	//Validate context
	err = alig_validate(context);
	if(err)
	{
		return err;
	}

	//Take reference only if a region is present
	if(context->region.valid)
	{
		if(array_list_size(context->filtered_list) == 0)
		{
			//Invalid region
			return ALIG_INVALID_REGION;
		}

		//Get first read position
		read = (bam1_t *)array_list_get(0, context->filtered_list);
		assert(read);
		ref_pos_begin = (read->core.pos - ALIG_REFERENCE_ADDITIONAL_OFFSET);
		if(ref_pos_begin < 0)
		{
			ref_pos_begin = 0;
		}

		//Get last read position + length
		read = (bam1_t *) array_list_get(array_list_size(context->filtered_list) - 1, context->filtered_list);
		assert(read);
		ref_pos_end = (read->core.pos + read->core.l_qseq) + ALIG_REFERENCE_ADDITIONAL_OFFSET;

		//Auxiliar
		aux_pos_begin = ref_pos_begin + 1;
		aux_pos_end = ref_pos_end + 1;
		//MANEJAR ALINEAMIENTOS CON MEMORIA!!!!!
		//Get reference
#ifdef __SSE2__
		ref_seq = (char *) _mm_malloc(((ref_pos_end - ref_pos_begin) + 2) * sizeof(char), MEM_ALIG_SSE_SIZE);
#else
		ref_seq = (char *) malloc(((ref_pos_end - ref_pos_begin) + 2) * sizeof(char));
#endif
		assert(ref_seq);
		flag = (uint32_t) read->core.flag;
		genome_read_sequence_by_chr_index(ref_seq, (flag & BAM_FREVERSE) ? 1 : 0, (unsigned int)read->core.tid, &aux_pos_begin, &aux_pos_end, context->genome);

		ref_length = (ref_pos_end - ref_pos_begin);

		//Set context reference
		if(context->reference.reference)
		{
			//Free previous reference
#ifdef __SSE2__
			free(context->reference.reference);
#else
			_mm_free(context->reference.reference);
#endif
		}
		context->reference.position = ref_pos_begin;
		context->reference.length = ref_length;
		context->reference.reference = ref_seq;
	}

	return NO_ERROR;
}

/**
 * Get haplotypes from present region.
 */
ERROR_CODE
alig_region_haplotype_process(alig_context_t *context)
{
	ERROR_CODE err;
	int i;

	//Read
	bam1_t* read;
	char *read_seq;
	char *comp_seq;
	size_t read_seq_l;
	uint32_t *read_cigar;
	size_t read_indels;
	size_t read_bases;
	int64_t read_disp_ref;

	//CIGARS
	uint32_t aux_cigar[MAX_CIGAR_LENGTH];
	size_t aux_cigar_l;

	//Reference
	char *ref_seq;
	uint32_t ref_miss;

	//Haplotype
	aux_indel_t *haplo;
	aux_indel_t *aux_haplo;
	size_t haplo_l;

	//List
	array_list_t *list;
	size_t list_l;

	assert(context);

	//Validate context
	err = alig_validate(context);
	if(err)
	{
		return err;
	}

	//Realign only if a region is present
	if(context->region.valid)
	{
		//Valid reference?
		if(!context->reference.reference)
		{
			LOG_ERROR("Warning: Trying to extract haplotypes with uninitialized reference\n");
			return ALIG_INVALID_REFERENCE;
		}

		//No need to filter

		//Get list
		list = context->filtered_list;

		//Get all haplotypes
		list_l = array_list_size(list);
		for(i = 0; i < list_l; i++)
		{
			//Get read
			read = array_list_get(i, list);
			assert(read);

			//Get read cigar
			read_cigar = bam1_cigar(read);
			if(read_cigar)	//Valid cigar?
			{
				//Convert read sequence and qualities to string
				read_seq_l = read->core.l_qseq;

#ifdef __SSE2__
				read_seq = (char *) _mm_malloc((read_seq_l + 1) * sizeof(char), MEM_ALIG_SSE_SIZE);
				comp_seq = (char *) _mm_malloc((read_seq_l + 1) * sizeof(char), MEM_ALIG_SSE_SIZE);
#else
				read_seq = (char *) malloc((read_seq_l + 1) * sizeof(char));
				comp_seq = (char *) malloc((read_seq_l + 1) * sizeof(char));
#endif

				assert(read);
				new_sequence_from_bam_ref(read, read_seq, read_seq_l + 1);

				//Get reference
				ref_seq = context->reference.reference;

				//Get relative displacement from reference
				read_disp_ref = read->core.pos - context->reference.position;
				if(read_disp_ref < 0)
					read_disp_ref = 0;

				//Get raw missmatch
				//nucleotide_compare(ref_seq + read_disp_ref, read_seq, read_seq_l, comp_seq, &ref_miss);
				ref_miss = 1;

				//If match reference perfectly, realign to reference and extract from process list
				if(ref_miss == 0)
				{
					//Count unclipped bases
					cigar32_count_nucleotides_not_clip(read_cigar, read->core.n_cigar, &read_bases);

					//Create new cigar with '$bases'M
					aux_cigar[0] = read_bases << BAM_CIGAR_SHIFT;	//ex: 108M
					aux_cigar_l = 1;

					//Reclip cigar
					cigar32_reclip(read_cigar, read->core.n_cigar, aux_cigar, aux_cigar_l, aux_cigar, &aux_cigar_l);

					//Replace read cigar
					cigar32_replace(read, aux_cigar, aux_cigar_l);

					//Logging
					//sprintf(log_msg, "%s \t%d:%d match reference perfectly\n", bam1_qname(read), read->core.tid + 1, read->core.pos + 1);
					//LOG_INFO(log_msg);
				}
				else
				{
					//Obtain haplotypes

					//Add this read to realign
					array_list_insert(read, context->realign_list);

					//This read have indels?
					cigar32_count_indels(read_cigar, read->core.n_cigar, &read_indels);
					if(read_indels == 1)	//Only one indel leftalign
					{
						//Left align cigar first
						if(context->flags & ALIG_LEFT_ALIGN)
						{
							//Leftmost cigar
							cigar32_leftmost(ref_seq + read_disp_ref, (context->reference.length) - read_disp_ref,
									read_seq, read->core.l_qseq,
									read_cigar, read->core.n_cigar,
									aux_cigar, &aux_cigar_l);
						} //Leftalign
						else
						{
							//Use original cigar
							memcpy(aux_cigar, read_cigar, read->core.n_cigar * sizeof(uint32_t));
							aux_cigar_l = read->core.n_cigar;
						}

						//Allocate haplotype
						haplo = (aux_indel_t *)malloc(sizeof(aux_indel_t) * read_indels);	//For now is only 1

						//Fill haplotype
						cigar32_get_indels(read->core.pos, aux_cigar, aux_cigar_l, haplo);

						//Check if haplotype is present in list
						aux_haplo = NULL;
						int h;
						haplo_l = array_list_size(context->haplo_list);
						for(h = 0; h < haplo_l; h++)
						{
							aux_haplo = array_list_get(h, context->haplo_list);
							assert(aux_haplo);
							if(aux_haplo->indel == haplo->indel && aux_haplo->ref_pos == haplo->ref_pos)
								break;
						}

						//Duplicate?
						if(h == haplo_l)
						{
							//Add to haplotype list (not duplicated)
							array_list_insert(haplo, context->haplo_list);
						}
						else
						{
							free(haplo);
						}
					} //Indel == 1 if

				}//Match reference if

				//Free
#ifdef __SSE2__
				_mm_free(read_seq);
				_mm_free(comp_seq);
#else
				free(read_seq);
				free(comp_seq);
#endif

			} //Valid cigar end if

		} //Haplotypes while

	} //Region if

	return NO_ERROR;
}

/**
 * Realign readings around indels using current haplotypes.
 */
ERROR_CODE
alig_region_indel_realignment(alig_context_t *context)
{
	ERROR_CODE err;
	int i;

	//Haplotype
	int alt_haplo_index;
	uint32_t alt_haplo_score;
	uint32_t ref_haplo_score;

	//Improvement
	double improvement;

	assert(context);

	//Validate context
	err = alig_validate(context);
	if(err)
	{
		return err;
	}

	//Realign only if a region is present and there is haplotypes
	if(!context->region.valid)
	{
		//LOG_ERROR("Warning: Trying to realign an invalid region\n");
		//return ALIG_INVALID_REGION;
		return NO_ERROR;
	}

	if(array_list_size(context->haplo_list) < 2)
	{
		LOG_INFO("Not enough haplotypes to realign\n");
		//return ALIG_INVALID_REGION;
		return NO_ERROR;
	}

	if(!context->reference.reference)
	{
		LOG_ERROR("Warning: Trying to indel realign with uninitialized reference\n");
		return ALIG_INVALID_REFERENCE;
	}

	//Get scores tables
	err = alig_get_scores(context);
	if(err)
	{
		sprintf(log_msg, "Error getting scores from region, error: %d\n", err);
		LOG_ERROR(log_msg);
		return ALIG_INVALID_SCORES;
	}

	alig_get_alternative_haplotype(context, &alt_haplo_index, &alt_haplo_score, &ref_haplo_score);

	//Improvement
	if(alt_haplo_score != 0)
	{
		improvement = (double)((double)ref_haplo_score / (double)alt_haplo_score) / 10.0;
	}
	else
	{
		improvement = ALIG_IMPROVEMENT_THREHOLD + 1;
	}

	//Enought improvement?
	if(improvement > ALIG_IMPROVEMENT_THREHOLD)
	{
		//Realign
		alig_indel_realign_from_haplo(context, alt_haplo_index);
	}

	return NO_ERROR;
}

/**
 * Clear context to be ready for next region.
 */
ERROR_CODE
alig_region_clear(alig_context_t *context)
{
	ERROR_CODE err;
	int i;

	//Lists
	size_t list_l;

	assert(context);

	//Clear region
	memset(&context->region, 0, sizeof(alig_region_t));

	//Clear bam lists
	array_list_clear(context->filtered_list, NULL);
	array_list_clear(context->realign_list, NULL);

	//Clear haplotype list
	array_list_clear(context->haplo_list, free);

	//Clear reference
	if(context->reference.reference)
	{
#ifdef __SSE2__
		_mm_free(context->reference.reference);
#else
		free(context->reference.reference);
#endif
	}
	memset(&context->reference, 0, sizeof(alig_reference_t));

	//Clear scores
	if(context->scores.m_total != 0)
	{
		if(context->scores.m_positions)
			free(context->scores.m_positions);
		if(context->scores.m_scores)
			free(context->scores.m_scores);
	}
	memset(&context->scores, 0, sizeof(alig_scores_t));

	return NO_ERROR;
}

/**
 * Indel realign one file.
 */
ERROR_CODE
alig_bam_file(char *bam_path, char *ref_name, char *ref_path, char *outbam)
{
	ERROR_CODE err;
	//	unsigned char outbam[30] = "output.bam";
	int i;

	//Files
	bam_file_t *bam_f;
	bam_file_t *out_bam_f;
	genome_t* ref;
	int bytes;

#ifdef D_TIME_DEBUG
	double init_time, end_time;
#endif

	//Read
	bam1_t *read;

	//Lists
	linked_list_t *in_list;
	size_t list_l;

	//Alignment context
	alig_context_t context;

	//aux
	int aux_count = 0;

	assert(bam_path);
	assert(ref_name);
	assert(ref_path);

#ifdef D_TIME_DEBUG
    init_time = omp_get_wtime();
#endif
	//Open bam
	{
		printf("Opening BAM from \"%s\" ...\n", bam_path);
		bam_f = bam_fopen(bam_path);
		assert(bam_f);
		printf("BAM opened!...\n");
	}

	//Open reference genome
	{
		printf("Opening reference genome from \"%s%s\" ...\n", ref_path, ref_name);
		ref = genome_new(ref_name, ref_path, BWT_MODE);
		assert(ref);
		printf("Reference opened!...\n");
	}

	//Create new bam
	{
		printf("Creating new bam file in \"%s\"...\n", outbam);
		//init_empty_bam_header(orig_bam_f->bam_header_p->n_targets, recal_bam_header);
		out_bam_f = bam_fopen_mode(outbam, bam_f->bam_header_p, "w");
		assert(out_bam_f);
		bam_fwrite_header(out_bam_f->bam_header_p, out_bam_f);
		out_bam_f->bam_header_p = NULL;
		printf("New BAM initialized!...\n");
	}

#ifdef D_TIME_DEBUG
	end_time = omp_get_wtime();
	time_add_time_slot(D_SLOT_INIT, TIME_GLOBAL_STATS, end_time - init_time);
#endif

	//Flush
	fflush(stdout);

	//Create input list
	in_list = linked_list_new(0);
	assert(in_list);

	//Fill input list
#ifdef D_TIME_DEBUG
    init_time = omp_get_wtime();
#endif
	bytes = 1;
	for(i = 0; i < ALIG_LIST_IN_SIZE && bytes > 0; i++)
	{
		//Get next read from disk
		read = bam_init1();
		bytes = bam_read1(bam_f->bam_fd, read);

		if(bytes > 0 && read)
		{
			linked_list_insert_last(read, in_list);
		}
		else
		{
			//Destroy empty last read
			bam_destroy1(read);
		}
	}
#ifdef D_TIME_DEBUG
	end_time = omp_get_wtime();
	list_l = linked_list_size(in_list);
	if(list_l > 0)
		time_add_time_slot(D_SLOT_READ, TIME_GLOBAL_STATS, (double)(end_time - init_time)/(double)list_l);
#endif

	//Create context
	//printf("Creating context...\n");
	fflush(stdout);
	alig_init(&context, in_list, ref, ALIG_LEFT_ALIGN | ALIG_REFERENCE_PRIORITY);

	//Validate context
	err = alig_validate(&context);
	if(err)
	{
		fprintf(stderr, "ERROR: Error creating alignment context, error code = %d\n", err);
		fflush(stdout);
		return err;
	}

	//Get next reads
#ifdef D_TIME_DEBUG
    init_time = omp_get_wtime();
#endif
	err = alig_region_next(&context);
	if(err)
	{
		fprintf(stderr, "ERROR: Cannot obtain next region, error code = %d\n", err);
		fflush(stdout);
		return err;
	}
#ifdef D_TIME_DEBUG
	end_time = omp_get_wtime();
	time_add_time_slot(D_SLOT_NEXT, TIME_GLOBAL_STATS, (double)(end_time - init_time)/(double)context.last_readed_count);
#endif

	//Iterate
	while(context.last_readed_count > 0)
	{
		//Load reference
#ifdef D_TIME_DEBUG
		init_time = omp_get_wtime();
#endif
		err = alig_region_load_reference(&context);
		if(err)
		{
			fprintf(stderr, "ERROR: Failed to load reference for a region, error code = %d\n", err);
			fflush(stdout);
			return err;
		}
#ifdef D_TIME_DEBUG
		end_time = omp_get_wtime();
		time_add_time_slot(D_SLOT_REFERENCE_LOAD, TIME_GLOBAL_STATS, (double)(end_time - init_time)/(double)context.last_readed_count);
#endif

		//Generate haplotype list
#ifdef D_TIME_DEBUG
		init_time = omp_get_wtime();
#endif
		err = alig_region_haplotype_process(&context);
		if(err)
		{
			fprintf(stderr, "ERROR: Failed to get haplotype list, error code = %d\n", err);
			fflush(stdout);
			return err;
		}
#ifdef D_TIME_DEBUG
		end_time = omp_get_wtime();
		list_l = array_list_size(context.filtered_list);
		time_add_time_slot(D_SLOT_HAPLO_GET, TIME_GLOBAL_STATS, (double)(end_time - init_time)/(double)list_l);
#endif

		//Realign
#ifdef D_TIME_DEBUG
		init_time = omp_get_wtime();
#endif
		err = alig_region_indel_realignment(&context);
		if(err)
		{
			fprintf(stderr, "ERROR: Cannot align region, error code = %d\n", err);
			fflush(stdout);
			return err;
		}
#ifdef D_TIME_DEBUG
		end_time = omp_get_wtime();
		list_l = array_list_size(context.realign_list);
		if(list_l > 0)
		{
			end_time = (double)(end_time - init_time)/(double)list_l;
			list_l = array_list_size(context.haplo_list);
			if(list_l > 0)
				time_add_time_slot(D_SLOT_REALIG_PER_HAPLO, TIME_GLOBAL_STATS, end_time/(double)list_l);
		}
#endif

		//Print progress
		if(context.read_count / 10000 > aux_count)
		{
			aux_count = context.read_count / 10000;
			printf("Total alignments readed: %d\r", (int)context.read_count);
			fflush(stdout);
		}

		//Write processed to disk
#ifdef D_TIME_DEBUG
		init_time = omp_get_wtime();
#endif
		list_l = context.last_readed_count;
		for(i = 0; i < list_l; i++)
		{
			//Get read
			read = linked_list_get_first(in_list);
			assert(read);

			//Write to disk
			bam_write1(out_bam_f->bam_fd, read);

			//Free read
			bam_destroy1(read);

			//Update input list
			linked_list_remove_first(in_list);
		}
#ifdef D_TIME_DEBUG
		end_time = omp_get_wtime();
		time_add_time_slot(D_SLOT_WRITE, TIME_GLOBAL_STATS, (double)(end_time - init_time)/(double)list_l);
#endif

		//Refill buffer
#ifdef D_TIME_DEBUG
		init_time = omp_get_wtime();
#endif
		list_l = linked_list_size(in_list);

		//Especial case: List is empty
		if(list_l == 0)
		{
			linked_list_clear(in_list, NULL);
		}

		//Fill list
		for(i = list_l; i < ALIG_LIST_IN_SIZE && bytes > 0; i++)
		{
			read = bam_init1();
			bytes = bam_read1(bam_f->bam_fd, read);

			if(bytes > 0 && read)
			{
				linked_list_insert_last(read, in_list);
			}
			else
			{
				//Destroy empty last read
				bam_destroy1(read);
			}
		}
#ifdef D_TIME_DEBUG
		end_time = omp_get_wtime();
		list_l = linked_list_size(in_list) - list_l;
		if(list_l > 0)
			time_add_time_slot(D_SLOT_READ, TIME_GLOBAL_STATS, (double)(end_time - init_time)/(double)list_l);
#endif

		//Clear context
		alig_region_clear(&context);

		//Get next reads
#ifdef D_TIME_DEBUG
		init_time = omp_get_wtime();
#endif
		err = alig_region_next(&context);
		if(err)
		{
			fprintf(stderr, "ERROR: Cannot obtain next region, error code = %d\n", err);
			fflush(stdout);
			return err;
		}

#ifdef D_TIME_DEBUG
		end_time = omp_get_wtime();
		if(context.last_readed_count > 0)
			time_add_time_slot(D_SLOT_NEXT, TIME_GLOBAL_STATS, (double)(end_time - init_time)/(double)context.last_readed_count);
#endif
	}

	//Free context
	printf("\nDestroying context...\n");
	fflush(stdout);
	alig_destroy(&context);

	printf("\nClosing BAM file...\n");
	bam_fclose(bam_f);
	printf("BAM closed.\n");

	printf("\nClosing reference file...\n");
	genome_free(ref);
	printf("Reference closed.\n");

	printf("Closing \"%s\" BAM file...\n", outbam);
	bam_fclose(out_bam_f);
	printf("BAM closed.\n");

	printf("Realignment around indels DONE.\n");

	return NO_ERROR;
}


//					H0			H1			H2			H3			H4
//				sco	pos
//	ERRXXXX		23	97		0	16		321	34
//	ERRYYYY
//	ERRZZZZ
//

/**
 * PRIVATE FUNCTION. Obtain score tables from present region.
 */
static ERROR_CODE
alig_get_scores(alig_context_t *context)
{
	ERROR_CODE err;
	int i, j;

	//Reads
	bam1_t *read;
	size_t read_l;
	size_t indels;
	size_t bases;
	char *read_seq;
	char *comp_seq;
	char *read_quals;
	uint32_t *read_cigar;
	uint32_t read_left_cigar[MAX_CIGAR_LENGTH];
	size_t read_left_cigar_l;
	size_t read_disp_ref = SIZE_MAX;
	size_t read_disp_clip;

	//Haplotypes
	aux_indel_t *haplo;
	size_t haplo_list_l;

	//Reference
	char *ref_seq;
	size_t ref_pos_begin;
	size_t ref_pos_end;

	//Compare
	char *read_seq_ref;
	size_t read_seq_ref_l;
	uint32_t aux_cigar[MAX_CIGAR_LENGTH];
	size_t aux_cigar_l;
	uint32_t misses;
	uint32_t misses_sum;
	size_t best_pos;
	size_t curr_pos;
	int64_t init_pos;
	int64_t end_pos;

	//Lists
	array_list_t *read_list;
	size_t read_list_l;

	//Matrix
	uint32_t *m_scores;
	size_t *m_positions;
	size_t m_total;
	size_t m_ldim;

	assert(context);

	//Validate context
	err = alig_validate(context);
	if(err)
	{
		return err;
	}

	//Get list ptr
	read_list = context->realign_list;

	//Lenghts
	haplo_list_l = array_list_size(context->haplo_list);
	read_list_l = array_list_size(read_list);

	//Reference
	ref_seq = context->reference.reference;
	ref_pos_begin = context->reference.position;
	ref_pos_end = ref_pos_begin + context->reference.length;

	//Free previous scores matrix
	if(context->scores.m_total != 0)
	{
		if(context->scores.m_positions)
			free(context->scores.m_positions);
		if(context->scores.m_positions)
			free(context->scores.m_positions);
	}
	memset(&context->scores, 0, sizeof(alig_scores_t));

	//Allocate new scores
	m_ldim = haplo_list_l + 1;
	m_total = read_list_l * m_ldim;
	m_positions = (size_t *)malloc(m_total * sizeof(size_t));
	m_scores = (uint32_t *)malloc(m_total * sizeof(uint32_t));
	context->scores.m_ldim = m_ldim;
	context->scores.m_total = m_total;
	context->scores.m_positions = m_positions;
	context->scores.m_scores = m_scores;

	assert(m_positions);
	assert(m_scores);

	//Init matrices
	for(i = 0; i < m_total; i++)
	{
		*((size_t*)m_positions + i) = SIZE_MAX;
		*((uint32_t*)m_scores + i) = UINT32_MAX;
	}

	//Iterate reads
	for(i = 0; i < read_list_l; i++)
	{
		//Get read
		read = array_list_get(i, read_list);
		assert(read);
		read_cigar = bam1_cigar(read);
		assert(read_cigar);

		//Convert read sequence to string
#ifdef __SSE2__
		read_seq = (char *) _mm_malloc((read->core.l_qseq + 1) * sizeof(char), MEM_ALIG_SSE_SIZE);
		comp_seq = (char *) _mm_malloc((read->core.l_qseq + 1) * sizeof(char), MEM_ALIG_SSE_SIZE);
		read_seq_ref = (char *) _mm_malloc((read->core.l_qseq + 1) * sizeof(char), MEM_ALIG_SSE_SIZE);
		read_quals = (char *) _mm_malloc((read->core.l_qseq + 1) * sizeof(char), MEM_ALIG_SSE_SIZE);
#else
		read_seq = (char *) malloc((read->core.l_qseq + 1) * sizeof(char));
		comp_seq = (char *) malloc((read->core.l_qseq + 1) * sizeof(char));
		read_seq_ref = (char *) malloc((read->core.l_qseq + 1) * sizeof(char));
		read_quals = (char *) malloc((read->core.l_qseq + 1) * sizeof(char));
#endif
		new_sequence_from_bam_ref(read, read_seq, read->core.l_qseq + 1);

		//Get qualities
		new_quality_from_bam_ref(read, 0, read_quals, read->core.l_qseq + 1);

		//Get initial clip displacement
		cigar32_count_clip_displacement(read_cigar, read->core.n_cigar, &read_disp_clip);
		cigar32_count_nucleotides_not_clip(read_cigar, read->core.n_cigar, &read_l);

		//Leftalign cigar
		read_disp_ref = read->core.pos - context->reference.position;
		cigar32_leftmost(context->reference.reference + read_disp_ref, context->reference.length - read_disp_ref,
				read_seq, read->core.l_qseq, bam1_cigar(read), read->core.n_cigar,
				read_left_cigar, &read_left_cigar_l);

		//Get raw score with reference
		nucleotide_miss_qual_sum(ref_seq + read_disp_ref, read_seq + read_disp_clip, read_quals + read_disp_clip, read_l, comp_seq, &misses, &misses_sum);
		m_scores[i * m_ldim] = misses_sum;
		m_positions[i * m_ldim] = read->core.pos;

		//Dont iterate haplotypes if perfect reference match
		if(misses != 0)
		{
			//Iterate haplotypes
			for(j = 0; j < haplo_list_l; j++)
			{
				//Get haplotype
				haplo = array_list_get(j, context->haplo_list);
				assert(haplo);

				//Get initial position to iterate
				read_disp_ref = SIZE_MAX;
				init_pos = haplo->ref_pos - read->core.l_qseq;

				//Initial position must be inside reference range
				if(init_pos < ref_pos_begin)
					init_pos = ref_pos_begin;

				//Get end pos
				end_pos = ref_pos_end - read->core.l_qseq;
				if(end_pos < init_pos)	//Check valid range
					end_pos = init_pos;

				//Iterate positions
				for(curr_pos = init_pos; curr_pos < end_pos; curr_pos++)
				{
					//Create cigar for this haplotype
					err = cigar32_from_haplo(read_cigar, read->core.n_cigar, haplo, curr_pos, aux_cigar, &aux_cigar_l);
					if(err)
					{
						if(err == CIGAR_INVALID_INDEL)
						{
							//Haplotype is invalid
							sprintf(log_msg, "Invalid haplotype %d\n", j + 1);
							LOG_ERROR(log_msg);
							break;
						}
						else
						{
							//Cant displace anymore
							break;
						}
					}

					//Get relative displacement from reference
					if(curr_pos < haplo->ref_pos)
					{
						read_disp_ref = curr_pos - ref_pos_begin;
					}
					else
					{
						read_disp_ref = haplo->ref_pos - ref_pos_begin;
					}

					//Get haplotype reference transform
					cigar32_create_ref(aux_cigar, aux_cigar_l,
							ref_seq + read_disp_ref, context->reference.length - read_disp_ref,
							read_seq, read->core.l_qseq,
							read_seq_ref, &read_seq_ref_l);

					if(read_seq_ref_l != read->core.l_qseq)
					{
						sprintf(log_msg, "Read-Ref: %d, Read: %d\n", read_seq_ref_l, read->core.l_qseq);
						LOG_ERROR(log_msg);
					}
					//assert(read_seq_ref_l == read->core.l_qseq);

					//Compare and miss with haplotype
					nucleotide_miss_qual_sum(read_seq, read_seq_ref, read_quals, read_seq_ref_l, comp_seq, &misses, &misses_sum);

					//Better?
					if(m_scores[(i * m_ldim) + (j+1)] > misses_sum)
					{
						//Update scores
						m_scores[(i * m_ldim) + (j+1)] = misses_sum;
						m_positions[(i * m_ldim) + (j+1)] = curr_pos;
					}

					//Perfect match?
					if(misses == 0)
					{
						break;
					}
				} //Iterate positions
			} //Iterate haplotypes
		} //Scores != 0 if

		//Free
#ifdef __SSE2__
		_mm_free(read_seq);
		_mm_free(comp_seq);
		_mm_free(read_seq_ref);
		_mm_free(read_quals);
#else
		free(read_seq);
		free(comp_seq);
		free(read_seq_ref);
		free(read_quals);
#endif

	}

	return NO_ERROR;
}

/**
 * PRIVATE FUNCTION. Obtain alternative haplotype from generated score tables.
 */
static inline ERROR_CODE
alig_get_alternative_haplotype(alig_context_t *context, int *out_haplo_index, uint32_t *out_haplo_score, uint32_t *out_ref_score)
{
	ERROR_CODE err;

	int i, j;
	uint32_t best_score;
	int best_haplo_index;
	size_t haplo_l;
	size_t bam_l;

	//Scores
	uint32_t *m_scores;
	size_t *m_positions;
	uint32_t aux_score;
	alig_scores_t *scores;

	assert(context);

	//Validate context
	err = alig_validate(context);
	if(err)
	{
		return err;
	}

	//Get scores ptr
	scores = &context->scores;
	assert(scores);

	//Valid scores?
	if(scores->m_total == 0)
		return ALIG_INVALID_SCORES;

	//Set lengths
	haplo_l = scores->m_ldim;
	bam_l = scores->m_total / scores->m_ldim;

	//Set matrix
	m_scores = scores->m_scores;
	m_positions = scores->m_positions;
	assert(m_scores);
	assert(m_positions);

	//Output ref score
	if(out_ref_score)
	{
		//Sum scores
		aux_score = 0;
		for(i = 0; i < bam_l; i++)
		{
			if(m_scores[(i * haplo_l)] != UINT32_MAX)	//If valid score
			{
				aux_score += m_scores[(i * haplo_l)];
			}
		}

		//Output
		*out_ref_score = aux_score;
	}


	//Find best haplotype
	best_score = UINT32_MAX;
	for(j = 1; j < haplo_l; j++)
	{
		//Sum scores
		aux_score = 0;
		for(i = 0; i < bam_l; i++)
		{
			if(m_scores[(i * haplo_l) + j] < m_scores[(i * haplo_l)])	//If haplotype score is better
			{
				//Haplotype is better so add its score
				aux_score += m_scores[(i * haplo_l) + j];
			}
			else
			{
				//Reference is better so add its score
				aux_score += m_scores[(i * haplo_l)];
			}
		}

		if(aux_score < best_score)
		{
			best_score = aux_score;
			best_haplo_index = j;
		}

		//ERASE
		/*{
			//Print haplotype score
			printf("H%d === total score = %d\n", j, aux_score);
		}*/

		//Found perfect match
		if(best_score == 0)
			break;
	}

	//Set output
	if(out_haplo_score)
		*out_haplo_score = best_score;

	if(out_haplo_index)
		*out_haplo_index = best_haplo_index;

	return NO_ERROR;
}

/**
 * PRIVATE FUNCTION. Realign around indels using an alternative haplotype.
 */
static ERROR_CODE
alig_indel_realign_from_haplo(alig_context_t *context, size_t alt_haplo_index)
{
	ERROR_CODE err;
	int i;

	//Read
	bam1_t *read;
	size_t read_indels;
	char *read_seq = NULL;
	char *read_seq_ref = NULL;
	char *read_quals = NULL;
	size_t read_seq_ref_l;
	size_t read_disp_ref = SIZE_MAX;

	//Reference
	char *ref_seq = NULL;
	size_t ref_length;

	//Comparation
	char *comp_aux;
	uint32_t misses;
	uint32_t misses_sum;

	//Scores
	uint32_t ref_score;
	uint32_t h_score;
	uint32_t min_score;

	//Matrix
	uint32_t *m_scores;
	size_t *m_positions;
	size_t m_ldim;

	//New cigar
	uint32_t best_cigar[MAX_CIGAR_LENGTH];
	size_t best_cigar_l;
	size_t bases;
	size_t best_pos;
	aux_indel_t *haplo;

	//List
	array_list_t *list;
	size_t list_l;

	assert(context);
	assert(alt_haplo_index > 0);

	//Validate context
	err = alig_validate(context);
	if(err)
	{
		return err;
	}

	//Set reference
	ref_seq = context->reference.reference;
	ref_length = context->reference.length;
	assert(ref_seq);

	//Set matrix
	m_scores = context->scores.m_scores;
	m_positions = context->scores.m_positions;
	m_ldim = context->scores.m_ldim;
	assert(m_scores);
	assert(m_positions);

	//Set list
	list = context->realign_list;

	//Set haplo
	haplo = array_list_get(alt_haplo_index - 1, context->haplo_list);
	assert(haplo);

	//Logging
	cigar32_to_string(&haplo->indel, 1, aux_msg);
	sprintf(log_msg, "Realign, alternative haplotype = %s:%d, %d reads\n", aux_msg, haplo->ref_pos + 1, array_list_size(list));
	LOG_INFO("************************************\n");
	LOG_INFO(log_msg);
	if(context->flags & ALIG_ORIGINAL_PRIORITY)
		LOG_INFO("Original cigar have priority\n");
	if(context->flags & ALIG_REFERENCE_PRIORITY)
		LOG_INFO("Reference haplotype have priority\n");

	//Iterate list to realign
	list_l = array_list_size(list);
	for(i = 0; i < list_l; i++)
	{
		//Get read
		read = array_list_get(i, list);
		assert(read);

		//Only if one or less indels
		cigar32_count_indels(bam1_cigar(read), read->core.n_cigar, &read_indels);
		if(read_indels > 1)
		{
			continue;
		}

		//Get scores
		ref_score = m_scores[(i * m_ldim)];
		h_score = m_scores[(i * m_ldim) + alt_haplo_index];

		//Get original score
		{
			//Convert read sequence to string
			read_seq = (char *) malloc((read->core.l_qseq + 1) * sizeof(char));
			read_seq_ref = (char *) malloc((read->core.l_qseq + 1) * sizeof(char));
			comp_aux = (char *) malloc((read->core.l_qseq + 1) * sizeof(char));
			new_sequence_from_bam_ref(read, read_seq, read->core.l_qseq + 1);

			//Get qualities
			read_quals = (char *) malloc((read->core.l_qseq + 1) * sizeof(char));
			new_quality_from_bam_ref(read, 0, read_quals, read->core.l_qseq + 1);

			read_disp_ref = read->core.pos - context->reference.position;

			//Get haplotype reference transform
			cigar32_create_ref(bam1_cigar(read), read->core.n_cigar, ref_seq + read_disp_ref, ref_length - read_disp_ref, read_seq, read->core.l_qseq, read_seq_ref, &read_seq_ref_l);
			assert(read_seq_ref_l == read->core.l_qseq);

			//Compare and miss with haplotype
			nucleotide_miss_qual_sum(read_seq_ref, read_seq, read_quals, read->core.l_qseq, comp_aux, &misses, &misses_sum);

			//Release resources
			free(read_seq);
			free(read_seq_ref);
			free(comp_aux);
			free(read_quals);

			//Is original cigar better?
			min_score = min(ref_score, h_score);
			if(((context->flags & ALIG_ORIGINAL_PRIORITY) && misses_sum <= min_score)	//Original cigar have priority
					|| misses_sum < min_score )	//Realigned cigar have priority
			{
				//Logging
				sprintf(log_msg, "NOT Realigned read: %s - %d (original is better, O:%d vs H0:%d vs HA:%d)\n",
						bam1_qname(read), read->core.pos, misses_sum, ref_score, h_score);
				LOG_INFO(log_msg);

				//Original is better, dont change!

				continue;
			}
		}

		//Is reference matching better?
		if( ((context->flags & ALIG_REFERENCE_PRIORITY) && ref_score <= h_score)	//Reference have priority
				|| ref_score < h_score )	//Haplotype have priority
		{
			//Unclip cigar
			//cigar32_unclip(bam1_cigar(read), read->core.n_cigar, best_cigar, &best_cigar_l);

			//Count unclipped bases
			cigar32_count_nucleotides_not_clip(bam1_cigar(read), read->core.n_cigar, &bases);

			//Create new cigar with '$bases'M
			best_cigar[0] = bases << BAM_CIGAR_SHIFT;	//ex: 108M
			best_cigar_l = 1;

			//Reclip cigar
			cigar32_reclip(bam1_cigar(read), read->core.n_cigar, best_cigar, best_cigar_l, best_cigar, &best_cigar_l);

			//Set best position
			best_pos = m_positions[(i * m_ldim)];

			//Logging
			sprintf(log_msg, "Realigned read: %s - %d (to reference O:%d vs H0:%d vs HA:%d)\n",
					bam1_qname(read), read->core.pos + 1, misses_sum, ref_score, h_score);
			LOG_INFO(log_msg);
		}
		else //Haplotype is better
		{
			//Set best position
			best_pos = m_positions[(i * m_ldim) + alt_haplo_index];

			//Create cigar for this haplotype
			err = cigar32_from_haplo(bam1_cigar(read), read->core.n_cigar, haplo, best_pos, best_cigar, &best_cigar_l);
			if(err)
			{
				fprintf(stderr, "Warning: Invalid cigar creation from haplotype, error code: %d\n", err);
			}

			//If haplotype position is minor then adjust
			if(best_pos > haplo->ref_pos)
			{
				best_pos = haplo->ref_pos;	//ex: 3D13M => 1D13M
			}

			//Logging
			sprintf(log_msg, "Realigned read: %s - %d (to alternative haplotype O:%d vs H0:%d vs HA:%d)\n",
					bam1_qname(read), read->core.pos + 1, misses_sum, ref_score, h_score);
			LOG_INFO(log_msg);

		}

		//Make realign effective
		{
			//Change CIGAR
			cigar32_replace(read, best_cigar, best_cigar_l);
			read->core.pos = best_pos;

			//Change MAPQ (add 10)
			if(read->core.qual != 255)
				read->core.qual = (uint8_t)min((int)read->core.qual + 10, 254);
		}

	}

	return NO_ERROR;
}


