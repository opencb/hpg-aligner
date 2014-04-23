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

typedef struct {
	bam1_t **buffer;
	size_t head_idx;
	size_t readed;
	size_t processed;
	size_t size;
	omp_lock_t lock;
	omp_lock_t writer;
	omp_lock_t reader;
	int end_condition;
} circular_buffer_t;

uint64_t cigar_changed = 0;

char log_msg[1024];
char aux_msg[512];

static inline ERROR_CODE alig_aux_write_to_disk(array_list_t *write_buffer, bam_file_t *output_bam_f, uint8_t force) __ATTR_HOT;
static inline ERROR_CODE alig_aux_read_from_disk(circular_buffer_t *read_buffer, bam_file_t *input_bam_f) __ATTR_HOT;
static inline ERROR_CODE alig_region_filter_read(bam1_t *read, array_list_t *list, alig_context_t *context) __ATTR_HOT __ATTR_INLINE;

/**
 * CONTEXT
 */

/**
 * Initialize empty realignment data structure.
 */
ERROR_CODE
alig_init(alig_context_t *context, genome_t *genome, uint8_t flags)
{
	ERROR_CODE err;

	assert(context);
	assert(genome);

	//Set all to 0 / NULL
	memset(context, 0, sizeof(alig_context_t));

	//Set fields
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
		array_list_free(context->haplo_list, NULL);
	}

	//Free reference
	if(context->reference.reference)
	{
		free(context->reference.reference);
	}

	//Free scores
	if(context->scores.m_positions)
		free(context->scores.m_positions);
	if(context->scores.m_scores)
		free(context->scores.m_scores);

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
alig_region_next(bam1_t **v_bams, size_t v_bams_l, int force_incomplete, alig_context_t *context)
{
	ERROR_CODE err, ret = NO_ERROR;
	int i;

	//List
	size_t list_l;

	//Read
	size_t reads;
	bam1_t* read = NULL;
	int32_t read_pos;
	size_t indels;
	size_t interval_read_begin;
	size_t interval_read_end;

	//Aux list
	static array_list_t *aux_list = NULL;

	assert(context);

	//Validate context
	err = alig_validate(context);
	if(err)
	{
		return err;
	}

	//Set last readed to 0
	context->last_readed_count = 0;

	//Reset interval
	memset(&context->region, 0, sizeof(alig_region_t));

	//Init aux list
	if(!aux_list)
	{
		aux_list = array_list_new(ALIG_LIST_NEXT_SIZE, 1.2f, 0);
	}
	array_list_clear(aux_list, NULL);

	//Read alignments until interval reached
	list_l = v_bams_l;
	reads = 0;
	for(i = 0; i < list_l; i++)
	{
		//Get next read
		read = v_bams[i];
		assert(read);

		//Filter read to obtain regions
		err = alig_region_filter_read(read, aux_list, context);
		if(err)
		{
			if(err == ALIG_PAST_INTERVAL)
			{
				//Interval obtained
				break;
			}
			else
			{
				//Error
				LOG_ERROR("Failed to filter read in context");
				return err;
			}
		}

		//Increment counter
		reads++;
	}

	//Interval is found?
	if(err != ALIG_PAST_INTERVAL)
	{
		if(force_incomplete == 0)
		{
			//Incomplete or inexistent interval
			return ALIG_INCOMPLETE_INTERVAL;
		}

		//Set return
		ret = ALIG_INCOMPLETE_INTERVAL;
	}

	//Update progress
	context->read_count += reads;
	context->last_readed_count = reads;

	//There is an valid interval?
	if(context->region.valid)
	{
		//Logging
		LOG_INFO("************************************\n");
		sprintf(log_msg, "INTERVAL %d - %d\n", context->region.start_pos + 1, context->region.end_pos + 1);
		LOG_INFO(log_msg);

		//Save interval reads in filtered list
		list_l = array_list_size(aux_list);
		for(i = 0; i < list_l; i++)
		{
			//Get read
			read = array_list_get(i, aux_list);
			assert(read);

			if(region_bam_overlap(read, &context->region))
			{
				//Logging
				sprintf(log_msg, "%s \t%d:%d\n", bam1_qname(read), read->core.tid + 1, read->core.pos + 1);
				LOG_INFO(log_msg);

				//Add to filtered read list
				array_list_insert(read, context->filtered_list);
			}
		} //For
	} //Valid region if

	return ret;
}

static inline ERROR_CODE
alig_region_filter_read(bam1_t *read, array_list_t *list, alig_context_t *context)
{
	int err;

	assert(read);
	assert(context);

	//Read
	int32_t read_pos;
	size_t interval_read_begin;
	size_t interval_read_end;

	//FILTERS: MAP QUALITY = 0, BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP, DIFFERENT MATE CHROM, 1 INDEL, CIGAR PRESENT
	if( read->core.qual != 0
		&& !(read->core.flag & BAM_DEF_MASK)
		&& read->core.mtid == read->core.tid
		&& read->core.n_cigar != 0
		)
	{
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
		else
		{
			//INTERVAL CASE

			//Is alignment inside interval?
			if (context->region.end_pos > read_pos && context->region.chrom == read->core.tid)
			{
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
			else //Alignment past interval
			{
				//Interval end!
				return ALIG_PAST_INTERVAL;
			}
		}//Cases if

		//Insert to list
		array_list_insert(read, list);

	}	// Filters if

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
		aux_pos_begin = ref_pos_begin + ALIG_REFERENCE_CORRECTION_OFFSET;
		aux_pos_end = ref_pos_end + ALIG_REFERENCE_CORRECTION_OFFSET;

		//Get reference
		ref_seq = (char *) malloc(((ref_pos_end - ref_pos_begin) + 2) * sizeof(char));

		assert(ref_seq);
		flag = (uint32_t) read->core.flag;
		genome_read_sequence_by_chr_index(ref_seq, (flag & BAM_FREVERSE) ? 1 : 0, (unsigned int)read->core.tid, &aux_pos_begin, &aux_pos_end, context->genome);

		ref_length = (ref_pos_end - ref_pos_begin);

		//Set context reference
		if(context->reference.reference)
		{
			//Free previous reference
			free(context->reference.reference);
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

				read_seq = (char *) malloc((read_seq_l + 1) * sizeof(char));
				comp_seq = (char *) malloc((read_seq_l + 1) * sizeof(char));

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
				//ref_miss = 1;

				//If match reference perfectly, realign to reference and extract from process list
				/*if(ref_miss == 0)
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
					sprintf(log_msg, "%s \t%d:%d match reference perfectly\n", bam1_qname(read), read->core.tid + 1, read->core.pos + 1);
					LOG_INFO(log_msg);
				}
				else*/
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
				free(read_seq);
				free(comp_seq);

			} //Valid cigar end if

		} //Haplotypes while

	} //Region if

	//Print info
	/*{
	  	printf("-----------\n");
		printf("Haplotypes in region: %d:%d\n", context->region.chrom, context->region.start_pos);
		static char str[20];
		int h;
		haplo_l = array_list_size(context->haplo_list);
		for(h = 0; h < haplo_l; h++)
		{
			aux_haplo = array_list_get(h, context->haplo_list);
			assert(aux_haplo);

			cigar32_to_string(&aux_haplo->indel, 1, str);
			printf("Pos:%d,  Indel:%s\n", aux_haplo->ref_pos, str);
		}
	}*/

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

	//Get best haplotype (alternative to H0)
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

	//Set last region
	context->last_region = context->region;

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
		free(context->reference.reference);
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

//Private function
static int compare_pos(const void *item1, const void *item2) {
	const bam1_t *read1 = *(const bam1_t **)item1;
	const bam1_t *read2 = *(const bam1_t **)item2;

	return (int) (read1->core.pos - read2->core.pos);
}

/**
 * Indel realign one file.
 */
ERROR_CODE
alig_bam_file(char *bam_path, char *ref_name, char *ref_path, char *outbam)
{
	ERROR_CODE err;
	//	unsigned char outbam[30] = "output.bam";
	int i, it;
	//int single_threaded = 0;

	//Times
	double init_time, end_time;
	double init_time_p, end_time_p;
	double init_time_p2, end_time_p2;
	double init_time_it, end_time_it;
	double aux_time, aux_time2;

	//Files
	bam_file_t *bam_f;
	bam_file_t *out_bam_f;
	genome_t* ref;
	int bytes;

	//Read
	bam1_t *read;
	size_t bam_processed = 0;
	bam1_t **v_reads;
	size_t v_reads_l;

	//Lists
	array_list_t *swap_buffer_ptr;
	array_list_t *write_buffer;
	array_list_t *proc_buffer;
	array_list_t *read_buffer;
	size_t list_l;

	//Input circular buffer
	circular_buffer_t in_buffer;

	//Alignment context
	alig_context_t context;

	//aux
	int aux_count = 0;
	int filled = 0;

	assert(bam_path);
	assert(ref_name);
	assert(ref_path);

	//Multhread info
	printf("Using %d threads\n", omp_get_max_threads());

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
		ref = genome_new(ref_name, ref_path);
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

	//Create buffer for write to disk
	write_buffer = array_list_new(ALIG_LIST_NEXT_SIZE, 1.2, 0);
	proc_buffer = array_list_new(ALIG_LIST_NEXT_SIZE, 1.2, 0);
	read_buffer = array_list_new(ALIG_LIST_NEXT_SIZE, 1.2, 0);
	assert(write_buffer);
	assert(proc_buffer);
	assert(read_buffer);

	//Init input buffer
	memset(&in_buffer, 0, sizeof(circular_buffer_t));
	in_buffer.buffer = (bam1_t **)malloc(ALIG_LIST_IN_SIZE * sizeof(bam1_t *));
	in_buffer.size = ALIG_LIST_IN_SIZE;
	omp_init_lock(&in_buffer.lock);
	omp_init_lock(&in_buffer.writer);
	omp_init_lock(&in_buffer.reader);
	omp_set_lock(&in_buffer.writer);

	//Create context
	//printf("Creating context...\n");
	fflush(stdout);
	alig_init(&context, ref, ALIG_LEFT_ALIGN | ALIG_REFERENCE_PRIORITY);

	//Validate context
	err = alig_validate(&context);
	if(err)
	{
		fprintf(stderr, "ERROR: Error creating alignment context, error code = %d\n", err);
		fflush(stdout);
	}

	//Init vars
	it = 0;
	init_time_it = omp_get_wtime();

	//Timing
#ifdef D_TIME_DEBUG
	init_time_it = 0;
#endif

	//Init multithreading
	omp_set_dynamic(1);
	omp_set_nested(1);
	if(omp_get_max_threads() < 3)
	{
		//Only 1 available thread
		//printf("Minimum threads required is 2, setting number of threads to 2\n");
		omp_set_num_threads(3);
		//single_threaded = 1;
	}
	#pragma omp parallel private(i, bytes, read, v_reads, v_reads_l, list_l, filled, err, aux_time)
	{
		//while(in_buffer.end_condition == 0 || in_buffer.readed > 0)
		{
			#pragma omp  sections
			{
				//Read section
				#pragma omp section
				{
					do
					{
#ifdef D_TIME_DEBUG
						aux_time = omp_get_wtime();
						aux_time2 = aux_time;
	#endif
						//Fill in_list
						alig_aux_read_from_disk(&in_buffer, bam_f);

#ifdef D_TIME_DEBUG
						aux_time2 = omp_get_wtime() - aux_time2;
						time_add_time_slot(D_SLOT_IT_READ, TIME_GLOBAL_STATS, aux_time2);
#endif
						//Clear context
						alig_region_clear(&context);

						//Load next reads until chrom change
						array_list_clear(read_buffer, NULL);
						filled = in_buffer.readed - in_buffer.processed;
						int32_t chrom;
						int chrom_changed = 0;
						int force_region = 0;
						chrom = in_buffer.buffer[(in_buffer.head_idx + in_buffer.processed) % in_buffer.size]->core.tid;
						for(i = 0; i < filled; i++)
						{
							//Get read
							read = in_buffer.buffer[(in_buffer.head_idx + in_buffer.processed + i) % in_buffer.size];

							//Check if chrom changed
							if(chrom == read->core.tid)
							{
								//Same chrom so continue
								array_list_insert(read, read_buffer);
							}
							else
							{
								//Chrom changed!
								chrom_changed = 1;
								break;
							}
						}

#ifdef D_TIME_DEBUG
						aux_time2 = omp_get_wtime();
#endif
						//Load next region
						v_reads = (bam1_t **)read_buffer->items;
						v_reads_l = read_buffer->size;
						err = alig_region_next(v_reads, v_reads_l, 1, &context);
						if(err)
						{
							if(err != ALIG_INCOMPLETE_INTERVAL)
							{
								LOG_ERROR("Cannot obtain next region\n");
								printf(stderr, "ERROR: Cannot obtain next region in buffer of size %d, error code = %d\n", err);
								fflush(stdout);
								omp_set_lock(&in_buffer.lock);
								in_buffer.end_condition = 1;
								omp_unset_lock(&in_buffer.lock);
								break;
							}
							else
							{
								//Provoked by chrom change or end contition?
								if(chrom_changed == 0 && in_buffer.end_condition == 0)
								{
									//Blocked
									printf("Forcing region extraction on whole input buffer, i need bigger input buffer!\n");
									LOG_WARN("Forcing region extraction on whole input buffer, i need bigger input buffer!\n");
								}
							}
						}
#ifdef D_TIME_DEBUG
						aux_time2 = omp_get_wtime() - aux_time2;
						if(context.last_readed_count > 0)
							time_add_time_slot(D_SLOT_NEXT, TIME_GLOBAL_STATS, aux_time2 / (double)context.last_readed_count);
#endif

						while(context.last_readed_count != 0)
						{
							//printf("Realigning %d\n", context.last_readed_count);
							//Load reference
#ifdef D_TIME_DEBUG
							init_time_p = omp_get_wtime();
#endif
							err = alig_region_load_reference(&context);
							if(!err)
							{
#ifdef D_TIME_DEBUG
								end_time_p = omp_get_wtime();
								time_add_time_slot(D_SLOT_REFERENCE_LOAD, TIME_GLOBAL_STATS, (double)(end_time_p - init_time_p)/(double)context.last_readed_count);
#endif

								//Generate haplotype list
#ifdef D_TIME_DEBUG
								init_time_p = omp_get_wtime();
#endif
								err = alig_region_haplotype_process(&context);
								if(!err)
								{
#ifdef D_TIME_DEBUG
									end_time_p = omp_get_wtime();
									time_add_time_slot(D_SLOT_HAPLO_GET, TIME_GLOBAL_STATS, (double)(end_time_p - init_time_p)/(double)context.last_readed_count);
#endif
									//Realign
#ifdef D_TIME_DEBUG
									init_time_p = omp_get_wtime();
#endif
									err = alig_region_indel_realignment(&context);
									if(!err)
									{
#ifdef D_TIME_DEBUG
										end_time_p = omp_get_wtime();
										if(context.last_readed_count > 0)
												time_add_time_slot(D_SLOT_REALIG_PER_HAPLO, TIME_GLOBAL_STATS, (double)(end_time_p - init_time_p)/(double)context.last_readed_count);
#endif
									}
									else
									{
										fprintf(stderr, "ERROR: Cannot align region, error code = %d\n", err);
										fflush(stdout);
									}
								}
								else
								{
									fprintf(stderr, "ERROR: Failed to get haplotype list, error code = %d\n", err);
									fflush(stdout);
								}
							}
							else
							{
								fprintf(stderr, "ERROR: Failed to load reference for a region, error code = %d\n", err);
								fflush(stdout);
							}

							//Update buffer
							in_buffer.processed += context.last_readed_count;

							//Update buffer pointer
							v_reads += context.last_readed_count;
							v_reads_l -= context.last_readed_count;

							//Clear context
							alig_region_clear(&context);

#ifdef D_TIME_DEBUG
							aux_time2 = omp_get_wtime();
#endif
							//Get next region
							err = alig_region_next(v_reads, v_reads_l, chrom_changed || in_buffer.end_condition, &context);
							if(err)
							{
								if(err != ALIG_INCOMPLETE_INTERVAL)
								{
									LOG_ERROR("Cannot obtain next region\n");
									fprintf(stderr, "ERROR: Cannot obtain next region in buffer of size %d, error code = %d\n", err);
									fflush(stdout);
									omp_set_lock(&in_buffer.lock);
									in_buffer.end_condition = 1;
									omp_unset_lock(&in_buffer.lock);
									break;
								}
							}
#ifdef D_TIME_DEBUG
							aux_time2 = omp_get_wtime() - aux_time2;
							if(context.last_readed_count > 0)
								time_add_time_slot(D_SLOT_NEXT, TIME_GLOBAL_STATS, aux_time2 / (double)context.last_readed_count);
#endif
						}//last_readed != 0 while

						//Get processed number
						list_l = in_buffer.processed;

						//Fill write buffer
						if(list_l > 0)
						{
							//Bound list
							for(i = 0; i < list_l; i++)
							{
								//Get next read
								//omp_set_lock(&in_buffer.lock);
								read = in_buffer.buffer[(in_buffer.head_idx + i) % in_buffer.size];
								//omp_unset_lock(&in_buffer.lock);
								assert(read);

								//Add to write list
								array_list_insert(read, proc_buffer);
							}

							//Update counters
							bam_processed += list_l;

							//Update buffer
							in_buffer.head_idx = (in_buffer.head_idx + list_l) % in_buffer.size;
							in_buffer.readed -= list_l;
							in_buffer.processed = 0;
						} //Fill write buffer

#ifdef D_TIME_DEBUG
						aux_time = omp_get_wtime() - aux_time;
						time_add_time_slot(D_SLOT_IT_PROCESS, TIME_GLOBAL_STATS, aux_time);
#endif
						omp_set_lock(&in_buffer.reader);
						swap_buffer_ptr = write_buffer;
						write_buffer = proc_buffer;
						proc_buffer = swap_buffer_ptr;
						omp_unset_lock(&in_buffer.writer);
					} while(/*single_threaded == 0 && */(in_buffer.end_condition == 0 || in_buffer.readed > 0));
				}//Read section

				//Write section
				#pragma omp section
				{
					do
					{
						omp_set_lock(&in_buffer.writer);
#ifdef D_TIME_DEBUG
						aux_time = omp_get_wtime();
#endif
						alig_aux_write_to_disk(write_buffer, out_bam_f, 1);

						//Print progress
						if(bam_processed / 10000 > aux_count)
						{
							aux_count = bam_processed / 10000;
							printf("Total alignments readed: %d\r", bam_processed);
							fflush(stdout);
						}
						it++;
#ifdef D_TIME_DEBUG
						aux_time = omp_get_wtime() - aux_time;
						time_add_time_slot(D_SLOT_IT_WRITE, TIME_GLOBAL_STATS, aux_time);
#endif
						//Accept more
						omp_unset_lock(&in_buffer.reader);
					}while(/*single_threaded == 0 &&*/ (in_buffer.end_condition == 0 || in_buffer.readed > 0));//Write section
				}

			}//OMP SECTIONS
		} //While

	}//OMP PARALLEL

	//Write last reads
	alig_aux_write_to_disk(write_buffer, out_bam_f, 1);

	//Info
	printf("\nReads processed from original BAM: %d\n", bam_processed);
	printf("Iterations: %d\n", it);

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

	//Free
	array_list_free(write_buffer, NULL);
	array_list_free(proc_buffer, NULL);
	array_list_free(read_buffer, NULL);
	omp_destroy_lock(&in_buffer.lock);
	free(in_buffer.buffer);

	if(err)
		return err;
	else
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
/*static*/ ERROR_CODE
alig_get_scores(alig_context_t *context)
{
	ERROR_CODE err;
	int i, j;

	//Reads
	bam1_t *read;
	size_t index;

	//Matrix
	uint32_t *m_scores;
	size_t *m_positions;
	size_t m_total;
	size_t m_ldim;

	//Lists
	array_list_t *read_list;
	size_t read_list_l;
	size_t haplo_list_l;

	//Score vector
	uint32_t *v_scores;
	size_t *v_positions;
	size_t v_total;

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

	//Free previous scores matrix
	if(context->scores.m_total != 0)
	{
		if(context->scores.m_positions)
			free(context->scores.m_positions);
		if(context->scores.m_scores)
			free(context->scores.m_scores);
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
	//#pragma omp parallel private(read, err, j, i, v_total, v_positions, v_scores, index) firstprivate(read_list_l)
	{
		//Allocate new scores
		v_total = context->scores.m_ldim;
		v_positions = (size_t *)malloc(v_total * sizeof(size_t));
		v_scores = (uint32_t *)malloc(v_total * sizeof(uint32_t));

		//#pragma omp for
		for(i = 0; i < read_list_l; i++)
		{
			//Get read
			read = array_list_get(i, read_list);
			assert(read);

			//Get scores
			err = alig_get_scores_from_read(read, context, v_scores, v_positions);

			//Merge scores
			index = i * v_total;
			memcpy((uint32_t*)m_scores + index, v_scores, v_total * sizeof(uint32_t));
			memcpy((size_t*)m_positions + index, v_positions, v_total * sizeof(size_t));

		}//FOR reads

		//Free
		free(v_scores);
		free(v_positions);
	}

	//Print table
	/*{
		int r, h;
		char indel[10];
		aux_indel_t *haplo;
		printf("------------------\n");
		printf("Table for region %d:%d:%d\n", context->region.chrom + 1, context->region.start_pos + 1, context->region.end_pos + 1);

		//Print haplos
		printf("%20s %10s ", " ", "H0");
		haplo_list_l = array_list_size(context->haplo_list);
		for(h = 0; h < haplo_list_l; h++)
		{
			haplo = array_list_get(h, context->haplo_list);
			assert(haplo);
			cigar32_to_string(&haplo->indel, 1, indel);
			printf("%10s ", indel);
		}
		printf("\n");

		//Print scores
		read_list_l = array_list_size(context->realign_list);
		for(r = 0; r < read_list_l; r++)
		{
			read = array_list_get(r, read_list);
			assert(read);
			printf("%20s ", bam1_qname(read));
			for(h = 0; h < haplo_list_l + 1; h++)
			{
				index = r * context->scores.m_ldim;
				printf("%10d ", m_scores[index + h]);
			}
			printf("\n");
		}

		//Print positions
		printf("---\n");
		read_list_l = array_list_size(context->realign_list);
		for(r = 0; r < read_list_l; r++)
		{
			read = array_list_get(r, read_list);
			assert(read);
			printf("%20s ", bam1_qname(read));
			for(h = 0; h < haplo_list_l + 1; h++)
			{
				index = r * context->scores.m_ldim;
				printf("%10d ", m_positions[index + h]);
			}
			printf("\n");
		}
	}*/

	return NO_ERROR;
}

/*static*/ ERROR_CODE
alig_get_scores_from_read(bam1_t *read, alig_context_t *context, uint32_t *v_scores, size_t *v_positions)
{
	int i, err;

	//Reads
	char *read_seq;
	char *comp_seq;
	char *quals_seq;
	size_t read_l;
	size_t indels;
	size_t bases;
	size_t read_disp_ref;
	size_t read_disp_clip;

	//Compare
	char *read_seq_ref;
	size_t read_seq_ref_l;
	uint32_t misses;
	uint32_t misses_sum;
	size_t best_pos;
	size_t curr_pos;
	int64_t init_pos;
	int64_t end_pos;

	//Cigars
	uint32_t *read_cigar;
	uint32_t *read_left_cigar;
	size_t read_left_cigar_l;
	uint32_t *aux_cigar;
	size_t aux_cigar_l;

	//Haplotype
	aux_indel_t *haplo;
	size_t haplo_list_l;

	//Reference
	alig_reference_t *reference;
	char *ref_seq;
	size_t ref_pos_begin;
	size_t ref_pos_end;

	//Score vector
	size_t v_total;

	assert(read);
	assert(context);

	//Lengths
	haplo_list_l = array_list_size(context->haplo_list);

	//Get read
	read_cigar = bam1_cigar(read);
	assert(read_cigar);

	//Get reference
	reference = &context->reference;
	ref_seq = reference->reference;
	ref_pos_begin = reference->position;
	ref_pos_end = ref_pos_begin + reference->length;

	//Allocate new read
	read_l = read->core.l_qseq;
	read_seq = (char *) malloc((read_l + 1) * sizeof(char));
	comp_seq = (char *) malloc((read_l + 1) * sizeof(char));
	read_seq_ref = (char *) malloc((read_l + 1) * sizeof(char));
	quals_seq = (char *) malloc((read_l + 1) * sizeof(char));

	//Allocate cigars
	aux_cigar = (uint32_t *)malloc(MAX_CIGAR_LENGTH * sizeof(uint32_t));
	read_left_cigar = (uint32_t *)malloc(MAX_CIGAR_LENGTH * sizeof(uint32_t));

	//Init scores
	v_total = context->scores.m_ldim;
	for(i = 0; i < v_total; i++)
	{
		*((size_t*)v_positions + i) = SIZE_MAX;
		*((uint32_t*)v_scores + i) = UINT32_MAX;
	}

	//Convert read sequence to string
	new_sequence_from_bam_ref(read, read_seq, read_l + 1);

	//Get qualities
	new_quality_from_bam_ref(read, 0, quals_seq, read->core.l_qseq + 1);

	//Get initial clip displacement
	cigar32_count_clip_displacement(read_cigar, read->core.n_cigar, &read_disp_clip);
	cigar32_count_nucleotides_not_clip(read_cigar, read->core.n_cigar, &read_l);

	//Get read displacement to reference
	read_disp_ref = read->core.pos  - reference->position;

	//Get clean read
	read_left_cigar[0] = (read_l << BAM_CIGAR_SHIFT) + BAM_CMATCH;
	cigar32_reclip(bam1_cigar(read), read->core.n_cigar, read_left_cigar, 1, aux_cigar, &aux_cigar_l);
	cigar32_create_ref(aux_cigar, aux_cigar_l,
			ref_seq + read_disp_ref, context->reference.length - read_disp_ref,
			read_seq, read->core.l_qseq,
			read_seq_ref, &read_seq_ref_l);

	//Get raw score with reference
	nucleotide_miss_qual_sum(read_seq_ref, read_seq, quals_seq, read_l, comp_seq, &misses, &misses_sum);
	v_scores[0] = misses_sum;
	v_positions[0] = read->core.pos;

	//Dont iterate haplotypes if perfect reference match
	//printf("%s - %d - %d\n", bam1_qname(read), misses, misses_sum);
	if(misses != 0)
	{
		//Leftalign cigar
		cigar32_leftmost(reference->reference + read_disp_ref, reference->length - read_disp_ref,
				read_seq, read->core.l_qseq, bam1_cigar(read), read->core.n_cigar,
				read_left_cigar, &read_left_cigar_l);

		//Iterate haplotypes
		for(i = 0; i < haplo_list_l; i++)
		{
			//Get haplotype
			haplo = array_list_get(i, context->haplo_list);
			assert(haplo);

			//Get initial position to iterate
			read_disp_ref = SIZE_MAX;
			init_pos = haplo->ref_pos - read->core.l_qseq;

			//Initial position must be inside reference range
			if(init_pos < ref_pos_begin)
				init_pos = ref_pos_begin;

			//Initial position must not overlap last region
			if(init_pos < context->last_region.end_pos)
				init_pos = context->last_region.end_pos;

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
						#pragma omp critical
						{
							sprintf(log_msg, "Invalid haplotype %d\n", i + 1);
							LOG_ERROR(log_msg);
						}
						break;
					}
					else
					{
						//Cant displace anymore
						break;
					}
				}
				else
				{
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
						#pragma omp critical
						{
							sprintf(log_msg, "Read-Ref: %d, Read: %d\n", read_seq_ref_l, read->core.l_qseq);
							LOG_ERROR(log_msg);
						}
					}
					//assert(read_seq_ref_l == read->core.l_qseq);

					//Compare and miss with haplotype
					nucleotide_miss_qual_sum(read_seq, read_seq_ref, quals_seq, read_seq_ref_l, comp_seq, &misses, &misses_sum);

					//Better?
					if(v_scores[i+1] > misses_sum)
					{
						//printf("%s - %d - %d\n", bam1_qname(read), misses, misses_sum);
						//Update scores
						v_scores[i+1] = misses_sum;
						v_positions[i+1] = curr_pos;
					}

					//Perfect match?
					if(misses == 0)
					{
						break;
					}
				}
			} //Iterate positions
		} //Iterate haplotypes
	} //Misses != 0 if

	//Free
	free(read_seq);
	free(comp_seq);
	free(read_seq_ref);
	free(quals_seq);
	free(read_left_cigar);
	free(aux_cigar);

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
	#pragma omp parallel for default(none) \
	shared(stderr) \
	firstprivate(list, list_l, m_ldim, haplo, ref_seq, ref_length, m_scores, m_positions, alt_haplo_index, context, log_level) \
	private(err, read, read_indels, ref_score, h_score, read_seq, read_seq_ref, read_seq_ref_l, comp_aux, read_quals, read_disp_ref, misses_sum, misses, min_score, log_msg, bases, best_cigar, best_cigar_l, best_pos)
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

static inline ERROR_CODE
alig_aux_write_to_disk(array_list_t *write_buffer, bam_file_t *output_bam_f, uint8_t force)
{
	int i;
	size_t write_buffer_l;
	bam1_t *read;
	static size_t last_index = 0;

	//Timing
	double time_w;

	//Get list length
	write_buffer_l = array_list_size(write_buffer);
	if(write_buffer_l == 0)
	{
		return NO_ERROR;
	}

	//Is empty?
	if(force || write_buffer_l > ALIG_LIST_COUNT_THRESHOLD_TO_WRITE)
	{
#ifdef D_TIME_DEBUG
		time_w = omp_get_wtime();
#endif
		//Sort buffer
		array_list_qsort(write_buffer, compare_pos);

		//Write processed to disk
		for(i = 0; i < write_buffer_l; i++)
		{
			//Get read
			read = array_list_get(i, write_buffer);
			assert(read);

			//Write to disk
			bam_write1(output_bam_f->bam_fd, read);

			//Free read
			bam_destroy1(read);
		}

		//Clear buffer
		array_list_clear(write_buffer, NULL);
		last_index = 0;

#ifdef D_TIME_DEBUG
		time_w = omp_get_wtime() - time_w;
		time_add_time_slot(D_SLOT_ALIG_WRITE, TIME_GLOBAL_STATS, time_w/(double)write_buffer_l);
#endif
	} //If write_buffer_l > 0

	return NO_ERROR;
}

static inline ERROR_CODE
alig_aux_read_from_disk(circular_buffer_t *read_buffer, bam_file_t *input_bam_f)
{
	int i, cond;
	ssize_t bytes;
	size_t filled;
	size_t buffer_l;
	bam1_t *bam;

	//Timing
	double time_r;

#ifdef D_TIME_DEBUG
	time_r = omp_get_wtime();
#endif

	//More reads?
	omp_set_lock(&read_buffer->lock);
	cond = read_buffer->end_condition;
	omp_unset_lock(&read_buffer->lock);
	if(cond)
		return NO_ERROR;

	//Get buffer free length
	omp_set_lock(&read_buffer->lock);
	buffer_l = read_buffer->size - read_buffer->readed;
	omp_unset_lock(&read_buffer->lock);

	//Fill buffer
	bytes = 1;
	filled = 0;
	for(i = 0; i < buffer_l; i++)
	{
		bam = bam_init1();
		bytes = bam_read1(input_bam_f->bam_fd, bam);

		//Valid read?
		if(bytes > 0)
		{
			assert(bam);

			//Add read to buffer
			omp_set_lock(&read_buffer->lock);
			read_buffer->buffer[(read_buffer->head_idx + read_buffer->readed) % read_buffer->size] = bam;
			read_buffer->readed++;
			omp_unset_lock(&read_buffer->lock);

			filled++;
		}
		else
		{
			//Destroy empty last read
			bam_destroy1(bam);

			//No more reads
			omp_set_lock(&read_buffer->lock);
			read_buffer->end_condition = 1;
			omp_unset_lock(&read_buffer->lock);

			break;
		}
	}

#ifdef D_TIME_DEBUG
	time_r = omp_get_wtime() - time_r;
	if(filled > 0)
		time_add_time_slot(D_SLOT_ALIG_READ, TIME_GLOBAL_STATS, time_r/(double)filled);
	//printf("RD:%f\n", (double)(end_time_r - init_time_r) * 1000000);
#endif

	return NO_ERROR;
}
