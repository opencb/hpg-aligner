/*
 * alig.c
 *
 *  Created on: Dec 12, 2013
 *      Author: rmoreno
 */

#include "alig.h"

uint64_t cigar_changed = 0;

/**
 * CONTEXT
 */

ERROR_CODE
alig_init(alig_context_t *context, linked_list_t *in_list, genome_t *genome)
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

	//Create process list
	context->process_list = array_list_new(ALIG_LIST_COUNT_THRESHOLD_TO_WRITE + 100, 1.2f, COLLECTION_MODE_SYNCHRONIZED);
	if(context->process_list == NULL)
	{
		return ALIG_INIT_FAIL;
	}

	//Create haplotype list
	context->haplo_list = array_list_new(100, 1.2f, COLLECTION_MODE_SYNCHRONIZED);
	if(context->haplo_list == NULL)
	{
		return ALIG_INIT_FAIL;
	}

	//Set flags
	context->flags = ALIG_LEFT_ALIGN;

	return NO_ERROR;
}

ERROR_CODE
alig_destroy(alig_context_t *context)
{
	int i;
	bam1_t *read;
	aux_indel_t *haplo;

	assert(context);

	//Free proccess list
	if(context->process_list)
	{
		//Free alignments
		for(i = 0; i < array_list_size(context->process_list); i++)
		{
			//Get read
			read = array_list_get(i, context->process_list);

			//Free read
			bam_destroy1(read);
		}

		//Free list
		array_list_free(context->process_list, NULL);
	}

	//Free haplotype list
	if(context->haplo_list)
	{
		//Free haplotypes
		/*for(i = 0; i < array_list_size(context->haplo_list); i++)
		{
			//Get haplo
			haplo = array_list_get(i, context->haplo_list);

			//Free haplo
			free(haplo);
		}*/

		//Free list
		array_list_free(context->haplo_list, free);
	}

	//Set all to NULL
	memset(context, 0, sizeof(alig_context_t));

	return NO_ERROR;
}

ERROR_CODE
alig_validate(alig_context_t *context)
{
	if(context == NULL)
		return ALIG_INVALID_CONTEXT;

	if(context->in_list == NULL)
		return INVALID_INPUT_BAM;
	if(context->genome == NULL)
		return INVALID_GENOME;

	//Have process list?
	if(context->process_list == NULL)
		return ALIG_INVALID_PROCESS_LIST;

	//Have process list?
	if(context->haplo_list == NULL)
		return ALIG_INVALID_HAPLO_LIST;

	return NO_ERROR;
}

/**
 * REGION OPERATIONS
 */

ERROR_CODE
alig_region_next(alig_context_t *context)
{
	ERROR_CODE err;

	//List
	linked_list_iterator_t *it;

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

			if(indels == 1)
			{
				//Get read position
				read_pos = read->core.pos;

				//Check if chrom is different
				if(last_read_chrom != read->core.tid)
				{
					//Interval end!

					//Return
					return NO_ERROR;
				}

				//Cases
				if(context->region.valid == 0)
				{
					//NO INTERVAL CASE

					//Get interval for this alignment
					region_get_from_bam1(read, &interval_read_begin, &interval_read_end);

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
				else if (context->region.end_pos < read_pos)
				{
					//ALIGNMENT PAST INTERVAL

					//ERASE
					{
						//printf("INTERVAL %d:%d-%d %d\n", last_read_chrom + 1, context->region.start_pos + 1, context->region.end_pos + 1, array_list_size(context->process_list));
						/*int i;
						for(i = 0; i < array_list_size(context->process_list); i++)
						{
							//Get read
							bam1_t *rread = array_list_get(i, context->process_list);

							printf("%s\n", bam1_qname(rread));
						}*/
					}

					//Interval end!

					//Return
					return NO_ERROR;

				}
			} //End indels

			//Save last read chrom
			last_read_chrom = read->core.tid;

		}	// Filters if

		//Increment progress
		context->read_count++;

		//Add read to process list
		array_list_insert(read, context->process_list);

		//Read next alignment
		read = (bam1_t *)linked_list_iterator_next(it);
	}

	linked_list_iterator_free(it);

	if(context->region.valid)
	{
		//There is an open interval
		//ERASE
		{
			//printf("LAST INTERVAL %d:%d-%d %d\n", last_read_chrom + 1, context->region.start_pos, context->region.end_pos, array_list_size(context->process_list));
		}
	}

	return NO_ERROR;
}

ERROR_CODE
alig_region_haplotype_process(alig_context_t *context)
{
	ERROR_CODE err;
	int i;

	//Read
	bam1_t* read;
	uint32_t *read_cigar;
	size_t read_indels;
	char *read_seq;
	size_t read_seq_l;
	int64_t read_disp_ref;

	//Reference
	char *ref_seq;
	uint32_t flag;
	int64_t ref_init_pos;
	int64_t ref_end_pos;
	size_t uinit_pos;
	size_t uend_pos;

	//Haplotype
	aux_indel_t *haplo;
	aux_indel_t *aux_haplo;
	size_t haplo_l;

	assert(context);

	//Validate context
	err = alig_validate(context);
	if(err)
	{
		return err;
	}

	//Get all haplotypes
	for(i = 0; i < array_list_size(context->process_list); i++)
	{
		//Get read
		read = array_list_get(i, context->process_list);
		assert(read);

		//FILTERS: MAP QUALITY = 0, BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP, DIFFERENT MATE CHROM, 1 INDEL, CIGAR PRESENT
		if( read->core.qual != 0
			&& !(read->core.flag & BAM_DEF_MASK)
			&& read->core.mtid == read->core.tid
			&& read->core.n_cigar != 0
			)
		{
			//Get read cigar
			read_cigar = bam1_cigar(read);
			if(read_cigar)
			{
				//Valid cigar

				//This read have indels?
				cigar32_count_indels(read_cigar, read->core.n_cigar, &read_indels);

				if(read_indels == 1)	//Only one indel support
				{
					char cigar_str[100];
					char cigar_str2[100];
					char cigar_str3[100];
					char ref_strx[200];
					size_t ref_strx_l;
					char ref_strx2[200];
					size_t ref_strx2_l;
					uint32_t clip_cigar[50];
					size_t clip_cigar_l;
					uint32_t unclip_cigar[50];
					size_t unclip_cigar_l;

					if(context->flags & ALIG_LEFT_ALIGN)
					{
						//Left align cigar first

						//Convert read sequence to string
						read_seq_l = read->core.l_qseq;
						read_seq = (char *) malloc((read_seq_l + 1) * sizeof(char));
						assert(read);
						new_sequence_from_bam_ref(read, read_seq, read_seq_l + 1);

						//Reference range
						ref_init_pos = read->core.pos - ALIG_REFERENCE_ADDITIONAL_OFFSET;
						ref_end_pos = read->core.pos + read_seq_l + ALIG_REFERENCE_ADDITIONAL_OFFSET;

						//Check limits
						if(ref_init_pos < 0)
						{
							ref_init_pos = 0;
						}
						if(ref_end_pos < 0)
						{
							ref_end_pos = INT64_MAX;
						}
						uinit_pos = ref_init_pos + 1;
						uend_pos = ref_end_pos + 1;

						//Get reference
						ref_seq = (char *)malloc((read_seq_l + (ALIG_REFERENCE_ADDITIONAL_OFFSET * 2) + 2) * sizeof(char));
						flag = (uint32_t) read->core.flag;
						genome_read_sequence_by_chr_index(ref_seq, (flag & BAM_FREVERSE) ? 1 : 0, (unsigned int)read->core.tid, &uinit_pos, &uend_pos, context->genome);

						//Get relative displacement from reference
						read_disp_ref = read->core.pos - ref_init_pos;
						if(read_disp_ref < 0)
							read_disp_ref = 0;

						//Leftmost cigar
						cigar32_leftmost(ref_seq + read_disp_ref, (ref_end_pos - ref_init_pos) - read_disp_ref, read_seq, read->core.l_qseq, read_cigar, read->core.n_cigar, clip_cigar, &clip_cigar_l);

						//ERASE
						cigar32_unclip(clip_cigar, clip_cigar_l, unclip_cigar, &unclip_cigar_l);
						cigar32_to_string(read_cigar, read->core.n_cigar, cigar_str);
						cigar32_to_string(unclip_cigar, unclip_cigar_l, cigar_str2);
						cigar32_to_string(clip_cigar, clip_cigar_l, cigar_str3);
						if(memcmp(cigar_str, cigar_str2, strlen(cigar_str)) && memcmp(cigar_str, cigar_str3, strlen(cigar_str)))
						{
							/*if(!cond)
							{
								cond = 1;
								printf("================================================\n");
								printf("************************************************\n");
								printf("INTERVAL: %d:%d-%d\n", read->core.tid, ref_pos_begin, ref_pos_end);
							}

							cigar32_create_ref(read_cigar, read->core.n_cigar, ref_seq + read_disp_ref, ref_length - read_disp_ref, read_seq, read->core.l_qseq, ref_strx, &ref_strx_l);
							cigar32_create_ref(clip_cigar, clip_cigar_l, ref_seq + read_disp_ref, ref_length - read_disp_ref, read_seq, read->core.l_qseq, ref_strx2, &ref_strx2_l);
							printf("************************************************\n");
							printf("CIGAR LEFTALIGNED: %s - %s - %s === %d:%d\n", cigar_str, cigar_str2, cigar_str3, read->core.tid, read->core.pos);
							printf("Read pos: %d === Disp: %d\n", read->core.pos, read_disp_ref);
							printf("Ref  : %s\n", ref_seq);
							printf("Read : ");
							for(int j = 0; j < read_disp_ref && j < 300; j++) printf(" ");
							printf("%s\n", read_seq);
							printf("Ref *: ");
							for(int j = 0; j < read_disp_ref && j < 300; j++) printf(" ");
							printf("%s\n", ref_strx);
							printf("Ref A: ");
							for(int j = 0; j < read_disp_ref && j < 300; j++) printf(" ");
							printf("%s\n", ref_strx2);*/

							print_log_message_with_format(LOG_INFO_LEVEL, "INFO", "alig.c", 378, "alig_bam_list",
									"CIGAR LEFTALIGNED: %s - %s - %s === %d:%d\n", cigar_str, cigar_str2, cigar_str3, read->core.tid + 1, read->core.pos + 1);
						}

						//Free
						free(ref_seq);
						free(read_seq);
					}

					//Allocate haplotype
					haplo = (aux_indel_t *)malloc(sizeof(aux_indel_t) * read_indels);	//For now is only 1

					//Fill haplotype
					cigar32_get_indels(read->core.pos, clip_cigar, clip_cigar_l, haplo);

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

				} //Indel count if
			} //Valid cigar end if
		}//Filters if
	}

	return NO_ERROR;
}

ERROR_CODE
alig_region_indel_realignment(alig_context_t *context)
{
	ERROR_CODE err;
	int i;

	//List
	array_list_t *list;
	size_t list_l;

	//Read
	bam1_t* read;
	int32_t read_pos;

	//Reference
	char *ref_seq = NULL;
	uint32_t flag;
	size_t ref_pos_begin = SIZE_MAX;
	size_t ref_pos_end = SIZE_MAX;
	size_t aux_pos_begin = SIZE_MAX;
	size_t aux_pos_end = SIZE_MAX;
	size_t ref_length = SIZE_MAX;

	//Haplotype
	aux_indel_t *haplo = NULL;
	aux_indel_t *aux_haplo = NULL;

	//Realign list
	array_list_t *alig_list;

	assert(context);

	//Validate context
	err = alig_validate(context);
	if(err)
	{
		return err;
	}

	//Get list ptr
	list = context->process_list;
	list_l = array_list_size(list);

	//Realign only if a region is present
	if(context->region.valid)
	{
		//printf("Alig..\n");
		//alig_bam_list(list, context->genome);

		//Get reference parameters
		/*{
			//Get first read position
			ref_pos_begin = ((bam1_t *)array_list_get(0, context->process_list))->core.pos;

			//Get last read position + length
			read = (bam1_t *) array_list_get(array_list_size(context->process_list) - 1, context->process_list);
			ref_pos_end = (read->core.pos + read->core.l_qseq) + ALIG_REFERENCE_ADDITIONAL_OFFSET;

			//Auxiliar
			aux_pos_begin = ref_pos_begin + 1;
			aux_pos_end = ref_pos_end + 1;

			//Get reference
			ref_seq = (char *) malloc(((ref_pos_end - ref_pos_begin) + 2) * sizeof(char));
			flag = (uint32_t) read->core.flag;
			genome_read_sequence_by_chr_index(ref_seq, (flag & BAM_FREVERSE) ? 1 : 0, (unsigned int)read->core.tid, &aux_pos_begin, &aux_pos_end, context->genome);

			ref_length = (ref_pos_end - ref_pos_begin);
		}*/

		//Allocate list
		alig_list = array_list_new(list_l, 1.2f, 0);

		//Iterate list to realign
		for(i = 0; i < list_l; i++)
		{
			//Get read
			read = array_list_get(i, list);
			assert(read);

			//FILTERS: MAP QUALITY = 0, BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP, DIFFERENT MATE CHROM, 1 INDEL, CIGAR PRESENT
			if( read->core.qual != 0
				&& !(read->core.flag & BAM_DEF_MASK)
				&& read->core.mtid == read->core.tid
				&& read->core.n_cigar != 0
				)
			{

				//Is in region?
				if(region_bam_overlap(read, &context->region))
				{
					//In region
					array_list_insert(read, alig_list);
				}

			}
		}

		//Align region reads
		//printf("Realigning %d reads\n", array_list_size(alig_list));
		alig_bam_list_realign(alig_list, context->haplo_list, context->genome);

		//Destroy list
		array_list_free(alig_list, NULL);
	}

	return NO_ERROR;
}


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

	//Clear process list
	array_list_clear(context->process_list, NULL);

	//Clear haplotype list
	array_list_clear(context->haplo_list, free);

	return NO_ERROR;
}


ERROR_CODE
alig_bam_file2(char *bam_path, char *ref_name, char *ref_path)
{
	ERROR_CODE err;
	unsigned char outbam[30] = "output.bam";
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

		if(bytes > 0)
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
	alig_init(&context, in_list, ref);

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
	list_l = array_list_size(context.process_list);
	time_add_time_slot(D_SLOT_NEXT, TIME_GLOBAL_STATS, (double)(end_time - init_time)/(double)list_l);
#endif

	//Iterate
	while(array_list_size(context.process_list) > 0)
	{
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
		list_l = array_list_size(context.process_list);
		time_add_time_slot(D_SLOT_HAPLO_GET, TIME_GLOBAL_STATS, (double)(end_time - init_time)/(double)list_l);
#endif

		//ERASE
		/*{
			//Printf haplotypes
			printf("Haplotypes:\n");
			int h;
			aux_indel_t *haplo;
			for(h = 0; h < array_list_size(context.haplo_list); h++)
			{
				haplo = array_list_get(h, context.haplo_list);
				assert(haplo);

				//Printf
				char cigar_str[50];
				cigar32_to_string(&haplo->indel, 1, cigar_str);
				printf("H%d: %s -- %d:%d\n", h, cigar_str, context.region.chrom + 1, haplo->ref_pos + 1);
			}
			printf("----\n");

		}*/

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
		list_l = array_list_size(context.process_list);
		end_time = (double)(end_time - init_time)/(double)list_l;
		list_l = array_list_size(context.haplo_list);
		time_add_time_slot(D_SLOT_REALIG_PER_HAPLO, TIME_GLOBAL_STATS, end_time/(double)list_l);
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
		list_l = array_list_size(context.process_list);
		for(i = 0; i < list_l; i++)
		{
			//Get read
			read = array_list_get(i, context.process_list);
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
		for(i = list_l; i < ALIG_LIST_IN_SIZE && bytes > 0; i++)
		{
			read = bam_init1();
			bytes = bam_read1(bam_f->bam_fd, read);

			if(bytes > 0)
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
		list_l = array_list_size(context.process_list);
		if(list_l > 0)
			time_add_time_slot(D_SLOT_NEXT, TIME_GLOBAL_STATS, (double)(end_time - init_time)/(double)list_l);
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

	return NO_ERROR;
}

//					H0			H1			H2			H3			H4
//				sco	pos
//	ERRXXXX		23	97		0	16		321	34
//	ERRYYYY
//	ERRZZZZ
//

ERROR_CODE
alig_bam_list_realign(array_list_t *bam_list, array_list_t *haplotype_list, genome_t* ref)
{
	int i, j;
	size_t bam_list_l;
	size_t h_list_l;

	//Reference
	char *ref_seq = NULL;
	uint32_t flag;
	size_t ref_pos_begin = SIZE_MAX;
	size_t ref_pos_end = SIZE_MAX;
	size_t aux_pos_begin = SIZE_MAX;
	size_t aux_pos_end = SIZE_MAX;
	size_t ref_length = SIZE_MAX;

	//Read
	bam1_t *read;
	char *read_seq = NULL;
	char *read_seq_ref = NULL;
	char *read_quals = NULL;
	size_t read_seq_ref_l;
	size_t read_disp_ref = SIZE_MAX;
	size_t indels;

	//Comparation
	char *comp_aux = NULL;
	uint32_t comp_cigar[MAX_CIGAR_LENGTH];
	uint32_t comp_cigar2[MAX_CIGAR_LENGTH];
	size_t comp_cigar_l;
	size_t comp_cigar_l2;
	size_t bases;
	uint32_t score = 0;
	uint32_t score2 = 0;
	uint32_t best_score;
	size_t best_haplo_index;
	uint32_t miss;
	uint32_t miss2;

	//Matrix
	uint32_t *m_score;
	size_t *m_pos;
	uint32_t *v_hscore;
	size_t m_total;
	size_t m_ldim;

	//Haplotype
	aux_indel_t *haplo;
	aux_indel_t best_haplo;

	//Alig
	uint32_t best_cigar[MAX_CIGAR_LENGTH];
	size_t best_cigar_l;
	uint32_t best_cigar_score;
	size_t best_pos;
	size_t curr_pos;

	assert(bam_list);
	assert(haplotype_list);
	assert(ref);

	//Empty lists?
	if(array_list_size(bam_list) == 0 || array_list_size(haplotype_list) == 0)
	{
		return NO_ERROR;
	}

	//Get reference
	{
		//Get first read position
		ref_pos_begin = ((bam1_t *)array_list_get(0, bam_list))->core.pos;

		//Get last read position + length
		read = (bam1_t *) array_list_get(array_list_size(bam_list) - 1, bam_list);
		ref_pos_end = (read->core.pos + read->core.l_qseq) + ALIG_REFERENCE_ADDITIONAL_OFFSET;

		//Auxiliar
		aux_pos_begin = ref_pos_begin + 1;
		aux_pos_end = ref_pos_end + 1;

		//Get reference
		ref_seq = (char *) malloc(((ref_pos_end - ref_pos_begin) + 2) * sizeof(char));
		read_seq_ref = (char *) malloc(((ref_pos_end - ref_pos_begin) + 2) * sizeof(char));
		comp_aux = (char *) malloc(((ref_pos_end - ref_pos_begin) + 2) * sizeof(char));
		flag = (uint32_t) read->core.flag;
		genome_read_sequence_by_chr_index(ref_seq, (flag & BAM_FREVERSE) ? 1 : 0, (unsigned int)read->core.tid, &aux_pos_begin, &aux_pos_end, ref);
		assert(ref);

		ref_length = (ref_pos_end - ref_pos_begin);
	}

	//Lengths
	bam_list_l = array_list_size(bam_list);
	h_list_l = array_list_size(haplotype_list);

	//Allocate miss arrays
	v_hscore = (uint32_t *)malloc((sizeof(size_t) * h_list_l));	//H0 is reference haplotype
	memset(v_hscore, 0, sizeof(size_t) * h_list_l);

	//Allocate score matrix
	m_ldim = h_list_l + 1;
	m_total = bam_list_l * m_ldim;
	m_pos = (size_t *)malloc(m_total * sizeof(size_t));
	m_score = (uint32_t *)malloc(m_total * sizeof(uint32_t));
	for(i = 0; i < m_total; i++)
	{
		*((size_t*)m_pos + i) = SIZE_MAX;
		*((uint32_t*)m_score + i) = UINT32_MAX;
	}

	//Iterate reads
	for(i = 0; i < bam_list_l; i++)
	{
		//Get read
		read = array_list_get(i, bam_list);

		//Only if one or less indels
		cigar32_count_indels(bam1_cigar(read), read->core.n_cigar, &indels);

		if(indels > 1)
		{
			continue;
		}

		//Convert read sequence to string
		read_seq = (char *) malloc((read->core.l_qseq + 1) * sizeof(char));
		new_sequence_from_bam_ref(read, read_seq, read->core.l_qseq + 1);

		//Get qualities
		read_quals = (char *) malloc((read->core.l_qseq + 1) * sizeof(char));
		new_quality_from_bam_ref(read, 0, read_quals, read->core.l_qseq + 1);

		//Get relative displacement from reference
		//read_disp_ref = read->core.pos - ref_pos_begin;

		//ERASE Print things
		/*{
			char erase_str[200];
			size_t disp;
			cigar32_to_string(bam1_cigar(read), read->core.n_cigar, erase_str);
			printf("---------------------\n");
			printf("Test %s - %d:%d - %s\n", bam1_qname(read), read->core.tid + 1, read->core.pos + 1, erase_str);
		}*/

		//Reference score
		{
			//Unclip cigar
			cigar32_unclip(bam1_cigar(read), read->core.n_cigar, comp_cigar, &comp_cigar_l);

			//Count unclipped bases
			cigar32_count_nucleotides_not_clip(comp_cigar, comp_cigar_l, &bases);

			//Create new cigar with '$bases'M
			comp_cigar[0] = bases << BAM_CIGAR_SHIFT;	//ex: 108M
			comp_cigar_l = 1;

			//Reclip cigar
			cigar32_reclip(bam1_cigar(read), read->core.n_cigar, comp_cigar, comp_cigar_l, comp_cigar, &comp_cigar_l);

			//Iterate positions
			for(curr_pos = read->core.pos; curr_pos < ref_pos_end - read->core.l_qseq; curr_pos++)
			{
				//Get relative displacement from reference
				read_disp_ref = curr_pos - ref_pos_begin;

				//Get reference transform for this read
				cigar32_create_ref(comp_cigar, comp_cigar_l, ref_seq + read_disp_ref, ref_length - read_disp_ref, read_seq, read->core.l_qseq, read_seq_ref, &read_seq_ref_l);

				//Test with reference
				nucleotide_compare(read_seq, read_seq_ref, read_seq_ref_l, comp_aux, &miss);

				//Calculate miss score
				score = 0;
				if(miss != 0)
				{
					int z;
					for(z = 0; z < read_seq_ref_l; z++)
					{
						if(comp_aux[z] == 0)
							score += read_quals[z];
					}
					//Better?
					if(m_score[i * m_ldim] > score)
					{
						m_score[i * m_ldim] = score;
						m_pos[i * m_ldim] = curr_pos;
					}
				}
				else
				{
					//Found perfect match
					m_score[i * m_ldim] = 0;
					m_pos[i * m_ldim] = curr_pos;
					break;
				}
			}

			//ERASE Print things
			/*{
				char cigar_str[50];
				//char erase_str[200];
				size_t disp;
				//cigar32_count_clip_displacement(comp_cigar, comp_cigar_l, &disp);
				cigar32_to_string(comp_cigar, comp_cigar_l, cigar_str);
				printf("H0 === MISS SCORE: %d:%d - Using CIGAR: %s\n", m_score[i * m_ldim], m_pos[i * m_ldim] + 1, cigar_str);
				//read_disp_ref = m_pos[i * m_ldim] - ref_pos_begin;
				//cigar32_create_ref(comp_cigar, comp_cigar_l, ref_seq + read_disp_ref, ref_length - read_disp_ref, read_seq, read->core.l_qseq, read_seq_ref, &read_seq_ref_l);
				//printf("%s\n%s\n", read_seq, read_seq_ref);

				//memcpy(erase_str, ref_seq + read_disp_ref, sizeof(char) * read->core.l_qseq);
				//erase_str[read->core.l_qseq] = '\0';
				//printf("Ref : ");
				//for(int z = 0; z < disp; z++) printf(" ");
				//printf("%s\n", erase_str);
				//memcpy(erase_str, read_seq, sizeof(char) * read->core.l_qseq);
				//erase_str[read->core.l_qseq] = '\0';
				//printf("Read: %s\n", erase_str);
				//printf("Ref*: %s\n", read_seq_ref);
			}*/
		}

		//Dont iterate haplotypes if perfect reference match
		if(m_score[i * m_ldim] != 0)
		{
			//Iterate haplotypes
			for(j = 0; j < h_list_l; j++)
			{
				//Get haplotype
				haplo = array_list_get(j, haplotype_list);
				assert(haplo);

				//IMPLEMENTAR FORWARD PRIMERO Y SI NO BACKWARD
				//TODO

				//Iterate positions
				read_disp_ref = SIZE_MAX;
				for(curr_pos = haplo->ref_pos - read->core.l_qseq; curr_pos < ref_pos_end - read->core.l_qseq; curr_pos++)
				{
					//Create cigar for this haplotype
					if(cigar32_from_haplo(bam1_cigar(read), read->core.n_cigar, haplo, curr_pos, comp_cigar, &comp_cigar_l))
					{
						//Cant displace anymore
						break;
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

					//ERASE Print things
					/*{
						char cigar_str[50];
						cigar32_to_string(comp_cigar, comp_cigar_l, cigar_str);
						printf("Testing %s\n", cigar_str);
					}*/

					//Get haplotype reference transform
					cigar32_create_ref(comp_cigar, comp_cigar_l, ref_seq + read_disp_ref, ref_length - read_disp_ref, read_seq, read->core.l_qseq, read_seq_ref, &read_seq_ref_l);

					//Compare and miss with haplotype
					nucleotide_compare(read_seq, read_seq_ref, read_seq_ref_l, comp_aux, &miss);

					//Calculate H miss score
					score = 0;
					if(miss != 0)
					{
						int z;
						for(z = 0; z < read_seq_ref_l; z++)
						{
							if(comp_aux[z] == 0)
								score += read_quals[z];
						}
						//Better?
						if(m_score[(i * m_ldim) + (j+1)] > score)
						{
							m_score[(i * m_ldim) + (j+1)] = score;
							m_pos[(i * m_ldim) + (j+1)] = curr_pos;
						}
					}
					else
					{
						//Found perfect match
						m_score[(i * m_ldim) + (j+1)] = score;
						m_pos[(i * m_ldim) + (j+1)] = curr_pos;
						break;
					}

					//ERASE
					/*{
						char cigar_str[50];
						//char erase_str[200];
						size_t disp;
						cigar32_from_haplo(bam1_cigar(read), read->core.n_cigar, haplo, curr_pos, comp_cigar, &comp_cigar_l);
						cigar32_count_clip_displacement(comp_cigar, comp_cigar_l, &disp);
						cigar32_to_string(comp_cigar, comp_cigar_l, cigar_str);
						printf("H%d === MISS SCORE: %d:%d - Using CIGAR: %s\n", j+1, score, curr_pos + 1, cigar_str);
						//memcpy(erase_str, ref_seq + read_disp_ref, sizeof(char) * read->core.l_qseq);
						//erase_str[read->core.l_qseq] = '\0';
						//printf("Ref : ");
						//for(int z = 0; z < disp; z++) printf(" ");
						//printf("%s\n", erase_str);
						//memcpy(erase_str, read_seq, sizeof(char) * read->core.l_qseq);
						//erase_str[read->core.l_qseq] = '\0';
						//printf("Read: %s\n", erase_str);
						//printf("Ref*: %s\n", read_seq_ref);
					}*/
				}

				//If reference maps better, not count haplotype score
				if(m_score[(i * m_ldim) + (j+1)] > m_score[(i * m_ldim)])
				{
					m_score[(i * m_ldim) + (j+1)] = UINT32_MAX;
					m_pos[(i * m_ldim) + (j+1)] = SIZE_MAX;
				}

				//ERASE Print things
				/*if(read_disp_ref != SIZE_MAX && m_pos[(i * m_ldim) + (j+1)] != SIZE_MAX)
				{
					char cigar_str[50];
					//char erase_str[200];
					size_t disp;
					cigar32_from_haplo(bam1_cigar(read), read->core.n_cigar, haplo, m_pos[(i * m_ldim) + (j+1)], comp_cigar, &comp_cigar_l);
					cigar32_count_clip_displacement(comp_cigar, comp_cigar_l, &disp);
					cigar32_to_string(comp_cigar, comp_cigar_l, cigar_str);
					printf("H%d === MISS SCORE: %d:%d - Using CIGAR: %s\n", j+1, m_score[(i * m_ldim) + (j+1)], m_pos[(i * m_ldim) + (j+1)] + 1, cigar_str);
					//memcpy(erase_str, ref_seq + read_disp_ref, sizeof(char) * read->core.l_qseq);
					//erase_str[read->core.l_qseq] = '\0';
					//printf("Ref : ");
					//for(int z = 0; z < disp; z++) printf(" ");
					//printf("%s\n", erase_str);
					//memcpy(erase_str, read_seq, sizeof(char) * read->core.l_qseq);
					//erase_str[read->core.l_qseq] = '\0';
					//printf("Read: %s\n", erase_str);
					//printf("Ref*: %s\n", read_seq_ref);
				}*/
			}
		}

		//Free
		free(read_seq);
		free(read_quals);
	}

	//ERASE
	/*{
		//Print table
		printf("---------------------\n");
		printf("*****Best scores*****\n");

		for(j = 0; j < m_ldim; j++)
		{
			printf("%8d ",j);
		}
		printf("\n");

		for(i = 0; i < bam_list_l; i++)
		{
			//Get read
			read = array_list_get(i, bam_list);

			for(j = 0; j < m_ldim; j++)
			{
				printf("%8d ", m_score[(i * m_ldim) + j]);
			}
			printf("%s\n", bam1_qname(read));
		}

		//Print table
		printf("*****Best positions*****\n");

		for(j = 0; j < m_ldim; j++)
		{
			printf("%8d ",j);
		}
		printf("\n");

		for(i = 0; i < bam_list_l; i++)
		{
			//Get read
			read = array_list_get(i, bam_list);

			for(j = 0; j < m_ldim; j++)
			{
				printf("%8d ", m_pos[(i * m_ldim) + j]);
			}
			printf("%s\n", bam1_qname(read));
		}
	}*/

	//Find best haplotype
	//printf("---------------------\n");
	best_score = UINT32_MAX;
	int count;
	for(j = 0; j < h_list_l; j++)
	{
		//Sum scores
		count = 0;
		for(i = 0; i < bam_list_l; i++)
		{
			if(m_score[(i * m_ldim) + (j+1)] != UINT32_MAX)	//If valid score
			{
				count++;
				v_hscore[j] += m_score[(i * m_ldim) + (j+1)];
			}
		}

		if(count != 0)
		{
			v_hscore[j] = v_hscore[j] / count;
		}

		if(v_hscore[j] < best_score)
		{
			best_score = v_hscore[j];
			best_haplo_index = j;
		}

		//ERASE
		/*{
			//Print haplotype score
			printf("H%d === total score = %d\n", j + 1, v_hscore[j]);
		}*/
	}

	//Get best haplotype
	haplo = array_list_get(best_haplo_index, haplotype_list);

	//Iterate bams to change its CIGARS
	for(i = 0; i < bam_list_l; i++)
	{
		//Get read
		read = array_list_get(i, bam_list);
		assert(read);

		//Only if one or less indels
		cigar32_count_indels(bam1_cigar(read), read->core.n_cigar, &indels);
		if(indels > 1)
		{
			continue;
		}

		//Dont change if original its a lot better
		{
			//Convert read sequence to string
			read_seq = (char *) malloc((read->core.l_qseq + 1) * sizeof(char));
			new_sequence_from_bam_ref(read, read_seq, read->core.l_qseq + 1);

			//Get qualities
			read_quals = (char *) malloc((read->core.l_qseq + 1) * sizeof(char));
			new_quality_from_bam_ref(read, 0, read_quals, read->core.l_qseq + 1);

			read_disp_ref = read->core.pos - ref_pos_begin;

			//Get haplotype reference transform
			cigar32_create_ref(bam1_cigar(read), read->core.n_cigar, ref_seq + read_disp_ref, ref_length - read_disp_ref, read_seq, read->core.l_qseq, read_seq_ref, &read_seq_ref_l);

			//Compare and miss with haplotype
			nucleotide_compare(read_seq, read_seq_ref, read_seq_ref_l, comp_aux, &miss);

			//Calculate H miss score
			score = 0;
			if(miss != 0)
			{
				int z;
				for(z = 0; z < read_seq_ref_l; z++)
				{
					if(comp_aux[z] == 0)
						score += read_quals[z];
				}
			}

			free(read_seq);
			free(read_quals);

			int mn = min(m_score[(i * m_ldim)], m_score[(i * m_ldim) + (best_haplo_index + 1)]);
			if(mn != 0)
			{
				float ratio = (float)score/(float)mn;
				if(ratio < 0.5f)
				{
					//Dont change reading
					continue;
				}
			}
		}

		//Is reference matching better?
		score = m_score[(i * m_ldim)];
		score2 = m_score[(i * m_ldim) + (best_haplo_index + 1)];
		if(score <= score2)
		{
			//Unclip cigar
			cigar32_unclip(bam1_cigar(read), read->core.n_cigar, best_cigar, &best_cigar_l);

			//Count unclipped bases
			cigar32_count_nucleotides_not_clip(best_cigar, best_cigar_l, &bases);

			//Create new cigar with '$bases'M
			best_cigar[0] = bases << BAM_CIGAR_SHIFT;	//ex: 108M
			best_cigar_l = 1;

			//Reclip cigar
			cigar32_reclip(bam1_cigar(read), read->core.n_cigar, best_cigar, best_cigar_l, best_cigar, &best_cigar_l);

			//Set best position
			best_pos = m_pos[(i * m_ldim)];
		}
		else //Haplotype is better
		{
			//Set best position
			best_pos = m_pos[(i * m_ldim) + (best_haplo_index + 1)];

			//Create cigar for this haplotype
			cigar32_from_haplo(bam1_cigar(read), read->core.n_cigar, haplo, best_pos, best_cigar, &best_cigar_l);

			//If haplotype position is minor then adjust
			if(best_pos > haplo->ref_pos)
			{
				best_pos = haplo->ref_pos;	//ex: 12M3D1M => 1D13M
			}
		}

		//ERASE
		/*{
			char cigar_str[50];
			char cigar_str2[50];
			cigar32_to_string(bam1_cigar(read), read->core.n_cigar, cigar_str);
			cigar32_to_string(best_cigar, best_cigar_l, cigar_str2);
			//Print new cigar
			printf("%s Cigar: %s => %s\n", bam1_qname(read), cigar_str, cigar_str2);
		}*/

		//Change CIGAR
		cigar32_replace(read, best_cigar, best_cigar_l);
		read->core.pos = best_pos;
	}

	//Free
	{
		free(ref_seq);
		free(read_seq_ref);
		free(v_hscore);
		free(m_score);
		free(m_pos);
		free(comp_aux);
		//free(comp_cigar);
	}

	return NO_ERROR;
}

