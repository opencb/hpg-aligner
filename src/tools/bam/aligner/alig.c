/*
 * alig.c
 *
 *  Created on: Dec 12, 2013
 *      Author: rmoreno
 */

#include "alig.h"

uint64_t cigar_changed = 0;

/**
 * BAM REALIGN
 */

ERROR_CODE
alig_bam_file(char *bam_path, char *ref_name, char *ref_path)
{
	bam_file_t *bam_f = NULL;
	bam_file_t *out_bam_f = NULL;
	bam_batch_t *batch = NULL;
	genome_t* ref = NULL;
	ERROR_CODE err;
	char outbam[30] = "output.bam";

	//BAM alignment
	bam1_t* bam_read = NULL;
	int bytes;
	uint64_t read_pos = 0;
	int32_t last_read_chrom = -1;
	alig_status status;

	//Indel interval
	size_t interval_begin = SIZE_MAX;
	size_t interval_end = SIZE_MAX;
	size_t interval_read_begin = SIZE_MAX;
	size_t interval_read_end = SIZE_MAX;

	//Lists
	array_list_t *write_list;
	array_list_t *process_list;

	size_t indels;
	uint64_t count = 0;
	int i;

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
		bam_fwrite_header(out_bam_f->bam_header_p, out_bam_f);
		out_bam_f->bam_header_p = NULL;
		printf("New BAM initialized!...\n");
	}

	//Create lists
	write_list = array_list_new(ALIG_LIST_COUNT_THRESHOLD_TO_WRITE, 1.2f, COLLECTION_MODE_SYNCHRONIZED);
	process_list = array_list_new(ALIG_LIST_COUNT_THRESHOLD_TO_WRITE, 1.2f, COLLECTION_MODE_SYNCHRONIZED);

	//Read first alignment
	bam_read = bam_init1();
	bytes = bam_read1(bam_f->bam_fd, bam_read);
	last_read_chrom = bam_read->core.tid;

	count = 1;
	status = NO_INTERVAL;
	while(bytes > 0 /*&& count < 100000*/)
	{
		//FILTERS: MAP QUALITY = 0, BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP, DIFFERENT MATE CHROM
		if(bam_read
				&& bam_read->core.qual != 0
				&& !(bam_read->core.flag & BAM_DEF_MASK)
				&& bam_read->core.mtid == bam_read->core.tid
				)
		{
			//Get read position
			read_pos = bam_read->core.pos;

			//Check if chrom is different
			if(last_read_chrom != bam_read->core.tid)
			{
				//Realign alignments in list
				alig_bam_list(process_list, ref);
				array_list_clear(process_list, NULL);

				//Reset interval
				interval_begin = SIZE_MAX;
				interval_end = SIZE_MAX;

				//Continue from this read
				last_read_chrom = bam_read->core.tid;
			}

			//Cases
			if(interval_begin == SIZE_MAX)
			{
				//NO INTERVAL CASE

				//Filter alignments with one indel
				cigar32_count_indels(bam1_cigar(bam_read), bam_read->core.n_cigar, &indels);
				if(indels == 1)
				{
					//Get interval for this alignment
					region_get_from_bam1(bam_read, &interval_read_begin, &interval_read_end);

					//This alignment have an interval?
					if(interval_read_begin != SIZE_MAX)
					{
						//Interval found
						//Set interval values
						interval_begin = interval_read_begin;
						interval_end = interval_read_end;
						assert(interval_begin <= interval_end);

						//Insert alignment in process list
						array_list_insert(bam_read, process_list);
					}
				}
			}
			else if (interval_end > read_pos)
			{
				//ALIGNMENT INSIDE INTERVAL CASE

				//Filter alignments with one indel
				cigar32_count_indels(bam1_cigar(bam_read), bam_read->core.n_cigar, &indels);
				if(indels == 1)
				{
					//Get interval for this alignment
					region_get_from_bam1(bam_read, &interval_read_begin, &interval_read_end);
				}
				else
				{
					interval_read_begin = SIZE_MAX;
				}

				//This alignment have an interval?
				if(interval_read_begin != SIZE_MAX)
				{
					//Interval found
					//Set interval begin
					if(interval_read_begin < interval_begin)
						interval_begin = interval_read_begin;

					//Set interval end
					if(interval_read_end > interval_end)
						interval_end = interval_read_end;

					//Check
					assert(interval_begin <= interval_end);
				}

				//Insert alignment in process list
				array_list_insert(bam_read, process_list);
			}
			else if (interval_end < read_pos)
			{
				//ALIGNMENT PAST INTERVAL

				//ERASE
				/*{
					printf("INTERVAL %d:%d-%d %d\n", last_read_chrom + 1, interval_begin, interval_end, array_list_size(process_list));
					for(i = 0; i < array_list_size(process_list); i++)
					{
						//Get read
						bam1_t *read = array_list_get(i, process_list);

						printf("%s\n", bam1_qname(read));
					}
				}*/

				//Realign alignments in process list
				alig_bam_list(process_list, ref);
				array_list_clear(process_list, NULL);

				//Reset interval
				interval_begin = SIZE_MAX;
				interval_end = SIZE_MAX;

				//Continue from this read
				last_read_chrom = bam_read->core.tid;
				continue;
			}

			//Save last read chrom
			last_read_chrom = bam_read->core.tid;

			/*if(interval_begin != SIZE_MAX)
			{
				//Insert alignment in process list
				array_list_insert(bam_read, process_list);
			}*/

		}	// Filters if

		//Put alignment to write list
		array_list_insert(bam_read, write_list);

		//Read next alignment
		bam_read = bam_init1();
		bytes = bam_read1(bam_f->bam_fd, bam_read);

		//If write list is big enough and process list is empty: write to disk
		if(array_list_size(write_list) > ALIG_LIST_COUNT_THRESHOLD_TO_WRITE
				&& array_list_size(process_list) == 0)
		{
			alig_bam_list_to_disk(write_list, out_bam_f);
		}

		//Show total progress
		count++;
		if(count % 100000 == 0)
		{
			printf("Total alignments readed: %d\r", count);
			fflush(stdout);
		}
	}

	//Destroy last empty read
	bam_destroy1(bam_read);

	//ERASE
	{
		printf("INTERVAL %d:%d-%d %d\n", last_read_chrom + 1, interval_begin, interval_end, array_list_size(process_list));
		for(i = 0; i < array_list_size(process_list); i++)
		{
			//Get read
			bam1_t *read = array_list_get(i, process_list);

			printf("%s\n", bam1_qname(read));
		}
	}

	//Realign lastest alignments
	alig_bam_list(process_list, ref);
	array_list_clear(process_list, NULL);

	//Write lastest reads
	alig_bam_list_to_disk(write_list, out_bam_f);
	printf("Total alignments readed: %d\r", count);
	fflush(stdout);

	//Print changed cigars
	//printf("\n");
	//printf("Total CIGAR aligned: %d\n", cigar_changed);
	//fflush(stdout);

	//Memory free
	{
		//Lists
		array_list_free(write_list, NULL);
		array_list_free(process_list, NULL);

		printf("\nClosing BAM file...\n");
		bam_fclose(bam_f);
		printf("BAM closed.\n");

		printf("\nClosing reference file...\n");
		genome_free(ref);
		printf("Reference closed.\n");

		printf("Closing \"%s\" BAM file...\n", outbam);
		bam_fclose(out_bam_f);
		printf("BAM closed.\n");
	}

	return NO_ERROR;
}

ERROR_CODE
alig_bam_list(array_list_t *bam_list, genome_t* ref)
{
	//Haplotypes list
	array_list_t *haplo_list = NULL;

	//Haplotype
	aux_indel_t *haplo = NULL;
	aux_indel_t *aux_haplo = NULL;

	//Read
	bam1_t *read = NULL;
	uint32_t *read_cigar = NULL;
	char *read_seq = NULL;
	size_t read_disp_ref = SIZE_MAX;
	size_t indels;

	//Reference
	char *ref_seq = NULL;
	uint32_t flag;
	size_t ref_pos_begin = SIZE_MAX;
	size_t ref_pos_end = SIZE_MAX;
	size_t aux_pos_begin = SIZE_MAX;
	size_t aux_pos_end = SIZE_MAX;
	size_t ref_length = SIZE_MAX;

	int i, cond = 0;

	assert(bam_list);
	assert(ref);

	//If empty list return
	if(array_list_size(bam_list) == 0)
		return NO_ERROR;

	//Get reference parameters
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
		flag = (uint32_t) read->core.flag;
		genome_read_sequence_by_chr_index(ref_seq, (flag & BAM_FREVERSE) ? 1 : 0, (unsigned int)read->core.tid, &aux_pos_begin, &aux_pos_end, ref);
		assert(ref);

		ref_length = (ref_pos_end - ref_pos_begin);
	}

	//Init haplotype list
	haplo_list = array_list_new(100, 1.2f, COLLECTION_MODE_SYNCHRONIZED);

	//Get all haplotypes
	for(i = 0; i < array_list_size(bam_list); i++)
	{
		//Get read
		read = array_list_get(i, bam_list);

		//Get read cigar
		read_cigar = bam1_cigar(read);

		//This read have indels?
		cigar32_count_indels(read_cigar, read->core.n_cigar, &indels);

		if(indels == 1)	//Only one indel support
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

			//Convert read sequence to string
			read_seq = (char *) malloc((read->core.l_qseq + 1) * sizeof(char));
			new_sequence_from_bam_ref(read, read_seq, read->core.l_qseq + 1);

			//Get relative displacement from reference
			read_disp_ref = read->core.pos - ref_pos_begin;

			//Leftmost cigar
			cigar32_leftmost(ref_seq + read_disp_ref, ref_length - read_disp_ref, read_seq, read->core.l_qseq, read_cigar, read->core.n_cigar, clip_cigar, &clip_cigar_l);

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

			//Allocate haplotype
			haplo = (aux_indel_t *)malloc(sizeof(aux_indel_t) * indels);	//For now is only 1

			//Fill haplotype
			cigar32_get_indels(read->core.pos, clip_cigar, clip_cigar_l, haplo);

			//Check if haplotype is present in list
			aux_haplo = NULL;
			int h;
			for(h = 0; h < array_list_size(haplo_list); h++)
			{
				aux_haplo = array_list_get(h, haplo_list);
				assert(aux_haplo);
				if(aux_haplo->indel == haplo->indel && aux_haplo->ref_pos == haplo->ref_pos)
					break;
			}

			//Duplicate?
			if(h == array_list_size(haplo_list))
			{
				//Add to haplotype list (not duplicated)
				array_list_insert(haplo, haplo_list);
			}

			//Free read
			free(read_seq);
		}
	}

	//ERASE
	if(array_list_size(haplo_list))
	{
		char cigar_str[20];
		aux_indel_t *indel_aux;

		printf("===========================\n");
		printf("Align %d reads\n", array_list_size(bam_list));
		printf("%d haplotypes\n", array_list_size(haplo_list));
		for(i = 0; i < array_list_size(haplo_list); i++)
		{
			indel_aux = array_list_get(i, haplo_list);
			assert(indel_aux);
			cigar32_to_string(&indel_aux->indel, 1, cigar_str);
			printf("H%d: %s === Pos -> %d:%d\n", i + 1, cigar_str, read->core.tid, indel_aux->ref_pos);
		}
	}

	//Indel local realignment
	alig_bam_list_realign(bam_list, haplo_list, ref);

	//ERASE
	/*if(array_list_size(haplo_list))
	{
		if(cond)
		{
			printf("Press to continue...\n");
			getchar();
		}
	}*/

	//Free memory
	{
		free(ref_seq);
		array_list_free(haplo_list, free);
	}

	return NO_ERROR;
}

ERROR_CODE
alig_bam_list_to_disk(array_list_t *bam_list, bam_file_t *bam_f)
{
	int i;
	bam1_t *read;

	assert(bam_list);
	assert(bam_f);

	//Iterate bams
	for(i = 0; i < array_list_size(bam_list); i++)
	{
		//Get read
		read = array_list_get(i, bam_list);

		//Write read to disk
		bam_write1(bam_f->bam_fd, read);

		//Free read
		bam_destroy1(read);
	}

	//Clear list
	array_list_clear(bam_list, NULL);

	return NO_ERROR;
}

//TODO
//ESCRIBIR NUEVA FUNCION QUE ALMACENE LOS SCORE DE TODAS LECTURAS CON TODOS LOS HAPLOTIPOS
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

		//Convert read sequence to string
		read_seq = (char *) malloc((read->core.l_qseq + 1) * sizeof(char));
		new_sequence_from_bam_ref(read, read_seq, read->core.l_qseq + 1);

		//Get qualities
		read_quals = (char *) malloc((read->core.l_qseq + 1) * sizeof(char));
		new_quality_from_bam_ref(read, 0, read_quals, read->core.l_qseq + 1);

		//Get relative displacement from reference
		//read_disp_ref = read->core.pos - ref_pos_begin;

		//ERASE Print things
		{
			char erase_str[200];
			size_t disp;
			cigar32_to_string(bam1_cigar(read), read->core.n_cigar, erase_str);
			printf("---------------------\n");
			printf("Test %s - %d:%d - %s\n", bam1_qname(read), read->core.tid + 1, read->core.pos + 1, erase_str);
		}

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
			{
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
			}
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
					{
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
					}
				}

				//If reference maps better, not count haplotype score
				if(m_score[(i * m_ldim) + (j+1)] > m_score[(i * m_ldim)])
				{
					m_score[(i * m_ldim) + (j+1)] = UINT32_MAX;
					m_pos[(i * m_ldim) + (j+1)] = SIZE_MAX;
				}

				//ERASE Print things
				/*if(read_disp_ref != SIZE_MAX)
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
	{
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
	}

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
			if(m_score[(i * m_ldim)] != UINT32_MAX)	//If match reference this score not count
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
		{
			//Print haplotype score
			printf("H%d === total score = %d\n", j + 1, v_hscore[j]);
		}
	}

	//Get best haplotype
	haplo = array_list_get(best_haplo_index, haplotype_list);

	//Iterate bams to change its CIGARS
	for(i = 0; i < bam_list_l; i++)
	{
		//Get read
		read = array_list_get(i, bam_list);
		assert(read);

		//Is reference matching better?
		score = m_score[(i * m_ldim)];
		score2 = m_score[(i * m_ldim) + (best_haplo_index + 1)];
		if(score < score2)
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
		{
			char cigar_str[50];
			char cigar_str2[50];
			cigar32_to_string(bam1_cigar(read), read->core.n_cigar, cigar_str);
			cigar32_to_string(best_cigar, best_cigar_l, cigar_str2);
			//Print new cigar
			printf("%s Cigar: %s => %s\n", bam1_qname(read), cigar_str, cigar_str2);
		}

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

