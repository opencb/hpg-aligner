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
	while(bytes > 0 /*&& count < 100000*/)
	{
		//Put alignment to write list
		array_list_insert(bam_read, write_list);

		//FILTERS: MAP QUALITY = 0, NOT PRIMARY ALIGNMENT, DIFFERENT MATE CHROM
		if(bam_read
				&& bam_read->core.qual != 0
				&& !(bam_read->core.flag & BAM_FSECONDARY)
				&& bam_read->core.mtid == bam_read->core.tid)
		{

			//Filter alignments with one indel
			cigar32_count_indels(bam1_cigar(bam_read), bam_read->core.n_cigar, &indels);
			if(indels == 1)
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

				//Get interval for this alignment
				if(region_get_from_bam1(bam_read, &interval_read_begin, &interval_read_end))
				{
					//ERROR HAPPENED
					char cigarro[100];
					cigar32_to_string(bam1_cigar(bam_read), bam_read->core.n_cigar, cigarro);
					//Invalid CIGAR
					printf("WARNING: Invalid CIGAR in %s - P:%d - l:%d - %s\n",
							bam1_qname(bam_read), bam_read->core.pos, bam_read->core.n_cigar, cigarro);

					//Read next alignment
					bam_destroy1(bam_read);
					bam_read = bam_init1();
					bytes = bam_read1(bam_f->bam_fd, bam_read);

					continue;
				}

				//Cases
				if(interval_begin == SIZE_MAX)
				{
					//NO INTERVAL CASE
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
				else if (interval_end > read_pos)
				{
					//ALIGNMENT INSIDE INTERVAL CASE

					//This alignment have an interval?
					if(interval_read_begin != SIZE_MAX)
					{
						//Interval found
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

					//Print interval
					//printf("%d:%d-%d\n", last_read_chrom + 1, interval_begin, interval_end);

					//Realign alignments in process list
					alig_bam_list(process_list, ref);
					array_list_clear(process_list, NULL);

					//Reset interval
					interval_begin = SIZE_MAX;
					interval_end = SIZE_MAX;
				}

				//Save last read chrom
				last_read_chrom = bam_read->core.tid;

			}	//Indels if
		}	// Filters if

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
				if(!cond)
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
				printf("%s\n", ref_strx2);
			}

			//Allocate haplotype
			haplo = (aux_indel_t *)malloc(sizeof(aux_indel_t) * indels);	//For now is only 1

			//Fill haplotype
			cigar32_get_indels(read->core.pos, clip_cigar, clip_cigar_l, haplo);

			//Add to haplotype list
			array_list_insert(haplo, haplo_list);

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
			printf("H%d: %s === Pos -> %d:%d\n", i, cigar_str, read->core.tid, indel_aux->ref_pos);
		}
	}

	//Indel local realignment
	alig_bam_list_realign(bam_list, haplo_list, ref);

	//ERASE
	if(array_list_size(haplo_list))
	{
		if(cond)
		{
			printf("Press to continue...\n");
			getchar();
		}
	}

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

ERROR_CODE
alig_bam_list_realign(array_list_t *bam_list, array_list_t *haplotype_list, genome_t* ref)
{
	int i, j;
	uint32_t *H_miss;
	uint32_t *ref_miss;

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
	size_t read_seq_ref_l;
	size_t read_disp_ref = SIZE_MAX;

	//Comparation
	char *comp_aux = NULL;
	uint32_t comp_cigar[MAX_CIGAR_LENGTH];
	size_t comp_cigar_l;
	size_t bases;

	//Haplotype
	aux_indel_t *haplo;

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

	//Allocate miss arrays
	H_miss = (uint32_t *)malloc(sizeof(size_t) * array_list_size(bam_list));
	ref_miss = (uint32_t *)malloc(sizeof(size_t) * array_list_size(bam_list));
	for(i = 0; i < array_list_size(bam_list); i++)
	{
		H_miss[i] = 0;
		ref_miss[i] = 0;
	}

	//Iterate bam list
	for(i = 0; i < array_list_size(bam_list); i++)
	{
		//Get read
		read = array_list_get(i, bam_list);

		//Convert read sequence to string
		read_seq = (char *) malloc((read->core.l_qseq + 1) * sizeof(char));
		new_sequence_from_bam_ref(read, read_seq, read->core.l_qseq + 1);

		//Get relative displacement from reference
		read_disp_ref = read->core.pos - ref_pos_begin;

		//Unclip cigar
		cigar32_unclip(bam1_cigar(read), read->core.n_cigar, comp_cigar, &comp_cigar_l);

		//Count unclipped bases
		cigar32_count_nucleotides_not_clip(comp_cigar, comp_cigar_l, &bases);

		//Create new cigar with '$bases'M
		comp_cigar[0] = bases << BAM_CIGAR_SHIFT;	//ex: 108M
		comp_cigar_l = 1;

		//Reclip cigar
		cigar32_reclip(bam1_cigar(read), read->core.n_cigar, comp_cigar, comp_cigar_l, comp_cigar, &comp_cigar_l);

		//Get reference transform for this read
		cigar32_create_ref(comp_cigar, comp_cigar_l, ref_seq + read_disp_ref, ref_length - read_disp_ref, read_seq, read->core.l_qseq, read_seq_ref, &read_seq_ref_l);

		//Test with reference
		nucleotide_compare(read_seq, read_seq_ref, read_seq_ref_l, comp_aux, ref_miss + i);

		//ERASE Print things
		{
			char erase_str[200];
			size_t disp;
			cigar32_to_string(bam1_cigar(read), read->core.n_cigar, erase_str);
			printf("---------\n", bam1_qname(read));
			printf("Test %s - %d:%d - %s\n", bam1_qname(read), read->core.tid + 1, read->core.pos + 1, erase_str);
			printf("---------\n");
			cigar32_count_clip_displacement(comp_cigar, comp_cigar_l, &disp);
			cigar32_to_string(comp_cigar, comp_cigar_l, erase_str);
			printf("Test with reference === MISS: %d - Using CIGAR: %s\n", ref_miss[i], erase_str);
			//memcpy(erase_str, ref_seq + read_disp_ref, sizeof(char) * read->core.l_qseq);
			//erase_str[read->core.l_qseq] = '\0';
			//printf("Ref : ");
			//for(int z = 0; z < disp; z++) printf(" ");
			//printf("%s\n", erase_str);
			//memcpy(erase_str, read_seq, sizeof(char) * read->core.l_qseq);
			//erase_str[read->core.l_qseq] = '\0';
			//printf("Read: %s\n", erase_str);
			//printf("Ref*: %s\n", read_seq_ref);
			//cigar32_to_string(&indel_aux->indel, 1, cigar_str);
		}

		//Dont match reference perfectly?
		if(ref_miss[i])
		{
			//Iterate haplotypes
			for(j = 0; j < array_list_size(haplotype_list); j++)
			{
				//Get haplotype
				haplo = array_list_get(j, haplotype_list);
				assert(haplo);

				//Create cigar for this haplotype
				cigar32_from_haplo(bam1_cigar(read), read->core.n_cigar, haplo, read->core.pos, comp_cigar, &comp_cigar_l);


				//printf("---------\n");
				//printf("Orig CIGAR: %s\n", cigar_str);

				//printf("Comp CIGAR: %s\n", cigar_str);

				//Get haplotype reference transform
				cigar32_create_ref(comp_cigar, comp_cigar_l, ref_seq + read_disp_ref, ref_length - read_disp_ref, read_seq, read->core.l_qseq, read_seq_ref, &read_seq_ref_l);

				//Compare and miss with haplotype
				nucleotide_compare(read_seq, read_seq_ref, read_seq_ref_l, comp_aux, H_miss + i);

				//ERASE Print things
				{
					char cigar_str[50];
					//char erase_str[200];
					size_t disp;
					cigar32_count_clip_displacement(comp_cigar, comp_cigar_l, &disp);
					cigar32_to_string(comp_cigar, comp_cigar_l, cigar_str);
					printf("Test with H%d === MISS: %d - Using CIGAR: %s\n", j+1, H_miss[i], cigar_str);
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
		}
	}

	//Free
	{
		free(ref_seq);
		free(read_seq_ref);
		free(H_miss);
		free(ref_miss);
		//free(comp_cigar);
	}

	return NO_ERROR;
}

