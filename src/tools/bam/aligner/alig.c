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

int
alig_bam_file(char *bam_path, char *ref_name, char *ref_path)
{
	bam_file_t *bam_f = NULL;
	bam_batch_t *batch = NULL;
	genome_t* ref = NULL;
	uint64_t count = 0;
	int i;

	//Open bam
	printf("Opening BAM from \"%s\" ...\n", bam_path);
	bam_f = bam_fopen(bam_path);
	printf("BAM opened!...\n");

	//Open reference genome
	printf("Opening reference genome from \"%s%s\" ...\n", ref_path, ref_name);
	ref = genome_new(ref_name, ref_path);
	assert(ref);
	printf("Reference opened!...\n");

	do
	{
		if(batch)
		{
			//Free batch
			bam_batch_free(batch, 1);
		}

		//Read batch
		batch = bam_batch_new(100000000, SINGLE_CHROM_BATCH);
		bam_fread_max_size(batch, 100000000, 1, bam_f);

		//Process batch
		alig_bam_batch(batch, ref);

		//Update read counter
		count += batch->num_alignments;

		//Show total progress
		printf("Total alignments readed: %d\r", count);
		fflush(stdout);

	}while(batch && batch->num_alignments != 0);

	if(batch)
	{
		//Free batch
		bam_batch_free(batch, 1);
	}

	//Print changed cigars
	printf("\n");
	printf("Total CIGAR aligned: %d\n", cigar_changed);
	fflush(stdout);

	//Memory free
	printf("\nClosing BAM file...\n");
	bam_fclose(bam_f);
	genome_free(ref);
	printf("BAM closed.\n");
}

int
alig_bam_batch(bam_batch_t* batch, genome_t* ref)
{
	int i;
	bam1_t *alig;

	//Positions
	size_t init_pos, end_pos;

	//Reference
	char *ref_seq;
	char *bam_seq;
	uint32_t flag;

	//CIGAR
	uint32_t new_cigar[60];
	size_t new_cigar_l;

	assert(batch);
	assert(ref);

	for(i = 0; i < batch->num_alignments; i++)
	{
		alig = batch->alignments_p[i];

		//FILTERS: MAP QUALITY = 0, NOT PRIMARY ALIGNMENT, DIFFERENT MATE CHROM
		if(alig->core.qual != 0 && !(alig->core.flag & BAM_FSECONDARY) && alig->core.mtid == alig->core.tid)
		{
			//Allocate
			ref_seq = (char *)malloc(sizeof(char) * alig->core.l_qseq + 32);
			bam_seq = (char *)malloc(sizeof(char) * alig->core.l_qseq + 2);

			//Get sequence
			new_sequence_from_bam_ref(alig, bam_seq, alig->core.l_qseq);

			//Get reference
			flag = (uint32_t) alig->core.flag;
			init_pos = alig->core.pos + 1;
			end_pos = alig->core.pos + alig->core.l_qseq + 30;
			genome_read_sequence_by_chr_index(ref_seq, (flag & BAM_FREVERSE) ? 1 : 0, (unsigned int)alig->core.tid, &init_pos, &end_pos, ref);

			//Leftmost CIGARS
			cigar_leftmost(ref_seq, bam_seq, alig->core.l_qseq, bam1_cigar(alig), alig->core.n_cigar, new_cigar, &new_cigar_l);

			//ERASE
			if(memcmp(bam1_cigar(alig), new_cigar, sizeof(uint32_t) * alig->core.n_cigar) != 0)
			{
				/*char str_cigar[200];
				char str_new_cigar[200];
				alig_aux_cigar32_to_string(bam1_cigar(alig), alig->core.n_cigar, str_cigar);
				alig_aux_cigar32_to_string(new_cigar, new_cigar_l, str_new_cigar);
				printf("NEW CIGAR: %s => %s, %d => %d\n", str_cigar, str_new_cigar, alig->core.n_cigar,new_cigar_l);
				printf("POSITION: %d:%d \n", alig->core.tid, alig->core.pos);
				printf("***********************\n");*/
				cigar_changed++;
			}

			//Free reference seq
			free(bam_seq);
			free(ref_seq);
		}
	}
}

ERROR_CODE
alig_bam_file2(char *bam_path, char *ref_name, char *ref_path)
{
	bam_file_t *bam_f = NULL;
	bam_file_t *out_bam_f = NULL;
	bam_batch_t *batch = NULL;
	genome_t* ref = NULL;
	ERROR_CODE err;
	char outbam[30] = "output.bam";

	//BAM alignment
	bam1_t* bam_read = NULL;
	size_t bytes;
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
	write_list = array_list_new(100, 1.5f, COLLECTION_MODE_SYNCHRONIZED);
	process_list = array_list_new(100, 1.5f, COLLECTION_MODE_SYNCHRONIZED);

	//Read first alignment
	bam_read = bam_init1();
	bytes = bam_read1(bam_f->bam_fd, bam_read);
	last_read_chrom = bam_read->core.tid;

	//Put alignment to write list
	array_list_insert(bam_read, write_list);

	count = 1;
	while(bytes > 0 /*&& count < 10000*/)
	{
		//FILTERS: MAP QUALITY = 0, NOT PRIMARY ALIGNMENT, DIFFERENT MATE CHROM
		if(bam_read
				&& bam_read->core.qual != 0
				&& !(bam_read->core.flag & BAM_FSECONDARY)
				&& bam_read->core.mtid == bam_read->core.tid)
		{
			//Get read position
			read_pos = bam_read->core.pos;

			//Check if chrom is different
			if(last_read_chrom != bam_read->core.tid)
			{
				//Realign alignments in list
				alig_bam_list(process_list);
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
				cigar_to_string(bam1_cigar(bam_read), bam_read->core.n_cigar, cigarro);
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
				alig_bam_list(process_list);
				array_list_clear(process_list, NULL);

				//Reset interval
				interval_begin = SIZE_MAX;
				interval_end = SIZE_MAX;
			}

			//Save last read chrom
			last_read_chrom = bam_read->core.tid;
		}

		//Read next alignment
		bam_read = bam_init1();
		bytes = bam_read1(bam_f->bam_fd, bam_read);

		//Put alignment to write list
		array_list_insert(bam_read, write_list);

		//If write list is big enought and process list is empty: write to disk
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

	//Write lastest reads
	alig_bam_list_to_disk(write_list, out_bam_f);

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
alig_bam_list(array_list_t *bam_list)
{
	assert(bam_list);

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

