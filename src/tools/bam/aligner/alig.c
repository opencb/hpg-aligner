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
