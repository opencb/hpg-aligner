/*
 * region.c
 *
 *  Created on: Dec 10, 2013
 *      Author: rmoreno
 */

#include "assert.h"
#include "aux_library.h"

int
region_get_from_read(const bam1_t *alig, int32_t *chrom, long int *r_pos, long int *r_end_pos)
{
	//CIGAR
	uint32_t *cigar;
	uint32_t cigar_l;

	//Region
	long int reg_pos = -1;
	long int reg_end_pos;

	//Read
	int read_pos;

	int i;

	//CHECK ARGUMENTS
	{
		//Check nulls
		assert(alig);
		assert(chrom);
		assert(r_pos);
	}

	//Get CIGAR
	cigar = bam1_cigar(alig);
	cigar_l = alig->core.n_cigar;

	//Get read position
	read_pos = alig->core.pos + 1;
	//printf("read pos: %d  ", read_pos);

	//Search for indels
	uint32_t elem, type;
	int disp = 0;
	for(i = 0; i < cigar_l; i++)
	{
		elem = cigar[i] >> BAM_CIGAR_SHIFT;	//Get number of bases from cigar
		type = cigar[i] & BAM_CIGAR_MASK;	//Get type from cigar

		switch(type)
		{
		case BAM_CINS:	//Insertion
			//printf("I");
			//Is first CIGAR elem?
			if(i == 0)
				break;

			//Set initial position
			if(reg_pos == -1)
			{
				reg_pos = read_pos + disp - 1;
			}

			//Set final position
			reg_end_pos = read_pos + disp;

			break;
		case BAM_CDEL:	//Deletion
			//printf("D");
			//Set initial position
			if(reg_pos == -1)
			{
				reg_pos = read_pos + disp;
			}

			//Set final position
			reg_end_pos = read_pos + disp + elem - 1;

			disp += elem;
			break;

		case BAM_CHARD_CLIP:
			if(i == 0)
			{
				break;
			}

		default:
			//printf("X");
			//Increment displacement
			disp += elem;
			break;
		}
	}
	//printf("\n");

	//Set output
	if(reg_pos != -1)
	{
		*chrom = alig->core.tid;
		*r_pos = reg_pos;
		*r_end_pos = reg_end_pos;
		return 0;
	}
	else
	{
		*chrom = -1;
		*r_pos = -1;
		return -1;
	}

}

int
region_get_from_batch(const bam_batch_t* batch)
{
	int i;
	bam1_t *alig;
	int32_t chrom;
	long int r_pos;
	long int r_end_pos;

	for(i = 0; i < batch->num_alignments; i++)
	{
		alig = batch->alignments_p[i];

		if(alig->core.qual != 0 && !(alig->core.flag & BAM_FSECONDARY) && alig->core.mtid == alig->core.tid)
		{
			region_get_from_read(alig, &chrom, &r_pos, &r_end_pos);
			if(r_pos != -1)
			{
				//Region found
				printf("Region = %d : %d-%d\n", chrom+1, r_pos, r_end_pos);
			}
		}
	}
}

int
region_get_from_file(const char *bam_path)
{
	bam_file_t *bam_f;
	bam_batch_t *batch;

	//Open bam
	printf("Opening BAM from \"%s\" ...\n", bam_path);
	bam_f = bam_fopen(bam_path);
	printf("BAM opened!...\n");

	//Read batch
	batch = bam_batch_new(10000000, SINGLE_CHROM_BATCH);
	bam_fread_max_size(batch, 10000000, 1, bam_f);

	//Process batch
	region_get_from_batch(batch);

	//Free batch
	bam_batch_free(batch, 1);

	//Memory free
	printf("\nClosing BAM file...\n");
	bam_fclose(bam_f);
	printf("BAM closed.\n");
}

