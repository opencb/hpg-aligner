/*
 * alig.c
 *
 *  Created on: Dec 12, 2013
 *      Author: rmoreno
 */

#include "alig.h"

/**
 * BAM REALIGN
 */

int
alig_bam_file(const char *bam_path, const char *ref_name, const char *ref_path)
{
	bam_file_t *bam_f = NULL;
	bam_batch_t *batch = NULL;
	genome_t* ref = NULL;
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

	//Read batch
	batch = bam_batch_new(10000000, SINGLE_CHROM_BATCH);
	bam_fread_max_size(batch, 10000000, 1, bam_f);

	//Process batch
	alig_bam_batch(batch, ref);

	//Free batch
	bam_batch_free(batch, 1);

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
			ref_seq = (char *)malloc(sizeof(char) * alig->core.l_qseq + 30);
			bam_seq = (char *)malloc(sizeof(char) * alig->core.l_qseq + 1);

			//Get sequence
			new_sequence_from_bam_ref(alig, bam_seq, alig->core.l_qseq);

			//Get reference
			flag = (uint32_t) alig->core.flag;
			init_pos = alig->core.pos + 1;
			end_pos = alig->core.pos + alig->core.l_qseq + 30;
			genome_read_sequence_by_chr_index(ref_seq, (flag & BAM_FREVERSE) ? 1 : 0, (unsigned int)alig->core.tid, &init_pos, &end_pos, ref);

			//Leftmost CIGARS
			alig_cigar_leftmost(ref_seq, bam_seq, alig->core.l_qseq, bam1_cigar(alig), alig->core.n_cigar, new_cigar, &new_cigar_l);

			//ERASE
			if(memcmp(bam1_cigar(alig), new_cigar, sizeof(uint32_t) * alig->core.n_cigar) != 0)
			{
				char str_cigar[200];
				char str_new_cigar[200];
				alig_aux_cigar32_to_string(bam1_cigar(alig), alig->core.n_cigar, str_cigar);
				alig_aux_cigar32_to_string(new_cigar, new_cigar_l, str_new_cigar);
				printf("NEW CIGAR: %s => %s, %d => %d\n", str_cigar, str_new_cigar, alig->core.n_cigar,new_cigar_l);
				printf("POSITION: %d \n", alig->core.pos);
				printf("***********************\n");
			}

			//Free reference seq
			free(bam_seq);
			free(ref_seq);
		}
	}
}

/**
 * CIGAR
 */

int
alig_cigar_leftmost(char *ref, char *read, size_t read_l, uint32_t *cigar, size_t cigar_l, uint32_t *new_cigar, size_t *new_cigar_l)
{
	//M blocks
	size_t blocks_c;

	//Indels
	size_t indels_c;
	size_t indel_index;
	size_t indel_l;

	//New cigar
	uint32_t unclip_cigar[30];
	uint32_t aux_cigar[30];
	size_t unclip_cigar_l;

	//References
	char *orig_ref;
	char *aux_ref;

	//Counters
	int count;

	//ERASE
	char str_cigar[200];
	char str_new_cigar[200];
	int printed = 0;

	assert(ref);
	assert(read);
	assert(read_l > 0);
	assert(cigar);
	assert(cigar_l < 30);
	assert(new_cigar);
	assert(new_cigar_l);

	//Unclip cigar
	alig_cigar_unclip(cigar, cigar_l, unclip_cigar, &unclip_cigar_l);

	//Count blocks and indels
	alig_cigar_count_all(unclip_cigar, unclip_cigar_l, &blocks_c, &indels_c, &indel_index);

	//Get indel length
	indel_l = unclip_cigar[indel_index] >> BAM_CIGAR_SHIFT;

	//By default return original CIGAR
	memcpy(new_cigar, cigar, cigar_l * sizeof(uint32_t));
	*new_cigar_l = cigar_l;

	//Only procceed if 1 indel (2 M blocks)
	if(blocks_c == 2 && indels_c == 1 && indel_l)
	{
		//Leftmost CIGAR
		{
			//Get reference for original CIGAR
			orig_ref = (char *)malloc(sizeof(char) * (read_l + 1));
			aux_ref = (char *)malloc(sizeof(char) * (read_l + 1));
			alig_aux_cigar32_create_ref(unclip_cigar, unclip_cigar_l, ref, read, read_l, orig_ref);

			//Shift left CIGAR
			alig_aux_cigar32_shift_left_indel(unclip_cigar, unclip_cigar_l, indel_index, aux_cigar);

			//Get new CIGAR ref
			alig_aux_cigar32_create_ref(aux_cigar, unclip_cigar_l, ref, read, read_l, aux_ref);

			//Is a valid ref?
			if(!memcmp(aux_ref, orig_ref, read_l * sizeof(char)))
			{
				//Equal so is a valid CIGAR
				memcpy(new_cigar, aux_cigar, sizeof(uint32_t) * unclip_cigar_l);
				*new_cigar_l = unclip_cigar_l;
				count = indel_l;

				//ERASE
				{
					read[read_l] = '\0';
					orig_ref[read_l] = '\0';
					alig_aux_cigar32_to_string(unclip_cigar, unclip_cigar_l, str_new_cigar);
					alig_aux_cigar32_to_string(cigar, cigar_l, str_cigar);
					printf("CIGAR = %s, Indel L: %d\n", str_cigar, indel_l);
					printf("CIGAR*= %s\n", str_new_cigar);
					printf("READ -> %s - L: %d\n", read, read_l);
					printf("REF  -> %s\n", ref);
					printf("REF* -> %s - %s\n", orig_ref, str_new_cigar);
					aux_ref[read_l] = '\0';
					alig_aux_cigar32_to_string(aux_cigar, unclip_cigar_l, str_new_cigar);
					printf("REF%d -> %s - %s ::: Retry - %d\n", 1, aux_ref, str_new_cigar, count);
					printed = 1;
				}
			}
			else
			{
				count = indel_l - 1;
			}

			if((aux_cigar[indel_index - 1] >> BAM_CIGAR_SHIFT) == 0)
				count = 0;

			int j = 1;
			while(count > 0)
			{
				count--;
				j++;

				//Shift left CIGAR
				alig_aux_cigar32_shift_left_indel(aux_cigar, unclip_cigar_l, indel_index, aux_cigar);

				//Get new CIGAR ref
				alig_aux_cigar32_create_ref(aux_cigar, unclip_cigar_l, ref, read, read_l, aux_ref);

				//Is a valid ref?
				if(!memcmp(aux_ref, orig_ref, read_l * sizeof(char)))
				{
					//Equal so is a valid CIGAR
					memcpy(new_cigar, aux_cigar, sizeof(uint32_t) * unclip_cigar_l);
					*new_cigar_l = unclip_cigar_l;
					count = indel_l;

					//ERASE
					{
						if(!printed)
						{
							read[read_l] = '\0';
							orig_ref[read_l] = '\0';
							alig_aux_cigar32_to_string(unclip_cigar, unclip_cigar_l, str_new_cigar);
							alig_aux_cigar32_to_string(cigar, cigar_l, str_cigar);
							printf("CIGAR = %s, Indel L: %d\n", str_cigar, indel_l);
							printf("CIGAR*= %s\n", str_new_cigar);
							printf("READ -> %s - L: %d\n", read, read_l);
							printf("REF  -> %s\n", ref);
							printf("REF* -> %s - %s\n", orig_ref, str_new_cigar);
							printed = 1;
						}
						aux_ref[read_l] = '\0';
						alig_aux_cigar32_to_string(aux_cigar, unclip_cigar_l, str_new_cigar);
						printf("REF%d -> %s - %s ::: Retry - %d\n", j, aux_ref, str_new_cigar, count);
					}
				}

				if((aux_cigar[indel_index - 1] >> BAM_CIGAR_SHIFT) == 0)
					count = 0;
			}

			//Free
			free(orig_ref);
			free(aux_ref);
		}
	}
	else
	{
		//Return original CIGAR
		memcpy(new_cigar, cigar, cigar_l * sizeof(uint32_t));
		*new_cigar_l = cigar_l;
	}

	return 0;
}


int
alig_cigar_unclip(uint32_t *cigar, size_t cigar_l, uint32_t *new_cigar, size_t *new_cigar_l)
{
	int i;
	int c_count;
	int c_type;
	int new_cigar_i;

	assert(cigar);
	assert(cigar_l > 0);
	assert(new_cigar);
	assert(new_cigar_l);

	//Iterate cigar elements
	new_cigar_i = 0;
	for(i = 0; i < cigar_l; i++)
	{
		c_count = cigar[i] >> BAM_CIGAR_SHIFT;	//Get number of bases from cigar
		c_type = cigar[i] & BAM_CIGAR_MASK;	//Get type from cigar

		//Skip if empty cigar
		if(c_count > 0)
		{
			//Is clip?
			if(c_type != BAM_CSOFT_CLIP
					&& c_type != BAM_CHARD_CLIP
					&& c_type != BAM_CPAD)
			{
				//Add to new cigar
				new_cigar[new_cigar_i] = cigar[i];
				new_cigar_i++;
			}
		}
	}

	//Set output cigar length
	*new_cigar_l = new_cigar_i;

	return 0;
}

int
alig_cigar_count_m_blocks(uint32_t *cigar, size_t cigar_l, size_t *blocks)
{
	int i;
	int c_count;
	int c_type;
	size_t block_c;
	char *str_cigar;

	assert(cigar);
	assert(cigar_l > 0);
	assert(blocks);

	//Iterate cigar elements
	block_c = 0;
	for(i = 0; i < cigar_l; i++)
	{
		c_count = cigar[i] >> BAM_CIGAR_SHIFT;	//Get number of bases from cigar
		c_type = cigar[i] & BAM_CIGAR_MASK;	//Get type from cigar

		switch(c_type)
		{

		case BAM_CMATCH:
		case BAM_CEQUAL:
		case BAM_CDIFF:
			block_c++;
			break;
		}
	}

	//Set output
	*blocks = block_c;

	return 0;
}

int
alig_cigar_count_indels(uint32_t *cigar, size_t cigar_l, size_t *indels)
{
	int i;
	int c_count;
	int c_type;
	size_t indel_c;
	char *str_cigar;

	assert(cigar);
	assert(cigar_l > 0);
	assert(indels);

	//Iterate cigar elements
	indel_c = 0;
	for(i = 0; i < cigar_l; i++)
	{
		c_count = cigar[i] >> BAM_CIGAR_SHIFT;	//Get number of bases from cigar
		c_type = cigar[i] & BAM_CIGAR_MASK;	//Get type from cigar

		switch(c_type)
		{
		case BAM_CINS:
		case BAM_CDEL:
			indel_c++;
			break;
		}
	}

	//Set output
	*indels = indel_c;

	return 0;
}

int
alig_cigar_count_all(uint32_t *cigar, size_t cigar_l, size_t *m_blocks, size_t *indels, size_t *first_indel_index)
{
	int i;
	int c_count;
	int c_type;
	size_t indel_c;
	size_t m_block_c;
	size_t indel_index;
	char *str_cigar;

	assert(cigar);
	assert(cigar_l > 0);
	//assert(indels);
	//assert(first_indel_offset);

	//Iterate cigar elements
	indel_c = 0;
	m_block_c = 0;
	indel_index = 0;
	for(i = 0; i < cigar_l; i++)
	{
		c_count = cigar[i] >> BAM_CIGAR_SHIFT;	//Get number of bases from cigar
		c_type = cigar[i] & BAM_CIGAR_MASK;	//Get type from cigar

		switch(c_type)
		{
		case BAM_CINS:
		case BAM_CDEL:
			indel_c++;
			if(!indel_index)
				indel_index = i;
			break;

		case BAM_CMATCH:
		case BAM_CEQUAL:
		case BAM_CDIFF:
			m_block_c++;
			break;
		}
	}

	//Set output
	if(indels)
		*indels = indel_c;
	if(m_blocks)
		*m_blocks = m_block_c;
	if(first_indel_index)
		*first_indel_index = indel_index;

	return 0;
}

