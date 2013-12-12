/*
 * alig.c
 *
 *  Created on: Dec 12, 2013
 *      Author: rmoreno
 */

#include "alig.h"

int
alig_cigar_leftmost(char *ref, char *read, size_t length, uint32_t *cigar, size_t cigar_l, uint32_t *new_cigar, size_t *new_cigar_l)
{
	//M blocks
	size_t blocks_c;

	//New cigar
	uint32_t aux_cigar[30];
	size_t aux_cigar_l;

	assert(ref);
	assert(read);
	assert(length > 0);
	assert(cigar_l < 30);
	assert(new_cigar);
	assert(new_cigar_l);

	//Count blocks
	alig_cigar_count_m_blocks(cigar, cigar_l, &blocks_c);

	//Only procceed if 1 indel (2 M blocks)
	if(blocks_c == 2)
	{
		//Unclip cigar
		alig_cigar_unclip(cigar, cigar_l, aux_cigar, &aux_cigar_l);

		//Leftmost CIGAR
		//...

		//Set output
		*new_cigar_l = 0;
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

