/*
 * aux_cigar.c
 *
 *  Created on: Dec 16, 2013
 *      Author: rmoreno
 */

#include "aux_cigar.h"

/**
 * CIGAR GENERATION
 */

/**
 * Leftalign one read cigar.
 */
ERROR_CODE
cigar32_leftmost(char *ref, size_t ref_l, char *read, size_t read_l, uint32_t *cigar, size_t cigar_l, uint32_t *out_cigar, size_t *out_cigar_l)
{
	//M blocks
	size_t blocks_c;

	//Indels
	size_t indels_c;
	size_t indel_index;
	size_t indel_l;

	//New cigar
	uint32_t unclip_cigar[MAX_CIGAR_LENGTH];
	uint32_t reclip_cigar[MAX_CIGAR_LENGTH];
	uint32_t aux_cigar[MAX_CIGAR_LENGTH];
	size_t unclip_cigar_l;
	size_t reclip_cigar_l;
	size_t aux_cigar_l;

	//References
	char *orig_ref;
	size_t orig_ref_l;
	char *aux_ref;
	size_t aux_ref_l;

	//Counters
	int retry;

	//ERASE
	char str_cigar[MAX_CIGAR_LENGTH*3];
	char str_new_cigar[MAX_CIGAR_LENGTH*3];
	int printed = 0;

	assert(ref);
	assert(ref_l > 0);
	assert(read);
	assert(read_l > 0);
	assert(cigar);
	//assert(cigar_l < 30);
	assert(out_cigar);
	assert(out_cigar_l);

	//Unclip cigar
	cigar32_unclip(cigar, cigar_l, unclip_cigar, &unclip_cigar_l);

	//Count blocks and indels
	cigar32_count_all(unclip_cigar, unclip_cigar_l, &blocks_c, &indels_c, &indel_index);

	//Get indel length
	indel_l = unclip_cigar[indel_index] >> BAM_CIGAR_SHIFT;

	//By default return original CIGAR
	memcpy(out_cigar, cigar, cigar_l * sizeof(uint32_t));
	*out_cigar_l = cigar_l;

	//Only procceed if 1 indel (2 M blocks)
	if(cigar_l < 30 && blocks_c == 2 && indels_c == 1 && indel_l)
	{
		//Leftmost CIGAR
		{
			//Get reference for original CIGAR
			orig_ref = (char *)malloc(sizeof(char) * (read_l + 1));
			aux_ref = (char *)malloc(sizeof(char) * (read_l + 1));
			cigar32_create_ref(cigar, cigar_l, ref, ref_l, read, read_l, orig_ref, &orig_ref_l);

			//Shift left CIGAR
			cigar32_shift_left_indel(unclip_cigar, unclip_cigar_l, indel_index, aux_cigar, &aux_cigar_l);

			//Reclip cigar
			cigar32_reclip(cigar, cigar_l, aux_cigar, aux_cigar_l, reclip_cigar, &reclip_cigar_l);

			//Get new CIGAR ref
			cigar32_create_ref(reclip_cigar, reclip_cigar_l, ref, ref_l, read, read_l, aux_ref, &aux_ref_l);

			//Is a valid ref?
			if(!memcmp(aux_ref, orig_ref, orig_ref_l * sizeof(char)))
			{
				//Equal so is a valid CIGAR
				memcpy(out_cigar, reclip_cigar, sizeof(uint32_t) * reclip_cigar_l);
				*out_cigar_l = reclip_cigar_l;
				retry = indel_l;

				//ERASE
				/*{
					//read[read_l] = '\0';
					//orig_ref[read_l] = '\0';
					cigar32_to_string(reclip_cigar, reclip_cigar_l, str_new_cigar);
					cigar32_to_string(cigar, cigar_l, str_cigar);
					printf("CIGAR = %s, Indel L: %d\n", str_cigar, indel_l);
					printf("CIGAR*= %s\n", str_new_cigar);
					printf("READ -> %s - L: %d\n", read, read_l);
					printf("REF  -> %s\n", ref);
					printf("REF* -> %s - %s\n", orig_ref, str_new_cigar);
					cigar32_to_string(aux_cigar, reclip_cigar_l, str_new_cigar);
					printf("REF%d -> %s - %s ::: Retry - %d\n", 1, aux_ref, str_new_cigar, retry);
					printed = 1;
				}*/
			}
			else
			{
				retry = indel_l - 1;
			}

			if((aux_cigar[indel_index - 1] >> BAM_CIGAR_SHIFT) == 0)
				retry = 0;

			int j = 1;
			while(retry > 0)
			{
				retry--;
				j++;

				//Shift left CIGAR
				cigar32_shift_left_indel(aux_cigar, unclip_cigar_l, indel_index, aux_cigar, &aux_cigar_l);

				//Reclip cigar
				cigar32_reclip(cigar, cigar_l, aux_cigar, aux_cigar_l, reclip_cigar, &reclip_cigar_l);

				//Get new CIGAR ref
				cigar32_create_ref(reclip_cigar, reclip_cigar_l, ref, ref_l, read, read_l, aux_ref, &aux_ref_l);

				//Is a valid ref?
				if(!memcmp(aux_ref, orig_ref, orig_ref_l * sizeof(char)))
				{
					//Equal so is a valid CIGAR
					memcpy(out_cigar, reclip_cigar, sizeof(uint32_t) * reclip_cigar_l);
					*out_cigar_l = reclip_cigar_l;
					retry = indel_l;

					//ERASE
					/*{
						if(!printed)
						{
							//read[read_l] = '\0';
							//orig_ref[read_l] = '\0';
							cigar32_to_string(unclip_cigar, unclip_cigar_l, str_new_cigar);
							cigar32_to_string(cigar, cigar_l, str_cigar);
							printf("CIGAR = %s, Indel L: %d\n", str_cigar, indel_l);
							printf("CIGAR*= %s\n", str_new_cigar);
							printf("READ -> %s - L: %d\n", read, read_l);
							printf("REF  -> %s\n", ref);
							printf("REF* -> %s - %s\n", orig_ref, str_new_cigar);
							printed = 1;
						}
						cigar32_to_string(aux_cigar, unclip_cigar_l, str_new_cigar);
						printf("REF%d -> %s - %s ::: Retry - %d\n", j, aux_ref, str_new_cigar, retry);
					}*/
				}

				if((aux_cigar[indel_index - 1] >> BAM_CIGAR_SHIFT) == 0)
					retry = 0;
			}

			//Free
			free(orig_ref);
			free(aux_ref);
		}
	}
	else
	{
		//Return original CIGAR
		memcpy(out_cigar, cigar, cigar_l * sizeof(uint32_t));
		*out_cigar_l = cigar_l;
	}

	return NO_ERROR;
}

/**
 * Remove clips from cigar.
 */
ERROR_CODE
cigar32_unclip(uint32_t *cigar, size_t cigar_l, uint32_t *out_cigar, size_t *out_cigar_l)
{
	int i;
	int c_count;
	int c_type;
	int new_cigar_i;

	assert(cigar);
	assert(cigar_l > 0);
	assert(out_cigar);
	assert(out_cigar_l);

	//Iterate cigar elements
	new_cigar_i = 0;
	for(i = 0; i < cigar_l && new_cigar_i < MAX_CIGAR_LENGTH; i++)
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
				out_cigar[new_cigar_i] = cigar[i];
				new_cigar_i++;
			}
		}
	}

	//Set output cigar length
	*out_cigar_l = new_cigar_i;

	return NO_ERROR;
}

/**
 * Add clips to cigar32 from other cigar32 (ideally the original one).
 */
ERROR_CODE
cigar32_reclip(uint32_t *clip_cigar, size_t clip_cigar_l, uint32_t *unclip_cigar, size_t unclip_cigar_l, uint32_t *out_cigar, size_t *out_cigar_l)
{
	int i;
	int c_count;
	int c_type;
	size_t new_l = 0;
	uint32_t aux_cigar[MAX_CIGAR_LENGTH];

	assert(clip_cigar);
	assert(clip_cigar_l > 0);
	assert(unclip_cigar);
	assert(unclip_cigar_l > 0);
	assert(out_cigar);
	assert(out_cigar_l);

	//Get first clips
	i = 0;
	while(i < clip_cigar_l && new_l < MAX_CIGAR_LENGTH)
	{
		c_type = clip_cigar[i] & BAM_CIGAR_MASK;	//Get type from cigar

		//Is clip?
		if(c_type == BAM_CSOFT_CLIP
				|| c_type == BAM_CHARD_CLIP
				|| c_type == BAM_CPAD)
		{
			//Set clip in new cigar
			aux_cigar[i] = clip_cigar[i];
			new_l++;
		}
		else
		{
			break;
		}

		//Increment index
		i++;
	}

	//Avoid buffer overflow
	if(new_l + unclip_cigar_l < MAX_CIGAR_LENGTH)
	{
		//Copy unclipped cigar
		memcpy(aux_cigar + new_l, unclip_cigar, unclip_cigar_l * sizeof(uint32_t));
		//i += unclip_cigar_l;
		new_l += unclip_cigar_l;

		//Search for last clips
		while(i < clip_cigar_l)
		{
			c_type = clip_cigar[i] & BAM_CIGAR_MASK;
			if(c_type == BAM_CSOFT_CLIP
							|| c_type == BAM_CHARD_CLIP
							|| c_type == BAM_CPAD)
			{
				break;
			}

			i++;
		}

		//Copy last clips
		while(i < clip_cigar_l && new_l < MAX_CIGAR_LENGTH)
		{
			c_type = clip_cigar[i] & BAM_CIGAR_MASK;	//Get type from cigar

			//Is clip?
			if(c_type == BAM_CSOFT_CLIP
					|| c_type == BAM_CHARD_CLIP
					|| c_type == BAM_CPAD)
			{
				//Set clip in new cigar
				aux_cigar[new_l] = clip_cigar[i];
				new_l++;

				//Increment index
				i++;
			}
			else
			{
				break;
			}
		}
	} //Overflow if

	//Set output cigar
	memcpy(out_cigar, aux_cigar, sizeof(uint32_t) * new_l);
	*out_cigar_l = new_l;

	return NO_ERROR;
}

/**
 * Shift a cigar32 indel one position in left direction.
 */
ERROR_CODE
cigar32_shift_left_indel(uint32_t *cigar, size_t cigar_l, size_t indel_index, uint32_t *out_cigar, size_t *out_cigar_l)
{
	int i, elem, type;

	//Index
	size_t actual_index;

	assert(cigar);
	assert(cigar_l > 0);
	assert(indel_index > 0);
	assert(out_cigar);
	assert(out_cigar_l);

	//printf("Shifting left: index-%d, length-%d\n", indel_index, cigar_l);

	//Avoid overflow
	if(indel_index > MAX_CIGAR_LENGTH)
	{
		memcpy(out_cigar, cigar, MAX_CIGAR_LENGTH);
	}

	//Copy first part of CIGAR
	if(indel_index > 1)
	{
		//printf("Copying first part of cigar: %d elements\n",indel_index - 1);
		memcpy(out_cigar, cigar, actual_index - 1);
	}

	//Previous element
	elem = cigar[indel_index - 1] >> BAM_CIGAR_SHIFT;	//Get number of bases from cigar
	type = cigar[indel_index - 1] & BAM_CIGAR_MASK;	//Get type from cigar
	out_cigar[indel_index - 1] = ((elem - 1) << BAM_CIGAR_SHIFT) + type;
	//printf("Copying previous to indel: %d-%d\n",elem - 1, type);

	//Indel element
	out_cigar[indel_index] = cigar[indel_index];
	//printf("Copying indel: %d-%d\n",cigar[indel_index] >> BAM_CIGAR_SHIFT, cigar[indel_index] & BAM_CIGAR_MASK);

	if(indel_index + 1 < MAX_CIGAR_LENGTH)
	{
		//Next element
		actual_index = indel_index + 1;
		elem = cigar[actual_index] >> BAM_CIGAR_SHIFT;	//Get number of bases from cigar
		type = cigar[actual_index] & BAM_CIGAR_MASK;	//Get type from cigar
		out_cigar[actual_index] = ((elem + 1) << BAM_CIGAR_SHIFT) + type;
		//printf("Copying post to indel: %d-%d\n",elem + 1, type);

		//Copy last part of CIGAR
		if(actual_index + 1 < cigar_l)
		{
			//printf("Copying last part of cigar: %d elements\n",cigar_l - (indel_index + 1));
			memcpy(out_cigar, cigar + actual_index, min(cigar_l - actual_index, MAX_CIGAR_LENGTH - actual_index));
		}
	}

	//Set output length
	*out_cigar_l = min(cigar_l, MAX_CIGAR_LENGTH);

	return NO_ERROR;
}

ERROR_CODE
cigar32_count_m_blocks(uint32_t *cigar, size_t cigar_l, size_t *blocks)
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

	return NO_ERROR;
}

ERROR_CODE
cigar32_count_indels(uint32_t *cigar, size_t cigar_l, size_t *indels)
{
	int i;
	int c_count;
	int c_type;
	size_t indel_c;
	char *str_cigar;

	assert(cigar);
	//assert(cigar_l > 0);
	assert(indels);

	if(cigar_l == 0)
		return 0;

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

	return NO_ERROR;
}

ERROR_CODE
cigar32_count_nucleotides(uint32_t *cigar, size_t cigar_l, size_t *bases)
{
	int i;
	int c_count;
	int c_type;
	size_t bases_c;

	assert(cigar);
	assert(cigar_l > 0);
	assert(bases);

	//Iterate cigar elements
	bases_c = 0;
	for(i = 0; i < cigar_l; i++)
	{
		c_count = cigar[i] >> BAM_CIGAR_SHIFT;	//Get number of bases from cigar
		c_type = cigar[i] & BAM_CIGAR_MASK;	//Get type from cigar

		switch(c_type)
		{
		case BAM_CMATCH:
		case BAM_CEQUAL:
		case BAM_CDIFF:
		case BAM_CINS:
		case BAM_CSOFT_CLIP:
			//Count bases
			bases_c += c_count;
			break;

		case BAM_CHARD_CLIP:
		case BAM_CDEL:
		case BAM_CREF_SKIP:
		case BAM_CPAD:
			//No bases
			break;
		}
	}

	//Set output
	*bases = bases_c;

	return NO_ERROR;
}
ERROR_CODE
cigar32_count_nucleotides_not_clip(uint32_t *cigar, size_t cigar_l, size_t *bases)
{
	int i;
	int c_count;
	int c_type;
	size_t bases_c;

	assert(cigar);
	assert(cigar_l > 0);
	assert(bases);

	//Iterate cigar elements
	bases_c = 0;
	for(i = 0; i < cigar_l; i++)
	{
		c_count = cigar[i] >> BAM_CIGAR_SHIFT;	//Get number of bases from cigar
		c_type = cigar[i] & BAM_CIGAR_MASK;	//Get type from cigar

		switch(c_type)
		{
		case BAM_CMATCH:
		case BAM_CEQUAL:
		case BAM_CDIFF:
		case BAM_CINS:
			//Count bases
			bases_c += c_count;
			break;

		case BAM_CHARD_CLIP:
		case BAM_CDEL:
		case BAM_CREF_SKIP:
		case BAM_CPAD:
		case BAM_CSOFT_CLIP:
			//No bases
			break;
		}
	}

	//Set output
	*bases = bases_c;

	return NO_ERROR;
}

ERROR_CODE
cigar32_count_all(uint32_t *cigar, size_t cigar_l, size_t *m_blocks, size_t *indels, size_t *first_indel_index)
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

	return NO_ERROR;
}

ERROR_CODE
cigar32_count_clip_displacement(uint32_t *cigar, size_t cigar_l, size_t *out_disp)
{
	size_t disp = 0;

	assert(cigar);
	assert(cigar_l > 0);

	int c_count;
	int c_type;

	c_count = cigar[0] >> BAM_CIGAR_SHIFT;	//Get number of bases from cigar
	c_type = cigar[0] & BAM_CIGAR_MASK;	//Get type from cigar

	if(c_type == BAM_CSOFT_CLIP)
	{
		disp = c_count;
	}

	//Set clip displacement
	*out_disp = disp;

	return NO_ERROR;
}

ERROR_CODE
cigar32_to_string(uint32_t *cigar, size_t cigar_l, char* str_cigar)
{
	int i, elem, type;

	//sprintf(str_cigar, "\0");
	str_cigar = "\0";

	for(i = 0; i < cigar_l; i++)
	{
		elem = cigar[i] >> BAM_CIGAR_SHIFT;	//Get number of bases from cigar
		type = cigar[i] & BAM_CIGAR_MASK;	//Get type from cigar

		switch(type)
		{
		case BAM_CMATCH:
			sprintf(str_cigar + strlen(str_cigar),"%dM", elem);
			break;
		case BAM_CHARD_CLIP:
			sprintf(str_cigar + strlen(str_cigar),"%dH", elem);
			break;
		case BAM_CINS:
			sprintf(str_cigar + strlen(str_cigar),"%dI", elem);
			break;
		case BAM_CDEL:
			sprintf(str_cigar + strlen(str_cigar),"%dD", elem);
			break;
		case BAM_CREF_SKIP:
			sprintf(str_cigar + strlen(str_cigar),"%dN", elem);
			break;
		case BAM_CSOFT_CLIP:
			sprintf(str_cigar + strlen(str_cigar),"%dS", elem);
			break;
		case BAM_CPAD:
			sprintf(str_cigar + strlen(str_cigar),"%dP", elem);
			break;
		case BAM_CEQUAL:
			sprintf(str_cigar + strlen(str_cigar),"%d=", elem);
			break;
		case BAM_CDIFF:
			sprintf(str_cigar + strlen(str_cigar),"%dX", elem);
			break;
		}
	}

	return NO_ERROR;
}

ERROR_CODE
cigar32_create_ref(uint32_t *cigar, size_t cigar_l, char *ref, size_t ref_l, char *read, size_t read_l, char *new_ref, size_t *new_ref_l)
{
	int i, elem, type, extra;
	char *aux_str;

	//Index
	size_t index_ref = 0;
	size_t index_read = 0;
	size_t index_aux = 0;
	size_t remain = 0;
	extra = 0;

	assert(cigar);
	assert(cigar_l > 0);
	assert(ref);
	assert(ref_l > 0);
	assert(read);
	assert(read_l > 0);
	assert(new_ref);
	assert(new_ref_l);

	//Allocate output
	aux_str = (char *)malloc(sizeof(char) * (read_l + 1));

	//Iterate CIGAR
	for(i = 0; i < cigar_l; i++)
	{
		elem = cigar[i] >> BAM_CIGAR_SHIFT;	//Get number of bases from cigar
		type = cigar[i] & BAM_CIGAR_MASK;	//Get type from cigar

		switch(type)
		{
		case BAM_CSOFT_CLIP:
		case BAM_CINS:	//Insertion
			if(index_read + elem > read_l)
				elem = read_l - index_read;
			//Copy insertion in aux
			memcpy(aux_str + index_aux, read + index_read, elem);
			//printf("I:%d-%d-%d\n", index_read, index_ref, index_aux);
			index_read += elem;
			index_aux += elem;
			break;
		case BAM_CDEL:	//Deletion
			if(index_ref + elem > ref_l)
				elem = ref_l - index_ref;
			//Increment index on reference to skip deletion
			//printf("D:%d-%d-%d\n", index_read, index_ref, index_aux);
			index_ref += elem;
			extra += elem;
			break;
		case BAM_CMATCH:
		case BAM_CDIFF:
		case BAM_CEQUAL:
			if(index_ref + elem > ref_l)
				elem = ref_l - index_ref;
			if(index_read + elem > read_l)
				elem = read_l - index_read;
			//Copy reference as it is
			memcpy(aux_str + index_aux, ref + index_ref, elem);
			//printf("M:%d-%d-%d\n", index_read, index_ref, index_aux);
			index_ref += elem;
			index_read += elem;
			index_aux += elem;
			break;

		/*case BAM_CSOFT_CLIP:
			//Clips
			if(index_read + elem > read_l)
				elem = read_l - index_read;
			index_read += elem;
			break;*/

		case BAM_CPAD:
		case BAM_CHARD_CLIP:
			break;

		default:
			fprintf(stderr, "WARNING: Unrecognised cigar N:%d T:%d\n", elem, type);
			fflush(stderr);
			abort();
		}
	}

	//Copy last nucleotides
	//remain = length - index_read;
	//if(remain)
	//	memcpy(aux_str + index_aux, ref + index_ref, remain);

	//Set output
	if(index_aux > read_l)
		index_aux = read_l;
	aux_str[index_aux] = '\0';
	memcpy(new_ref, aux_str, sizeof(char) * index_aux + 1);

	//Set output length
	*new_ref_l = index_aux;

	//Free
	free(aux_str);

	return NO_ERROR;
}

ERROR_CODE
cigar32_get_indels(size_t ref_pos, uint32_t *cigar, size_t cigar_l, aux_indel_t *out_indels)
{
	int i, elem, type, indel_index;
	size_t current_pos;

	assert(ref_pos != SIZE_MAX);
	assert(cigar);
	assert(cigar_l > 0);
	assert(out_indels);

	//Iterate CIGAR
	indel_index = 0;
	current_pos = ref_pos;
	for(i = 0; i < cigar_l; i++)
	{
		elem = cigar[i] >> BAM_CIGAR_SHIFT;	//Get number of bases from cigar
		type = cigar[i] & BAM_CIGAR_MASK;	//Get type from cigar

		switch(type)
		{
		case BAM_CINS:	//Insertion
		case BAM_CDEL:	//Deletion
			//Create indel
			out_indels[indel_index].indel = cigar[i];
			out_indels[indel_index].ref_pos = current_pos;
			indel_index++;
			current_pos += elem;
			break;
		case BAM_CMATCH:
		case BAM_CDIFF:
		case BAM_CEQUAL:
			current_pos += elem;
			break;

		case BAM_CPAD:
		case BAM_CHARD_CLIP:
		case BAM_CSOFT_CLIP:
			break;

		default:
			fprintf(stderr, "WARNING: Unrecognised cigar N:%d T:%d\n", elem, type);
			fflush(stderr);
			abort();
		}
	}

	return NO_ERROR;
}

ERROR_CODE
cigar32_from_haplo(uint32_t *cigar, size_t cigar_l, aux_indel_t *haplo, size_t read_pos, uint32_t *new_cigar, size_t *new_cigar_l)
{
	int i, elem, type;
	size_t current_pos;
	size_t current_cigar_elem;
	size_t bases;

	//Read
	size_t aux_read_pos;
	int bases_left;

	//Indel
	int indel_size;
	int indel_type;
	size_t disp_to_indel;

	//Generated cigar
	uint32_t gen_cigar[MAX_CIGAR_LENGTH];
	size_t gen_cigar_l;

	assert(cigar);
	assert(cigar_l > 0);
	assert(haplo);
	assert(read_pos != SIZE_MAX);
	assert(new_cigar);
	assert(new_cigar_l);

	//Indel params
	indel_size = haplo->indel >> BAM_CIGAR_SHIFT;
	indel_type = haplo->indel & BAM_CIGAR_MASK;

	//Haplotype position must be posterior to read position
	if(haplo->ref_pos < read_pos - indel_size)
	{
		memcpy(new_cigar, cigar, cigar_l * sizeof(uint32_t));
		*new_cigar_l = cigar_l;
		return INVALID_INPUT_PARAMS;
	}

	//Get read displacement from haplotype
	cigar32_count_clip_displacement(cigar, cigar_l, &disp_to_indel);
	aux_read_pos = read_pos; //+ disp_to_indel;
	disp_to_indel = haplo->ref_pos - aux_read_pos;

	//Count unclipped bases
	cigar32_count_nucleotides_not_clip(cigar, cigar_l, &bases);

	//Indel type cases
	switch(indel_type)
	{
	//Indel is an insertion
	case BAM_CINS:
		//Where is the indel
		if(disp_to_indel > 0)
		{
			//Insertion is not in in read beginning
			bases_left = (int)bases - disp_to_indel;

			//Insertion inside read
			if(bases_left - indel_size > 0)
			{
				//Remaining bases including insertion
				bases_left -= indel_size;

				//Create cigar
				gen_cigar[0] = (disp_to_indel << BAM_CIGAR_SHIFT) + BAM_CMATCH;
				gen_cigar[1] = (indel_size << BAM_CIGAR_SHIFT) + indel_type;
				gen_cigar[2] = (bases_left << BAM_CIGAR_SHIFT) + BAM_CMATCH;
				gen_cigar_l = 3;
			}
			//Insertion at read end
			else
			{
				//Set indel size
				indel_size = bases_left;

				//If
				if(indel_size > 0)
				{
					//Last element size != 0
					gen_cigar[0] = ((bases - indel_size) << BAM_CIGAR_SHIFT) + BAM_CMATCH;
					gen_cigar[1] = (indel_size << BAM_CIGAR_SHIFT) + (haplo->indel & BAM_CIGAR_MASK);
					gen_cigar_l = 2;
				}
				else
				{
					//Last element size is 0
					gen_cigar[0] = (bases << BAM_CIGAR_SHIFT) + BAM_CMATCH;
					gen_cigar_l = 1;
				}
			}
		}
		else
		{
			//Insertion is in read beginning

			//In this case, negative displacement (read_pos > indel_pos) mean a shorter insertion
			indel_size += disp_to_indel;
			if(indel_size > 0)
			{
				//Set length of final cigar element
				bases -= indel_size;

				//Create cigar
				gen_cigar[0] = (indel_size << BAM_CIGAR_SHIFT) + indel_type;
				gen_cigar[1] = (bases << BAM_CIGAR_SHIFT) + BAM_CMATCH;
				gen_cigar_l = 2;
			}
			else
			{
				//Read is too far from insertion
				gen_cigar[0] = (bases << BAM_CIGAR_SHIFT) + BAM_CMATCH;
				gen_cigar_l = 1;
			}
		}

		break;

	//Indel is a deletion
	case BAM_CDEL:

		//Where is the indel
		if(disp_to_indel > 0)
		{
			//Deletion is not in in read beginning
			bases_left = (int)bases - disp_to_indel;

			//Deletion inside read
			if(bases_left > 0)
			{
				//Create cigar
				gen_cigar[0] = (disp_to_indel << BAM_CIGAR_SHIFT) + BAM_CMATCH;
				gen_cigar[1] = (indel_size << BAM_CIGAR_SHIFT) + indel_type;
				gen_cigar[2] = (bases_left << BAM_CIGAR_SHIFT) + BAM_CMATCH;
				gen_cigar_l = 3;
			}
			//Deletion at read end
			else
			{
				//Deletion is not useful at read end
				gen_cigar[0] = (bases << BAM_CIGAR_SHIFT) + BAM_CMATCH;
				gen_cigar_l = 1;
			}
		}
		else
		{
			//Deletion is in read beginning
			//Deletion is not useful at read beginning
			gen_cigar[0] = (bases << BAM_CIGAR_SHIFT) + BAM_CMATCH;
			gen_cigar_l = 1;
		}

		break;

	default:
		//Unrecognised indel????
		return CIGAR_INVALID_INDEL;
	}

	//Set output
	cigar32_reclip(cigar, cigar_l, gen_cigar, gen_cigar_l, new_cigar, new_cigar_l);

	return NO_ERROR;
}

ERROR_CODE
cigar32_replace(bam1_t *read, uint32_t *cigar, size_t cigar_l)
{
	size_t new_data_l;
	size_t new_data_ml;
	uint8_t *new_data;
	size_t n_offset;
	size_t r_offset;

	assert(read);
	assert(cigar);

	//Do nothing
	if(cigar_l == 0 || cigar_l >= MAX_CIGAR_LENGTH)
	{
		return NO_ERROR;
	}

	//Check if same length cigar
	if(cigar_l == read->core.n_cigar)
	{
		//Same length so memcpy
		memcpy(bam1_cigar(read), cigar, read->core.n_cigar * sizeof(uint32_t));
	}
	else
	{
		//Get new length for bam data
		new_data_l = (read->data_len - (read->core.n_cigar * sizeof(uint32_t))) + (cigar_l * sizeof(uint32_t));

		//Get allocate size
		if (read->m_data < new_data_l) {
			new_data_ml = new_data_l;
			kroundup32(new_data_ml);
		}
		else
		{
			new_data_ml = read->m_data;
		}

		//Different length so reallocate
		new_data = (uint8_t *)malloc(new_data_ml * sizeof(uint8_t));

		//Copy contents
		memcpy(new_data, bam1_qname(read), read->core.l_qname);
		r_offset = read->core.l_qname;
		n_offset = read->core.l_qname;
		memcpy(new_data + n_offset, cigar, cigar_l * sizeof(uint32_t));
		r_offset += read->core.n_cigar * sizeof(uint32_t);
		n_offset += cigar_l * sizeof(uint32_t);
		memcpy(new_data + n_offset, bam1_seq(read), read->data_len - r_offset);

		//Replace data
		free(read->data);
		read->data = new_data;
		read->data_len = new_data_l;
		read->m_data = new_data_ml;
		read->core.n_cigar = cigar_l;
	}

	return NO_ERROR;
}

