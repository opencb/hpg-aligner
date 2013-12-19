/*
 * aux_cigar.c
 *
 *  Created on: Dec 16, 2013
 *      Author: rmoreno
 */

#include "aux_cigar.h"

ERROR_CODE
cigar32_leftmost(char *ref, char *read, size_t read_l, uint32_t *cigar, size_t cigar_l, uint32_t *new_cigar, size_t *new_cigar_l)
{
	//M blocks
	size_t blocks_c;

	//Indels
	size_t indels_c;
	size_t indel_index;
	size_t indel_l;

	//New cigar
	uint32_t unclip_cigar[30];
	uint32_t reclip_cigar[30];
	uint32_t aux_cigar[30];
	size_t unclip_cigar_l;
	size_t reclip_cigar_l;

	//References
	char *orig_ref;
	size_t orig_ref_l;
	char *aux_ref;
	size_t aux_ref_l;

	//Counters
	int retry;

	//ERASE
	char str_cigar[200];
	char str_new_cigar[200];
	int printed = 0;

	assert(ref);
	assert(read);
	assert(read_l > 0);
	assert(cigar);
	//assert(cigar_l < 30);
	assert(new_cigar);
	assert(new_cigar_l);

	//Unclip cigar
	cigar32_unclip(cigar, cigar_l, unclip_cigar, &unclip_cigar_l);

	//Count blocks and indels
	cigar32_count_all(unclip_cigar, unclip_cigar_l, &blocks_c, &indels_c, &indel_index);

	//Get indel length
	indel_l = unclip_cigar[indel_index] >> BAM_CIGAR_SHIFT;

	//By default return original CIGAR
	memcpy(new_cigar, cigar, cigar_l * sizeof(uint32_t));
	*new_cigar_l = cigar_l;

	//Only procceed if 1 indel (2 M blocks)
	if(cigar_l < 30 && blocks_c == 2 && indels_c == 1 && indel_l)
	{
		//Leftmost CIGAR
		{
			//Get reference for original CIGAR
			orig_ref = (char *)malloc(sizeof(char) * (read_l + 1));
			aux_ref = (char *)malloc(sizeof(char) * (read_l + 1));
			cigar32_create_ref(cigar, cigar_l, ref, read, read_l, orig_ref, &orig_ref_l);

			//Shift left CIGAR
			cigar32_shift_left_indel(unclip_cigar, unclip_cigar_l, indel_index, aux_cigar);

			//Reclip cigar
			cigar32_reclip(cigar, cigar_l, aux_cigar, unclip_cigar_l, reclip_cigar, &reclip_cigar_l);

			//Get new CIGAR ref
			cigar32_create_ref(reclip_cigar, reclip_cigar_l, ref, read, read_l, aux_ref, &aux_ref_l);

			//Is a valid ref?
			if(!memcmp(aux_ref, orig_ref, orig_ref_l * sizeof(char)))
			{
				//Equal so is a valid CIGAR
				memcpy(new_cigar, aux_cigar, sizeof(uint32_t) * reclip_cigar_l);
				*new_cigar_l = reclip_cigar_l;
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
				cigar32_shift_left_indel(aux_cigar, unclip_cigar_l, indel_index, aux_cigar);

				//Reclip cigar
				cigar32_reclip(cigar, cigar_l, aux_cigar, unclip_cigar_l, reclip_cigar, &reclip_cigar_l);

				//Get new CIGAR ref
				cigar32_create_ref(reclip_cigar, reclip_cigar_l, ref, read, read_l, aux_ref, &aux_ref_l);

				//Is a valid ref?
				if(!memcmp(aux_ref, orig_ref, orig_ref_l * sizeof(char)))
				{
					//Equal so is a valid CIGAR
					memcpy(new_cigar, reclip_cigar, sizeof(uint32_t) * reclip_cigar_l);
					*new_cigar_l = reclip_cigar_l;
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
		memcpy(new_cigar, cigar, cigar_l * sizeof(uint32_t));
		*new_cigar_l = cigar_l;
	}

	return NO_ERROR;
}


ERROR_CODE
cigar32_unclip(uint32_t *cigar, size_t cigar_l, uint32_t *new_cigar, size_t *new_cigar_l)
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

	return NO_ERROR;
}

ERROR_CODE
cigar32_reclip(uint32_t *clip_cigar, size_t clip_cigar_l, uint32_t *unclip_cigar, size_t unclip_cigar_l, uint32_t *new_cigar, size_t *new_cigar_l)
{
	int i;
	int c_count;
	int c_type;
	size_t new_l = 0;

	assert(clip_cigar);
	assert(clip_cigar_l > 0);
	assert(unclip_cigar);
	assert(unclip_cigar_l > 0);
	assert(new_cigar);
	assert(new_cigar_l);

	//Get first clips
	i = 0;
	while(i < clip_cigar_l)
	{
		c_type = clip_cigar[i] & BAM_CIGAR_MASK;	//Get type from cigar

		//Is clip?
		if(c_type == BAM_CSOFT_CLIP
				|| c_type == BAM_CHARD_CLIP
				|| c_type == BAM_CPAD)
		{
			//Set clip in new cigar
			new_cigar[i] = clip_cigar[i];
			new_l++;
		}
		else
		{
			break;
		}

		//Increment index
		i++;
	}

	//Copy unclipped cigar
	memcpy(new_cigar + i, unclip_cigar, unclip_cigar_l * sizeof(uint32_t));
	i += unclip_cigar_l;
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
	while(i < clip_cigar_l)
	{
		c_type = clip_cigar[i] & BAM_CIGAR_MASK;	//Get type from cigar

		//Is clip?
		if(c_type == BAM_CSOFT_CLIP
				|| c_type == BAM_CHARD_CLIP
				|| c_type == BAM_CPAD)
		{
			//Set clip in new cigar
			new_cigar[new_l] = clip_cigar[i];
			new_l++;

			//Increment index
			i++;
		}
		else
		{
			break;
		}
	}

	//Set output cigar length
	*new_cigar_l = new_l;
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
cigar32_to_string(uint32_t *cigar, size_t cigar_l, char* str_cigar)
{
	int i, elem, type;

	sprintf(str_cigar, "");
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
cigar32_create_ref(uint32_t *cigar, size_t cigar_l, char *ref, char *read, size_t length, char *new_ref, size_t *new_ref_l)
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
	assert(read);
	assert(length);
	assert(new_ref);
	assert(new_ref_l);

	//Allocate output
	aux_str = (char *)malloc(sizeof(char) * length);

	//Iterate CIGAR
	for(i = 0; i < cigar_l; i++)
	{
		elem = cigar[i] >> BAM_CIGAR_SHIFT;	//Get number of bases from cigar
		type = cigar[i] & BAM_CIGAR_MASK;	//Get type from cigar

		switch(type)
		{
		case BAM_CINS:	//Insertion
			//Copy insertion in aux
			memcpy(aux_str + index_aux, read + index_read, elem);
			//printf("I:%d-%d-%d\n", index_read, index_ref, index_aux);
			index_read += elem;
			index_aux += elem;
			break;
		case BAM_CDEL:	//Deletion
			//Increment index on reference to skip deletion
			//printf("D:%d-%d-%d\n", index_read, index_ref, index_aux);
			index_ref += elem;
			extra += elem;
			break;
		case BAM_CMATCH:
		case BAM_CDIFF:
		case BAM_CEQUAL:
			//Copy reference as it is
			memcpy(aux_str + index_aux, ref + index_ref, elem);
			//printf("M:%d-%d-%d\n", index_read, index_ref, index_aux);
			index_ref += elem;
			index_read += elem;
			index_aux += elem;
			break;

		case BAM_CSOFT_CLIP:
			//Clips
			index_read += elem;
			break;

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
	if(index_aux > length)
		index_aux = length;
	aux_str[index_aux] = '\0';
	memcpy(new_ref, aux_str, sizeof(char) * index_aux + 1);

	//Set output length
	*new_ref_l = index_aux;

	//Free
	free(aux_str);

	return NO_ERROR;
}

ERROR_CODE
cigar32_shift_left_indel(uint32_t *cigar, size_t cigar_l, size_t indel_index, uint32_t *new_cigar)
{
	int i, elem, type;

	//Index
	size_t actual_index;

	assert(cigar);
	assert(cigar_l > 0);
	assert(indel_index > 0);
	assert(new_cigar);

	//printf("Shifting left: index-%d, length-%d\n", indel_index, cigar_l);

	//Copy first part of CIGAR
	if(indel_index > 1)
	{
		//printf("Copying first part of cigar: %d elements\n",indel_index - 1);
		memcpy(new_cigar, cigar, indel_index - 1);
	}

	//Previous element
	elem = cigar[indel_index - 1] >> BAM_CIGAR_SHIFT;	//Get number of bases from cigar
	type = cigar[indel_index - 1] & BAM_CIGAR_MASK;	//Get type from cigar
	new_cigar[indel_index - 1] = ((elem - 1) << BAM_CIGAR_SHIFT) + type;
	//printf("Copying previous to indel: %d-%d\n",elem - 1, type);

	//Indel element
	new_cigar[indel_index] = cigar[indel_index];
	//printf("Copying indel: %d-%d\n",cigar[indel_index] >> BAM_CIGAR_SHIFT, cigar[indel_index] & BAM_CIGAR_MASK);

	//Next element
	elem = cigar[indel_index + 1] >> BAM_CIGAR_SHIFT;	//Get number of bases from cigar
	type = cigar[indel_index + 1] & BAM_CIGAR_MASK;	//Get type from cigar
	new_cigar[indel_index + 1] = ((elem + 1) << BAM_CIGAR_SHIFT) + type;
	//printf("Copying post to indel: %d-%d\n",elem + 1, type);

	//Copy last part of CIGAR
	if(indel_index + 2 < cigar_l)
	{
		//printf("Copying last part of cigar: %d elements\n",cigar_l - (indel_index + 1));
		memcpy(new_cigar, cigar + indel_index + 1, cigar_l - (indel_index + 1));
	}

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

