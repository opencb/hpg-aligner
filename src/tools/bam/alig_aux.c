/*
 * alig_aux.c
 *
 *  Created on: Dec 12, 2013
 *      Author: rmoreno
 */

#include "alig_aux.h"


int
alig_aux_cigar32_to_string(uint32_t *cigar, size_t cigar_l, char* str_cigar)
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
}

int
alig_aux_cigar32_create_ref(uint32_t *cigar, size_t cigar_l, char *ref, char *read, size_t length, char *new_ref)
{
	int i, elem, type;
	char *aux_str;

	//Index
	size_t index_ref = 0;
	size_t index_read = 0;
	size_t index_aux = 0;
	size_t remain = 0;

	assert(cigar);
	assert(cigar_l > 0);
	assert(ref);
	assert(read);
	assert(length);
	assert(new_ref);

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
			break;
		case BAM_CMATCH:
		case BAM_CHARD_CLIP:
		case BAM_CSOFT_CLIP:
		case BAM_CPAD:
		case BAM_CEQUAL:
			//Copy reference as it is
			memcpy(aux_str + index_aux, ref + index_ref, elem);
			//printf("M:%d-%d-%d\n", index_read, index_ref, index_aux);
			index_ref += elem;
			index_read += elem;
			index_aux += elem;
			break;

		default:
			fprintf(stderr, "WARNING: Unrecognised cigar N:%d T:%d\n", elem, type);
			fflush(stderr);
			abort();
		}
	}

	//Copy last nucleotides
	remain = length - index_aux;
	if(remain)
		memcpy(aux_str + index_aux, ref + index_ref, remain);

	//Set output
	memcpy(new_ref, aux_str, sizeof(char) * length);

	return 0;
}

int
alig_aux_cigar32_shift_left_indel(uint32_t *cigar, size_t cigar_l, size_t indel_index, uint32_t *new_cigar)
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

	return 0;
}

