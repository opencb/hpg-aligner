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
