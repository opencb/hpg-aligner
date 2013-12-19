#ifndef AUX_NUCLEOTIDE_H_
#define AUX_NUCLEOTIDE_H_

#include "aux_library.h"

/***************************
 * NUCLEOTIDE OPERATIONS
 **************************/

EXTERNC ERROR_CODE nucleotide_compare(char *ref_seq, char *bam_seq, size_t bam_seq_l, char *comp_res, size_t *miss_count);


#endif /* AUX_NUCLEOTIDE_H_ */
