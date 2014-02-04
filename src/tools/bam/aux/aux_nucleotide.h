#ifndef AUX_NUCLEOTIDE_H_
#define AUX_NUCLEOTIDE_H_

#include "aux_library.h"

/***************************
 * NUCLEOTIDE OPERATIONS
 **************************/

EXTERNC ERROR_CODE nucleotide_compare(char *ref_seq, char *bam_seq, size_t bam_seq_l, char *comp_seq, uint32_t *miss_count);
EXTERNC ERROR_CODE nucleotide_miss_qual_sum(char *ref_seq, char *bam_seq, char *bam_qual, size_t bam_seq_l, char *comp_seq, uint32_t *out_miss_count, uint32_t *out_sum_quals);


#endif /* AUX_NUCLEOTIDE_H_ */
