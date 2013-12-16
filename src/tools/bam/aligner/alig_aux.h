/*
 * alig_aux.h
 *
 *  Created on: Dec 12, 2013
 *      Author: rmoreno
 */

#ifndef ALIG_AUX_H_
#define ALIG_AUX_H_

#include <assert.h>

#include <stdlib.h>
#include <stdint.h>

#include "bioformats/bam/samtools/bam.h"

int alig_aux_cigar32_to_string(uint32_t *cigar, size_t cigar_l, char* str_cigar);

int alig_aux_cigar32_create_ref(uint32_t *cigar, size_t cigar_l, char *ref, char *read, size_t length, char *new_ref);

int alig_aux_cigar32_shift_left_indel(uint32_t *cigar, size_t cigar_l, size_t indel_index, uint32_t *new_cigar);

#endif /* ALIG_AUX_H_ */
