/*
 * aux_cigar.h
 *
 *  Created on: Dec 16, 2013
 *      Author: rmoreno
 */

#ifndef AUX_CIGAR_H_
#define AUX_CIGAR_H_

#include "aux_library.h"

#define MAX_CIGAR_LENGTH 64

/***************************
 * CIGAR OPERATIONS
 **************************/

typedef struct {
	uint32_t indel;
	size_t ref_pos;
} aux_indel_t;

/**
 * CIGAR GENERATION
 */

EXTERNC ERROR_CODE cigar32_leftmost(char *ref, size_t ref_l, char *read, size_t read_l, uint32_t *cigar, size_t cigar_l, uint32_t *out_cigar, size_t *out_cigar_l);
EXTERNC ERROR_CODE cigar32_unclip(uint32_t *cigar, size_t cigar_l, uint32_t *out_cigar, size_t *out_cigar_l);
EXTERNC ERROR_CODE cigar32_reclip(uint32_t *clip_cigar, size_t clip_cigar_l, uint32_t *unclip_cigar, size_t unclip_cigar_l, uint32_t *out_cigar, size_t *out_cigar_l);
EXTERNC ERROR_CODE cigar32_shift_left_indel(uint32_t *cigar, size_t cigar_l, size_t indel_index, uint32_t *out_cigar, size_t *out_cigar_l);

/**
 * CIGAR COUNT
 */

EXTERNC ERROR_CODE cigar32_count_m_blocks(uint32_t *cigar, size_t cigar_l, size_t *blocks);
EXTERNC ERROR_CODE cigar32_count_indels(uint32_t *cigar, size_t cigar_l, size_t *indels);
EXTERNC ERROR_CODE cigar32_count_nucleotides(uint32_t *cigar, size_t cigar_l, size_t *bases);
EXTERNC ERROR_CODE cigar32_count_nucleotides_not_clip(uint32_t *cigar, size_t cigar_l, size_t *bases);
EXTERNC ERROR_CODE cigar32_count_all(uint32_t *cigar, size_t cigar_l, size_t *m_blocks, size_t *indels, size_t *first_indel_index);
EXTERNC ERROR_CODE cigar32_count_clip_displacement(uint32_t *cigar, size_t cigar_l, size_t *out_disp);

/**
 * CIGAR FUNCTIONS
 */
EXTERNC ERROR_CODE cigar32_to_string(uint32_t *cigar, size_t cigar_l, char* str_cigar);
EXTERNC ERROR_CODE cigar32_create_ref(uint32_t *cigar, size_t cigar_l, char *ref, size_t ref_l, char *read, size_t read_l, char *new_ref, size_t *new_ref_l);
EXTERNC ERROR_CODE cigar32_get_indels(size_t ref_pos, uint32_t *cigar, size_t cigar_l, aux_indel_t *out_indels);
EXTERNC ERROR_CODE cigar32_from_haplo(uint32_t *cigar, size_t cigar_l, aux_indel_t *haplo, size_t read_pos, uint32_t *new_cigar, size_t *new_cigar_l);
EXTERNC ERROR_CODE cigar32_replace(bam1_t *read, uint32_t *cigar, size_t cigar_l);

#endif /* AUX_CIGAR_H_ */
