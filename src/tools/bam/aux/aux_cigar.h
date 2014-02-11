/**
* Copyright (C) 2013 Raúl Moreno Galdón
*
* This program is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public License
* along with this program.  If not, see <http://www.gnu.org/licenses/>.
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

/**
 * Leftalign one read cigar.
 * \param[in] ref Reference sequence to match.
 * \param[in] ref_l Reference sequence length.
 * \param[in] read Read sequence.
 * \param[in] read_l Read sequence length.
 * \param[in] cigar 32bits cigar from read.
 * \param[in] cigar_l Length of 32bits cigar.
 * \param[out] out_cigar Output 32bits leftaligned cigar.
 * \param[out] out_cigar_l Output 32bits leftaligned cigar length.
 */
EXTERNC ERROR_CODE cigar32_leftmost(char *ref, size_t ref_l, char *read, size_t read_l, uint32_t *cigar, size_t cigar_l, uint32_t *out_cigar, size_t *out_cigar_l);

/**
 * Remove clips from cigar.
 * \param[in] cigar Cigar32 to unclip.
 * \param[in] cigar_l Lenght of cigar32 to unclip.
 * \param[out] out_cigar Unclipped cigar32.
 * \param[out] out_cigar_l Unclipped cigar32 length.
 */
EXTERNC ERROR_CODE cigar32_unclip(uint32_t *cigar, size_t cigar_l, uint32_t *out_cigar, size_t *out_cigar_l);

/**
 * Add clips to cigar32 from other cigar32 (ideally the original one).
 * \param[in] clip_cigar Cigar32 with clips.
 * \param[in] clip_cigar_l Cigar32 with clips length.
 * \param[in] unclip_cigar Unclipped cigar32.
 * \param[in] unclip_cigar_l Unclipped cigar32 length.
 * \param[out] out_cigar Output clipped cigar32.
 * \param[out] out_cigar_l Output clipped cigar32 length.
 */
EXTERNC ERROR_CODE cigar32_reclip(uint32_t *clip_cigar, size_t clip_cigar_l, uint32_t *unclip_cigar, size_t unclip_cigar_l, uint32_t *out_cigar, size_t *out_cigar_l);

/**
 * Shift a cigar32 indel one position in left direction.
 * \param[in] cigar Cigar32 to shift.
 * \param[in] cigar_l Length of cigar32 to shift.
 * \param[in] indel_index Indel index position to shift in cigar32.
 * \param[out] out_cigar Output shifted cigar32.
 * \param[out] out_cigar_l Output shifted cigar32 length.
 */
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
