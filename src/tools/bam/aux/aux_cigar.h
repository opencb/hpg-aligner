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

/**
 * Count number of M nucleotides blocks. For example, cigar '23M4I56M' have 2 Mblocks.
 * \param[in] cigar Cigar32 to M block count.
 * \param[in] cigar_l Length of input cigar32.
 * \param[out] blocks Number of M blocks.
 */
EXTERNC ERROR_CODE cigar32_count_m_blocks(uint32_t *cigar, size_t cigar_l, size_t *blocks);

/**
 * Count number of indels present on cigar32. 
 * \param[in] cigar Count indels in this cigar32.
 * \param[in] cigar_l Length of input cigar32.
 * \param[out] indels Output indels count. 
 */
EXTERNC ERROR_CODE cigar32_count_indels(uint32_t *cigar, size_t cigar_l, size_t *indels);

/**
 * Count number of bases must have a read from one cigar32. 
 * \param[in] cigar Cigar32 to count bases. 
 * \param[in] cigar_l Length of input cigar32.
 * \param[out] bases Number of bases. 
 */
EXTERNC ERROR_CODE cigar32_count_nucleotides(uint32_t *cigar, size_t cigar_l, size_t *bases);

/**
 * Count number of unclipped bases must have a read from one cigar32. 
 * \param[in] cigar Cigar32 to count bases. 
 * \param[in] cigar_l Length of input cigar32.
 * \param[out] bases Number of unclipped bases. 
 */
EXTERNC ERROR_CODE cigar32_count_nucleotides_not_clip(uint32_t *cigar, size_t cigar_l, size_t *bases);

/**
 * Count M blocks, number of indels and get first indel index from cigar32. 
 * \param[in] cigar Cigar32 to process. 
 * \param[in] cigar_l Length of input cigar32.
 * \param[out] m_blocks Number of M blocks. Optional, pass NULL in case. 
 * \param[out] indels Number of indels. Options, pass NULL in case. 
 * \param[out] first_indel_index Index of first indel (If no indel, output is 0). Optional, pass NULL in case. 
 */
EXTERNC ERROR_CODE cigar32_count_all(uint32_t *cigar, size_t cigar_l, size_t *m_blocks, size_t *indels, size_t *first_indel_index);

/**
 * Count initial clip displacement. For example, '23S34M' have a clip displacement of 23.
 * \param[in] cigar Cigar32 to count displacement. 
 * \param[in] cigar_l Length of input cigar32.
 * \param[out] out_disp Initial clip displacement. 
 */
EXTERNC ERROR_CODE cigar32_count_clip_displacement(uint32_t *cigar, size_t cigar_l, size_t *out_disp);

/**
 * CIGAR FUNCTIONS
 */

/**
 * Convert cigar32 to char string. 
 * \param[in] cigar Cigar32 to convert. 
 * \param[in] cigar_l Length of input cigar32.
 * \param[out] str_cigar Preallocated output string, size will be same as cigar_l.
 */
EXTERNC ERROR_CODE cigar32_to_string(uint32_t *cigar, size_t cigar_l, char* str_cigar);

/**
 * Create a processed read sequence that can be compared directly to reference. This function is usefull to obtain missmatches between read sequences and reference ones. 
 * \param[in] cigar Cigar32 of input read sequence. 
 * \param[in] cigar_l Length of input cigar32.
 * \param[in] ref Reference sequence. 
 * \param[in] ref_l Length of reference sequence. 
 * \param[in] read Read sequence. 
 * \param[in] read_l Length of read sequence. 
 * \param[out] new_ref Output processed sequence. 
 * \param[out] new_ref_l Output processed sequence length. 
 */
EXTERNC ERROR_CODE cigar32_create_ref(uint32_t *cigar, size_t cigar_l, char *ref, size_t ref_l, char *read, size_t read_l, char *new_ref, size_t *new_ref_l, char *mask);

/**
 * Get an array of aux_indel_t representing indels present on a cigar32.  
 * \param[in] ref_pos Cigar32 position on reference sequence. 
 * \param[in] cigar Cigar32 to obtain indels. 
 * \param[in] cigar_l Length of input cigar32.
 * \param[out] out_indels Preallocated aux_indel_t array to store indels. WARNING: Size must be at least equal to indel count in cigar. 
 */
EXTERNC ERROR_CODE cigar32_get_indels(size_t ref_pos, uint32_t *cigar, size_t cigar_l, aux_indel_t *out_indels);

/**
 * Create cigar32 from haplotype store as aux_indel_t. 
 * \param[in] cigar Cigar32 in which haplotype will be inserted.
 * \param[in] cigar_l Length of input cigar32.
 * \param[in] read_pos Read sequence position in reference. 
 * \param[out] new_cigar Output cigar32 with haplotype inserted. 
 * \param[out] new_cigar_l Output cigar32 length. 
 */
EXTERNC ERROR_CODE cigar32_from_haplo(uint32_t *cigar, size_t cigar_l, aux_indel_t *haplo, size_t read_pos, uint32_t *new_cigar, size_t *new_cigar_l);

/**
 * Replace cigar32 with new one. 
 * \param[in] read Target read in bam1_t format to reemplace cigar. 
 * \param[in] cigar New cigar32. 
 * \param[in] cigar_l Length of new cigar32.
 */
EXTERNC ERROR_CODE cigar32_replace(bam1_t *read, uint32_t *cigar, size_t cigar_l);

#endif /* AUX_CIGAR_H_ */
