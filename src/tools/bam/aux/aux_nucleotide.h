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

#ifndef AUX_NUCLEOTIDE_H_
#define AUX_NUCLEOTIDE_H_

#include "aux_library.h"

/***************************
 * NUCLEOTIDE OPERATIONS
 **************************/

/**
 * Compare two sequences and obtain missmatches. 
 * \param[in] ref_seq First sequence.
 * \param[in] bam_seq Second sequence. 
 * \param[in] bam_seq_l Length of input sequences. 
 * \param[out] comp_seq Output sequence to store individual missmatch value. 0xFF means equal, 0x00 means differ. 
 * \param[out] miss_count Number of missmatches in comparation. Optional, set to NULL in case.
 */
EXTERNC ERROR_CODE nucleotide_compare(char *ref_seq, char *bam_seq, size_t bam_seq_l, char *comp_seq, uint32_t *miss_count);

/**
 * Compare two sequences and obtain missmatches and missmatches qualities summatory. 
 * \param[in] ref_seq First sequence.
 * \param[in] bam_seq Second sequence. 
 * \param[in] bam_qual Sequence qualities. 
 * \param[in] bam_seq_l Length of input sequences. 
 * \param[out] comp_seq Output sequence to store individual missmatch value. 0xFF means equal, 0x00 means differ. 
 * \param[out] out_miss_count Number of missmatches in comparation. Optional, set to NULL in case.
 * \param[out] out_sum_quals Summatory of missmatch qualities. 
 */
EXTERNC ERROR_CODE nucleotide_miss_qual_sum(char *ref_seq, char *bam_seq, char *bam_qual, size_t bam_seq_l, char *comp_seq, uint32_t *out_miss_count, uint32_t *out_sum_quals);


#endif /* AUX_NUCLEOTIDE_H_ */
