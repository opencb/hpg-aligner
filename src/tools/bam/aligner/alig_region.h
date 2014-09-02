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

#ifndef ALIG_REGION_H_
#define ALIG_REGION_H_

#include <stdlib.h>
#include <stdint.h>
#include <limits.h>

#include "assert.h"
#include "aux/aux_library.h"

#include "containers/khash.h"
#include "containers/linked_list.h"
#include "bioformats/bam/samtools/bam.h"

/**
 * Structure representing genome region.
 */
typedef struct {
	size_t start_pos;
	size_t end_pos;
	int32_t chrom;
	uint8_t valid;
} alig_region_t;

/**
 * Hash table to store regions.
 */
KHASH_MAP_INIT_INT(32, linked_list_t*);
typedef struct {
	//Hash to store chrom regions alig_region_t
	khash_t(32) *chroms;
	//Others?
} alig_region_table_t;

/**
 * REGION STRUCT FUNCTIONS
 */

int region_create(alig_region_t *region);
int region_destroy(alig_region_t *region);

/**
 * REGION TABLE BY CHROM
 */

int region_table_create(alig_region_table_t *table);
int region_table_destroy(alig_region_table_t *table);

int region_table_insert(alig_region_table_t *table, alig_region_t *region);

/**
 * REGION DISCOVER
 */

/**
 * \brief Get a region from string cigar.
 *
 * \param[in] cigar String representation of input cigar.
 * \param[in] cigar_len Number of input cigar elements.
 * \param[in] pos Read position in reference genome.
 * \param[out] r_pos Region start position in reference genome.
 * \param[out] r_end_pos Region end position in reference genome.
 */
int region_get_from_cigar(char *cigar, size_t cigar_len, size_t pos, size_t *r_pos, size_t *r_end_pos);

/**
 * \brief Get a region from bam1_t.
 *
 * \param[in] alig Input BAM.
 * \param[out] r_pos Region start position in reference genome.
 * \param[out] r_end_pos Region end position in reference genome.
 */
int region_get_from_bam1(const bam1_t *alig, size_t *r_pos, size_t *r_end_pos);

/**
 * \brief Get a region from cigar32.
 *
 * \param[in] cigar Input cigar32.
 * \param[in] cigar_l Number of input cigar elements.
 * \param[in] pos Read position in reference genome.
 * \param[out] r_pos Region start position in reference genome.
 * \param[out] r_end_pos Region end position in reference genome.
 */
int region_get(uint32_t *cigar, uint32_t cigar_l, size_t pos, size_t *r_pos, size_t *r_end_pos);

/**
 * \brief Get region table from BAM batch.
 *
 * \param[in] batch Input BAM batch.
 * \param[in/out] region_table Table to store regions.
 */
int region_get_from_batch(const bam_batch_t* batch, alig_region_table_t *region_table);

/**
 * \brief Get region table from BAM file.
 *
 * \param[in] bam_path Path to input BAM file.
 */
int region_get_from_file(char *bam_path);

/**
 * \brief Get if BAM read overlaps a region.
 *
 * \param[in] read BAM read to check.
 * \param[in] region Region to check.
 * \return 1 if overlap, 0 if not.
 */
int region_bam_overlap(bam1_t *read, alig_region_t *region);

#endif /* ALIG_REGION_H_ */
