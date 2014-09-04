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
static inline int region_get_from_bam1(const bam1_t *alig, size_t *r_pos, size_t *r_end_pos) __ATTR_HOT __ATTR_INLINE;

/**
 * \brief Get a region from cigar32.
 *
 * \param[in] cigar Input cigar32.
 * \param[in] cigar_l Number of input cigar elements.
 * \param[in] pos Read position in reference genome.
 * \param[out] r_pos Region start position in reference genome.
 * \param[out] r_end_pos Region end position in reference genome.
 */
static inline int region_get(uint32_t *cigar, uint32_t cigar_l, size_t pos, size_t *r_pos, size_t *r_end_pos) __ATTR_HOT __ATTR_INLINE;

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
static inline int region_bam_overlap(bam1_t *read, alig_region_t *region) __ATTR_HOT __ATTR_INLINE;


/**
 * Get a region from cigar32.
 */
static inline int
region_get(uint32_t *cigar, uint32_t cigar_l, size_t pos, size_t *r_pos, size_t *r_end_pos)
{
	//Region
	size_t reg_pos = SIZE_MAX;
	size_t reg_end_pos;

	int i;

	//CHECK ARGUMENTS
	{
		//Check nulls
		assert(cigar);
		//assert(cigar_l > 0);
		assert(pos >= 0);
		assert(r_pos);
		assert(r_end_pos);
	}

	//If cigar length == 0 then return error
	if(cigar_l == 0)
		return -1;

	//Search for indels
	uint32_t elem, type;
	int disp = 0;
	for(i = 0; i < cigar_l; i++)
	{
		elem = cigar[i] >> BAM_CIGAR_SHIFT;	//Get number of bases from cigar
		type = cigar[i] & BAM_CIGAR_MASK;	//Get type from cigar

		switch(type)
		{
		case BAM_CINS:	//Insertion
			//printf("I");
			//Is first CIGAR elem?
			if(i == 0)
				break;

			//Set initial position
			if(reg_pos == SIZE_MAX)
			{
				reg_pos = pos + disp - 1;
			}

			//Set final position
			reg_end_pos = pos + disp;

			break;
		case BAM_CDEL:	//Deletion
			//printf("D");
			//Set initial position
			if(reg_pos == SIZE_MAX)
			{
				reg_pos = pos + disp;
			}

			//Set final position
			reg_end_pos = pos + disp + elem - 1;

			disp += elem;
			break;

		case BAM_CMATCH:
		case BAM_CDIFF:
		case BAM_CEQUAL:
			//Increment displacement
			disp += elem;
			break;

		case BAM_CHARD_CLIP:
		case BAM_CSOFT_CLIP:
			/*if(i == 0)
			{
				break;
			}*/

		default:
			//printf("X");

			break;
		}
	}
	//printf("\n");

	//Set output
	if(reg_pos != SIZE_MAX)
	{
		*r_pos = reg_pos;
		*r_end_pos = reg_end_pos;
		return 0;
	}
	else
	{
		*r_pos = SIZE_MAX;
		return 0;
	}

	return 0;
}

/**
 * Get a region from bam1_t.
 */
static inline int
region_get_from_bam1(const bam1_t *alig, size_t *r_pos, size_t *r_end_pos)
{
	//CIGAR
	uint32_t *cigar;
	uint32_t cigar_l;

	//Read pos
	size_t read_pos;

	//CHECK ARGUMENTS
	{
		assert(alig);
		assert(r_pos);
		assert(r_end_pos);
	}

	//Get CIGAR
	cigar = bam1_cigar(alig);
	cigar_l = alig->core.n_cigar;

	//Get read position
	read_pos = alig->core.pos + 1;

	//Call region get
	return region_get(cigar, cigar_l, read_pos, r_pos, r_end_pos);
}

/**
 * Get if BAM read overlaps a region.
 */
int
region_bam_overlap(bam1_t *read, alig_region_t *region)
{
	size_t read_pos;
	size_t read_end;

	assert(read);
	assert(region);

	//Same chrom?
	if(read->core.tid != region->chrom)
	{
		//Dont overlap
		return 0;
	}

	//Get read position
	read_pos = read->core.pos;
	read_end = read_pos + read->core.l_qseq;

	//Is in region?
	if(read_end < region->start_pos || read_pos > region->end_pos)
	{
		//Dont overlap
		return 0;
	}

	//Overlap
	return 1;
}

#endif /* ALIG_REGION_H_ */
