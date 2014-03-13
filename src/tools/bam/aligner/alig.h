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

#ifndef ALIG_H_
#define ALIG_H_

#include <assert.h>

#include <stdlib.h>
#include <stdint.h>
#include <limits.h>
#include <stddef.h>
#include <unistd.h>

#include <omp.h>

#include <bioformats/bam/samtools/bam.h>
#include <bioformats/bam/bam_file.h>
#include <aligners/bwt/genome.h>
#include "aux/aux_common.h"
#include "aux/aux_library.h"
#include "aux/timestats.h"
#include "alig_region.h"

//VERSION
#define ALIG_VER_CURRENT		"0"
#define ALIG_VER_REVISION		"4"
#define ALIG_VER_AGE			"0"
#define ALIG_VER 			ALIG_VER_CURRENT"."ALIG_VER_REVISION"."ALIG_VER_AGE

//OPTIONS
#define ALIG_LIST_IN_SIZE	10000
#define ALIG_LIST_NEXT_SIZE 1000
#define ALIG_LIST_COUNT_THRESHOLD_TO_WRITE 500

#define ALIG_REFERENCE_ADDITIONAL_OFFSET 100
#define ALIG_REFERENCE_CORRECTION_OFFSET 4

#define ALIG_IMPROVEMENT_THREHOLD 0.0

//FLAGS
#define ALIG_LEFT_ALIGN 0x01	//Left align cigars?
#define ALIG_REFERENCE_PRIORITY 0x02	//Reference haplotype have priority over alternative? (if equals case)
#define ALIG_ORIGINAL_PRIORITY 0x04	//Original cigar have preference over realigned one? (if equals case)

/**
 * Time measures
 */
//#define D_TIME_DEBUG
#ifdef D_TIME_DEBUG
	enum alig_slots {
		//GENERAL
		D_SLOT_TOTAL,
		D_SLOT_INIT,

		//ITERATION
		D_SLOT_IT_PROCESS,
		D_SLOT_IT_READ,
		D_SLOT_IT_WRITE,

		//BAM I/0
		D_SLOT_ALIG_READ,
		D_SLOT_ALIG_WRITE,

		//REALIGN
		D_SLOT_NEXT,
		D_SLOT_REFERENCE_LOAD,
		D_SLOT_HAPLO_GET,
		D_SLOT_REALIG_PER_HAPLO,
		D_SLOT_NUCLEO_CMP
	};

#endif


/**
 * REALIGNMENT REFERENCE SEQUENCE
 */
typedef struct {
	char *reference;
	size_t position;
	size_t length;
} alig_reference_t;

/**
 * REALIGNMENT SCORE TABLE
 */
typedef struct {
	uint32_t *m_scores;
	size_t *m_positions;
	size_t m_total;
	size_t m_ldim;
} alig_scores_t;
/**
 * 					H0			H1			H2			H3			H4
 * 					sco	pos
 * 		ERRXXXX		23	97		0	16		321	34
 * 		ERRYYYY
 * 		ERRZZZZ
 */

/**
 * REALIGNER CONTEXT
 */
typedef struct {
	//Reference genome
	genome_t *genome;

	//BAM lists
	array_list_t *filtered_list;
	array_list_t *realign_list;

	//Alignments readed
	size_t read_count;
	size_t last_readed_count;

	//Haplotypes
	array_list_t *haplo_list;

	//Current region
	alig_region_t region;

	//Last region
	alig_region_t last_region;

	//Current reference
	alig_reference_t reference;

	//Scores
	alig_scores_t scores;

	//Left align
	uint8_t flags;

} alig_context_t;


/**
 * INTERVAL STATUS
 */
typedef enum {
	NO_INTERVAL,
	INTERVAL
} alig_status;

/**
 * CONTEXT
 */

/**
 * \brief Initialize empty realignment data structure.
 *
 * \param[in] context Context to initialize.
 * \param[in] in_list Input list from which realigner take readings.
 * \param[in] genome Context containing reference genome.
 * \param[in] flags Flags to configure realigner behavior. Can be ALIG_LEFT_ALIGN, ALIG_REFERENCE_PRIORITY and ALIG_ORIGINAL_PRIORITY.
 */
EXTERNC ERROR_CODE alig_init(alig_context_t *context, genome_t *genome, uint8_t flags);

/**
 * \brief Free resources from realigner.
 *
 * \param[in] context Context to destroy.
 */
EXTERNC ERROR_CODE alig_destroy(alig_context_t *context);

/**
 * \brief Check if a realigner context is valid.
 *
 * \param[in] context Context to validate.
 */
EXTERNC ERROR_CODE alig_validate(alig_context_t *context);

/**
 * REGION OPERATIONS
 */

/**
 * \brief Get next region of reads to process.
 *
 * \param[in] context Context to process.
 */
EXTERNC ERROR_CODE alig_region_next(bam1_t **v_bams, size_t v_bams_l, int force_incomplete, alig_context_t *context);

/**
 * \brief Load reference sequence for present region in context.
 *
 * \param[in] context Context to process.
 */
EXTERNC ERROR_CODE alig_region_load_reference(alig_context_t *context);

/**
 * \brief Get haplotypes from present region.
 *
 * \param[in] context Context to process.
 */
EXTERNC ERROR_CODE alig_region_haplotype_process(alig_context_t *context);

/**
 * \brief Realign readings around indels using current haplotypes.
 *
 * \param[in] context Context to process.
 */
EXTERNC ERROR_CODE alig_region_indel_realignment(alig_context_t *context);

/**
 * \brief Clear context to be ready for next region.
 *
 * \param[in] context Context to process.
 */
EXTERNC ERROR_CODE alig_region_clear(alig_context_t *context);

/**
 * \brief Indel realign one file.
 *
 * \param[in] bam_path Path to the input BAM file.
 * \param[in] ref_name Reference file name (not including path).
 * \param[in] ref_path Path to reference file (not including name).
 * \param[in] outbam Path to output BAM file.
 */
EXTERNC ERROR_CODE alig_bam_file(char *bam_path, char *ref_name, char *ref_path, char *outbam);


/**
 * PRIVATE FUNCTIONS
 */

/**
 * \brief PRIVATE FUNCTION. Obtain score tables from present region.
 *
 * \param[in] context Context to process.
 */
/*static*/ ERROR_CODE alig_get_scores(alig_context_t *context);
/*static*/ ERROR_CODE alig_get_scores_from_read(bam1_t *read, alig_context_t *context, uint32_t *v_scores, size_t *v_positions);

/**
 * \brief PRIVATE FUNCTION. Obtain alternative haplotype from generated score tables.
 *
 * \param[in] context Context to process.
 * \param[out] out_haplo_index Index in haplotype list of alternative haplotype.
 * \param[out] out_haplo_score Alternative haplotype score.
 * \param[out] out_ref_score Reference haplotype (H0) score.
 */
static ERROR_CODE alig_get_alternative_haplotype(alig_context_t *context, int *out_haplo_index, uint32_t *out_haplo_score, uint32_t *out_ref_score);

/**
 * \brief PRIVATE FUNCTION. Realign around indels using an alternative haplotype.
 *
 * \param[in] context Context to process.
 * \param[in] alt_haplo_index Index of alternative haplotype in haplotypes list.
 */
static ERROR_CODE alig_indel_realign_from_haplo(alig_context_t *context, size_t alt_haplo_index);


#endif /* ALIG_H_ */
