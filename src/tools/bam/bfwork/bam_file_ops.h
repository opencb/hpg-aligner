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

#ifndef BAM_FILE_OPS_H_
#define BAM_FILE_OPS_H_

#include <assert.h>
#include <omp.h>

#include "aux/aux_common.h"
#include "aux/aux_library.h"
#include "aux/timestats.h"
#include "recalibrate/recal_config.h"
#include "recalibrate/recal_structs.h"
#include "recalibrate/bam_recal_library.h"
#include "aligner/alig.h"

/***********************************************
 * FRAMEWORK REALIGNER
 **********************************************/

/**
 * \brief Realign around indels BAM file
 *
 * \param[in] bam_path Input BAM file.
 * \param[in] ref_path Path to reference 'dna_compression.bin'.
 * \param[in] outbam Output BAM path. OPTIONAL, if NULL no output will be written.
 */
EXTERNC ERROR_CODE alig_bam_file(const char *bam_path, const char *ref_path, const char *outbam, const char *stats_path);


/***********************************************
 * FRAMEWORK RECALIBRATOR
 **********************************************/

/**
 * \brief Recalibrate BAM file
 *
 * \param[in] flags Flags for recalibration. Can be RECALIBRATE_COLLECT and RECALIBRATE_RECALIBRATE.
 * \param[in] bam_path Input BAM file.
 * \param[in] ref Path to reference 'dna_compression.bin'. OPTIONAL if not RECALIBRATE_COLLECT.
 * \param[in/out] data_file Path to data file. OPTIONAL, can be input if RECALIBRATE_RECALIBRATE only or output if not.
 * \param[in] info_file Path to text format of data file. OPTIONAL, if NULL no info file will be written.
 * \param[in] outbam Output BAM path. OPTIONAL, if NULL no output will be written.
 * \param[in] cycles Max cycles to recalibrate.
 * \param[in] stats_path Path to timing stats output folder. OPTIONAL, if NULL no timing output will be written.
 */
EXTERNC ERROR_CODE recal_bam_file(uint8_t flags, const char *bam_path, const char *ref_path, const char *data_file, const char *info_file, const char *outbam, int cycles, const char *stats_path);

/**
 * \brief Realign and recalibrate BAM file
 *
 * \param[in] bam_path Input BAM file.
 * \param[in] ref Path to reference 'dna_compression.bin'. OPTIONAL if not RECALIBRATE_COLLECT.
 * \param[in/out] data_file Path to data file. OPTIONAL, can be input if RECALIBRATE_RECALIBRATE only or output if not.
 * \param[in] info_file Path to text format of data file. OPTIONAL, if NULL no info file will be written.
 * \param[in] outbam Output BAM path. OPTIONAL, if NULL no output will be written.
 * \param[in] cycles Max cycles to recalibrate.
 * \param[in] stats_path Path to timing stats output folder. OPTIONAL, if NULL no timing output will be written.
 */
EXTERNC ERROR_CODE alig_recal_bam_file(const char *bam_path, const char *ref_path, const char *data_file, const char *info_file, const char *outbam, int cycles, const char *stats_path);

#endif /* BAM_FILE_OPS_H_ */
