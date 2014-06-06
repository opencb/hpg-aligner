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

#ifndef BAM_RECAL_LIBRARY_H_
#define BAM_RECAL_LIBRARY_H_

#include <bioformats/bam/samtools/bam.h>
#include <bioformats/bam/bam_file.h>
#include <bioformats/bam/alignment.h>
#include <aligners/bwt/genome.h>
#include "aux/aux_common.h"

/**
 * \brief Recalibration data storage.
 *
 * This struct hold all data necessary to perform recalibrations.
 */
struct recal_info {
	U_QUALS min_qual;	//Set minor stored quality
	U_QUALS num_quals;	//Set range of stored qualities
	U_CYCLES num_cycles;	//Set maximum number of cycles stored
	U_DINUC num_dinuc;	//Set number of maximum dinucleotides

	double total_miss;				//Total misses
	U_BASES total_bases;				//Total bases
	double total_delta;			//Global delta
	double total_estimated_Q;	//Global estimated Quality

	double* qual_miss;				//Misses per quality
	U_BASES* qual_bases;				//Bases per quality
    double* qual_delta;			//Delta per quality

    double* qual_cycle_miss;		//Misses per quality-cycle pair
    U_BASES* qual_cycle_bases;		//Bases per quality-cycle pair
    double* qual_cycle_delta;		//Deltas per quality-cycle pair

    double* qual_dinuc_miss;		//Misses per quality-dinuc pair
    U_BASES* qual_dinuc_bases;		//Bases per quality-dinuc pair
    double* qual_dinuc_delta;		//Deltas per quality-dinuc pair
};

/**
 * \brief Data collect environment storage.
 *
 * This struct hold all environment necessary to perform data collection.
 */
struct data_collect_env;

/**
 * \brief Recalibration environment storage.
 *
 * This struct hold all environment necessary to perform recalibration.
 */
struct recalibration_env;
typedef struct recal_info recal_info_t;
typedef struct data_collect_env recal_data_collect_env_t;
typedef struct recalibration_env recal_recalibration_env_t;

/**
 * \brief Dinucleotide enumeration.
 *
 * Represents all possible dinucleotide cases.
 */
typedef enum DINUC
{
	dAA = 0,
	dAG = 1,
	dAC = 2,
	dAT = 3,
	dGA = 4,
	dGG = 5,
	dGC = 6,
	dGT = 7,
	dCA = 8,
	dCG = 9,
	dCC = 10,
	dCT = 11,
	dTA = 12,
	dTG = 13,
	dTC = 14,
	dTT = 15,
	d_X = 16	/* Not a dinucleotide. For example "NT" or "NN".*/
} DINUCLEOTIDE;


/***********************************************
 * DATA MANAGEMENT
 **********************************************/

/**
 * \brief Initialize empty recalibration data struct.
 *
 * \param cycles Number of maximum cycles to stat.
 * \param out_info Previously allocated info struct to initialize.
 */
EXTERNC ERROR_CODE recal_init_info(const U_CYCLES cycles, recal_info_t *out_data);

/**
 * \brief Free all resources of recalibration.
 *
 * Free all data struct attributes including itself at the end.
 *
 * \param data Data struct to free
 */
EXTERNC ERROR_CODE recal_destroy_info(recal_info_t *data);

/**
 * \brief Reduce data.
 *
 * Reduce all data into one struct.
 *
 * \param dst_data Data struct to contain reduction
 * \param src_data Data struct to reduce with dst_data
 */
EXTERNC ERROR_CODE recal_reduce_info(recal_info_t *dst_data, recal_info_t *src_data);

/**
 * \brief Initialize empty data collection environment struct.
 *
 * \param cycles Number of maximum cycles to handle.
 * \param collect_env Previously allocated info struct to initialize.
 */
EXTERNC ERROR_CODE recal_get_data_init_env(const U_CYCLES cycles, recal_data_collect_env_t *collect_env);

/**
 * \brief Free all resources of data collect environment.
 *
 * Free all data struct attributes including itself at the end.
 *
 * \param data Data collect struct to free
 */
EXTERNC ERROR_CODE recal_get_data_destroy_env(recal_data_collect_env_t *collect_env);

/**
 * \brief Initialize empty recalibration environment struct.
 *
 * \param cycles Number of maximum cycles to handle.
 * \param recalibration_env Previously allocated info struct to initialize.
 */
EXTERNC ERROR_CODE recal_recalibration_init_env(const U_CYCLES cycles, recal_recalibration_env_t *recalibration_env);

/**
 * \brief Free all resources of recalibration environment.
 *
 * Free all data struct attributes including itself at the end.
 *
 * \param data Recalibration environment struct to free
 */
EXTERNC ERROR_CODE recal_recalibration_destroy_env(recal_recalibration_env_t *recalibration_env);

/**
 * \brief Add recalibration data from one base.
 *
 * \param data Data struct to add stats.
 * \param qual Quality to add.
 * \param cycle Cycle to add.
 * \param dinuc Dinucleotide to add.
 * \param match Indicate if match(!=0) or not (0)
 */
EXTERNC ERROR_CODE recal_add_base(recal_info_t *data, const char qual, const U_CYCLES cycle, const char dinuc, const double match) __ATTR_HOT;

/**
 * \brief Add recalibration data from vector of bases.
 *
 * Function checks if sequence nucleotides are valid.
 *
 * \param data Vector of data structs to add stats.
 * \param seq Sequence vector.
 * \param quals Qualities vector to add.
 * \param init_cycle Init cycle to add.
 * \param num_cycles Number of cycles to add.
 * \param dinuc Vector of dinucleotides to add.
 * \param matches Vector of match(!=0) or not (0)
 */
EXTERNC ERROR_CODE recal_add_base_v(recal_info_t *data, const char *seq, const char *quals, const U_CYCLES init_cycle, const U_CYCLES num_cycles, const char *dinuc, const char *matches) __ATTR_HOT;

/**
 * \brief Compute deltas from bases and misses.
 *
 * \param data Data to compute deltas.
 */
EXTERNC ERROR_CODE recal_calc_deltas(recal_info_t* data);


/***********************************************
 * DINUC ENUM OPERATIONS
 **********************************************/

/**
 * \brief Return dinucleotide enumeration from two bases.
 *
 * \param A First base.
 * \param B Second base.
 * \param out_dinuc Pointer to output dinucleotide.
 */
EXTERNC ERROR_CODE recal_get_dinuc(const char A, const char B, U_DINUC *out_dinuc) __ATTR_HOT;


/***********************************************
 * BAM RECALIBRATION PHASE 1 - DATA COLLECT
 **********************************************/

/**
 * \brief Get recalibration data from BAM path.
 *
 * \param bam_path Path to BAM.
 * \param ref_name String with the name of the reference genome.
 * \param ref_path Path to reference.
 * \param out_info Data struct to fill.
 */
EXTERNC ERROR_CODE recal_get_data_from_file(const char *bam_path, const char *ref_name, const char *ref_path, recal_info_t *out_info);

/**
 * \brief Get recalibration data from BAM file.
 *
 * \param bam BAM file struct to process.
 * \param ref Reference genome struct.
 * \param out_info Data struct to fill.
 */
EXTERNC ERROR_CODE recal_get_data_from_bam(const bam_file_t *bam, const genome_t* ref, recal_info_t* output_data);

/**
 * \brief Get recalibration data from BAM batch of alignments.
 *
 * \param batch BAM batch struct to process.
 * \param ref Reference genome struct.
 * \param out_info Data struct to fill.
 */
EXTERNC ERROR_CODE recal_get_data_from_bam_batch(const bam_batch_t* batch, const genome_t* ref, recal_info_t* output_data);

/**
 * \brief Get recalibration data from alignment.
 *
 * \param batch BAM alignment struct to process.
 * \param ref Reference genome struct.
 * \param out_info Data struct to fill.
 * \param collect_env Enviroment struct neccessary for data collection.
 */
EXTERNC ERROR_CODE recal_get_data_from_bam_alignment(const bam1_t* alig, const genome_t* ref, recal_info_t* output_data, recal_data_collect_env_t *collect_env) __ATTR_HOT;

/***********************************************
 * BAM RECALIBRATION PHASE 2 - RECALIBRATION
 **********************************************/

/**
 * \brief Recalibrate BAM file from path and store in file.
 *
 * \param orig_bam_path Path to BAM which will be recalibrated.
 * \param bam_info Data struct with recalibration info.
 * \param recal_bam_path Path to output BAM.
 */
EXTERNC ERROR_CODE recal_recalibrate_bam_file(const char *orig_bam_path, const recal_info_t *bam_info, const char *recal_bam_path);

/**
 * \brief Recalibrate BAM file and store in file.
 *
 * \param orig_bam_f BAM file struct to recalibrate.
 * \param bam_info Data struct with recalibration info.
 * \param recal_bam_f Recalibrated BAM output file struct.
 */
EXTERNC ERROR_CODE recal_recalibrate_bam(const bam_file_t *orig_bam_f, const recal_info_t *bam_info, bam_file_t *recal_bam_f);

/**
 * \brief Recalibrate BAM batch of alignments and store in file.
 *
 * \param batch Batch struct to recalibrate.
 * \param bam_info Data struct with recalibration info.
 * \param recal_bam_f Recalibrated BAM output file struct.
 */
EXTERNC ERROR_CODE recal_recalibrate_batch(const bam_batch_t* batch, const recal_info_t *bam_info);

/**
 * \brief Recalibrate alignment and store in file.
 *
 * \param alig BAM alignment struct to recalibrate.
 * \param bam_info Data struct with recalibration info.
 * \param recal_bam_f Recalibrated BAM output file struct.
 * \param recalibration_env Enviroment struct neccessary for data collection.
 */
EXTERNC ERROR_CODE recal_recalibrate_alignment(const bam1_t* alig, const recal_info_t *bam_info, recal_recalibration_env_t *recalibration_env) __ATTR_HOT;


/***********************************************
 * FILE OPERATIONS
 **********************************************/

/**
 * \brief Print to file data from recalibration.
 *
 * \param data Data struct with recalibration info.
 * \param path File path to print the info.
 */
EXTERNC ERROR_CODE recal_fprint_info(const recal_info_t *data, const char *path);

/**
 * \brief Save to file recalibration data.
 *
 * \param data Data struct with recalibration info.
 * \param path File path to save data.
 */
EXTERNC ERROR_CODE recal_save_recal_info(const recal_info_t *data, const char *path);

/**
 * \brief Load from file recalibration data.
 *
 * \param path File path to load data.
 * \param data Data struct to store info.
 */
EXTERNC ERROR_CODE recal_load_recal_info(const char *path, recal_info_t *data);

#endif /* BAM_RECAL_LIBRARY_H_ */
