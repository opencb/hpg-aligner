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

#ifndef AUX_LIBRARY_H_
#define AUX_LIBRARY_H_

#include <bioformats/bam/samtools/bam.h>
#include <bioformats/bam/bam_file.h>

#include "aux_common.h"

/***************************
 * MATH OPERATIONS
 **************************/

/**
 * Returns quality from probability.
 * \param p Error probability [0,1]
 */
EXTERNC double Qvalue(double p);

/**
 * Return probability from quality.
 * \param q Error quality. Default -> [P_SANGER_MIN,P_SANGER_MAX]
 */
EXTERNC double Pvalue(double q);

/**
 * Return Solexa quality from probability.
 * \param p Error probability [0,1]
 */
EXTERNC double Qsolexa(double p);

/**
 * Return probability from Solexa quality.
 * \param q Error quality [P_SOLEXA_MIN,P_SOLEXA_MAX]
 */
EXTERNC double Psolexa(double q);

/**
 * Return Sanger quality from probability.
 * \param p Error probability [0,1]
 */
EXTERNC double Qsanger(double p);

/**
 * Return probability from Sanger quality.
 * \param q Error quality [P_SANGER_MIN,P_SANGER_MAX]
 */
EXTERNC double Psanger(double q);

/**
 * Return evaluation of gaussian function.
 * g: a + b * exp(-pow(value-c,2)/(2*d*d))
 * \param value Value to evaluate.
 * \param a A constant
 * \param b B constant
 * \param c C constant
 * \param d D constant
 */
EXTERNC double gaussian_function(double value, double a, double b, double c, double d);

/**
 * Return approximation for log10(n!)
 * \param n Number to apply factorial and logarithm.
 */
EXTERNC double log10_gamma(uint64_t n);

/***************************
 * QUALITY OPERATIONS
 **************************/

/**
 * Obtains estimated global quality from a vector of bases.
 * \param v_bases 	Vector containing number of bases for each quality.
 * 					Each component represents a count per each quality, beginning from start_quality value.
 * 					For i = 0 -> number of bases for quality start_quality.
 * 					For i = 1 -> number of bases for quality start_quality + 1 ...
 * \param count	Number of vector components.
 * \param start_quality Indicates Quality value for first component of v_bases.
 * \param estimated_Q Pointer to store result of estimated global quality.
 */
EXTERNC ERROR_CODE recal_get_estimated_Q(U_BASES *v_bases, size_t count, U_QUALS start_quality, double *estimated_Q);

/**
 * Obtains empirical quality from a count of bases and misses and a theorical quality.
 * \param miss Counts how many errors are in bases.
 * \param bases Counts how many bases.
 * \param initial_quality Theorical quality.
 * \param emp_qual Pointer to store result of empirical quality.
 */
EXTERNC ERROR_CODE recal_get_empirical_Q(double miss, U_BASES bases, double initial_quality, double *emp_qual);

/**
 * Obtains how much one empirical quality approximates one theorical quality, expressed by logarithm.
 * \param Qemp Empirical quality to be evaluated.
 * \param Qreported Theorical quality to compare.
 * \param log Pointer to store logarithm of comparation result.
 */
EXTERNC ERROR_CODE log10_Qemp_Reported(double Qemp, double Qreported, double *log);

/**
 * Obtains how much one empirical quality approximates one theorical quality expressed by a count of observations and errors.
 * Result is expressed by a logarithm.
 * \param Qemp Empirical quality to be evaluated.
 * \param obs Count of observations.
 * \param err Count of errors in this observations.
 * \param log Pointer to store logarithm of comparation result.
 */
EXTERNC ERROR_CODE log10_Qemp_likelihood(double Qemp, U_BASES obs, U_BASES err, double *log);

/***************************
 * BAM OPERATIONS
 **************************/

/**
 * Init bam_header_t struct to work with recalibration.
 * \param num_chroms Number of target chromosomes
 * \param header Pointer to previously allocated empty header struct
 */
EXTERNC ERROR_CODE init_empty_bam_header(const unsigned int num_chroms, bam_header_t *header);

EXTERNC ERROR_CODE compare_bams_qual(const char* bamPath0, const char* bamPath1, const int cycles);

/**
 * Get string containing bam1 sequence nucleotides.
 * \param bam1 Pointer to bam1 struct containing the sequence we want to extract.
 * \return Allocated string containing bam1 sequence nucleotides.
 */
EXTERNC char * new_sequence_from_bam(bam1_t *bam1);

/**
 * Get string containing bam1 sequence nucleotides.
 * \param bam1 Pointer to bam1 struct containing the sequence we want to extract.
 * \param seq Pointer to allocated character vector for sequence storage.
 * \param max_l Maximum length can handle seq vector.
 */
EXTERNC ERROR_CODE new_sequence_from_bam_ref(bam1_t *bam1, char *seq, U_CYCLES max_l);

/**
 * Get string containing bam1 quality for every nucleotide.
 * \param bam1 Pointer to bam1 struct containing the qualities we want to extract.
 * \param base_quality What value must have quality representing probability zero.
 * \return Allocated string containing bam1 nucleotides qualities.
 */
EXTERNC char * new_quality_from_bam(bam1_t *bam1, int base_quality);

/**
 * Get string containing bam1 quality for every nucleotide.
 * \param bam1 Pointer to bam1 struct containing the qualities we want to extract.
 * \param base_quality What value must have quality representing probability zero.
 * \param qual Pointer to allocated character vector for quality storage.
 * \param max_l Maximum length can handle qual vector.
 */
EXTERNC ERROR_CODE new_quality_from_bam_ref(bam1_t *bam1, U_QUALS base_quality, char *qual, U_CYCLES max_l);

/**
 * Decomponse a string cigar in two vectors, containing number of ocurrences for every cigar operation (I, D, =, M, etc..)
 * \param cigar String containing cigar to be decomposed.
 * \param cigar_l Length of cigar string.
 * \param n_elem Output vector containing ocurrences of every cigar operation.
 * \param type Output vector containing every operation in cigar.
 * \param types_l Lenght of output vectors.
 * \param max_types_length Maximum components output vectors can handle.
 */
EXTERNC ERROR_CODE decompose_cigar(char *cigar, uint8_t cigar_l, char *n_elem, char *type, uint8_t *types_l, uint8_t max_types_length);

EXTERNC ERROR_CODE supress_indels(char *seq, U_CYCLES seq_l, char *cigar_elem, char *cigar_type, uint8_t cigar_type_l, char *seq_res, U_CYCLES *seq_res_l) __attribute__((deprecated));

/**
 * Get sequence and quality vector supressing insertions and adding deletions.
 * This is used to compare directly with reference sequence.
 * \param seq Vector containing original nucleotide sequence.
 * \param qual Vector containing original nucleotide qualities.
 * \param cigar Cigar string to apply.
 * \param cigar_l Cigar string length.
 * \param seq_res Pointer to store sequence vector.
 * \param qual_res Pointer to store qualities vector.
 * \param seq_res_l Pointer to store length of result vectors.
 * \param max_res_l Maximum number of components result vectors can handle.
 */
EXTERNC ERROR_CODE supress_indels_from_32_cigar(char *seq, char *qual, int32_t seq_l, uint32_t *cigar, uint16_t cigar_l, char *seq_res, char *qual_res, uint32_t *seq_res_l, uint32_t max_res_l);

EXTERNC ERROR_CODE batch_count_chroms(bam_batch_t *batch, size_t *chrom_l);

EXTERNC ERROR_CODE batch_split_by_chrom(bam_batch_t *batch, bam_batch_t *v_batchs, size_t *res_batch_l, size_t max_res_batchs);

/***************************
 * VECTOR OPERATIONS
 **************************/

/**
 * Initializes vector with initial values.
 * \param vector Vector to initialize
 * \param size Size of vector
 * \param value Value to initialize vector
 */
EXTERNC ERROR_CODE initialize_vector(uint32_t *vector, const size_t size, const uint32_t value);

/**
 * Return vector of integers with initial values.
 * \param size Size of the new vector.
 * \param value Initial value in all elements of the vector
 * \param out_vector Pointer to pointer which stores the vector (allocated in function)
 */
EXTERNC ERROR_CODE new_vector_uint32(const size_t size, const uint32_t value, uint32_t **out_vector);

/**
 * Return vector of double with initial values.
 * \param size Size of the new vector.
 * \param value Initial value in all elements of the vector
 * \param out_vector Pointer to pointer which stores the vector (allocated in function)
 */
EXTERNC ERROR_CODE new_vector_double(const size_t size, const double value, double **out_vector);

/**
 * Get max value of vector.
 * \param vector Vector to evaluate.
 * \param size Size of vector.
 * \param max Pointer to store maximum value.
 */
EXTERNC ERROR_CODE max_value(double *vector, size_t size, double *max);

/**
 * Get component index of vector max value.
 * \param vector Vector to evaluate.
 * \param size Size of vector.
 * \param max_i Pointer to store maximum value component index.
 */
EXTERNC ERROR_CODE max_index(double *vector, size_t size, uint16_t *max_i);

/***************************
 * MISCELANEA OPERATIONS
 **************************/

EXTERNC void printf_proc_features();

EXTERNC void print_binary(unsigned int num);

#endif /* AUX_LIBRARY_H_ */
