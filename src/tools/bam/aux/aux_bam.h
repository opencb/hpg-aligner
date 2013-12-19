/*
 * aux_bam.h
 *
 *  Created on: Jun 25, 2013
 *      Author: rmoreno
 */

#ifndef AUX_BAM_H_
#define AUX_BAM_H_

#include "aux_library.h"

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

#endif /* AUX_BAM_H_ */
