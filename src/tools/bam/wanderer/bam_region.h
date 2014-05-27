/*
 * bam_region.h
 *
 *  Created on: 24/04/2014
 *      Author: rmoreno
 */

#ifndef BAM_REGION_H_
#define BAM_REGION_H_

#include <assert.h>

#include "aux/aux_common.h"
#include "aux/aux_library.h"
#include "aux/timestats.h"

#define BAM_REGION_DEFAULT_SIZE 10000

#define FILTER_ZERO_QUAL 1
#define FILTER_DIFF_MATE_CHROM 2
#define FILTER_NO_CIGAR 4
#define FILTER_DEF_MASK 8

#define EMPTY_CHROM -1


/**
 *	REGION STRUCT
 */
typedef struct {
	bam1_t **reads;
	size_t size;
	size_t max_size;
	size_t processed;
	bam1_t *next_read;

	//Locus
	size_t init_pos;
	size_t end_pos;
	int chrom;
} bam_region_t;

typedef struct {
	bam_region_t *region;
	bam1_t **filter_reads;
	size_t size;
	uint8_t filter_flags;

	//Locus
	size_t init_pos;
	size_t end_pos;
} bam_region_window_t;

/**
 * REGION OPERATIONS
 */
EXTERNC void breg_init(bam_region_t *region);
EXTERNC void breg_destroy(bam_region_t *region, int free_bam);

EXTERNC void breg_fill(bam_region_t *region, bam_file_t *input_file);
EXTERNC void breg_write_processed(bam_region_t *region, bam_file_t *output_file);

EXTERNC void breg_load_window(bam_region_t *region, size_t init_pos, size_t end_pos, uint8_t filters, bam_region_window_t *window);

/**
 * WINDOW OPERATIONS
 */
EXTERNC void breg_window_init(bam_region_window_t *window);
EXTERNC void breg_window_destroy(bam_region_window_t *window);
EXTERNC void breg_window_filter(bam_region_window_t *window, uint8_t filters);

#endif /* BAM_REGION_H_ */
