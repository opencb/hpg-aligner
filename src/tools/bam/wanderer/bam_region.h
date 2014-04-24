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

#define BAM_REGION_DEFAULT_SIZE 100000

/**
 *	REGION STRUCT
 */
typedef struct {
	bam1_t **reads;
	size_t init_pos;
	size_t end_pos;
	int chrom;
	size_t size;
	size_t max_size;

	size_t processed;

	bam1_t *next_read;
} bam_region_t;

typedef struct {
	bam_region_t *region;
	size_t init_pos;
	size_t end_pos;
	size_t init_index;
	size_t size;
} bam_region_window_t;

/**
 * REGION OPERATIONS
 */
EXTERNC void breg_init(bam_region_t *region);
EXTERNC void breg_destroy(bam_region_t *region, int free_bam);

//EXTERNC void breg_load(bam_region_t *region, bam_file *file, size_t init_pos, size_t end_pos, int chrom);
EXTERNC void breg_fill(bam_region_t *region, bam_file_t *input_file);
EXTERNC void breg_write_processed(bam_region_t *region, bam_file_t *output_file);

#endif /* BAM_REGION_H_ */
