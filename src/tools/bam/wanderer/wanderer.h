/*
 * wanderer.h
 *
 *  Created on: 24/04/2014
 *      Author: rmoreno
 */

#ifndef WANDERER_H_
#define WANDERER_H_

#include <assert.h>

#include "aux/aux_common.h"
#include "aux/aux_library.h"
#include "aux/timestats.h"
#include "bam_region.h"

#define WANDERER_WINDOWS_MAX 1000

#define WANDERER_SUCCESS 0
#define WANDERER_IN_PROGRESS 1
#define WANDERER_ERROR -1

/**
 * WANDERING FUNCTION DEFINITION
 */
typedef int (*wanderer_function)(void *) ;
typedef int (*processor_function)(void *, void *) ;

/**
 * BAM WANDERER STRUCT
 */
typedef struct {
	//I/O
	bam_file_t *input_file;
	bam_file_t *output_file;

	//Loaded reads
	bam_region_t *current_region;
	size_t processed;

	//Windows to process
	bam_region_window_t *windows;
	size_t windows_l;

	//Function to execute
	wanderer_function wander_f;
	processor_function processing_f;
} bam_wanderer_t;

/**
 * BAM WANDERER OPERATIONS
 */
EXTERNC void bwander_init(bam_wanderer_t *wanderer);
EXTERNC void bwander_destroy(bam_wanderer_t *wanderer);

EXTERNC void bwander_configure(bam_wanderer_t *wanderer, bam_file_t *in_file, bam_file_t *out_file, wanderer_function wf, processor_function pf);

EXTERNC void bwander_run(bam_wanderer_t *wanderer);

/**
 * WINDOWS
 */
EXTERNC int bwander_window_register(bam_wanderer_t *wanderer, bam_region_window_t *window);
EXTERNC void bwander_window_clear(bam_wanderer_t *wanderer);

#endif /* WANDERER_H_ */
