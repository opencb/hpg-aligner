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

#define WANDERER_SUCCESS 0
#define WANDERER_IN_PROGRESS 1
#define WANDERER_ERROR -1

/**
 * WANDERING FUNCTION DEFINITION
 */
typedef int (*wandering_function)(void *) ;

/**
 * BAM WANDERER STRUCT
 */
typedef struct {
	//I/O
	bam_file_t *input_file;
	bam_file_t *output_file;

	//Loaded reads
	bam_region_t *current_region;

	//Function to execute
	wandering_function wander_f;
} bam_wanderer_t;

/**
 * BAM WANDERER OPERATIONS
 */
EXTERNC void bwander_init(bam_wanderer_t *wanderer);
EXTERNC void bwander_destroy(bam_wanderer_t *wanderer);

EXTERNC void bwander_configure(bam_wanderer_t *wanderer, bam_file_t *in_file, bam_file_t *out_file, wandering_function f);

EXTERNC void bwander_run(bam_wanderer_t *wanderer);

/**
 * REGISTER WINDOW
 */
EXTERNC void bwander_window_register(bam_wanderer_t *wanderer, bam_region_window_t *window);

#endif /* WANDERER_H_ */
