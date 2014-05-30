/*
 * wanderer.h
 *
 *  Created on: 24/04/2014
 *      Author: rmoreno
 */

#ifndef WANDERER_H_
#define WANDERER_H_

#include <assert.h>
#include <omp.h>

#include "aux/aux_common.h"
#include "aux/aux_library.h"
#include "aux/timestats.h"
#include "bam_region.h"

#define WANDERER_REGIONS_MAX 2000

#define WANDERER_SUCCESS 0
#define WANDERER_IN_PROGRESS 1
#define WANDERER_ERROR -1

/**
 * WANDERING FUNCTION DEFINITION
 */
typedef int (*wanderer_function)(void *, bam_region_t *, bam1_t *) ;
typedef int (*processor_function)(void *, bam_region_t *) ;

/**
 * BAM WANDERER STRUCT
 */
typedef struct {
	//I/O
	bam_file_t *input_file;
	bam_file_t *output_file;
	omp_lock_t output_file_lock;

	//Regions
	omp_lock_t regions_lock;
	bam_region_t **regions;
	size_t regions_l;

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

EXTERNC int bwander_run(bam_wanderer_t *wanderer);


/**
 * WINDOWS
 */
/*EXTERNC int bwander_window_register(bam_wanderer_t *wanderer, bam_region_window_t *window);
EXTERNC void bwander_window_clear(bam_wanderer_t *wanderer);*/

int filter_read(bam1_t *read, uint8_t filters);

int bwander_region_insert(bam_wanderer_t *wanderer, bam_region_t *region);

#endif /* WANDERER_H_ */
