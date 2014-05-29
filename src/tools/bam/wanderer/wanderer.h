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

	//Regions
	bam_region_t **regions;
	size_t regions_l;

	//Function to execute
	wanderer_function wander_f;
	processor_function processing_f;
} bam_wanderer_t;

/**
 * BAM WANDERER OPERATIONS
 */
static INLINE int filter_read(bam1_t *read, uint8_t filters);
EXTERNC void bwander_init(bam_wanderer_t *wanderer);
EXTERNC void bwander_destroy(bam_wanderer_t *wanderer);

EXTERNC void bwander_configure(bam_wanderer_t *wanderer, bam_file_t *in_file, bam_file_t *out_file, wanderer_function wf, processor_function pf);

EXTERNC int bwander_run(bam_wanderer_t *wanderer);

static INLINE int bwander_region_insert(bam_wanderer_t *wanderer, bam_region_t *region);

/**
 * WINDOWS
 */
/*EXTERNC int bwander_window_register(bam_wanderer_t *wanderer, bam_region_window_t *window);
EXTERNC void bwander_window_clear(bam_wanderer_t *wanderer);*/

static INLINE int
filter_read(bam1_t *read, uint8_t filters)
{
	assert(read);

	//Filter read
	if(filters != 0)
	{
		if(filters & FILTER_ZERO_QUAL)
		{
			if(read->core.qual == 0)
				return 1;
		}

		if(filters & FILTER_DIFF_MATE_CHROM)
		{
			if(read->core.tid != read->core.mtid)
				return 1;
		}

		if(filters & FILTER_NO_CIGAR)
		{
			if(read->core.n_cigar == 0)
				return 1;
		}

		if(filters & FILTER_DEF_MASK)
		{
			if(read->core.flag & BAM_DEF_MASK)
				return 1;
		}
	}

	return 0;
}

static INLINE int
bwander_region_insert(bam_wanderer_t *wanderer, bam_region_t *region)
{
	assert(wanderer);
	assert(region);

	if(wanderer->regions_l >= WANDERER_REGIONS_MAX)
	{
		LOG_FATAL("Not enough region slots\n");
	}

	wanderer->regions[wanderer->regions_l] = region;
	wanderer->regions_l++;

	LOG_INFO_F("Registering region %d:%d-%d with %d reads\n",
				region->chrom + 1, region->init_pos + 1,
				region->end_pos + 1, region->size);
	LOG_INFO_F("Regions to process %d\n", wanderer->regions_l);

	return NO_ERROR;
}

#endif /* WANDERER_H_ */
