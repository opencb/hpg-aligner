/*
 * wanderer.h
 *
 *  Created on: 24/04/2014
 *      Author: rmoreno
 */

#ifndef BFWORK_H_
#define BFWORK_H_

#include <assert.h>
#include <omp.h>
#include <libgen.h>

#include "containers/linked_list.h"
#include "aux/aux_common.h"
#include "aux/aux_library.h"
#include "aux/timestats.h"
#include "bam_region.h"

//VERSION
#define FWORK_VER_CURRENT		"1"
#define FWORK_VER_REVISION		"0"
#define FWORK_VER_AGE			"0"
#define FWORK_VER 			FWORK_VER_CURRENT"."FWORK_VER_REVISION"."FWORK_VER_AGE

//FIXED SIZES
#define FWORK_REGIONS_MAX 1000
#define FWORK_CONTEXT_MAX 	16
#define FWORK_PROC_FUNC_MAX 	16

/**
 * WANDERING FUNCTION DEFINITION
 */
typedef int (*wanderer_function)(void *, bam_region_t *, bam1_t *) ;
typedef int (*processor_function)(void *, bam_region_t *) ;

/**
 * WANDERING CONTEXT
 */
typedef struct {
	//Wandering function
	wanderer_function wander_f;

	//Processing functions
	processor_function processing_f[FWORK_PROC_FUNC_MAX];
	size_t processing_f_l;

	//User data
	void *user_data;
	omp_lock_t user_data_lock;
	void **local_user_data;
} bfwork_context_t;

/**
 * BAM WANDERER STRUCT
 */
typedef struct {
	//I/O
	char *input_file_str;
	char *output_file_str;
	char *reference_str;
	bam_file_t *input_file;
	bam_file_t *output_file;
	genome_t *reference;
	omp_lock_t output_file_lock;
	omp_lock_t reference_lock;

	//Regions
	omp_lock_t regions_lock;
	linked_list_t *regions_list;
	omp_lock_t free_slots;
	//size_t regions_l;

	//Current wandering context
	bfwork_context_t *context;

	//Contexts
	bfwork_context_t *v_context[FWORK_CONTEXT_MAX];
	size_t v_context_l;

	//Timing
	p_timestats time_stats;
} bam_fwork_t;

/**
 * BAM WANDERER OPERATIONS
 */
EXTERNC void bfwork_init(bam_fwork_t *fwork);
EXTERNC void bfwork_destroy(bam_fwork_t *fwork);

EXTERNC int bfwork_configure(bam_fwork_t *fwork, const char *in_file, const char *out_file, const char *reference, bfwork_context_t *context);

EXTERNC int bfwork_run(bam_fwork_t *fwork);

/**
 * WANDERING CONTEXT OPERATIONS
 */
EXTERNC void bfwork_context_init(bfwork_context_t *context, wanderer_function wf, processor_function pf);
EXTERNC void bfwork_context_destroy(bfwork_context_t *context);

EXTERNC int bfwork_context_add_proc(bfwork_context_t *context, processor_function pf);

/**
 * USER DATA
 */
EXTERNC int bfwork_context_set_user_data(bfwork_context_t *context, void *user_data);

static int bfwork_lock_user_data(bam_fwork_t *fwork, void **user_data);
static int bfwork_unlock_user_data(bam_fwork_t *fwork);
static int bfwork_local_user_data(bam_fwork_t *fwork, void **user_data);
static int bfwork_local_user_data_set(bam_fwork_t *fwork, void *user_data);

static int bfwork_context_local_user_data_reduce(bfwork_context_t *context, void *reduced, void (*cb_reduce)(void *, void *));
static int bfwork_context_local_user_data_free(bfwork_context_t *context, void (*cb_free)(void *));

/**
 * TIMING
 */
EXTERNC int bfwork_init_timing(bam_fwork_t *fwork, const char *tag);
EXTERNC void bfwork_destroy_timing(bam_fwork_t *fwork);
EXTERNC int bfwork_print_times(bam_fwork_t *fwork);

/**
 * FILTER
 */
static int filter_read(bam1_t *read, uint8_t filters);

/**
 * INLINE FUNCTIONS
 */

static inline int
bfwork_lock_user_data(bam_fwork_t *fwork, void **user_data)
{
	assert(fwork);

	//Take the lock
	omp_set_lock(&fwork->context->user_data_lock);

	//No user data?
	if(fwork->context->user_data == NULL)
	{
		LOG_WARN("Getting uninitialized user data\n");
	}

	//Set return
	*user_data = fwork->context->user_data;

	return NO_ERROR;
}
static inline int
bfwork_unlock_user_data(bam_fwork_t *fwork)
{
	assert(fwork);

	//Remove the lock
	omp_unset_lock(&fwork->context->user_data_lock);

	return NO_ERROR;
}

static inline int
bfwork_local_user_data(bam_fwork_t *fwork, void **user_data)
{
	int thread_id;

	assert(fwork);
	assert(user_data);

	//Get thread id
	thread_id = omp_get_thread_num();
	assert(thread_id >= 0);

	//Get user data
	*user_data = fwork->context->local_user_data[thread_id];

	return NO_ERROR;
}

static inline int
bfwork_local_user_data_set(bam_fwork_t *fwork, void *user_data)
{
	int thread_id;

	assert(fwork);
	assert(user_data);

	//Get thread id
	thread_id = omp_get_thread_num();
	assert(thread_id >= 0);

	//Get user data
	fwork->context->local_user_data[thread_id] = user_data;

	return NO_ERROR;
}

static inline int
bfwork_context_local_user_data_reduce(bfwork_context_t *context, void *reduced, void (*cb_reduce)(void *, void *))
{
	int i, threads;
	void *data;

	assert(context);
	assert(reduced);
	assert(cb_reduce);

	//Get num threads
	threads = omp_get_max_threads();

	//Iterate data
	for(i = 0; i < threads; i++)
	{
		//Get next data
		data = context->local_user_data[i];

		//Free if data exists
		if(data != NULL)
		{
			//Callback
			if(cb_reduce != NULL)
			{
				cb_reduce(data, reduced);
			}
		}
	}

	return NO_ERROR;
}

static inline int
bfwork_context_local_user_data_free(bfwork_context_t *context, void (*cb_free)(void *))
{
	int i, threads;
	void *data;

	assert(context);

	//Get num threads
	threads = omp_get_max_threads();

	//Iterate data
	for(i = 0; i < threads; i++)
	{
		//Get next data
		data = context->local_user_data[i];

		//Free if data exists
		if(data != NULL)
		{
			//Callback
			if(cb_free != NULL)
			{
				cb_free(data);
			}

			//Free memory
			free(data);
			context->local_user_data[i] = NULL;
		}
	}

	return NO_ERROR;
}

static inline int
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

	return NO_ERROR;
}

#endif /* BFWORK_H_ */
