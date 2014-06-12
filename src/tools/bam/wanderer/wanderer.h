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

#include "containers/linked_list.h"
#include "aux/aux_common.h"
#include "aux/aux_library.h"
#include "aux/timestats.h"
#include "bam_region.h"

#define WANDERER_REGIONS_MAX 1000
#define WANDERER_PROC_FUNC_MAX 	16

#define WANDERER_SUCCESS 0
#define WANDERER_IN_PROGRESS 1
#define WANDERER_ERROR -1

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
	processor_function processing_f[WANDERER_PROC_FUNC_MAX];
	size_t processing_f_l;

	//User data
	void *user_data;
	omp_lock_t user_data_lock;
	void **local_user_data;
} bwander_context_t;

/**
 * BAM WANDERER STRUCT
 */
typedef struct {
	//I/O
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

	//Wandering context
	bwander_context_t *context;
} bam_wanderer_t;

/**
 * BAM WANDERER OPERATIONS
 */
EXTERNC void bwander_init(bam_wanderer_t *wanderer);
EXTERNC void bwander_destroy(bam_wanderer_t *wanderer);

EXTERNC int bwander_configure(bam_wanderer_t *wanderer, bam_file_t *in_file, bam_file_t *out_file, genome_t *reference, bwander_context_t *context);

EXTERNC int bwander_run(bam_wanderer_t *wanderer);

/**
 * WANDERING CONTEXT OPERATIONS
 */
EXTERNC void bwander_context_init(bwander_context_t *context, wanderer_function wf, processor_function pf);
EXTERNC void bwander_context_destroy(bwander_context_t *context);

EXTERNC int bwander_context_add_proc(bwander_context_t *context, processor_function pf);

/**
 * USER DATA
 */
EXTERNC int bwander_context_set_user_data(bwander_context_t *context, void *user_data);

static int bwander_lock_user_data(bam_wanderer_t *wanderer, void **user_data);
static int bwander_unlock_user_data(bam_wanderer_t *wanderer);
static int bwander_local_user_data(bam_wanderer_t *wanderer, void **user_data);
static int bwander_local_user_data_set(bam_wanderer_t *wanderer, void *user_data);

static int bwander_context_local_user_data_reduce(bwander_context_t *context, void *reduced, void (*cb_reduce)(void *, void *));
static int bwander_context_local_user_data_free(bwander_context_t *context, void (*cb_free)(void *));

/**
 * FILTER
 */
static int filter_read(bam1_t *read, uint8_t filters);

/**
 * INLINE FUNCTIONS
 */

static inline int
bwander_lock_user_data(bam_wanderer_t *wanderer, void **user_data)
{
	assert(wanderer);

	//Take the lock
	omp_set_lock(&wanderer->context->user_data_lock);

	//No user data?
	if(wanderer->context->user_data == NULL)
	{
		LOG_WARN("Getting uninitialized user data\n");
	}

	//Set return
	*user_data = wanderer->context->user_data;

	return NO_ERROR;
}
static inline int
bwander_unlock_user_data(bam_wanderer_t *wanderer)
{
	assert(wanderer);

	//Remove the lock
	omp_unset_lock(&wanderer->context->user_data_lock);

	return NO_ERROR;
}

static inline int
bwander_local_user_data(bam_wanderer_t *wanderer, void **user_data)
{
	int thread_id;

	assert(wanderer);
	assert(user_data);

	//Get thread id
	thread_id = omp_get_thread_num();
	assert(thread_id >= 0);

	//Get user data
	*user_data = wanderer->context->local_user_data[thread_id];

	return NO_ERROR;
}

static inline int
bwander_local_user_data_set(bam_wanderer_t *wanderer, void *user_data)
{
	int thread_id;

	assert(wanderer);
	assert(user_data);

	//Get thread id
	thread_id = omp_get_thread_num();
	assert(thread_id >= 0);

	//Get user data
	wanderer->context->local_user_data[thread_id] = user_data;

	return NO_ERROR;
}

static inline int
bwander_context_local_user_data_reduce(bwander_context_t *context, void *reduced, void (*cb_reduce)(void *, void *))
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
bwander_context_local_user_data_free(bwander_context_t *context, void (*cb_free)(void *))
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

#endif /* WANDERER_H_ */
