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

//CONTEXT EXECUTION
#define FWORK_CONTEXT_QUEUE_SEQUENTIAL	0x01
#define FWORK_CONTEXT_QUEUE_PARALLEL 0x02

//FIXED SIZES
#define FWORK_REGIONS_MAX 1000
#define FWORK_CONTEXT_MAX 	16
#define FWORK_PROC_FUNC_MAX 	16

//ALIGNMENTS FILTERS
#define FILTER_ZERO_QUAL 1
#define FILTER_DIFF_MATE_CHROM 2
#define FILTER_NO_CIGAR 4
#define FILTER_DEF_MASK 8

/**
 * WANDERING FUNCTION DEFINITION
 */
typedef int (*wanderer_function)(void *, bam_region_t *, bam1_t *) ;
typedef int (*processor_function)(void *, bam_region_t *) ;

/**
 * CONTEXT REDUCTION FUNCTION AFTER RUN
 */
typedef int (*reducer_function)(void *, void *) ;

/**
 * WANDERING CONTEXT
 */
typedef struct {
	//Wandering function
	wanderer_function wander_f;

	//Processing functions
	processor_function processing_f[FWORK_PROC_FUNC_MAX];
	size_t processing_f_l;

	//Callback function
	reducer_function reduce;
	void *reduce_dest;

	//User data
	void *user_data;
	omp_lock_t user_data_lock;
	void **local_user_data;

	//Timing
	p_timestats time_stats;
	char *tag;
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
} bam_fwork_t;

/**
 * BAM WANDERER OPERATIONS
 */

/**
 * \brief Initialize empty BAM framework data structure.
 *
 * \param[in] fwork Framework to initialize.
 */
EXTERNC void bfwork_init(bam_fwork_t *fwork);

/**
 * \brief Destroy BAM framework data structure.
 *
 * \param[in] fwork Framework to destroy.
 */
EXTERNC void bfwork_destroy(bam_fwork_t *fwork);

/**
 * \brief Configure BAM framework data structure for wandering.
 *
 * \param[in] fwork Framework to configure.
 * \param[in] in_file Path to input BAM file.
 * \param[in] out_file Path to output BAM file. OPTIONAL, set NULL for no output.
 * \param[in] reference Path to reference file 'dna_compression.bin'.
 * \param[in] context Initial framework context to use.
 */
EXTERNC int bfwork_configure(bam_fwork_t *fwork, const char *in_file, const char *out_file, const char *reference, bfwork_context_t *context);

/**
 * \brief Add additional context to execute in framework.
 *
 * \param[in] fwork Framework to add context.
 * \param[in] context Context to be used in framework.
 * \param[in] flags Specify execution queue for this context. Can be FWORK_CONTEXT_QUEUE_SEQUENTIAL or FWORK_CONTEXT_QUEUE_PARALLEL.
 */
EXTERNC int bfwork_add_context(bam_fwork_t *fwork, bfwork_context_t *context, uint8_t flags);

/**
 * \brief Run framework contexts.
 *
 * \param[in] fwork Framework to run.
 */
EXTERNC int bfwork_run(bam_fwork_t *fwork);

/**
 * WANDERING CONTEXT OPERATIONS
 */

/**
 * \brief Initialize empty BAM context data structure.
 *
 * \param[in] context Context to configure.
 * \param[in] wf Wandering function to use with this context.
 * \param[in] pf Processing function to use with this context.
 */
EXTERNC void bfwork_context_init(bfwork_context_t *context, wanderer_function wf, processor_function pf, reducer_function rf, void *reduce_dest);

/**
 * \brief Destroy BAM context data structure.
 *
 * \param[in] context Context to destroy.
 */
EXTERNC void bfwork_context_destroy(bfwork_context_t *context);

/**
 * \brief Add additional processing function to this context.
 *
 * \param[in] context Target context.
 * \param[in] pf Processing function to add.
 */
EXTERNC int bfwork_context_add_proc(bfwork_context_t *context, processor_function pf);

/**
 * USER DATA
 */

/**
 * \brief Set a pointer to be shared as a global user data among all threads in a context.
 *
 * \param[in] context Target context.
 * \param[in] user_data Pointer to be shared.
 */
EXTERNC int bfwork_context_set_user_data(bfwork_context_t *context, void *user_data);

/**
 * \brief Lock global user data to be used in current context of the framework.
 * This data is in mutual exclusion so unlock must be used with this function to avoid deadlocks.
 * This function must be only used in processing functions.
 *
 * \param[in] fwork Target framework which contains current context.
 * \param[out] user_data Global user data returned.
 */
static int bfwork_lock_user_data(bam_fwork_t *fwork, void **user_data);

/**
 * \brief Unlock global user data previously locked by 'bfwork_lock_user_data'.
 * This function must be only used in processing functions.
 *
 * \param[in] fwork Target framework which contains current context.
 */
static int bfwork_unlock_user_data(bam_fwork_t *fwork);

/**
 * \brief Get local user data to be used in current context of the framework.
 * This data is only visible for owner thread so must not be used to share data among threads.
 * A thread can get its local data using this function only inside processing function.
 *
 * \param[in] fwork Target framework which contains current context.
 * \param[out] user_data Local user data returned.
 */
static int bfwork_local_user_data(bam_fwork_t *fwork, void **user_data);

/**
 * \brief Set local user data to be used in current context of the framework.
 * This data is only visible for owner thread so must not be used to share data among threads.
 * A thread can set its local data using this function only inside processing function.
 *
 * \param[in] fwork Target framework which contains current context.
 * \param[in] user_data Local user data saved to context.
 */
static int bfwork_local_user_data_set(bam_fwork_t *fwork, void *user_data);

/**
 * \brief Reduce all local data pointers of every threads of a context using custom function.
 * 'cb_reduce' must take as a first argument the destination data and as a second the local data of a thread.
 *
 * \param[in] context Target context to reduce its threads local data.
 * \param[out] reduced Destination data where reduction must be done.
 * \param[in] reducer_function Reduction function.
 */
static int bfwork_context_local_user_data_reduce(bfwork_context_t *context, void *reduced, reducer_function rf);

/**
 * \brief Free all threads local data using custom function.
 * 'cb_free' must take as argument a thread local data that must be free.
 *
 * \param[in] context Target context to destroy its threads local data.
 * \param[in] cb_free Callback free function. OPTIONAL, set NULL to only use default free.
 */
static int bfwork_context_local_user_data_free(bfwork_context_t *context, void (*cb_free)(void *));

/**
 * TIMING
 */

/**
 * \brief Initialize and activate timing of a context.
 * If tag is set, will be used to prepend in output file names.
 *
 * \param[in] context Target context to activate timing.
 * \param[in] tag Tag to be used in output file. OPTIONAL.
 * \param[in] path_folder Folder where store stats. OPTIONAL.
 */
EXTERNC int bfwork_context_init_timing(bfwork_context_t *context, const char *tag, const char *path_folder);

/**
 * \brief Destroy timing of a context.
 *
 * \param[in] context Target context to destroy timing.
 */
EXTERNC void bfwork_context_destroy_timing(bfwork_context_t *context);

/**
 * \brief Output timings in standard output.
 *
 * \param[in] context Target context to output timing.
 */
EXTERNC int bfwork_context_print_times(bfwork_context_t *context);

/**
 * FILTER
 */

/**
 * \brief Filter a read using input filters.
 *
 * \param[in] read Alignment to be filtered.
 * \param[in] filters Filters flags.
 * \return 0 if pass the filters, if not, returns the filter that fails.
 */
static int filter_read(bam1_t *read, uint8_t filters);

/**
 * INLINE FUNCTIONS
 */

/**
 * Lock global user data to be used in current context of the framework.
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

/**
 * Unlock global user data previously locked by 'bfwork_lock_user_data'.
 */
static inline int
bfwork_unlock_user_data(bam_fwork_t *fwork)
{
	assert(fwork);

	//Remove the lock
	omp_unset_lock(&fwork->context->user_data_lock);

	return NO_ERROR;
}

/**
 * Get local user data to be used in current context of the framework.
 */
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

/**
 * Set local user data to be used in current context of the framework.
 */
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

/**
 * Reduce all local data pointers of every threads of a context using custom function.
 */
static inline int
bfwork_context_local_user_data_reduce(bfwork_context_t *context, void *reduced, reducer_function rf)
{
	int i, threads;
	void *data;

	assert(context);
	assert(reduced);
	assert(rf);

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
			rf(reduced, data);
		}
	}

	return NO_ERROR;
}

/**
 * Free all threads local data using custom function.
 */
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

/**
 * Filter a read using input filters.
 */
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
				return FILTER_ZERO_QUAL;
		}

		if(filters & FILTER_DIFF_MATE_CHROM)
		{
			if(read->core.tid != read->core.mtid)
				return FILTER_DIFF_MATE_CHROM;
		}

		if(filters & FILTER_NO_CIGAR)
		{
			if(read->core.n_cigar == 0)
				return FILTER_NO_CIGAR;
		}

		if(filters & FILTER_DEF_MASK)
		{
			if(read->core.flag & BAM_DEF_MASK)
				return FILTER_DEF_MASK;
		}
	}

	return NO_ERROR;
}

#endif /* BFWORK_H_ */
