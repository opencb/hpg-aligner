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

#include "bam_file_ops.h"

/***********************************************
 * WANDERERS DECLARATIONS
 **********************************************/

/**
 * Realign wanderer
 */
int realigner_wanderer(bam_fwork_t *fwork, bam_region_t *region, bam1_t *read);
int realigner_processor(bam_fwork_t *fwork, bam_region_t *region);

/**
 * Recalibrate wanderer
 */
int recalibrate_wanderer(bam_fwork_t *fwork, bam_region_t *region, bam1_t *read);
int recalibrate_collect_processor(bam_fwork_t *fwork, bam_region_t *region);
int recalibrate_recalibrate_processor(bam_fwork_t *fwork, bam_region_t *region);
void recalibrate_reduce_data(void *dest, void *data);
void recalibrate_destroy_data(void *data);

/***********************************************
 * FRAMEWORK REALIGNER
 **********************************************/

/**
 * Realign around indels BAM file
 */
ERROR_CODE
alig_bam_file(const char *bam_path, const char *ref_path, const char *outbam, const char *stats_path)
{
	//Times
	double times;

	//Wanderer
	bam_fwork_t fwork;
	bfwork_context_t context;

	assert(bam_path);
	assert(ref_path);

	//Init wandering
	bfwork_init(&fwork);

	//Create realigner context
	bfwork_context_init(&context,
			(int (*)(void *, bam_region_t *, bam1_t *))realigner_wanderer,
			(int (*)(void *, bam_region_t *))realigner_processor,
			NULL,	//No reduction needed
			NULL
	);

#ifdef D_TIME_DEBUG
	//Init timing
	bfwork_context_init_timing(&context, "realigner", stats_path);
#endif

	//Configure wanderer
	bfwork_configure(&fwork, bam_path, outbam, ref_path, &context);

	//Run wander
	bfwork_run(&fwork);

#ifdef D_TIME_DEBUG
	//Print times
	bfwork_context_print_times(&context);

	//Destroy wanderer time
	bfwork_context_destroy_timing(&context);
#endif

	//Destroy wanderer
	bfwork_destroy(&fwork);

	//Destroy context
	bfwork_context_destroy(&context);

	return NO_ERROR;
}

/***********************************************
 * FRAMEWORK RECALIBRATOR
 **********************************************/

/**
 * Recalibrate BAM file
 */
ERROR_CODE
recal_bam_file(uint8_t flags, const char *bam_path, const char *ref, const char *data_file, const char *info_file, const char *outbam, int cycles, const char *stats_path)
{
	int bytes;

	//Data
	recal_info_t info;
	U_CYCLES aux_cycles;
	int cycles_param;

	//Times
	double times;

	//Wanderer
	bam_fwork_t fwork;
	bfwork_context_t collect_context;
	bfwork_context_t recal_context;

	assert(bam_path);

	//Get phases
	if((flags & RECALIBRATE_RECALIBRATE) && !(flags & RECALIBRATE_COLLECT))
	{
		//Second phase only, needs data file
		assert(data_file);
	}

	//Create new data
	aux_cycles = cycles;
	recal_init_info(aux_cycles, &info);

	//Init wandering
	bfwork_init(&fwork);

	//Configure framework
	bfwork_configure(&fwork, bam_path, outbam, ref, NULL);

	//Collection is needed?
	if(flags & RECALIBRATE_COLLECT)
	{
		assert(ref);

		//Create data collection context
		bfwork_context_init(&collect_context,
				(int (*)(void *, bam_region_t *, bam1_t *))recalibrate_wanderer,
				(int (*)(void *, bam_region_t *))recalibrate_collect_processor,
				(int (*)(void *, void *))recalibrate_reduce_data,
				&info
		);

#ifdef D_TIME_DEBUG
		//Init timing
		bfwork_context_init_timing(&collect_context, "collect", stats_path);
#endif

		//Set user data
		cycles_param = cycles;
		bfwork_context_set_user_data(&collect_context, &cycles_param);

		//Add context for data collection
		bfwork_add_context(&fwork, &collect_context, FWORK_CONTEXT_SEQUENTIAL);

		printf("Cycles: %d\n",cycles);
	}

	//Recalibration is needed?
	if(flags & RECALIBRATE_RECALIBRATE)
	{
		//Create recalibration context
		bfwork_context_init(&recal_context,
						(int (*)(void *, bam_region_t *, bam1_t *))recalibrate_wanderer,
						(int (*)(void *, bam_region_t *))recalibrate_recalibrate_processor,
						NULL,	//No reduction needed
						NULL
		);

#ifdef D_TIME_DEBUG
		//Init timing
		bfwork_context_init_timing(&recal_context, "recalibrate", stats_path);
#endif

		//Set context user data
		bfwork_context_set_user_data(&recal_context, &info);

		//Previous collect?
		if(!(flags & RECALIBRATE_COLLECT))
		{
			//No previous collect, load data file
			if(data_file)
			{
				//Load data from disk
				recal_load_recal_info(data_file, &info);
			}
			else
			{
				//If only recalibrate input data must be present
				LOG_FATAL("Only recalibration phase 2 requested but no input data specified\n");
			}
		}

		//Add context for recalibration
		bfwork_add_context(&fwork, &recal_context, FWORK_CONTEXT_SEQUENTIAL);
	}

	//Run wander
	bfwork_run(&fwork);

	//Collection last operations
	if(flags & RECALIBRATE_COLLECT)
	{
		//Save data file
		if(data_file)
		{
			recal_save_recal_info(&info, data_file);
		}

		//Save info file
		if(info_file)
		{
			recal_fprint_info(&info, info_file);
		}

		bfwork_context_local_user_data_free(&collect_context, recalibrate_destroy_data);

#ifdef D_TIME_DEBUG
		//Print times
		bfwork_context_print_times(&collect_context);

		//Destroy wanderer time
		bfwork_context_destroy_timing(&collect_context);
#endif

		//Destroy context
		bfwork_context_destroy(&collect_context);
	}

	//Recalibration last operations
	if(flags & RECALIBRATE_RECALIBRATE)
	{
		//Free local data
		bfwork_context_local_user_data_free(&recal_context, recalibrate_destroy_data);

#ifdef D_TIME_DEBUG
		//Print times
		bfwork_context_print_times(&recal_context);

		//Destroy wanderer time
		bfwork_context_destroy_timing(&recal_context);
#endif

		//Destroy context
		bfwork_context_destroy(&recal_context);
	}

	//Destroy wanderer
	bfwork_destroy(&fwork);

	//Free data memory
	recal_destroy_info(&info);

	return NO_ERROR;
}

/**
 * Realign and recalibrate BAM file
 */
ERROR_CODE
alig_recal_bam_file(const char *bam_path, const char *ref_path, const char *data_file, const char *info_file, const char *outbam, int cycles, const char *stats_path)
{
	int bytes;

	//Data
	recal_info_t info;
	U_CYCLES aux_cycles;
	int cycles_param;

	//Times
	double times;

	//Wanderer
	bam_fwork_t fwork;
	bfwork_context_t realign_context;
	bfwork_context_t recal_context;

	assert(bam_path);
	assert(ref_path);

	//Create new data
	aux_cycles = cycles;
	recal_init_info(aux_cycles, &info);

	//Init wandering
	bfwork_init(&fwork);

	//Configure framework
	bfwork_configure(&fwork, bam_path, outbam, ref_path, NULL);

	//Create data realign and collection context
	bfwork_context_init(&realign_context,
			(int (*)(void *, bam_region_t *, bam1_t *))realigner_wanderer,
			(int (*)(void *, bam_region_t *))realigner_processor,
			(int (*)(void *, void *))recalibrate_reduce_data,
			&info
	);
	bfwork_context_add_proc(&realign_context, (int (*)(void *, bam_region_t *))recalibrate_collect_processor);
	bfwork_context_set_output(&realign_context, NULL);

	//Create recalibration context
	bfwork_context_init(&recal_context,
					(int (*)(void *, bam_region_t *, bam1_t *))recalibrate_wanderer,
					(int (*)(void *, bam_region_t *))recalibrate_recalibrate_processor,
					NULL,	//No reduction needed
					NULL
	);

#ifdef D_TIME_DEBUG
	//Init timing
	bfwork_context_init_timing(&realign_context, "aligcollect", stats_path);
	bfwork_context_init_timing(&recal_context, "recalibrate", stats_path);
#endif

	//Set user data
	cycles_param = cycles;
	bfwork_context_set_user_data(&realign_context, &cycles_param);
	bfwork_context_set_user_data(&recal_context, &info);
	printf("Cycles: %d\n",cycles);

	//Add context for recalibration
	bfwork_add_context(&fwork, &realign_context, FWORK_CONTEXT_SEQUENTIAL);
	bfwork_add_context(&fwork, &recal_context, FWORK_CONTEXT_SEQUENTIAL);

	//Run wander
	bfwork_run(&fwork);

	//Save data file
	if(data_file)
	{
		recal_save_recal_info(&info, data_file);
	}

	//Save info file
	if(info_file)
	{
		recal_fprint_info(&info, info_file);
	}

	//Free local data
	bfwork_context_local_user_data_free(&realign_context, recalibrate_destroy_data);
	bfwork_context_local_user_data_free(&recal_context, recalibrate_destroy_data);

#ifdef D_TIME_DEBUG
	//Print times
	bfwork_context_print_times(&realign_context);
	bfwork_context_print_times(&recal_context);

	//Destroy wanderer time
	bfwork_context_destroy_timing(&realign_context);
	bfwork_context_destroy_timing(&recal_context);
#endif

	//Destroy context
	bfwork_context_destroy(&realign_context);
	bfwork_context_destroy(&recal_context);

	//Destroy wanderer
	bfwork_destroy(&fwork);

	//Free data memory
	recal_destroy_info(&info);

	return NO_ERROR;
}

/***********************************************
 * WANDERERS
 **********************************************/

/**
 * Realign wanderer
 */

int
realigner_wanderer(bam_fwork_t *fwork, bam_region_t *region, bam1_t *read)
{
	int i, err;

	//Current region
	size_t aux_init_pos;
	size_t aux_end_pos;
	size_t read_pos;

	assert(fwork);
	assert(region);
	assert(read);

	//Filter read
	if(filter_read(read, FILTER_ZERO_QUAL | FILTER_DIFF_MATE_CHROM | FILTER_NO_CIGAR | FILTER_DEF_MASK))
	{
		//Read is not valid for process
		return WANDER_READ_FILTERED;
	}

	//Get read position
	read_pos = read->core.pos;

	//Inside this region?
	if(region->end_pos != SIZE_MAX)
	{
		if(	region->chrom != read->core.tid
				|| region->end_pos < read->core.pos)
		{
			//Not in window region
			return WANDER_REGION_CHANGED;
		}
	}

	//Get interval for this alignment
	err = region_get_from_bam1(read, &aux_init_pos, &aux_end_pos);
	if(err)
	{
		LOG_ERROR_F("Trying to get region from invalid read: %s\n", bam1_qname(read));
		return INVALID_INPUT_BAM;
	}

	//This alignment have an interval?
	if(aux_init_pos != SIZE_MAX && aux_end_pos != SIZE_MAX)
	{
		//Interval found

		//Update region chrom
		region->chrom = read->core.tid;

		//Update region start position
		if(region->init_pos == SIZE_MAX || region->init_pos > aux_init_pos)
		{
			region->init_pos = aux_init_pos;
		}

		//Update region end position
		if(region->end_pos == SIZE_MAX || region->end_pos < aux_end_pos)
		{
			region->end_pos = aux_end_pos;
		}
	}

	return NO_ERROR;
}

int
realigner_processor(bam_fwork_t *fwork, bam_region_t *region)
{
	int err;
	alig_context_t context;
	bam1_t **v_reads;
	size_t v_reads_l;

	//Create contexts
	omp_set_lock(&fwork->reference_lock);
	alig_init(&context, fwork->reference, ALIG_LEFT_ALIGN | ALIG_REFERENCE_PRIORITY);
	omp_unset_lock(&fwork->reference_lock);

	//Load region reads in aligner context
	v_reads = region->reads;
	v_reads_l = region->size;
	err = alig_region_next(v_reads, v_reads_l, 1, &context);
	if(err && err != ALIG_INCOMPLETE_INTERVAL)
	{
		LOG_ERROR_F("Cannot obtain next region in aligner, error code = %d\n", err);
		return err;
	}

	//Load reference for this region
	omp_set_lock(&fwork->reference_lock);
	err = alig_region_load_reference(&context);
	omp_unset_lock(&fwork->reference_lock);
	if(err)
	{
		LOG_ERROR_F("Loading reference sequence, error code = %d\n", err);
		return err;
	}

	//Obtain haplotypes
	err = alig_region_haplotype_process(&context);
	if(err)
	{
		LOG_ERROR_F("Obtaining haplotypes, error code = %d\n", err);
		return err;
	}

	//Realign
	err = alig_region_indel_realignment(&context);
	if(err)
	{
		LOG_ERROR_F("Realigning, error code = %d\n", err);
		return err;
	}

	//Destroy region
	alig_destroy(&context);

	return NO_ERROR;
}

/**
 * Recalibrate wanderer
 */

int
recalibrate_wanderer(bam_fwork_t *fwork, bam_region_t *region, bam1_t *read)
{
	assert(fwork);
	assert(region);
	assert(read);

	//Filter read
	if(filter_read(read, FILTER_ZERO_QUAL | FILTER_DIFF_MATE_CHROM | FILTER_NO_CIGAR | FILTER_DEF_MASK))
	{
		//Read is not valid for process
		return WANDER_READ_FILTERED;
	}

	//Update region bounds
	if(region->init_pos > read->core.pos)
	{
		region->init_pos = read->core.pos;
		region->chrom = read->core.tid;
	}
	if(region->end_pos < read->core.pos)
	{
		region->end_pos = read->core.pos;
	}

	return NO_ERROR;
}

int
recalibrate_collect_processor(bam_fwork_t *fwork, bam_region_t *region)
{
	int err, i;
	recal_info_t *data;
	recal_data_collect_env_t *collect_env;
	bam1_t *read;
	size_t *cycles;

	//Get data
	bfwork_local_user_data(fwork, (void **)&data);
	if(data == NULL)
	{
		//Local data is not initialized
		data = (recal_info_t *)malloc(sizeof(recal_info_t));

		//Lock cycles
		bfwork_lock_user_data(fwork, (void **)&cycles);
		recal_init_info(*cycles, data);
		bfwork_unlock_user_data(fwork);

		//Set local data
		bfwork_local_user_data_set(fwork, data);
	}

	//Initialize get data environment
	collect_env = (recal_data_collect_env_t *) malloc(sizeof(recal_data_collect_env_t));
	recal_get_data_init_env(data->num_cycles, collect_env);

	//Obtain data from all reads in region
	for(i = 0; i < region->size; i++)
	{
		//Get next read
		read = region->reads[i];
		assert(read);

		//Get data
		omp_set_lock(&fwork->reference_lock);
		recal_get_data_from_bam_alignment(read, fwork->reference, data, collect_env);
		omp_unset_lock(&fwork->reference_lock);
	}

	//Destroy environment
	recal_get_data_destroy_env(collect_env);

	return NO_ERROR;
}

int
recalibrate_recalibrate_processor(bam_fwork_t *fwork, bam_region_t *region)
{
	int err, i;
	recal_info_t *gdata;
	recal_info_t *data;
	recal_recalibration_env_t *recal_env;
	bam1_t *read;

	//Get data
	bfwork_local_user_data(fwork, (void **)&data);
	if(data == NULL)
	{
		//Local data is not initialized
		data = (recal_info_t *)malloc(sizeof(recal_info_t));

		//Lock data
		bfwork_lock_user_data(fwork, (void **)&gdata);

		//Init struct
		recal_init_info(gdata->num_cycles, data);

		//Clone global data
		if(recal_reduce_info(data, gdata))	//Copy data into local
		{
			LOG_FATAL_F("In local recalibration data copy\nLOCAL: Cycles: %d, Min Q: %d, Num Q: %d, Dinuc: %d\nGLOBAL: Cycles: %d, Min Q: %d, Num Q: %d, Dinuc: %d\n",
					data->num_cycles, data->min_qual, data->num_quals, data->num_dinuc, gdata->num_cycles, gdata->min_qual, gdata->num_quals, gdata->num_dinuc);
		}

		//Recalculate deltas (reduce only merge bases and misses)
		recal_calc_deltas(data);

		bfwork_unlock_user_data(fwork);

		//Set local data
		bfwork_local_user_data_set(fwork, data);
	}

	//Initialize get data environment
	recal_env = (recal_recalibration_env_t *) malloc(sizeof(recal_recalibration_env_t));
	recal_recalibration_init_env(data->num_cycles, recal_env);

	//Recalibrate region
	for(i = 0; i < region->size; i++)
	{
		//Get next read
		read = region->reads[i];
		assert(read);

		//Recalibrate read
		recal_recalibrate_alignment(read, data, recal_env);
	}

	//Destroy environment
	recal_recalibration_destroy_env(recal_env);

	return NO_ERROR;
}

void
recalibrate_reduce_data(void *dest, void *data)
{
	recal_info_t *data_ptr, *dest_ptr;
	assert(data);
	assert(dest);

	//Cast pointers
	data_ptr = (recal_info_t *)data;
	dest_ptr = (recal_info_t *)dest;

	//Combine
	recal_reduce_info(dest_ptr, data_ptr);

	//Delta processing
	recal_calc_deltas(dest_ptr);
	//printf("Estimated %.2f \tEmpirical %.2f \t TotalDelta %.2f\n", dest_ptr->total_estimated_Q, dest_ptr->total_delta + dest_ptr->total_estimated_Q, dest_ptr->total_delta);
}

void
recalibrate_destroy_data(void *data)
{
	recal_info_t *aux;

	assert(data);

	aux = (recal_info_t *)data;
	recal_destroy_info(aux);
}

