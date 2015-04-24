/*
 * dummy_wander.c
 *
 *  Created on: 20/06/2014
 *      Author: rmoreno
 */

#include "dummy_wander.h"

int
dummy_wanderer(bam_fwork_t *fwork, bam_region_t *region, bam1_t *read)
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
dummy_processor(bam_fwork_t *fwork, bam_region_t *region)
{
	//Do nothing...

	return NO_ERROR;
}

ERROR_CODE
dummy_bam_file(const char *bam_path, const char *ref_path, const char *outbam, const char *stats_path)
{
	//Times


	//Wanderer
	bam_fwork_t fwork;
	bfwork_context_t context1;
	bfwork_context_t context2;
	bfwork_context_t context3;
	bfwork_context_t context4;

	assert(bam_path);
	assert(ref_path);

	//Init wandering
	bfwork_init(&fwork);

	//Create context 1
	{
		bfwork_context_init(&context1,
				(int (*)(void *, bam_region_t *, bam1_t *))dummy_wanderer,
				(int (*)(void *, bam_region_t *))dummy_processor,
				NULL,	//No reduction needed
				NULL
		);

		//bfwork_context_set_output(&context1, NULL);

		//Configure wanderer
		bfwork_configure(&fwork, bam_path, outbam, ref_path, &context1);
	}

	//Create context 2
	{
		bfwork_context_init(&context2,
				(int (*)(void *, bam_region_t *, bam1_t *))dummy_wanderer,
				(int (*)(void *, bam_region_t *))dummy_processor,
				NULL,	//No reduction needed
				NULL
		);

		bfwork_context_set_output(&context2, NULL);

		//Configure wanderer
		bfwork_add_context(&fwork, &context2, FWORK_CONTEXT_SEQUENTIAL);
	}

	//Create context 3
	{
		bfwork_context_init(&context3,
				(int (*)(void *, bam_region_t *, bam1_t *))dummy_wanderer,
				(int (*)(void *, bam_region_t *))dummy_processor,
				NULL,	//No reduction needed
				NULL
		);

		//bfwork_context_set_output(&context3, "test.bam");

		//Configure wanderer
		bfwork_add_context(&fwork, &context3, FWORK_CONTEXT_SEQUENTIAL);
	}

	//Create context 4
	{
		bfwork_context_init(&context4,
				(int (*)(void *, bam_region_t *, bam1_t *))dummy_wanderer,
				(int (*)(void *, bam_region_t *))dummy_processor,
				NULL,	//No reduction needed
				NULL
		);

		//This should not be effective because is the last context
		bfwork_context_set_output(&context4, "phantom.bam");

		//Configure wanderer
		bfwork_add_context(&fwork, &context4, FWORK_CONTEXT_SEQUENTIAL);
	}

	//Run wander
	bfwork_run(&fwork);

	//Destroy wanderer
	bfwork_destroy(&fwork);

	//Destroy context
	bfwork_context_destroy(&context1);
	bfwork_context_destroy(&context2);
	bfwork_context_destroy(&context3);
	bfwork_context_destroy(&context4);

	return NO_ERROR;
}
