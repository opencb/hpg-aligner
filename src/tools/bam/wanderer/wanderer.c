/*
 * wanderer.c
 *
 *  Created on: 24/04/2014
 *      Author: rmoreno
 */

#include "wanderer.h"

void
bwander_init(bam_wanderer_t *wanderer)
{
	assert(wanderer);

	//Set all to zero
	memset(wanderer, 0, sizeof(bam_wanderer_t));

	//Create region struct
	wanderer->current_region = (bam_region_t *)malloc(sizeof(bam_region_t));
	breg_init(wanderer->current_region);
}

void
bwander_destroy(bam_wanderer_t *wanderer)
{
	assert(wanderer);

	//Region exists?
	if(wanderer->current_region)
	{
		breg_destroy(wanderer->current_region, 0);
		free(wanderer->current_region);
	}
}

void 
bwander_configure(bam_wanderer_t *wanderer, bam_file_t *in_file, bam_file_t *out_file, wandering_function f)
{
	assert(wanderer);
	assert(in_file);
	assert(f);

	//Set I/O
	wanderer->input_file = in_file;
	wanderer->output_file = out_file;		

	//Set wandering function
	wanderer->wander_f = f;

	//Logging
	LOG_INFO("Wanderer configured\n");
}

void
bwander_run(bam_wanderer_t *wanderer)
{
	int err;
	char str[100];

	assert(wanderer);
	assert(wanderer->input_file);
	assert(wanderer->current_region);

	//Logging
	LOG_INFO("Wanderer is now running\n");

	//Load first region
	breg_fill(wanderer->current_region, wanderer->input_file);
	breg_filter(wanderer->current_region, FILTER_ZERO_QUAL | FILTER_DIFF_MATE_CHROM | FILTER_NO_CIGAR | FILTER_DEF_MASK);
	
	sprintf(str, "Wandering over region %d:%d-%d with %d reads\n",
			wanderer->current_region->chrom + 1, wanderer->current_region->init_pos + 1,
			wanderer->current_region->end_pos + 1, wanderer->current_region->size);
	printf(str);

	//Run loop
	err = wanderer->wander_f(wanderer->current_region);
	while(err != 0)
	{
		//Check error
		if(err <= WANDERER_ERROR)
		{
			LOG_ERROR("In wandering function WANDERER_ERROR\n");
		}

		//Write processed
		breg_write_processed(wanderer->current_region, wanderer->output_file);

		//Load next region
		breg_fill(wanderer->current_region, wanderer->output_file);
		breg_filter(wanderer->current_region, FILTER_ZERO_QUAL | FILTER_DIFF_MATE_CHROM | FILTER_NO_CIGAR | FILTER_DEF_MASK);

		//Logging
		sprintf(str, "Wandering over region %d:%d-%d with %d reads\n",
				wanderer->current_region->chrom, wanderer->current_region->init_pos,
				wanderer->current_region->end_pos, wanderer->current_region->size);
		printf(str);
		printf("%d:%d-%d - Proc: %d - Size: %d\n", wanderer->current_region->chrom,
				 wanderer->current_region->init_pos,  wanderer->current_region->end_pos,
				 wanderer->current_region->processed,  wanderer->current_region->size);

		//Execute wandering
		//err = wanderer->wander_f(wanderer->current_region);
		err = WANDERER_SUCCESS;
	}

	//Logging
	LOG_INFO("Wanderer end execution\n");
}

