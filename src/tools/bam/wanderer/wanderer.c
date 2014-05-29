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
	int i;
	bam_region_window_t *window;

	assert(wanderer);

	//Set all to zero
	memset(wanderer, 0, sizeof(bam_wanderer_t));

	//Create region struct
	wanderer->current_region = (bam_region_t *)malloc(sizeof(bam_region_t));
	breg_init(wanderer->current_region);

	//Create windows
	wanderer->windows = (bam_region_window_t *)malloc(WANDERER_WINDOWS_MAX * sizeof(bam_region_window_t));
	for(i = 0; i < WANDERER_WINDOWS_MAX; i++)
	{
		window = &wanderer->windows[i];
		breg_window_init(window);
	}
}

void
bwander_destroy(bam_wanderer_t *wanderer)
{
	int i;
	bam_region_window_t *window;

	assert(wanderer);

	//Region exists?
	if(wanderer->current_region)
	{
		breg_destroy(wanderer->current_region, 0);
		free(wanderer->current_region);
	}

	//Free windows
	for(i = 0; i < WANDERER_WINDOWS_MAX; i++)
	{
		window = &wanderer->windows[i];
		breg_window_destroy(window);
	}
	free(wanderer->windows);
}

void 
bwander_configure(bam_wanderer_t *wanderer, bam_file_t *in_file, bam_file_t *out_file, wanderer_function wf, processor_function pf)
{
	assert(wanderer);
	assert(in_file);
	assert(wf);
	assert(pf);

	//Set I/O
	wanderer->input_file = in_file;
	wanderer->output_file = out_file;		

	//Set functions
	wanderer->wander_f = wf;
	wanderer->processing_f = pf;

	//Logging
	LOG_INFO("Wanderer configured\n");
}

void
bwander_run(bam_wanderer_t *wanderer)
{
	int i, err;
	bam_region_window_t *window;
	size_t total;

	assert(wanderer);
	assert(wanderer->input_file);
	assert(wanderer->current_region);

	//Logging
	LOG_INFO("Wanderer is now running\n");

	//Load first region
	breg_fill(wanderer->current_region, wanderer->input_file);
	
	LOG_INFO_F("\n============== Wandering over region %d:%d-%d with %d reads ============== \n",
			wanderer->current_region->chrom + 1, wanderer->current_region->init_pos + 1,
			wanderer->current_region->end_pos + 1, wanderer->current_region->size);

	//Run loop
	err = wanderer->wander_f(wanderer);
	total = 0;
	while(err != 0)
	{
		//Check error
		if(err <= WANDERER_ERROR)
		{
			LOG_ERROR("In wandering function WANDERER_ERROR\n");
			abort();
		}

		//Run process function
		//TODO for every window
		LOG_INFO_F("Processing %d windows\n", wanderer->windows_l);
		for(i = 0; i < wanderer->windows_l; i++)
		{
			//Get next window
			window = &wanderer->windows[i];

			//Process window
			err = wanderer->processing_f(wanderer, window);
			if(err <= WANDERER_ERROR)
			{
				LOG_ERROR("In process function WANDERER_ERROR\n");
				abort();
			}
		}

		//Clear windows
		bwander_window_clear(wanderer);

		//Counter
		total += wanderer->processed;

		//Write processed
		breg_write_n(wanderer->current_region, wanderer->processed, wanderer->output_file);
		wanderer->processed = 0;

		//Output
		printf("\rProcessed reads: %d", total);
		fflush(stdout);
		//LOG_INFO(str);

		//Load next region
		err = breg_fill(wanderer->current_region, wanderer->input_file);
		if(err)
		{
			if(err == WANDER_READ_EOF)
			{
				//End execution
				break;
			}
			else
			{
				LOG_ERROR("Filling region\n");
				abort();
			}
		}

		//Logging
		LOG_INFO_F("\n============== Wandering over region %d:%d-%d with %d reads ============== \n",
				wanderer->current_region->chrom + 1, wanderer->current_region->init_pos + 1,  wanderer->current_region->end_pos + 1, wanderer->current_region->size);

		//Execute wandering
		err = wanderer->wander_f(wanderer);
		//err = WANDERER_SUCCESS;
	}
	printf("\n");

	//Logging
	LOG_INFO("Wanderer SUCCESS!\n");
}

/**
 * REGISTER WINDOW
 */
int
bwander_window_register(bam_wanderer_t *wanderer, bam_region_window_t *window)
{
	bam_region_window_t *dest_window;

	assert(wanderer);
	assert(wanderer->current_region);
	assert(window);

	//Logging
	LOG_INFO_F("Registering window %d:%d-%d with %d reads\n",
			window->region->chrom + 1, window->init_pos + 1,
			window->end_pos + 1, window->size);

	//Slots free?
	if(wanderer->windows_l >= WANDERER_WINDOWS_MAX)
	{
		LOG_WARN("Wanderer bam_region_window_t buffer is full\n");
		return WANDER_WINDOW_BUFFER_FULL;
	}

	//Duplicate and register
	dest_window = &wanderer->windows[wanderer->windows_l];
	dest_window->init_pos = window->init_pos;
	dest_window->end_pos = window->end_pos;
	dest_window->filter_flags = window->filter_flags;
	dest_window->size = window->size;
	dest_window->region = window->region;
	memcpy(dest_window->filter_reads, window->filter_reads, window->size * sizeof(bam1_t *));

	//Increase windows counter
	wanderer->windows_l++;

	return NO_ERROR;
}

void
bwander_window_clear(bam_wanderer_t *wanderer)
{
	int i;
	bam_region_window_t *window;

	assert(wanderer);
	assert(wanderer->current_region);

	//Logging
	LOG_INFO_F("Clearing %d windows\n", wanderer->windows_l);

	//Iterate windows
	for(i = 0; i < wanderer->windows_l; i++)
	{
		//Get next window
		window = &wanderer->windows[i];

		//Clear window
		breg_window_clear(window);
	}

	//Set size to zero
	wanderer->windows_l = 0;
}

