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
	bam_region_t *region;

	assert(wanderer);

	//Set all to zero
	memset(wanderer, 0, sizeof(bam_wanderer_t));

	//Create regions
	wanderer->regions = (bam_region_t **)malloc(WANDERER_REGIONS_MAX * sizeof(bam_region_t *));
	/*for(i = 0; i < WANDERER_REGIONS_MAX; i++)
	{
		region = wanderer->regions[i];
		breg_init(region);
	}*/
}

void
bwander_destroy(bam_wanderer_t *wanderer)
{
	int i;
	bam_region_t *region;

	assert(wanderer);

	//Regions exists?
	if(wanderer->regions)
	{
		for(i = 0; i < wanderer->regions_l; i++)
		{
			//Get region
			region = wanderer->regions[i];
			breg_destroy(region, 0);
			free(region);
		}
		free(wanderer->regions);
	}
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

int
bwander_run(bam_wanderer_t *wanderer)
{
	int i, err;
	bam_region_t *region;
	bam_region_t *current_region;
	size_t total;
	int bytes;
	int while_cond = 1;

	bam1_t *read;

	assert(wanderer);
	assert(wanderer->input_file);
	assert(wanderer->regions);

	//Logging
	LOG_INFO("Wanderer is now running\n");

	//Allocate current region
	current_region = (bam_region_t *)malloc(sizeof(bam_region_t));
	breg_init(current_region);

	//Get first read from file
	read = bam_init1();
	assert(read);
	bytes = bam_read1(wanderer->input_file->bam_fd, read);

	while(bytes > 0)
	{
		//Region read
		while(wanderer->regions_l < 6)
		{
			//Wander this read
			err = wanderer->wander_f(wanderer, current_region, read);
			switch(err)
			{
			case WANDER_READ_FILTERED:
				//This read dont pass the filters
			case NO_ERROR:
				//Add read to region
				current_region->reads[current_region->size]  = read;
				current_region->size++;

				//Get next read from file
				read = bam_init1();
				assert(read);
				bytes = bam_read1(wanderer->input_file->bam_fd, read);
				break;

			case WANDER_REGION_CHANGED:
				//The region have changed

				//Add region to wanderer regions
				bwander_region_insert(wanderer, current_region);

				//Create new current region
				current_region = (bam_region_t *)malloc(sizeof(bam_region_t));
				breg_init(current_region);

				//Reprocess this read with new region
				break;

			default:
				//Unknown error
				LOG_ERROR_F("Wanderer fails with error code: %d\n", err);
				return err;
			}
		}

		//Region process
		{
			//Iterate regions
			for(i = 0; i < wanderer->regions_l; i++)
			{
				//Get next region
				region = wanderer->regions[i];

				//Process region
				wanderer->processing_f(wanderer, region);
			}

		} //End region process

		//Region writing (in region order)
		{
			//Iterate regions
			for(i = 0; i < wanderer->regions_l; i++)
			{
				//Get next region
				region = wanderer->regions[i];

				//This region is processed?
				//TODO if()
				{
					//Write region to disk
					breg_write_n(region, region->size, wanderer->output_file);
					wanderer->regions_l--;
				}
				//else break;
			}

		} //End region writing
	}

	//Check read error
	if(bytes <= 0)
	{
		//Destroy bam
		bam_destroy1(read);

		//End of file
		if(bytes == -1)
		{
			LOG_INFO("EOF\n");
			return WANDER_READ_EOF;
		}
		else
		{
			LOG_INFO("TRUNCATED\n");
			return WANDER_READ_TRUNCATED;
		}
	}

	//Logging
	LOG_INFO("Wanderer SUCCESS!\n");

	return NO_ERROR;
}

/**
 * REGISTER WINDOW
 */
/*int
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
}*/

/*void
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
}*/

