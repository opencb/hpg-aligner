/*
 * bam_region.c
 *
 *  Created on: 24/04/2014
 *      Author: rmoreno
 */

#include "bam_region.h"

void
breg_init(bam_region_t *region)
{
	assert(region);

	//Set all to zero
	memset(region, 0, sizeof(bam_region_t));

	//Allocate reads array
	region->reads = (bam1_t **)malloc(BAM_REGION_DEFAULT_SIZE * sizeof(bam1_t *));
	assert(region->reads);
	region->max_size = BAM_REGION_DEFAULT_SIZE;
}

void
breg_destroy(bam_region_t *region, int free_bam)
{
	int i;
	bam1_t *read;

	assert(region);

	//Free reads
	if(region->reads)
	{
		if(free_bam)
		{
			//Free reads
			for(i = 0; i < region->size; i++)
			{
				read = region->reads[i];
				bam_destroy1(read);
			}
		}

		//Free memory
		free(region->reads);
		region->reads = NULL;
	}
}

void
breg_fill(bam_region_t *region, bam_file_t *input_file)
{
	int free_slots, i;
	bam1_t **ptr;
	bam1_t *read;
	size_t bytes;
	size_t added;
	size_t init_pos = 0;
	size_t end_pos = 0;
	int chrom = 0;

	assert(region);
	assert(input_file);

	//Get first read
	if(region->size > 0)
	{
		read = region->reads[0];
	}
	else
	{
		//Get first read from file
		read = bam_init1();
		assert(read);
		bytes = bam_read1(input_file->bam_fd, read);

		//Check
		if(bytes == 0)
		{
			//End of file
			bam_destroy1(read);
			return;
		}

		//Add read to region
		region->reads[0] = read;
		region->size = 1;
	}

	//Get free slots
	free_slots = region->max_size - region->size;
	if(free_slots == 0)
	{
		LOG_ERROR("NOT ENOUGHT FREE SLOTS, CANT FILL BUFFER. Aborting...\n");
		abort();
	}
	printf("FREE SLOTS: %d\n", free_slots);

	//Get region chrom and initial position
	init_pos = read->core.pos;
	chrom = read->core.tid;

	//Fill remaining slots
	ptr = region->reads + region->size;
	added = 0;
	for(i = 0; i < free_slots; i++)
	{
		//Read next bam read
		read = bam_init1();
		assert(read);
		bytes = bam_read1(input_file->bam_fd, read);

		//Valid read?
		if(bytes > 0)
		{
			//Same chrom?
			if(read->core.tid != chrom)
			{
				//Chrom changed
				printf("Chrom changed!! %d->%d\n", chrom, read->core.tid);
				break;
			}

			//Add read to region
			ptr[i] = read;
			added++;
		}
		else
		{
			LOG_INFO("End of file\n");
			bam_destroy1(read);
			read = NULL;
			break;
		}
	}
	//Save next read
	region->next_read = read;

	//Update region size
	region->size += added;

	printf("ADDED: %d\n", added);

	//Get end position of region
	end_pos = region->reads[region->size - 1]->core.pos;

	//Update region attributes
	region->init_pos = init_pos;
	region->end_pos = end_pos;
	region->chrom = chrom;
}

//Private compare function
static int compare_pos(const void *item1, const void *item2) {
	const bam1_t *read1 = *(const bam1_t **)item1;
	const bam1_t *read2 = *(const bam1_t **)item2;

	return (int) (read1->core.pos - read2->core.pos);
}

void
breg_write_processed(bam_region_t *region, bam_file_t *output_file)
{
	int i;
	bam1_t *read;

	assert(region);
	assert(output_file);

	//Sort reads
	qsort(region->reads, region->processed, sizeof(void *), compare_pos);

	//Iterate reads
	for(i = 0; i < region->processed; i++)
	{
		//Get read
		read = region->reads[i];
		assert(read);

		//Write to disk
		bam_write1(output_file->bam_fd, read);
		bam_destroy1(read);
	}

	//Update array
	region->size -= region->processed;
	if(region->size > 0)
	{
		memmove(region->reads, region->reads + region->processed, region->size);
	}
	region->processed = 0;

	//Write next read
	if(region->next_read)
	{
		region->reads[region->size] = region->next_read;
		region->next_read = NULL;
		region->size++;
	}
}
