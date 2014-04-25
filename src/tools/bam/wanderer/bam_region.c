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

	assert(region);
	assert(input_file);

	//Add next read
	if(region->next_read != NULL)
	{
		region->reads[region->size] = region->next_read;
		region->next_read = NULL;
		region->size++;
	}

	//Get first read
	if(region->size > 0)
	{
		printf("First read\n");
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
	assert(read);

	//Get free slots
	free_slots = region->max_size - region->size;
	if(free_slots == 0)
	{
		LOG_ERROR("NOT ENOUGHT FREE SLOTS, CANT FILL BUFFER. Aborting...\n");
		abort();
	}
	printf("FREE SLOTS: %d\n", free_slots);

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

	printf("ADDED: %d - SIZE: %d\n", added, region->size);
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
		region->reads[i] = NULL;
	}

	//Update array
	region->size -= region->processed;
	if(region->size > 0)
	{
		printf("Moving %d from index %d\n", region->size, region->processed);
		memmove(region->reads, region->reads + region->processed, region->size);
		memset(region->reads + region->size, 0, region->max_size - region->size);
	}
	region->processed = 0;
}

void
breg_load_window(bam_region_t *region, int chrom, size_t init_pos, size_t end_pos, uint8_t filters, bam_region_window_t *window)
{
	assert(region);
	assert(window);

	//Setup window
	window->region = region;
	window->chrom = chrom;
	window->init_pos = init_pos;
	window->end_pos = end_pos;

	//Obtain filtered reads
	breg_window_filter(window, filters);

	printf("Window loaded, %d reads\n", window->size);
	chrom != EMPTY_CHROM ? printf("Region %d:%d - %d\n", chrom, init_pos, end_pos) : printf("Whole region\n");
}

/**
 * WINDOW OPERATIONS
 */
void
breg_window_init(bam_region_window_t *window)
{
	assert(window);

	//Set all to zero
	memset(window, 0, sizeof(bam_region_window_t));

	//Invalid chrom
	window->chrom = EMPTY_CHROM;
}

void
breg_window_destroy(bam_region_window_t *window)
{
	assert(window);

	//Free filtered reads
	if(window->filter_reads)
	{
		free(window->filter_reads);
	}
}

void
breg_window_filter(bam_region_window_t *window, uint8_t filters)
{
	int i;
	bam1_t *read;
	size_t reads_l;
	bam_region_t *region;

	assert(window);
	assert(window->region);

	//Clean filter
	window->size = 0;

	//Allocate if not
	if(!window->filter_reads)
	{
		window->filter_reads = (bam1_t **)malloc(BAM_REGION_DEFAULT_SIZE * sizeof(bam1_t *));
		assert(window->filter_reads);
	}

	//Iterate reads
	region = window->region;
	reads_l = region->size;
	for(i = 0; i < reads_l; i++)
	{
		//Get next read
		read = region->reads[i];
		assert(read);

		//Filter read
		{
			if(filters & FILTER_ZERO_QUAL)
			{
				if(read->core.qual == 0)
					continue;
			}

			if(filters & FILTER_DIFF_MATE_CHROM)
			{
				if(read->core.tid != read->core.mtid)
					continue;
			}

			if(filters & FILTER_NO_CIGAR)
			{
				if(read->core.n_cigar == 0)
					continue;
			}

			if(filters & FILTER_DEF_MASK)
			{
				if(read->core.flag & BAM_DEF_MASK)
					continue;
			}
		}

		//Check window region
		if(window->chrom != EMPTY_CHROM)
		{
			//Is in window region?
			if(	window->chrom != read->core.tid
					|| window->init_pos > read->core.pos + read->core.l_qseq
					|| window->end_pos < read->core.pos)
			{
				//Not in window region
				continue;
			}
		}

		//Read is valid and inside region
		window->filter_reads[window->size] = read;
		window->size++;
	}
}
