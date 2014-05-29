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
	memset(region->reads, 0, BAM_REGION_DEFAULT_SIZE * sizeof(bam1_t *));
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

int
breg_fill(bam_region_t *region, bam_file_t *input_file)
{
	int free_slots, i, err;
	bam1_t *read;
	int bytes;
	size_t added;
	int64_t end_pos = -1;

	assert(region);
	assert(input_file);

	added = region->size;

	//Get free slots
	free_slots = region->max_size - region->size;
	if(free_slots == 0)
	{
		LOG_FATAL("NOT ENOUGH FREE SLOTS, CANT FILL BUFFER. Aborting...\n");
	}
	//printf("FREE SLOTS: %d\n", free_slots);

	//Add next read
	if(region->next_read != NULL)
	{
		//LOG_INFO_F("NEXT at %d\n", region->size);
		region->reads[region->size] = region->next_read;
		region->next_read = NULL;
		region->size++;
	}

	//Get first read
	if(region->size > 0)
	{
		//LOG_INFO_F("First read of %d reads\n", region->size);
		read = region->reads[0];
	}
	else
	{
		//Get first read from file
		read = bam_init1();
		assert(read);
		bytes = bam_read1(input_file->bam_fd, read);

		//Check
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

		//Add read to region
		region->reads[0] = read;
		region->size++;
	}
	assert(read);

	//Set region init locus and chrom
	region->init_pos = read->core.pos;
	region->chrom = read->core.tid;

	//Fill remaining slots
	err = NO_ERROR;
	for(i =  region->size; i < region->max_size; i++)
	{
		//Read next bam read
		read = bam_init1();
		assert(read);
		bytes = bam_read1(input_file->bam_fd, read);

		//Valid read?
		if(bytes > 0)
		{
			//Different chrom?
			if(read->core.tid != region->chrom)
			{
				LOG_INFO_F("Different chrom - Reg:%d - Read:%d\n", region->chrom+1, read->core.tid+1);

				//Save next read
				region->next_read = read;

				break;
			}

			//Add read to region
			region->reads[i] = read;
			region->size++;
		}
		else
		{
			//Destroy bam
			bam_destroy1(read);
			read = NULL;

			if(bytes == -1)
			{
				LOG_INFO("EOF\n");
				err = WANDER_READ_EOF;
			}
			else
			{
				LOG_INFO("TRUNCATED\n");
				err = WANDER_READ_TRUNCATED;
			}

			break;
		}
	}

	//Update region size
	added = region->size - added;

	//printf("ADDED: %d - SIZE: %d\n", added, region->size);

	//Set region end locus
	region->end_pos = region->reads[region->size - 1]->core.pos;

	return err;
}

//Private compare function
static int compare_pos(const void *item1, const void *item2) {
	const bam1_t *read1 = *(const bam1_t **)item1;
	const bam1_t *read2 = *(const bam1_t **)item2;

	return (int) (read1->core.pos - read2->core.pos);
}

void
breg_write_n(bam_region_t *region, size_t n, bam_file_t *output_file)
{
	int i;
	bam1_t *read;

	assert(region);
	assert(output_file);

	LOG_INFO_F("Writing %d processed reads\n", n);

	//Sort reads
	qsort(region->reads, n, sizeof(void *), compare_pos);

	//Iterate reads
	for(i = 0; i < n; i++)
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
	region->size -= n;
	if(region->size > 0)
	{
		//LOG_INFO_F("Moving %d from index %d, new size: %d\n", region->size, region->processed, region->size);

		memmove(region->reads, region->reads + n, region->size * sizeof(bam1_t *));
		memset(region->reads + region->size, 0, (region->max_size - region->size) * sizeof(bam1_t *));
	}
}

void
breg_load_window(bam_region_t *region, size_t init_pos, size_t end_pos, uint8_t filters, bam_region_window_t *window)
{
	assert(region);
	assert(window);

	//Setup window
	window->region = region;
	window->init_pos = init_pos;
	window->end_pos = end_pos;

	//Obtain filtered reads
	breg_window_filter(window, filters);
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

	//Allocate windows array
	window->filter_reads = (bam1_t **)malloc(BAM_REGION_DEFAULT_SIZE * sizeof(bam1_t *));
	assert(window->filter_reads);
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
	assert(window->filter_reads);

	//Clean filter
	window->size = 0;

	//Iterate reads
	region = window->region;
	reads_l = region->size;
	for(i = 0; i < reads_l; i++)
	{
		//Check size
		if(window->size >= BAM_REGION_DEFAULT_SIZE)
		{
			LOG_WARN("Window bam1_t buffer is full\n");
			break;
		}

		//Get next read
		read = region->reads[i];
		if(read == NULL)
		{
			printf("ERROR en lectura: %d\n", i);
			continue;
		}

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
		if(window->region->chrom != EMPTY_CHROM)
		{
			//Is in window region?
			if(	window->region->chrom != read->core.tid
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

void
breg_window_clear(bam_region_window_t *window)
{
	assert(window);

	//Set all to zero
	window->size = 0;
	window->filter_flags = 0;
	window->init_pos = 0;
	window->end_pos = 0;
	window->region = NULL;
}
