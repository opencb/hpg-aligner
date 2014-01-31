/*
 * region.c
 *
 *  Created on: Dec 10, 2013
 *      Author: rmoreno
 */

#include "alig_region.h"


/**
 * REGION STRUCT FUNCTIONS
 */

int
region_create(alig_region_t *region)
{
	assert(region);

	region->chrom = -1;
	region->end_pos = 0;
	region->start_pos = 0;
}

int
region_destroy(alig_region_t *region)
{
	assert(region);

	free(region);
}

/**
 * REGION TABLE BY CHROM
 */

int
region_table_create(alig_region_table_t *table)
{
	int i, ret;
	khash_t(32) *h;
	//linked_list_t *aux_list;
	//khiter_t iter;

	assert(table);
	h = table->chroms;
	assert(h);
	h = kh_init(32);
	table->chroms = h;

	//Create chrom lists
	/*for(i = 0; i < 32; i++)
	{
		//Create list
		aux_list = linked_list_new(0);

		//Insert list in hash
		iter = kh_put(32, table->chroms, i, &ret);
		assert(ret);
		kh_value(table->chroms, iter) = aux_list;
	}*/
}

int
region_table_destroy(alig_region_table_t *table)
{
	int i;
	linked_list_t *aux_list;
	khash_t(32) *h = table->chroms;
	//khiter_t iter;

	assert(table);
	assert(h);

	//printf("Iterating...\n");
	//Free table items
	for (i = kh_begin(h); i != kh_end(h); ++i)
	{
		if (kh_exist(h, i))
		{
			//printf("Iterator: %d\n", i);
			//Get list from hash for i chrom
			aux_list = kh_value(h, i);

			//Free list calling region_destroy in its regions
			linked_list_free(aux_list,(void *)region_destroy);
		}
	}

	//Free hash table
	kh_destroy(32, h);
}

int
region_table_insert(alig_region_table_t *table, alig_region_t *region)
{
	khiter_t iter;
	linked_list_t *chrom_list = NULL;
	khash_t(32) *h = table->chroms;

	alig_region_t *list_region;
	linked_list_iterator_t *region_it = NULL;
	int ret, found;

	assert(table);
	assert(region);
	assert(h);

	//Get list for region chromosome
	/*iter = kh_get(32, table->chroms, region->chrom);
	chrom_list = kh_value(table->chroms, iter);
	assert(chrom_list);*/

	//Get chrom list
	iter = kh_get(32, h, region->chrom);
	if(iter == kh_end(h))
	{
		//Create list
		chrom_list = linked_list_new(0);
		iter =kh_put(32, h, region->chrom, &ret);
		//assert(ret == 0);
		kh_value(h, iter) = chrom_list;
	}
	else
	{
		//Get chrom list from hash
		chrom_list = kh_value(h, iter);
	}

	assert(chrom_list);

	//Insert region in list
	{
		//Search for existing overlap region
		region_it = linked_list_iterator_new(chrom_list);
		assert(region_it);
		found = 0;
		list_region = (alig_region_t *)linked_list_iterator_last(region_it);

		//Iterate regions for this chrom
		while(list_region)
		{
			//This region overlaps the new region?
			if(list_region->chrom == region->chrom)	//Unneccesary?
			{
				//This region begins inside other region
				if(region->start_pos > list_region->start_pos && region->start_pos < list_region->end_pos)
				{
					//Merge
					list_region->end_pos = region->end_pos;
					found = 1;
					break;
				}
				else if(region->end_pos > list_region->start_pos && region->end_pos < list_region->end_pos)
				{
					//Merge
					list_region->start_pos = region->start_pos;
					found = 1;
					break;
				}
				else if(region->start_pos > list_region->start_pos && region->end_pos < list_region->end_pos)
				{
					//Noth-1ing to do
					found = 1;
					break;
				}
				else if(region->start_pos < list_region->start_pos && region->end_pos > list_region->end_pos)
				{
					//Merge
					list_region->start_pos = region->start_pos;
					list_region->end_pos = region->end_pos;
					found = 1;
					break;
				}
			}

			//Get next region
			list_region = (alig_region_t *)linked_list_iterator_prev(region_it);
		}

		//Free iterator
		linked_list_iterator_free(region_it);

		//Insert new region in list
		if(!found)
			linked_list_insert(region, chrom_list);
		//printf("Table insertion: %d:%u-%u\n", region->chrom, region->start_pos, region->start_pos + region->length);
	}
}

/**
 * REGION DISCOVER
 */

int
region_get_from_cigar(char *cigar, size_t cigar_len, size_t pos, size_t *r_pos, size_t *r_end_pos)
{
	//CIGAR
	uint32_t *cigar32;
	uint32_t cigar_l;

	//CHECK ARGUMENTS
	{
		assert(cigar);
		assert(cigar_len > 0);
		assert(pos >= 0);
		assert(r_pos);
		assert(r_end_pos);
	}

	//Get CIGAR 32 bits
	convert_to_cigar_uint32_t((uint8_t *)cigar32, cigar, cigar_len);

	//Call region get
	region_get(cigar32, cigar_l, pos, r_pos, r_end_pos);

	free(cigar32);
}

int
region_get_from_bam1(const bam1_t *alig, size_t *r_pos, size_t *r_end_pos)
{
	//CIGAR
	uint32_t *cigar;
	uint32_t cigar_l;

	//Read pos
	size_t read_pos;

	//CHECK ARGUMENTS
	{
		assert(alig);
		assert(r_pos);
		assert(r_end_pos);
	}

	//Get CIGAR
	cigar = bam1_cigar(alig);
	cigar_l = alig->core.n_cigar;

	//Get read position
	read_pos = alig->core.pos + 1;

	//Call region get
	return region_get(cigar, cigar_l, read_pos, r_pos, r_end_pos);
}

int
region_get(uint32_t *cigar, uint32_t cigar_l, size_t pos, size_t *r_pos, size_t *r_end_pos)
{
	//Region
	size_t reg_pos = SIZE_MAX;
	size_t reg_end_pos;

	int i;

	//CHECK ARGUMENTS
	{
		//Check nulls
		assert(cigar);
		//assert(cigar_l > 0);
		assert(pos >= 0);
		assert(r_pos);
		assert(r_end_pos);
	}

	//If cigar length == 0 then return error
	if(cigar_l == 0)
		return -1;

	//Search for indels
	uint32_t elem, type;
	int disp = 0;
	for(i = 0; i < cigar_l; i++)
	{
		elem = cigar[i] >> BAM_CIGAR_SHIFT;	//Get number of bases from cigar
		type = cigar[i] & BAM_CIGAR_MASK;	//Get type from cigar

		switch(type)
		{
		case BAM_CINS:	//Insertion
			//printf("I");
			//Is first CIGAR elem?
			if(i == 0)
				break;

			//Set initial position
			if(reg_pos == SIZE_MAX)
			{
				reg_pos = pos + disp - 1;
			}

			//Set final position
			reg_end_pos = pos + disp;

			break;
		case BAM_CDEL:	//Deletion
			//printf("D");
			//Set initial position
			if(reg_pos == SIZE_MAX)
			{
				reg_pos = pos + disp;
			}

			//Set final position
			reg_end_pos = pos + disp + elem - 1;

			disp += elem;
			break;

		case BAM_CMATCH:
		case BAM_CDIFF:
		case BAM_CEQUAL:
			//Increment displacement
			disp += elem;
			break;

		case BAM_CHARD_CLIP:
		case BAM_CSOFT_CLIP:
		if(i == 0)
		{
			break;
		}

		default:
			//printf("X");

			break;
		}
	}
	//printf("\n");

	//Set output
	if(reg_pos != SIZE_MAX)
	{
		*r_pos = reg_pos;
		*r_end_pos = reg_end_pos;
		return 0;
	}
	else
	{
		*r_pos = SIZE_MAX;
		return 0;
	}

}

int
region_get_from_batch(const bam_batch_t* batch, alig_region_table_t *region_table)
{
	int i;
	bam1_t *alig;
	size_t r_pos;
	size_t r_end_pos;

	assert(region_table);

	//Region
	alig_region_t *region;

	for(i = 0; i < batch->num_alignments; i++)
	{
		alig = batch->alignments_p[i];

		//FILTERS: MAP QUALITY = 0, NOT PRIMARY ALIGNMENT, DIFFERENT MATE CHROM
		if(alig->core.qual != 0 && !(alig->core.flag & BAM_FSECONDARY) && alig->core.mtid == alig->core.tid)
		{
			region_get_from_bam1(alig, &r_pos, &r_end_pos);
			if(r_pos != -1)
			{
				//Region found
				/*if(r_pos != r_end_pos)
					printf("Region = %d : %d-%d\n", alig->core.tid+1, r_pos, r_end_pos);
				else
					printf("Region = %d : %d\n", alig->core.tid+1, r_pos);*/

				//Create region
				region = (alig_region_t *)malloc(sizeof(alig_region_t));
				region_create(region);

				//Fill region
				region->chrom = alig->core.tid;
				region->start_pos = r_pos;
				region->end_pos = r_end_pos;

				//Insert in hash table
				region_table_insert(region_table, region);
			}
		}
	}
}

int
region_get_from_file(const char *bam_path)
{
	bam_file_t *bam_f = NULL;
	bam_batch_t *batch = NULL;
	linked_list_t *chrom_list = NULL;
	linked_list_iterator_t *region_it = NULL;
	khiter_t iter;
	alig_region_table_t region_table;
	alig_region_t *region;
	int i;

	//Open bam
	printf("Opening BAM from \"%s\" ...\n", bam_path);
	bam_f = bam_fopen(bam_path);
	printf("BAM opened!...\n");

	//Read batch
	batch = bam_batch_new(10000000, SINGLE_CHROM_BATCH);
	bam_fread_max_size(batch, 10000000, 1, bam_f);

	//Create region table
	region_table_create(&region_table);

	//Process batch
	region_get_from_batch(batch, &region_table);

	//Print regions to screen
	for (i = kh_begin(region_table.chroms); i != kh_end(region_table.chroms); ++i)
	{
		//Get list
		iter = kh_get(32, region_table.chroms, i);
		if (iter != kh_end(region_table.chroms) && kh_exist(region_table.chroms, iter))
		{
			chrom_list = kh_value(region_table.chroms, iter);
			assert(chrom_list);

			//Get iterator
			region_it = linked_list_iterator_new(chrom_list);
			assert(region_it);

			region = (alig_region_t *)linked_list_iterator_last(region_it);
			while(region)
			{
				//Print region
				if(region->start_pos != region->end_pos)
					printf("%d:%d-%d\n", region->chrom + 1,  region->start_pos, region->end_pos);
				else
					printf("%d:%d\n", region->chrom + 1, region->start_pos);

				//Get next region
				region = (alig_region_t *)linked_list_iterator_prev(region_it);
			}
			//Free iterator
			linked_list_iterator_free(region_it);
		}
	}

	//Free region
	region_table_destroy(&region_table);

	//Free batch
	bam_batch_free(batch, 1);

	//Memory free
	printf("\nClosing BAM file...\n");
	bam_fclose(bam_f);
	printf("BAM closed.\n");
}

int
region_bam_overlap(bam1_t *read, alig_region_t *region)
{
	size_t read_pos;
	size_t read_end;

	assert(read);
	assert(region);

	//Same chrom?
	if(read->core.tid != region->chrom)
	{
		//Dont overlap
		return 0;
	}

	//Get read position
	read_pos = read->core.pos;
	read_end = read_pos + read->core.l_qseq;

	//Is in region?
	if(read_end < region->start_pos || read_pos > region->end_pos)
	{
		//Dont overlap
		return 0;
	}

	//Overlap
	return 1;
}

