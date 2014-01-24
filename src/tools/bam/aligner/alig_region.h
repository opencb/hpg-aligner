/*
 * region.h
 *
 *  Created on: Dec 11, 2013
 *      Author: rmoreno
 */

#ifndef ALIG_REGION_H_
#define ALIG_REGION_H_

#include <stdlib.h>
#include <stdint.h>
#include <limits.h>

#include "assert.h"
#include "aux/aux_library.h"

#include "containers/khash.h"
#include "containers/linked_list.h"
#include "bioformats/bam/samtools/bam.h"

typedef struct {
	size_t start_pos;
	size_t end_pos;
	int32_t chrom;
	uint8_t valid;
} alig_region_t;

KHASH_MAP_INIT_INT(32, linked_list_t*);
typedef struct {
	//Hash to store chrom regions alig_region_t
	khash_t(32) *chroms;
	//Others?
} alig_region_table_t;

/**
 * REGION STRUCT FUNCTIONS
 */

int region_create(alig_region_t *region);
int region_destroy(alig_region_t *region);

/**
 * REGION TABLE BY CHROM
 */

int region_table_create(alig_region_table_t *table);
int region_table_destroy(alig_region_table_t *table);

int region_table_insert(alig_region_table_t *table, alig_region_t *region);

/**
 * REGION DISCOVER
 */

int region_get_from_cigar(char *cigar, size_t cigar_len, size_t pos, size_t *r_pos, size_t *r_end_pos);
int region_get_from_bam1(const bam1_t *alig, size_t *r_pos, size_t *r_end_pos);
int region_get(uint32_t *cigar, uint32_t cigar_l, size_t pos, size_t *r_pos, size_t *r_end_pos);
int region_get_from_batch(const bam_batch_t* batch, alig_region_table_t *region_table);
int region_get_from_file(const char *bam_path);

int region_bam_overlap(bam1_t *read, alig_region_t *region);

#endif /* ALIG_REGION_H_ */
