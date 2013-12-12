/*
 * alig.h
 *
 *  Created on: Dec 12, 2013
 *      Author: rmoreno
 */

#ifndef ALIG_H_
#define ALIG_H_

#include <assert.h>

#include <stdlib.h>
#include <stdint.h>

#include "bioformats/bam/samtools/bam.h"


/**
 * CIGAR
 */

int alig_cigar_leftmost(char *ref, char *read, size_t length, uint32_t *cigar, size_t cigar_l, uint32_t *new_cigar, size_t *new_cigar_l);
int alig_cigar_unclip(uint32_t *cigar, size_t cigar_l, uint32_t *new_cigar, size_t *new_cigar_l);
int alig_cigar_count_m_blocks(uint32_t *cigar, size_t cigar_l, size_t *blocks);


#endif /* ALIG_H_ */
