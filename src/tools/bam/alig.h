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
#include <bioformats/bam/bam_file.h>
#include <aligners/bwt/genome.h>
#include "aux_library.h"
#include "alig_aux.h"

/**
 * BAM REALIGN
 */

int alig_bam_file(char *bam_path, char *ref_name, char *ref_path);
int alig_bam_batch(bam_batch_t* batch, genome_t* ref);


/**
 * CIGAR
 */

int alig_cigar_leftmost(char *ref, char *read, size_t read_l, uint32_t *cigar, size_t cigar_l, uint32_t *new_cigar, size_t *new_cigar_l);

int alig_cigar_unclip(uint32_t *cigar, size_t cigar_l, uint32_t *new_cigar, size_t *new_cigar_l);
int alig_cigar_count_m_blocks(uint32_t *cigar, size_t cigar_l, size_t *blocks);
int alig_cigar_count_indels(uint32_t *cigar, size_t cigar_l, size_t *indels);
int alig_cigar_count_all(uint32_t *cigar, size_t cigar_l, size_t *m_blocks, size_t *indels, size_t *first_indel_index);


#endif /* ALIG_H_ */
