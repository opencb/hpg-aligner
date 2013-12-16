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

#include <bioformats/bam/samtools/bam.h>
#include <bioformats/bam/bam_file.h>
#include <aligners/bwt/genome.h>
#include "aux/aux_library.h"

/**
 * BAM REALIGN
 */

int alig_bam_file(char *bam_path, char *ref_name, char *ref_path);
int alig_bam_batch(bam_batch_t* batch, genome_t* ref);


#endif /* ALIG_H_ */
