#ifndef UTILS_H
#define UTILS_H

/*
 * utils.h
 *
 *  Created on: Oct 30, 2014
 *      Author: sgallego
 */

#include "bioformats/bam/bam_file.h"
#include "samtools/bam.h"

//------------------------------------------------------------------------

void revcomp_seq(char* seq);

//------------------------------------------------------------------------

int is_pair(char *bam_filename);
int max_quality(char *bam_filename);
char *bam1_get_sequence(bam1_t *bam1, char *sequence_string);
char *bam1_get_quality(bam1_t *bam1, char *quality_string);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif // end of UTILS_H
