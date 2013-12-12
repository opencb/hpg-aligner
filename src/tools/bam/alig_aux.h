/*
 * alig_aux.h
 *
 *  Created on: Dec 12, 2013
 *      Author: rmoreno
 */

#ifndef ALIG_AUX_H_
#define ALIG_AUX_H_

#include <stdlib.h>
#include <stdint.h>

#include "bioformats/bam/samtools/bam.h"

int alig_aux_cigar32_to_string(uint32_t *cigar, size_t cigar_l, char* str_cigar);


#endif /* ALIG_AUX_H_ */
