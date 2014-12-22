#ifndef COMMONS_FASTQ_H
#define COMMONS_FASTQ_H

/*
 * fastq_commons.h
 *
 *  Created on: May 22, 2013
 *      Author: jtarraga
 *              
 */

#include <stdlib.h>
#include <string.h>

#include "argtable/argtable2.h"

#include "bioformats/fastq/fastq_filter.h"

//------------------------------------------------------------------------

#define NO_VALUE       -1
#define MIN_VALUE       0
#define MAX_VALUE  100000

//------------------------------------------------------------------------

int parse_range(int *min, int *max, char *range, char *msg);

void free_argtable(int num_filter_options, void **argtable);
void usage_argtable(char *exec_name, char *command_name, void **argtable);

#endif	/*  COMMONS_FASTQ_H  */

//------------------------------------------------------------------------
//------------------------------------------------------------------------
