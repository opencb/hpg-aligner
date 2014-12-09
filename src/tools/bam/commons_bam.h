#ifndef COMMONS_BAM_H
#define COMMONS_BAM_H

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

#include "bioformats/features/region/region_table.h"
#include "bioformats/features/region/region_table_utils.h"
#include "bioformats/bam/bam_file.h"

//------------------------------------------------------------------------

#define NO_VALUE       -1
#define MIN_VALUE       0
#define MAX_VALUE  100000

//------------------------------------------------------------------------

extern int read_progress;

//------------------------------------------------------------------------

void free_argtable(int num_filter_options, void **argtable);
void usage_argtable(char *exec_name, char *command_name, void **argtable);

int parse_range(int *min, int *max, char *range, char *msg);

region_table_t *build_region_table(char *bam_filename, char *by_string, 
				   char *by_gff_file);

#endif	/*  COMMONS_BAM_H  */

//------------------------------------------------------------------------
//------------------------------------------------------------------------
