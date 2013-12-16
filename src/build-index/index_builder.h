#ifndef INDEX_BUILDER_H
#define INDEX_BUILDER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aligners/bwt/bwt.h"
#include "aligners/bwt/genome.h"
#include "commons/file_utils.h"

void run_index_builder(char *genome_filename, char *bwt_dirname, 
		       int bwt_ratio, bool duplicate_strand, char *nucleotides);

void help_index_builder();

//void run_index_builder_bs(char *genome_filename, char *bwt_dirname, int bwt_ratio, char *bases);


#endif // INDEX_BUILDER_H
