#ifndef INDEX_BUILDER_H
#define INDEX_BUILDER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aligners/bwt/bwt.h"
#include "aligners/bwt/genome.h"
#include "commons/file_utils.h"
#include "argtable2.h"
#include "sa/sa_index3.h"

#define SA_INDEX  0
#define BWT_INDEX 1

#define NUM_INDEX_OPTIONS 3
#define NUM_INDEX_BWT_OPTIONS 2

typedef struct index_options {
  int index_ratio;
  int bs_index;
  int help;
  char *ref_genome;
  char *index_filename;  
} index_options_t;

index_options_t *index_options_new();
void index_options_free(index_options_t *options);

void run_index_builder_bwt(char *genome_filename, char *bwt_dirname, 
			   int bwt_ratio, bool duplicate_strand, 
			   char *nucleotides);

void run_index_builder_sa(char *genome_filename, uint k_value, 
			  char *bwt_dirname);

void help_index_builder();

//void run_index_builder_bs(char *genome_filename, char *bwt_dirname, int bwt_ratio, char *bases);


#endif // INDEX_BUILDER_H
