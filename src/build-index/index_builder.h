#ifndef INDEX_BUILDER_H
#define INDEX_BUILDER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "aligners/bwt/bwt.h"
#include "aligners/bwt/genome.h"
#include "commons/file_utils.h"
#include "options.h"
#include "argtable2.h"
#include "sa/sa_index3.h"

#define SA_INDEX  0
#define BWT_INDEX 1

#define NUM_INDEX_OPTIONS     5
#define NUM_INDEX_BWT_OPTIONS 0

#define BWT_RATIO_DEFAULT  8

typedef struct index_options {
  int mode;
  int version;
  int index_ratio;
  int help;
  char *decoy_genome;
  char *ref_genome;
  char *index_filename;  
  char *cmdline;
} index_options_t;

index_options_t *index_options_new();
void index_options_free(index_options_t *options);

void run_index_builder_bwt(char *genome_filename, char *bwt_dirname, 
			   int bwt_ratio, bool duplicate_strand, 
			   char *nucleotides);

void help_index_builder();
void index_options_display(index_options_t *options);

void run_index_builder(int argc, char **argv, char *mode_str);

#endif // INDEX_BUILDER_H
