#ifndef PREPROCESS_RNA_H
#define PREPROCESS_RNA_H

#include "aligners/bwt/bwt.h"
#include "aligners/bwt/genome.h"
#include "buffers.h"
#include "cal_seeker.h"

struct preprocess_rna_input {
  size_t max_intron_length;
  size_t flank_length;
  size_t seeds_max_distance;
  size_t seed_size;
  genome_t *genome;
};

typedef struct preprocess_data {
  size_t data_number;
  int *num_cal_targets;
  unsigned char *negative_strand;
  unsigned char **associate_cals;
  int **cal_targets;
} preprocess_data_t;

preprocess_data_t *preprocess_data_new(size_t data_number, preprocess_data_t *preprocess_data);
void preprocess_data_free(preprocess_data_t *preprocess_data);

void preprocess_rna_input_init(size_t max_intron_length, size_t flank_length, 
			       size_t seeds_max_distance, size_t seed_size, 
			       genome_t *genome, preprocess_rna_input_t *preprocess_rna);


int apply_preprocess_rna(preprocess_rna_input_t* input, batch_t *batch);


#endif
