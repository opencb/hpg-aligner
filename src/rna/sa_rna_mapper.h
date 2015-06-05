#ifndef _SA_RNA_MAPPER_H
#define _SA_RNA_MAPPER_H

#include <stdio.h>
#include <stdlib.h>

#include "workflow_functions.h"
#include "batch_writer.h"
#include "sa/sa_search.h"
#include "dna/sa_dna_commons.h"
#include "dna/doscadfun.h"
#include "dna/sa/sa_mapper_stage.h"
#include "rna/rna_server.h"
#include "breakpoint.h"
#include "pair_server.h"

//#include "aligners/bwt/bwt.h"

//--------------------------------------------------------------------
//Fastq Reader & Writer
//--------------------------------------------------------------------

int write_to_file(void *data);
void *sa_alignments_reader_rna(void *input);
void *sa_fq_reader_rna(void *input);
int sa_sam_writer_rna(void *data);
int sa_bam_writer_rna(void *data);

//--------------------------------------------------------------------
// sa_batch_t struct
//--------------------------------------------------------------------

typedef struct sa_rna_input {
  cal_optarg_t *cal_optarg;
  genome_t *genome;
  avls_list_t *avls_list;
  metaexons_t *metaexons;
  sw_optarg_t *sw_optarg;
  FILE *file1;
  FILE *file2;
  pair_server_input_t *pair_input;
  int min_score;
  int max_alig;
} sa_rna_input_t;

typedef struct sa_batch {
  size_t num_reads;
  array_list_t *fq_reads;
  array_list_t **mapping_lists;
} sa_batch_t;


sa_batch_t *sa_batch_new(array_list_t *fq_reads);
sa_batch_t *sa_batch_simple_new(array_list_t *fq_reads);
void sa_batch_free(sa_batch_t *p);


//--------------------------------------------------------------------
// sa mapper
//--------------------------------------------------------------------

int sa_rna_mapper(void *data);
int sa_rna_mapper_last(void *data);

//--------------------------------------------------------------------
//--------------------------------------------------------------------

#endif
