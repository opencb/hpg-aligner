#ifndef _SA_RNA_MAPPER_H
#define _SA_RNA_MAPPER_H

#include <stdio.h>
#include <stdlib.h>

#include "batch_writer.h"
#include "sa/sa_search.h"
#include "dna/sa_dna_commons.h"
#include "dna/doscadfun.h"
#include "dna/sa_mapper_stage.h"
#include "rna/rna_server.h"
#include "breakpoint.h"

//#include "aligners/bwt/bwt.h"

//--------------------------------------------------------------------
//Fastq Reader & Writer
//--------------------------------------------------------------------

void *sa_alignments_reader_rna(void *input);
void *sa_fq_reader_rna(void *input);
int sa_sam_writer_rna(void *data);

//--------------------------------------------------------------------
// sa_batch_t struct
//--------------------------------------------------------------------

typedef struct sa_rna_input {
  genome_t *genome;
  avls_list_t *avls_list;
  metaexons_t *metaexons;
  sw_optarg_t *sw_optarg;
  FILE *file1;
  FILE *file2;
} sa_rna_input_t;

typedef struct sa_batch {
  size_t num_reads;
  array_list_t *fq_reads;
  array_list_t **mapping_lists;
} sa_batch_t;


inline sa_batch_t *sa_batch_new(array_list_t *fq_reads) {

  fastq_read_t *read;
  size_t num_reads = array_list_size(fq_reads);

  sa_batch_t *p = (sa_batch_t *) malloc(sizeof(sa_batch_t));
  p->num_reads = num_reads;
  p->fq_reads = fq_reads;
  p->mapping_lists = (array_list_t **) malloc(num_reads * sizeof(array_list_t *));

  for (size_t i = 0; i < num_reads; i++) {
    p->mapping_lists[i] = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  }

  return p;

}

inline sa_batch_t *sa_batch_simple_new(array_list_t *fq_reads) {

  fastq_read_t *read;
  size_t num_reads = array_list_size(fq_reads);

  sa_batch_t *p = (sa_batch_t *) malloc(sizeof(sa_batch_t));
  p->num_reads = num_reads;
  p->fq_reads = fq_reads;

  return p;

}


inline void sa_batch_free(sa_batch_t *p) {
  if (p) {
    if (p->fq_reads) { 
      array_list_free(p->fq_reads, (void *) fastq_read_free);
    }
    if (p->mapping_lists) { 
      free(p->mapping_lists); 
    }
    free(p);
  }
}


//--------------------------------------------------------------------
// sa mapper
//--------------------------------------------------------------------

int sa_rna_mapper(void *data);
int sa_rna_mapper_last(void *data);

//--------------------------------------------------------------------
//--------------------------------------------------------------------

#endif
