#include "bioformats/fastq/fastq_batch_reader.h"
#include "buffers.h"

/*
#include "dna/dna_aligner.h"
#include "rna/rna_aligner.h"
#include "build-index/index_builder.h"

#include "sa/sa_builder.h"
#include "sa/sa_index2.h"

#include "aligners/bwt/bwt.h"
*/

//--------------------------------------------------------------------
// structs
//--------------------------------------------------------------------

typedef struct wf_batch {
  int num_threads;
  void *index;
  void *optarg;
  batch_writer_input_t *writer_input;
  mapping_batch_t *mapping_batch;  
} wf_batch_t;

wf_batch_t *wf_batch_new(void *index, void *optarg,
			 batch_writer_input_t *writer_input,
			 mapping_batch_t *mapping_batch);

void wf_batch_free(wf_batch_t *p);

//--------------------------------------------------------------------

typedef struct wf_in {
  int num_threads;
  fastq_batch_reader_input_t *fq_reader_input;
  wf_batch_t *wf_batch;
} wf_in_t;

wf_in_t *wf_in_new(fastq_batch_reader_input_t *fq_reader_input,
		   wf_batch_t *wf_batch);

void wf_in_free(wf_in_t *p);

//--------------------------------------------------------------------
// global variables
//--------------------------------------------------------------------

extern size_t total_num_mappings;

//--------------------------------------------------------------------
// function prototypes
//--------------------------------------------------------------------

void *fq_reader1(void *data);
int bam_writer1(void *data);
int bam_writer2(void *data);

int sa_mapper_stage(void *data);
int sa_binary_mapper_stage(void *data);

int bwt_mapper_stage(void *data);

int binary_search_stage(void *data);
int pre_binary_search_stage(void *data);

//--------------------------------------------------------------------

void read_fastq(char *fastq_filename, array_list_t *reads);
int cmp_read(void const *a, void const *b);

//--------------------------------------------------------------------
//--------------------------------------------------------------------
