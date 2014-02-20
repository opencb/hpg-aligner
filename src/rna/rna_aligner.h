#ifndef RNA_ALIGNER_H
#define RNA_ALIGNER_H

#include <sys/syscall.h>
#include <sched.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include <pthread.h>
#include <stdlib.h>
#include <getopt.h>
#include <unistd.h>

#include "containers/cprops/trie.h"

#include "commons/log.h"
#include "commons/file_utils.h"

#include "commons/workflow_scheduler.h"

#include "bioformats/fastq/fastq_batch_reader.h"

#include "batch_writer.h"
#include "error.h"
#include "timing.h"
#include "buffers.h"
#include "bwt_server.h"
#include "region_seeker.h"
#include "cal_seeker.h"
#include "sw_server.h"
#include "pair_server.h"
#include "preprocess_rna.h"
#include "options.h"
#include "statistics.h"
#include "workflow_functions.h"
#include "rna_server.h"

//For SA Mapping
#include "sa/sa_index3.h"
#include "sa/sa_search.h"
#include "rna/sa_rna_mapper.h"

#include "dna/sa_dna_commons.h"
#include "dna/doscadfun.h"
#include "dna/sa_io_stages.h"
#include "dna/sa_mapper_stage.h"

typedef struct buffer_reader_input {
  fastq_batch_reader_input_t *reader_input;
  linked_list_t *buffer;
} buffer_reader_input_t;

void buffer_reader_input_init(fastq_batch_reader_input_t *reader_input, 
			      linked_list_t *buffer, 
			      buffer_reader_input_t *buffer_reader_input);

typedef struct extra_stage {
  workflow_t *workflow;
  linked_list_t *align_list;
  pair_mng_t *pair_mng;
  char *intron_filename;
} extra_stage_t;
/*
typedef struct buffer_item {
  fastq_read_t *read;
  array_list_t *alignments_list;
}buffer_item_t;
*/
typedef struct buffer_pair_item {
  fastq_read_t *read_1;
  array_list_t *alignments_list_1;
  fastq_read_t *read_2;
  array_list_t *alignments_list_2;
}buffer_pair_item_t;

/*
void run_rna_aligner(genome_t *genome, bwt_index_t *bwt_index, 
		     pair_mng_t *pair_mng,
		     bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
		     report_optarg_t *report_optarg, metaexons_t *metaexons, 
		     options_t *options);
*/
void rna_aligner(options_t *options);

typedef struct exon {
  int chr;
  int strand;
  int start;
  int end;
  int exon_number;
  char *chr_name;
  char *gene_id;
  char *transcript_id;
  char *exon_id;
} exon_t;

static inline exon_t *exon_new(int chr, int strand, int start, int end, int exon_number,
			     char *chr_name,  char *gene_id, char *transcript_id,
			     char *exon_id) {
  exon_t *p = (exon_t *) malloc(sizeof(exon_t));

  p->chr = chr;
  p->strand = strand;
  p->start = start;
  p->end = end;
  p->exon_number = exon_number;
  p->chr_name = chr_name;
  p->gene_id = gene_id;
  p->transcript_id = transcript_id;
  p->exon_id = exon_id;

  return p;
}

static inline void exon_free(exon_t *p) {
  if (p) {
    if (p->chr_name) free(p->chr_name);
    if (p->gene_id) free(p->gene_id);
    if (p->transcript_id) free(p->transcript_id);
    if (p->exon_id) free(p->exon_id);

    free(p);
  }
}

static inline void exon_display(exon_t *p) {
  if (p) {
    printf("chr %s (%i): strand %c: %i - %i, gene_id = %s, transcript_id = %s, exon_number = %i, exon_id = %s\n",
	   p->chr_name, p->chr, (p->strand ? '-' : '+'), p->start, p->end, p->gene_id,
	   p->transcript_id, p->exon_number, p->exon_id);
  }
}

void load_transcriptome(char *filename, genome_t *genome, 
			avls_list_t *avls_list, metaexons_t *metaexons);

#endif
