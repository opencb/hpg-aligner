#ifndef WORKFLOW_FUNCTIONS_H
#define WORKFLOW_FUNCTIONS_H

#include <stdio.h>

#include "commons/log.h"
#include "commons/file_utils.h"

#include "commons/workflow_scheduler.h"

#include "bioformats/fastq/fastq_batch_reader.h"
#include "bioformats/bam/bam_file.h"

#include "rna/rna_aligner.h"
#include "rna/rna_server.h"

#include "buffers.h"
#include "bwt_server.h"
#include "region_seeker.h"
#include "cal_seeker.h"
#include "pair_server.h"
#include "sw_server.h"
#include "batch_writer.h"

#define MAX_READS_RNA 200

extern int global_status;

//--------------------------------------------------------------------
// workflow producers
//--------------------------------------------------------------------

void *fastq_reader(void *input);
void *buffer_reader(void *input);
void *file_reader(void *input);
void *file_reader_2(void *input);

//--------------------------------------------------------------------
// workflow consumer
//--------------------------------------------------------------------

int bam_writer(void *data);

//--------------------------------------------------------------------
// stage functions
//--------------------------------------------------------------------

int bwt_stage(void *data);
int bwt_stage_bs(void *data);

int seeding_stage(void *data);
int seeding_stage_bs(void *data);

int cal_stage(void *data);
int cal_stage_bs(void *data);

int pre_pair_stage(void *data);
int pre_pair_stage_bs(void *data);

int rna_preprocess_stage(void *data);

int sw_stage(void *data);
int sw_stage_bs(void *data);

int post_pair_stage(void *data);
int post_pair_stage_bs(void *data);

int bs_status_stage(void *data);

int rna_last_stage(void *data);

int rna_last_hc_stage(void *data);

//---------------------------------------------------------------------
//workflow input
//--------------------------------------------------------------------- 

typedef struct wf_input {
  fastq_batch_reader_input_t *fq_reader_input;
  batch_t *batch;
} wf_input_t;

typedef struct wf_input_buffer {
  linked_list_t *buffer;
  batch_t *batch;
} wf_input_buffer_t;

typedef struct wf_input_file {
  FILE *file;
  batch_t *batch;
} wf_input_file_t;

wf_input_t *wf_input_new(fastq_batch_reader_input_t *fq_reader_input,
                         batch_t *batch);

wf_input_buffer_t *wf_input_buffer_new(linked_list_t *buffer,
				       batch_t *batch);

wf_input_file_t *wf_input_file_new(FILE *fd,
				   batch_t *batch);

void wf_input_free(wf_input_t *wfi);

void wf_input_file_free(wf_input_file_t *wfi);

//------------------------------------------------------------------------------------

extern pthread_cond_t cond_sp;
extern pthread_mutex_t mutex_sp;

#endif
