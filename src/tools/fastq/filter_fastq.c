/*
 * filter_fastq.c
 *
 *  Created on: May 22, 2013
 *      Author: jtarraga
 */

#include "filter_fastq.h"

//====================================================================
// W O R K F L O W     F O R      S T A T I S T I C S
//====================================================================

#define CONSUMER_STAGE   -1

//--------------------------------------------------------------------
// workflow input
//--------------------------------------------------------------------

typedef struct fastq_filter_wf_input {
  size_t num_passed;
  size_t num_failed;
  filter_options_t *options;
  fastq_file_t *in_file;
  fastq_file_t *passed_file;
  fastq_file_t *failed_file;
} fastq_filter_wf_input_t;

fastq_filter_wf_input_t *fastq_filter_wf_input_new(filter_options_t *opts,
						  fastq_file_t *in_file,
						  fastq_file_t *passed_file,
						  fastq_file_t *failed_file) {
  
  fastq_filter_wf_input_t *b = (fastq_filter_wf_input_t *) calloc(1, sizeof(fastq_filter_wf_input_t));

  b->num_passed = 0;
  b->num_failed = 0;
  b->options = opts;
  b->in_file = in_file;
  b->passed_file = passed_file;
  b->failed_file = failed_file;

  return b;
}

void fastq_filter_wf_input_free(fastq_filter_wf_input_t *b) {
  if (b) {
    free(b);
  }
}

//--------------------------------------------------------------------
// batch structure between the different workflow stages
//--------------------------------------------------------------------

typedef struct fastq_filter_wf_batch {
  size_t *num_passed;
  size_t *num_failed;
  filter_options_t *options;
  array_list_t *fq_reads;
  array_list_t *passed_reads;
  array_list_t *failed_reads;
  fastq_file_t *passed_file;
  fastq_file_t *failed_file;
} fastq_filter_wf_batch_t;

fastq_filter_wf_batch_t *fastq_filter_wf_batch_new(size_t *num_passed,
						   size_t *num_failed,
						   filter_options_t *opts,
						   array_list_t *fq_reads,
						   fastq_file_t *passed_file,
						   fastq_file_t *failed_file) {
  
  fastq_filter_wf_batch_t *b = (fastq_filter_wf_batch_t *) calloc(1, 
								  sizeof(fastq_filter_wf_batch_t));
  
  size_t num_reads = array_list_size(fq_reads) / 2;
  if (num_reads <= 0) num_reads = 100;

  b->num_passed = num_passed;
  b->num_failed = num_failed;
  b->options = opts;
  b->fq_reads = fq_reads;
  b->passed_reads = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);  
  b->failed_reads = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);  
  b->passed_file = passed_file;
  b->failed_file = failed_file;
  
  return b;
}

void fastq_filter_wf_batch_free(fastq_filter_wf_batch_t *b) {
  if (b) {
    array_list_free(b->fq_reads, (void *) fastq_read_free);
    array_list_free(b->passed_reads, NULL);
    array_list_free(b->failed_reads, NULL);

    free(b);
  }
}

//--------------------------------------------------------------------
// workflow producer : read FastQ file
//--------------------------------------------------------------------

void *fastq_filter_producer(void *input) {
  
  fastq_filter_wf_input_t *wf_input = (fastq_filter_wf_input_t *) input;
  fastq_filter_wf_batch_t *wf_batch = NULL;
  size_t max_num_reads = wf_input->options->batch_size;

  array_list_t *fq_reads = array_list_new(max_num_reads,
					  1.25f, COLLECTION_MODE_ASYNCHRONIZED);  
  size_t num_reads = fastq_fread_se(fq_reads, max_num_reads, wf_input->in_file);

  if (num_reads) {
    wf_batch = fastq_filter_wf_batch_new(&wf_input->num_passed,
					 &wf_input->num_failed,
					 wf_input->options,
					 fq_reads,
					 wf_input->passed_file,
					 wf_input->failed_file);				      
  } else {
    array_list_free(fq_reads, NULL);
  }
    
  return wf_batch;
}

//--------------------------------------------------------------------
// workflow worker : filter reads
//--------------------------------------------------------------------

int fastq_filter_worker(void *data) {
  fastq_filter_wf_batch_t *batch = (fastq_filter_wf_batch_t *) data;

  filter_options_t *o = batch->options;
  
  fastq_filter_options_t *opts;
  opts = fastq_filter_options_new(o->min_read_length, o->max_read_length,
				  o->min_read_quality, o->max_read_quality,
				  o->max_out_of_quality, o->left_length,
				  o->min_left_quality, o->max_left_quality,
				  o->right_length, o->min_right_quality,
				  o->max_right_quality, o->max_N);


  fastq_filter(batch->fq_reads, batch->passed_reads, 
	       batch->failed_reads, opts);

  // free memory
  fastq_filter_options_free(opts);

  return CONSUMER_STAGE;
}

//--------------------------------------------------------------------
// workflow consumer : write reads into separated files
//--------------------------------------------------------------------

int fastq_filter_consumer(void *data) {
  fastq_filter_wf_batch_t *batch = (fastq_filter_wf_batch_t *) data;

  // write passed reads
  fastq_fwrite(batch->passed_reads, batch->passed_file);
  (*batch->num_passed) += array_list_size(batch->passed_reads);

  // write failed reads
  fastq_fwrite(batch->failed_reads, batch->failed_file);
  (*batch->num_failed) += array_list_size(batch->failed_reads);

  // free memory
  fastq_filter_wf_batch_free(batch);
}

//--------------------------------------------------------------------
// workflow description
//--------------------------------------------------------------------

void filter_fastq(filter_options_t *opts) {

  int name_length = strlen(opts->out_dirname) + 100;

  char passed_filename[name_length];
  sprintf(passed_filename, "%s/passed.fq", opts->out_dirname);
  fastq_file_t *passed_file = fastq_fopen_mode(passed_filename, "w");

  char failed_filename[name_length];
  sprintf(failed_filename, "%s/failed.fq", opts->out_dirname);
  fastq_file_t *failed_file = fastq_fopen_mode(failed_filename, "w");

  fastq_file_t *in_file = fastq_fopen(opts->in_filename);

  // update opts
  if (opts->min_read_length == NO_VALUE) opts->min_read_length = MIN_VALUE;
  if (opts->max_read_length == NO_VALUE) opts->max_read_length = MAX_VALUE;
  if (opts->min_read_quality == NO_VALUE) opts->min_read_quality = MIN_VALUE;
  if (opts->max_read_quality == NO_VALUE) opts->max_read_quality = MAX_VALUE;
  if (opts->max_out_of_quality == NO_VALUE) opts->max_out_of_quality = MAX_VALUE;
  if (opts->left_length == NO_VALUE) opts->left_length = MIN_VALUE;
  if (opts->min_left_quality == NO_VALUE) opts->min_left_quality = MIN_VALUE;
  if (opts->max_left_quality == NO_VALUE) opts->max_left_quality = MAX_VALUE;
  if (opts->right_length == NO_VALUE) opts->right_length = MIN_VALUE;
  if (opts->min_right_quality == NO_VALUE) opts->min_right_quality = MIN_VALUE;
  if (opts->max_right_quality == NO_VALUE) opts->max_right_quality = MAX_VALUE;
  if (opts->max_N == NO_VALUE) opts->max_N = MAX_VALUE;

  //------------------------------------------------------------------
  // workflow management
  //
  fastq_filter_wf_input_t *wf_input = fastq_filter_wf_input_new(opts,
								in_file,
								passed_file,
								failed_file);
  
  // create and initialize workflow
  workflow_t *wf = workflow_new();
  
  workflow_stage_function_t stage_functions[] = {fastq_filter_worker};
  char *stage_labels[] = {"FASTQ filter worker"};
  workflow_set_stages(1, &stage_functions, stage_labels, wf);
  
  // optional producer and consumer functions
  workflow_set_producer(fastq_filter_producer, "FASTQ filter producer", wf);
  workflow_set_consumer(fastq_filter_consumer, "FASTQ filter consumer", wf);
  
  workflow_run_with(opts->num_threads, wf_input, wf);

  printf("\n\n");
  printf("RESULTS\n");
  printf("=================================================\n");
  printf("Num. passed reads: %lu (%s)\n", wf_input->num_passed, passed_filename);
  printf("Num. failed reads: %lu (%s)\n", wf_input->num_failed, failed_filename);
  printf("=================================================\n");

  // free memory
  workflow_free(wf);
  fastq_filter_wf_input_free(wf_input);
  //
  // end of workflow management
  //------------------------------------------------------------------

  // free memory and close files
  fastq_fclose(in_file);
  fastq_fclose(passed_file);
  fastq_fclose(failed_file);
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
