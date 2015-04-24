/*
 * filter_bam.c
 *
 *  Created on: May 22, 2013
 *      Author: jtarraga
 */

#include "filter_bam.h"

//====================================================================
// W O R K F L O W     F O R      S T A T I S T I C S
//====================================================================

#define CONSUMER_STAGE   -1

//--------------------------------------------------------------------
// workflow input
//--------------------------------------------------------------------

typedef struct bam_filter_wf_input {
  size_t num_passed;
  size_t num_failed;
  filter_options_t *options;
  bam_file_t *in_file;
  bam_file_t *passed_file;
  bam_file_t *failed_file;
} bam_filter_wf_input_t;

bam_filter_wf_input_t *bam_filter_wf_input_new(filter_options_t *opts,
					       bam_file_t *in_file,
					       bam_file_t *passed_file,
					       bam_file_t *failed_file) {
  
  bam_filter_wf_input_t *p = (bam_filter_wf_input_t *) calloc(1, sizeof(bam_filter_wf_input_t));

  p->num_passed = 0;
  p->num_failed = 0;
  p->options = opts;
  p->in_file = in_file;
  p->passed_file = passed_file;
  p->failed_file = failed_file;

  return p;
}

void bam_filter_wf_input_free(bam_filter_wf_input_t *p) {
  if (p) {
    free(p);
  }
}

//--------------------------------------------------------------------
// batch structure between the different workflow stages
//--------------------------------------------------------------------

typedef struct bam_filter_wf_batch {
  size_t *num_passed;
  size_t *num_failed;
  filter_options_t *options;
  array_list_t *bam1s;
  array_list_t *passed_bam1s;
  array_list_t *failed_bam1s;
  bam_file_t *passed_file;
  bam_file_t *failed_file;
} bam_filter_wf_batch_t;

bam_filter_wf_batch_t *bam_filter_wf_batch_new(size_t *num_passed,
					       size_t *num_failed,
					       filter_options_t *opts,
					       array_list_t *bam1s,
					       bam_file_t *passed_file,
					       bam_file_t *failed_file) {
  
  bam_filter_wf_batch_t *p = (bam_filter_wf_batch_t *) calloc(1, 
							      sizeof(bam_filter_wf_batch_t));
  
  size_t num_items = array_list_size(bam1s) / 2;
  if (num_items <= 0) num_items = 100;

  p->num_passed = num_passed;
  p->num_failed = num_failed;
  p->options = opts;
  p->bam1s = bam1s;
  p->passed_bam1s = array_list_new(num_items, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);  
  p->failed_bam1s = array_list_new(num_items, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);  
  p->passed_file = passed_file;
  p->failed_file = failed_file;
  
  return p;
}

void bam_filter_wf_batch_free(bam_filter_wf_batch_t *p) {
  if (p) {
    if (p->bam1s) {
      bam1_t *bam1;
      size_t num_items = array_list_size(p->bam1s);
      for (size_t i = 0; i < num_items; i++) {
	bam1 = array_list_get(i, p->bam1s);
	bam_destroy1(bam1);
      }
      array_list_free(p->bam1s, NULL);
    }
    if (p->passed_bam1s) array_list_free(p->passed_bam1s, NULL);
    if (p->failed_bam1s) array_list_free(p->failed_bam1s, NULL);
    
    free(p);
  }
}

//--------------------------------------------------------------------
// workflow producer : read BAM file
//--------------------------------------------------------------------

void *bam_filter_producer(void *input) {
  
  bam_filter_wf_input_t *wf_input = (bam_filter_wf_input_t *) input;
  bam_filter_wf_batch_t *new_batch = NULL;
  int max_num_bam1s = wf_input->options->batch_size;

  bamFile bam_file = wf_input->in_file->bam_fd;
  bam1_t *bam1;

  array_list_t *bam1_list = array_list_new(max_num_bam1s, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);  
  for (int i = 0; i < max_num_bam1s; i++) {
    bam1 = bam_init1();    
    if (bam_read1(bam_file, bam1) > 0) {

      if (++read_progress % 500000 == 0) {
	LOG_INFO_F("%i reads extracting from disk...\n", read_progress);
      }

      array_list_insert(bam1, bam1_list);
    } else {
      bam_destroy1(bam1);
      break;
    }
  }

  size_t num_items = array_list_size(bam1_list);
  
  if (num_items == 0) {
    array_list_free(bam1_list, NULL);
  } else {
    new_batch = bam_filter_wf_batch_new(&wf_input->num_passed,
					&wf_input->num_failed,
					wf_input->options,
					bam1_list,
					wf_input->passed_file,
					wf_input->failed_file);				      
  }
    
  return new_batch;
}

//--------------------------------------------------------------------
// workflow worker : filter reads
//--------------------------------------------------------------------

int bam_filter_worker(void *data) {
  bam_filter_wf_batch_t *batch = (bam_filter_wf_batch_t *) data;

  filter_options_t *o = batch->options;

  bam_filter_options_t *opts;
  opts = bam_filter_options_new(o->unique, o->proper_pairs,
				o->min_length, o->max_length,
				o->min_quality, o->max_quality,
				o->min_num_errors, o->max_num_errors,
				o->region_table);
  
  bam_filter(batch->bam1s, batch->passed_bam1s, 
	     batch->failed_bam1s, opts);

  // free memory
  bam_filter_options_free(opts);

  return CONSUMER_STAGE;
}

//--------------------------------------------------------------------
// workflow consumer : write reads into separated files
//--------------------------------------------------------------------

int bam_filter_consumer(void *data) {
  bam_filter_wf_batch_t *batch = (bam_filter_wf_batch_t *) data;

  // write passed reads
  bam_fwrite_bam1s(batch->passed_bam1s, batch->passed_file);
  (*batch->num_passed) += array_list_size(batch->passed_bam1s);

  // write failed reads
  bam_fwrite_bam1s(batch->failed_bam1s, batch->failed_file);
  (*batch->num_failed) += array_list_size(batch->failed_bam1s);

  // free memory
  bam_filter_wf_batch_free(batch);

  return 0;

}

//--------------------------------------------------------------------
// workflow description
//--------------------------------------------------------------------

void filter_bam(filter_options_t *opts) {

  int name_length = strlen(opts->out_dirname) + 100;

  bam_file_t *in_file = bam_fopen(opts->in_filename);

  char passed_filename[name_length];
  sprintf(passed_filename, "%s/passed.bam", opts->out_dirname);
  bam_file_t *passed_file = bam_fopen_mode(passed_filename, in_file->bam_header_p, "w");
  bam_fwrite_header(in_file->bam_header_p, passed_file);

  char failed_filename[name_length];
  sprintf(failed_filename, "%s/failed.bam", opts->out_dirname);
  bam_file_t *failed_file = bam_fopen_mode(failed_filename, in_file->bam_header_p, "w");
  bam_fwrite_header(in_file->bam_header_p, failed_file);

  // update opts
  if (opts->min_length == NO_VALUE) opts->min_length = MIN_VALUE;
  if (opts->max_length == NO_VALUE) opts->max_length = MAX_VALUE;
  if (opts->min_quality == NO_VALUE) opts->min_quality = MIN_VALUE;
  if (opts->max_quality == NO_VALUE) opts->max_quality = MAX_VALUE;
  if (opts->min_num_errors == NO_VALUE) opts->min_num_errors = MIN_VALUE;
  if (opts->max_num_errors == NO_VALUE) opts->max_num_errors = MAX_VALUE;

  //------------------------------------------------------------------
  // workflow management
  //
  bam_filter_wf_input_t *wf_input = bam_filter_wf_input_new(opts,
							    in_file,
							    passed_file,
							    failed_file);
  
  // create and initialize workflow
  workflow_t *wf = workflow_new();
  
  workflow_stage_function_t stage_functions[] = {bam_filter_worker};
  char *stage_labels[] = {"BAM filter worker"};
  workflow_set_stages(1, stage_functions, stage_labels, wf);
  
  // optional producer and consumer functions
  workflow_set_producer((workflow_producer_function_t *)bam_filter_producer, "BAM filter producer", wf);
  workflow_set_consumer((workflow_consumer_function_t *)bam_filter_consumer, "BAM filter consumer", wf);
  
  workflow_run_with(opts->num_threads, wf_input, wf);

  printf("\n\n");
  printf("RESULTS\n");
  printf("======================================================\n");
  printf("Num. passed alignments: %lu (%s)\n", wf_input->num_passed, passed_filename);
  printf("Num. failed alignments: %lu (%s)\n", wf_input->num_failed, failed_filename);
  printf("======================================================\n");

  // free memory
  workflow_free(wf);
  bam_filter_wf_input_free(wf_input);
  //
  // end of workflow management
  //------------------------------------------------------------------

  // free memory and close files
  in_file->bam_header_p = NULL;
  bam_fclose(in_file);     
  passed_file->bam_header_p = NULL;
  bam_fclose(passed_file); 
  bam_fclose(failed_file);
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
