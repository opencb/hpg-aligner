/*
 * edit_fastq.c
 *
 *  Created on: May 22, 2013
 *      Author: jtarraga
 */

#include "edit_fastq.h"

//====================================================================
// W O R K F L O W     F O R      S T A T I S T I C S
//====================================================================

#define CONSUMER_STAGE   -1

//--------------------------------------------------------------------
// workflow input
//--------------------------------------------------------------------

typedef struct fastq_edit_wf_input {
  size_t num_edited;
  size_t num_passed;
  size_t num_failed;
  edit_options_t *options;
  fastq_file_t *in_file;
  fastq_file_t *edit_file;
  fastq_file_t *failed_file;
} fastq_edit_wf_input_t;

fastq_edit_wf_input_t *fastq_edit_wf_input_new(edit_options_t *opts,
					       fastq_file_t *in_file,
					       fastq_file_t *edit_file,
					       fastq_file_t *failed_file) {
  
  fastq_edit_wf_input_t *b = (fastq_edit_wf_input_t *) calloc(1, sizeof(fastq_edit_wf_input_t));

  b->num_edited = 0;
  b->num_passed = 0;
  b->num_failed = 0;
  b->options = opts;
  b->in_file = in_file;
  b->edit_file = edit_file;
  b->failed_file = failed_file;

  return b;
}

void fastq_edit_wf_input_free(fastq_edit_wf_input_t *b) {
  if (b) {
    free(b);
  }
}

//--------------------------------------------------------------------
// batch structure between the different workflow stages
//--------------------------------------------------------------------

typedef struct fastq_edit_wf_batch {
  size_t num_tmp_edited;
  size_t *num_edited;
  size_t *num_passed;
  size_t *num_failed;
  edit_options_t *options;
  array_list_t *fq_reads;
  array_list_t *passed_reads;
  array_list_t *failed_reads;
  fastq_file_t *edit_file;
  fastq_file_t *failed_file;
} fastq_edit_wf_batch_t;

fastq_edit_wf_batch_t *fastq_edit_wf_batch_new(size_t *num_edited,
					       size_t *num_passed,
					       size_t *num_failed,
					       edit_options_t *opts,
					       array_list_t *fq_reads,
					       fastq_file_t *edit_file,
					       fastq_file_t *failed_file) {
  
  fastq_edit_wf_batch_t *b = (fastq_edit_wf_batch_t *) calloc(1, 
							      sizeof(fastq_edit_wf_batch_t));
  
  size_t num_reads = array_list_size(fq_reads) / 2;
  if (num_reads <= 0) num_reads = 100;

  b->num_tmp_edited = 0;
  b->num_edited = num_edited;
  b->num_passed = num_passed;
  b->num_failed = num_failed;
  b->options = opts;
  b->fq_reads = fq_reads;
  b->passed_reads = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);  
  b->failed_reads = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);  
  b->edit_file = edit_file;
  b->failed_file = failed_file;
  
  return b;
}

void fastq_edit_wf_batch_free(fastq_edit_wf_batch_t *b) {
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

void *fastq_edit_producer(void *input) {
  
  fastq_edit_wf_input_t *wf_input = (fastq_edit_wf_input_t *) input;
  fastq_edit_wf_batch_t *wf_batch = NULL;
  size_t max_num_reads = wf_input->options->batch_size;

  array_list_t *fq_reads = array_list_new(max_num_reads,
					  1.25f, COLLECTION_MODE_ASYNCHRONIZED);  
  size_t num_reads = fastq_fread_se(fq_reads, max_num_reads, wf_input->in_file);

  if (num_reads) {
    wf_batch = fastq_edit_wf_batch_new(&wf_input->num_edited,
				       &wf_input->num_passed,
				       &wf_input->num_failed,
				       wf_input->options,
				       fq_reads,
				       wf_input->edit_file,
				       wf_input->failed_file);				      
  } else {
    array_list_free(fq_reads, NULL);
  }
    
  return wf_batch;
}

//--------------------------------------------------------------------
// workflow worker : edit reads
//--------------------------------------------------------------------

int fastq_edit_worker(void *data) {
  fastq_edit_wf_batch_t *batch = (fastq_edit_wf_batch_t *) data;

  edit_options_t *o = batch->options;

  fastq_edit_options_t *opts;
  opts = fastq_edit_options_new(o->left_length, o->min_left_quality, 
				o->max_left_quality, o->right_length, 
				o->min_right_quality, o->max_right_quality, 
				0, 0, 0);


  batch->num_tmp_edited = fastq_edit(batch->fq_reads, opts);

  if (o->filter_on) {
    fastq_filter_options_t *filter_opts;

    filter_opts = fastq_filter_options_new(o->min_read_length, o->max_read_length,
					   o->min_read_quality, o->max_read_quality,
					   o->max_out_of_quality, 
					   MIN_VALUE, MIN_VALUE, MAX_VALUE,
					   MIN_VALUE, MIN_VALUE, MAX_VALUE,
					   o->max_N);
    
    fastq_filter(batch->fq_reads, batch->passed_reads, 
		 batch->failed_reads, filter_opts);
    
    // free memory
    fastq_filter_options_free(filter_opts);
  }

  // free memory
  fastq_edit_options_free(opts);


  return CONSUMER_STAGE;
}

//--------------------------------------------------------------------
// workflow consumer : write reads into separated files
//--------------------------------------------------------------------

int fastq_edit_consumer(void *data) {
  fastq_edit_wf_batch_t *batch = (fastq_edit_wf_batch_t *) data;

  edit_options_t *o = batch->options;
  if (o->filter_on) {
    // write passed reads
    fastq_fwrite(batch->passed_reads, batch->edit_file);
    (*batch->num_passed) += array_list_size(batch->passed_reads);
    
    // write failed reads
    fastq_fwrite(batch->failed_reads, batch->failed_file);
    (*batch->num_failed) += array_list_size(batch->failed_reads);
  } else {
    // write passed reads
    fastq_fwrite(batch->fq_reads, batch->edit_file);
  }
  (*batch->num_edited) += batch->num_tmp_edited;

  /*
  */
  // free memory
  fastq_edit_wf_batch_free(batch);
}

//--------------------------------------------------------------------
// workflow description
//--------------------------------------------------------------------

void edit_fastq(edit_options_t *opts) {

  int name_length = strlen(opts->out_dirname) + 100;

  char edit_filename[name_length];
  sprintf(edit_filename, "%s/edit.fq", opts->out_dirname);
  fastq_file_t *edit_file = fastq_fopen_mode(edit_filename, "w");
  
  char failed_filename[name_length];
  fastq_file_t *failed_file = NULL;
  if (opts->filter_on) {
    sprintf(failed_filename, "%s/failed.fq", opts->out_dirname);
    failed_file = fastq_fopen_mode(failed_filename, "w");
  }

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
  fastq_edit_wf_input_t *wf_input = fastq_edit_wf_input_new(opts,
							    in_file,
							    edit_file,
							    failed_file);
  // create and initialize workflow
  workflow_t *wf = workflow_new();
  
  workflow_stage_function_t stage_functions[] = {fastq_edit_worker};
  char *stage_labels[] = {"FASTQ edit worker"};
  workflow_set_stages(1, &stage_functions, stage_labels, wf);
  
  // optional producer and consumer functions
  workflow_set_producer(fastq_edit_producer, "FASTQ edit producer", wf);
  workflow_set_consumer(fastq_edit_consumer, "FASTQ edit consumer", wf);
  
  workflow_run_with(opts->num_threads, wf_input, wf);

  printf("\n\n");
  printf("RESULTS\n");
  printf("=================================================\n");
  printf("Num. edited reads : %lu\n", wf_input->num_edited);
  printf("Output file       : %s\n", edit_filename);
  if (opts->filter_on) {
    printf("\nFiltering : Enabled\n");
    printf("\tNum. passed reads : %lu (%s)\n", wf_input->num_passed, edit_filename);
    printf("\tNum. failed reads : %lu (%s)\n", wf_input->num_failed, failed_filename);
  }
  printf("=================================================\n");

  // free memory
  workflow_free(wf);
  fastq_edit_wf_input_free(wf_input);
  //
  // end of workflow management
  //------------------------------------------------------------------

  // free memory and close files
  fastq_fclose(in_file);
  fastq_fclose(edit_file);
  if (failed_file) fastq_fclose(failed_file);
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
