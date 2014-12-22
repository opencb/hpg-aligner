/*
 * stats_fastq.c
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

#include "stats_fastq.h"

//--------------------------------------------------------------------
// stats counters
//--------------------------------------------------------------------

stats_counters_t *stats_counters_new() {
  stats_counters_t *sc = calloc(1, sizeof(stats_counters_t));

  sc->num_reads = 0;
  sc->phred = QUALITY_PHRED33_VALUE;

  sc->filter_on = 0;
  sc->num_passed = 0;
  sc->num_failed = 0;

  sc->min_length = 100000;
  sc->max_length = 0;
  sc->acc_length = 0;
  sc->mean_length = 0.0f;

  sc->acc_quality = 0;
  sc->mean_quality = 0.0f;

  sc->kmers_on = 0;
  for(int i = 0; i < NUM_KMERS; i++) {
    sc->kmers[i].counter = 0;
    sc->kmers[i].counter_by_pos_size = 0;
  }
  //  memset(sc->kmers, 0, sizeof(NUM_KMERS * sizeof(kmer_t)));

  sc->kh_length_histogram = kh_init(32);
  sc->kh_quality_histogram = kh_init(32);
  sc->kh_gc_histogram = kh_init(32);

  sc->kh_count_quality_per_nt = kh_init(32);
  sc->kh_acc_quality_per_nt = kh_init(32);
  
  sc->kh_num_As_per_nt = kh_init(32);
  sc->kh_num_Cs_per_nt = kh_init(32);
  sc->kh_num_Ts_per_nt = kh_init(32);
  sc->kh_num_Gs_per_nt = kh_init(32);
  sc->kh_num_Ns_per_nt = kh_init(32);

  return sc;
}

void stats_counters_free(stats_counters_t *sc) {
  if (sc) {


    if (sc->kmers_on) {
      kmer_t *kmer;
      for (int i = 0; i < NUM_KMERS; i++) {
	kmer = &sc->kmers[i];
	if (kmer->counter_by_pos_size > 0) {
	  free(kmer->counter_by_pos);
	}
      }
    }

    kh_destroy(32, sc->kh_length_histogram);
    kh_destroy(32, sc->kh_quality_histogram);
    kh_destroy(32, sc->kh_gc_histogram);

    kh_destroy(32, sc->kh_count_quality_per_nt);
    kh_destroy(32, sc->kh_acc_quality_per_nt);
    
    kh_destroy(32, sc->kh_num_As_per_nt);
    kh_destroy(32, sc->kh_num_Cs_per_nt);
    kh_destroy(32, sc->kh_num_Ts_per_nt);
    kh_destroy(32, sc->kh_num_Gs_per_nt);
    kh_destroy(32, sc->kh_num_Ns_per_nt);
    
    free(sc);
  }
}

//====================================================================
// W O R K F L O W     F O R      S T A T I S T I C S
//====================================================================

#define CONSUMER_STAGE   -1

//--------------------------------------------------------------------
// workflow input
//--------------------------------------------------------------------

typedef struct fastq_stats_wf_input {
  stats_options_t *options;
  fastq_file_t *in_file;
  fastq_read_stats_options_t *read_stats_options;
  stats_counters_t *stats_counters;
} fastq_stats_wf_input_t;

fastq_stats_wf_input_t *fastq_stats_wf_input_new(stats_options_t *opts,
						 fastq_file_t *in_file,
						 stats_counters_t *stats_counters) {
  
  fastq_stats_wf_input_t *p = (fastq_stats_wf_input_t *) calloc(1, sizeof(fastq_stats_wf_input_t));

  p->in_file = in_file;
  p->options = opts;
  p->stats_counters = stats_counters;
  p->read_stats_options = fastq_read_stats_options_new(opts->kmers_on, 1);

  return p;
}

void fastq_stats_wf_input_free(fastq_stats_wf_input_t *p) {
  if (p) {
    if (p->read_stats_options) fastq_read_stats_options_free(p->read_stats_options);

    free(p);
  }
}

//--------------------------------------------------------------------
// structure between the different workflow stages
//--------------------------------------------------------------------

typedef struct fastq_stats_wf_batch {
  stats_options_t *options;
  fastq_read_stats_options_t *read_stats_options;
  array_list_t *fq_reads;
  array_list_t *fq_stats;
  stats_counters_t *stats_counters;
  array_list_t *passed_reads;
  array_list_t *failed_reads;
} fastq_stats_wf_batch_t;

fastq_stats_wf_batch_t *fastq_stats_wf_batch_new(stats_options_t *opts,
						 fastq_read_stats_options_t *read_stats_options,
						 array_list_t *fq_reads,
						 array_list_t *fq_stats,
						 stats_counters_t *stats_counters) {
  
  fastq_stats_wf_batch_t *p = (fastq_stats_wf_batch_t *) calloc(1, 
								sizeof(fastq_stats_wf_batch_t));
  
  p->options = opts;
  p->read_stats_options = read_stats_options;
  p->fq_reads = fq_reads;
  p->fq_stats = fq_stats;
  p->stats_counters = stats_counters;
  p->passed_reads = NULL;
  p->failed_reads = NULL;
  
  return p;
}

void fastq_stats_wf_batch_free(fastq_stats_wf_batch_t *p) {
  if (p) {
    if (p->fq_reads) array_list_free(p->fq_reads, (void *) fastq_read_free);
    if (p->fq_stats) fastq_reads_stats_free(p->fq_stats);
    if (p->passed_reads) array_list_free(p->passed_reads, NULL);
    if (p->failed_reads) array_list_free(p->failed_reads, NULL);    

    free(p);
  }
}

//--------------------------------------------------------------------
// workflow producer : read FastQ file
//--------------------------------------------------------------------

void *fastq_stats_producer(void *input) {
  
  fastq_stats_wf_input_t *wf_input = (fastq_stats_wf_input_t *) input;
  fastq_stats_wf_batch_t *wf_batch = NULL;
  size_t max_num_reads = wf_input->options->batch_size;

  
  array_list_t *fq_reads = array_list_new(max_num_reads,
					  1.25f, COLLECTION_MODE_ASYNCHRONIZED);  
  size_t num_reads = fastq_fread_se(fq_reads, max_num_reads, wf_input->in_file);

  if (num_reads) {
    wf_batch = fastq_stats_wf_batch_new(wf_input->options,
					wf_input->read_stats_options,
					fq_reads,
					NULL, //fq_stats,
					wf_input->stats_counters);				      
  } else {
    array_list_free(fq_reads, NULL);
  }
    
  return wf_batch;
}

//--------------------------------------------------------------------
// workflow worker : compute statistics per read
//--------------------------------------------------------------------

int fastq_stats_worker(void *data) {
  fastq_stats_wf_batch_t *batch = (fastq_stats_wf_batch_t *) data;

  if (batch->stats_counters->filter_on) {
    // prepare filter options
    stats_options_t *o = batch->options;
    fastq_filter_options_t *opts;
    opts = fastq_filter_options_new(o->min_read_length, o->max_read_length,
				    o->min_read_quality, o->max_read_quality,
				    o->max_out_of_quality, o->left_length,
				    o->min_left_quality, o-> max_left_quality,
				    o->right_length, o->min_right_quality,
				    o->max_right_quality, o->max_N);
    
    // prepare filter lists
    size_t num_reads = array_list_size(batch->fq_reads) / 2;
    if (num_reads <= 0) num_reads = 100;

    batch->passed_reads = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);  
    batch->failed_reads = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);  

    // and finally, filter
    fastq_filter(batch->fq_reads, batch->passed_reads, batch->failed_reads, opts);

    if ((num_reads = array_list_size(batch->passed_reads)) > 0) {
      batch->fq_stats = array_list_new(num_reads,
				       1.25f, COLLECTION_MODE_ASYNCHRONIZED);  
      // compute statistics for those reads that passed the filters
      fastq_reads_stats(batch->passed_reads,
			batch->read_stats_options,
			batch->fq_stats);
    }

    // free memory
    fastq_filter_options_free(opts);

  } else {
    // no filter
    size_t num_reads = array_list_size(batch->fq_reads);
    batch->fq_stats = array_list_new(num_reads,
				     1.25f, COLLECTION_MODE_ASYNCHRONIZED);  

    fastq_reads_stats(batch->fq_reads,
		      batch->read_stats_options,
		      batch->fq_stats);

  }

  return CONSUMER_STAGE;
}

//--------------------------------------------------------------------
// workflow consumer : merge statistics
//--------------------------------------------------------------------

int fastq_stats_consumer(void *data) {
  fastq_stats_wf_batch_t *batch = (fastq_stats_wf_batch_t *) data;

  khiter_t k;
  khash_t(32) *h;

  fastq_read_t *fq_read;
  array_list_t *fq_reads = batch->fq_reads;
  fastq_read_stats_t *read_stats;
  stats_counters_t *counters = batch->stats_counters;

  if (counters->filter_on) {
    counters->num_passed += array_list_size(batch->passed_reads);
    counters->num_failed += array_list_size(batch->failed_reads);
    fq_reads = batch->passed_reads;
  }

  if (array_list_size(batch->fq_stats) > 0) {
    int ret;
    size_t value, read_length;
    size_t num_reads = array_list_size(fq_reads);
    
    kmer_t *kmer_in, *kmer_out;
    int len, kmers_on = batch->read_stats_options->kmers_on;
    
    for (size_t i = 0 ; i < num_reads; i++) {
      fq_read = array_list_get(i, fq_reads);
      read_stats = array_list_get(i, batch->fq_stats);
      if (read_stats) {
	counters->num_reads++;
	
	// length, min, max, acc
	read_length = read_stats->length;
	counters->acc_length += read_length;
	if (counters->min_length > read_length) 
	  counters->min_length = read_length;
	if (counters->max_length < read_length) 
	  counters->max_length = read_length;
	
	// quality
	counters->acc_quality += read_stats->quality_average;
	
	// A, C, T, G, N
	counters->num_As += read_stats->num_A;
	counters->num_Cs += read_stats->num_C;
	counters->num_Ts += read_stats->num_T;
	counters->num_Gs += read_stats->num_G;
	counters->num_Ns += read_stats->num_N;
	
	// length histogram
	value = read_length;
	h = counters->kh_length_histogram;
	k = kh_put(32, h, value, &ret);
	if (ret) {
	  kh_value(h, k) = 1;
	} else {
	  kh_value(h, k) =  kh_value(h, k) + 1;
	}
	
	// quality histogram
	value = round(read_stats->quality_average);
	h = counters->kh_quality_histogram;
	k = kh_put(32, h, value, &ret);
	if (ret) {
	  kh_value(h, k) = 1;
	} else {
	  kh_value(h, k) =  kh_value(h, k) + 1;
	}
	
	// GC histogram
	value = 100 * (read_stats->num_G + read_stats->num_C) / read_length;
	h = counters->kh_gc_histogram;
	k = kh_put(32, h, value, &ret);
	if (ret) {
	  kh_value(h, k) = 1;
	} else {
	  kh_value(h, k) =  kh_value(h, k) + 1;
	}
	
	
	// stats per nt position
	for (size_t j = 0; j < read_length; j++) {
	  
	  // quality counter
	  h = counters->kh_count_quality_per_nt;
	  k = kh_put(32, h, j, &ret);
	  if (ret) {
	    kh_value(h, k) = 1;
	  } else {
	    kh_value(h, k) =  kh_value(h, k) + 1;
	  }
	  
	  // quality accumulator 
	  h = counters->kh_acc_quality_per_nt;
	  k = kh_put(32, h, j, &ret);
	  if (ret) {
	    kh_value(h, k) = fq_read->quality[j];
	  } else {
	    kh_value(h, k) =  kh_value(h, k) + fq_read->quality[j];
	  }
	  
	  // A, T, C, G, N
	  h = NULL;
	  switch (fq_read->sequence[j]) {
	  case 'A': 	h = counters->kh_num_As_per_nt;
	    break;
	  case 'T':	h = counters->kh_num_Ts_per_nt;
	    break;
	  case 'C':	h = counters->kh_num_Cs_per_nt;
	    break;
	  case 'G':	h = counters->kh_num_Gs_per_nt;
	    break;
	  case 'N':	h = counters->kh_num_Ns_per_nt;
	    break;
	  default:	break;
	  }
	  
	  if (h) {
	    k = kh_put(32, h, j, &ret);
	    if (ret) {
	      kh_value(h, k) = 1;
	    } else {
	      kh_value(h, k) =  kh_value(h, k) + 1;
	    }
	  }
	}
	
	// merge kmers results
	if (kmers_on) {
	  for (int j = 0; j < NUM_KMERS; j++) {
	    kmer_in = &read_stats->kmers[j];
	    //	  printf("kmer_in: %i %s, counter = %lu, len = %i\n", kmer_in->id, kmer_in->string, kmer_in->counter, kmer_in->counter_by_pos_size);
	    if (kmer_in->counter > 0) {
	      len = kmer_in->counter_by_pos_size;
	      kmer_out = &counters->kmers[j];
	      kmer_out->counter += kmer_in->counter;
	      //	    printf("kmer %i, len = %i, (out len = %i)\n", j, len, kmer_out->counter_by_pos_size);
	      if (kmer_out->counter_by_pos_size == 0) {
		//	      printf("\tcalloc...\n");
		kmer_out->counter_by_pos_size = len;
		kmer_out->counter_by_pos = (size_t *) calloc(len, sizeof(size_t));
		//	      printf("\t\tcalloc...done\n");
	      } else if (kmer_out->counter_by_pos_size < len) {
		//	      printf("\trealloc...\n");
		kmer_out->counter_by_pos_size = len;
		kmer_out->counter_by_pos = realloc(kmer_out->counter_by_pos, len);
		//	      printf("\t\trealloc...done\n");
	      }
	      for (int k = 0; k < len; k++) {
		kmer_out->counter_by_pos[k] += kmer_in->counter_by_pos[k];
	      }
	    }
	  }
	}
      }
    }
  }
  
  // free memory
  fastq_stats_wf_batch_free(batch);
}

//--------------------------------------------------------------------
// workflow description
//--------------------------------------------------------------------

void stats_fastq(stats_options_t *opts) {

  fastq_file_t *fq_file = fastq_fopen(opts->in_filename);
  stats_counters_t *counters = stats_counters_new();
  counters->phred = opts->quality_encoding_value;
  counters->kmers_on = opts->kmers_on;
  
  // update opts
  if (opts->filter_on) {
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
  }
  counters->filter_on = opts->filter_on;
  
  //------------------------------------------------------------------
  // workflow management
  //
  fastq_stats_wf_input_t *wf_input = fastq_stats_wf_input_new(opts,
							      fq_file,
							      counters);
  
  // create and initialize workflow
  workflow_t *wf = workflow_new();
  
  workflow_stage_function_t stage_functions[] = {fastq_stats_worker};
  char *stage_labels[] = {"FASTQ stats worker"};
  workflow_set_stages(1, &stage_functions, stage_labels, wf);
  
  // optional producer and consumer functions
  workflow_set_producer(fastq_stats_producer, "FASTQ stats producer", wf);
  workflow_set_consumer(fastq_stats_consumer, "FASTQ stats consumer", wf);
  
  workflow_run_with(opts->num_threads, wf_input, wf);
  
  // free memory
  workflow_free(wf);
  fastq_stats_wf_input_free(wf_input);
  //
  // end of workflow management
  //------------------------------------------------------------------

  // update and sort kmers according to number of counts
  kmer_t *kmer;
  for (int i = 0; i < NUM_KMERS; i++) {
    kmer = &counters->kmers[i];
    kmer->id = i;
    kmers_string(i, kmer->string);
  }
  qsort(counters->kmers, NUM_KMERS, sizeof(kmer_t), kmers_sort);
    
  // and report statistics
  stats_report(counters, opts);

  printf("\n\n");
  printf("RESULTS\n");
  printf("=================================================\n");
  printf("Report files and images were stored in '%s' directory\n", opts->out_dirname);
  if (counters->filter_on) {
    printf("\nFiltering: enabled\n");
    printf("\tSo, statistics were computed for %lu of %lu reads.\n",
	   counters->num_passed, counters->num_passed + counters->num_failed);
  } else {
    printf("\nFiltering: disabled\n");
    printf("\tSo, statistics were computed for the whole input file.\n");
  }
  printf("=================================================\n");


  // free memory
  stats_counters_free(counters);
  fastq_fclose(fq_file);
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
