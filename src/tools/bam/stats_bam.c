/*
 * stats_bam.c
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

#include "stats_bam.h"

//====================================================================
// W O R K F L O W     F O R      S T A T I S T I C S
//====================================================================

#define CONSUMER_STAGE   -1

workflow_t *workflow;

//--------------------------------------------------------------------
// workflow input
//--------------------------------------------------------------------

typedef struct bam_stats_wf_input {
  stats_options_t *options;
  bam_file_t *in_file;
  bam_stats_options_t *bam_stats_options;
  stats_counters_t *counters;
} bam_stats_wf_input_t;

//--------------------------------------------------------------------

bam_stats_wf_input_t *bam_stats_wf_input_new(stats_options_t *options,
					     bam_file_t *in_file,
					     bam_stats_options_t *bam_stats_options,
					     stats_counters_t *counters) {
  
  bam_stats_wf_input_t *p = (bam_stats_wf_input_t *) calloc(1, sizeof(bam_stats_wf_input_t));

  p->in_file = in_file;
  p->options = options;
  p->bam_stats_options = bam_stats_options;
  p->counters = counters;

  return p;
}

//--------------------------------------------------------------------

void bam_stats_wf_input_free(bam_stats_wf_input_t *p) {
  if (p) {
    if (p->bam_stats_options) bam_stats_options_free(p->bam_stats_options);

    free(p);
  }
}

//--------------------------------------------------------------------
// structure between the different workflow stages
//--------------------------------------------------------------------

typedef struct bam_stats_wf_batch {
  stats_options_t *options;
  bam_stats_options_t *bam_stats_options;
  array_list_t *bam1s;
  array_list_t *bam_stats;
  stats_counters_t *counters;
  array_list_t *passed_bam1s;
  array_list_t *failed_bam1s;
} bam_stats_wf_batch_t;

//--------------------------------------------------------------------

bam_stats_wf_batch_t *bam_stats_wf_batch_new(stats_options_t *options,
					     bam_stats_options_t *bam_stats_options,
					     array_list_t *bam1s,
					     array_list_t *bam_stats,
					     stats_counters_t *counters) {
  bam_stats_wf_batch_t *p = (bam_stats_wf_batch_t *) calloc(1, sizeof(bam_stats_wf_batch_t));
  
  p->options = options;
  p->bam_stats_options = bam_stats_options;
  p->bam1s = bam1s;
  p->bam_stats = bam_stats;
  p->counters = counters;
  p->passed_bam1s = NULL;
  p->failed_bam1s = NULL;
  
  return p;
}

//--------------------------------------------------------------------

void bam_stats_wf_batch_free(bam_stats_wf_batch_t *p) {
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

    if (p->bam_stats) array_list_free(p->bam_stats, (void *) bam_stats_free);
    if (p->passed_bam1s) array_list_free(p->passed_bam1s, NULL);
    if (p->failed_bam1s) array_list_free(p->failed_bam1s, NULL);    

    free(p);
  }
}

//--------------------------------------------------------------------
// workflow producer
//--------------------------------------------------------------------

void *bam_stats_producer(void *input) {
  bam_stats_wf_input_t *wf_input = (bam_stats_wf_input_t *) input;
  bam_stats_wf_batch_t *new_batch = NULL;
  int max_num_bam1s = wf_input->options->batch_size;

  bamFile bam_file = wf_input->in_file->bam_fd;
  bam1_t *bam1;

  //  printf("producer: active items = %i of %i\n", workflow_get_num_items(workflow), workflow->max_num_work_items);

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
    new_batch = bam_stats_wf_batch_new(wf_input->options,
				       wf_input->bam_stats_options,
				       bam1_list,  
				       NULL,
				       wf_input->counters);				      
  }

  return new_batch;
}

//--------------------------------------------------------------------
// workflow worker
//--------------------------------------------------------------------

int bam_stats_worker(void *data) {
  bam_stats_wf_batch_t *batch = (bam_stats_wf_batch_t *) data;

  //  printf("worker: active items = %i of %i\n", workflow_get_num_items(workflow), workflow->max_num_work_items);

  // filter ?
  if (batch->options->region_table) {
    // prepare filter options
    bam_filter_options_t *opts;
    opts = bam_filter_options_new(0, 0,
				  MIN_VALUE, MAX_VALUE,
				  MIN_VALUE, MAX_VALUE,
				  MIN_VALUE, MAX_VALUE,
				  batch->options->region_table);
    
    // prepare filter lists
    size_t num_items = array_list_size(batch->bam1s) / 2;
    if (num_items <= 0) num_items = 100;

    batch->passed_bam1s = array_list_new(num_items, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);  
    batch->failed_bam1s = array_list_new(num_items, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);  

    // and finally, filter
    bam_filter(batch->bam1s, batch->passed_bam1s, 
	       batch->failed_bam1s, opts);

    if ((num_items = array_list_size(batch->passed_bam1s)) > 0) {
      // compute statistics for those reads that passed the filters
      batch->bam_stats = array_list_new(num_items,
					1.25f, COLLECTION_MODE_ASYNCHRONIZED);  

      bam1s_stats(batch->passed_bam1s,
		  batch->bam_stats_options,
		  batch->bam_stats);
    } else {
      batch->bam_stats = array_list_new(1, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);  
    }

    // free memory
    bam_filter_options_free(opts);

  } else {
    size_t num_items = array_list_size(batch->bam1s);
    batch->bam_stats = array_list_new(num_items,
				      1.25f, COLLECTION_MODE_ASYNCHRONIZED);  
    
    bam1s_stats(batch->bam1s,
		batch->bam_stats_options,
		batch->bam_stats);
    
  }

  return CONSUMER_STAGE;
}

//--------------------------------------------------------------------
// workflow consumer
//--------------------------------------------------------------------
int bam_progress = 0;

int bam_stats_consumer(void *data) {

  //  printf("consumer: active items = %i of %i\n", workflow_get_num_items(workflow), workflow->max_num_work_items);

  bam1_t *bam1;
  bam_stats_wf_batch_t *batch = (bam_stats_wf_batch_t *) data;
  
  bam_stats_t *stats;
  stats_counters_t *counters = batch->counters;

  int strand, seq_id, start, end, len, quality, gc;
  size_t num_items = array_list_size(batch->bam_stats);

  array_list_t *bam1s = batch->bam1s;

  if (batch->options->region_table) {
    counters->num_passed += array_list_size(batch->passed_bam1s);
    counters->num_failed += array_list_size(batch->failed_bam1s);
    bam1s = batch->passed_bam1s;
  }

  // variables for storint stats in db
  int db_on = batch->options->db_on;
  sqlite3 *db;
  sqlite3_stmt* stmt;
  bam_query_fields_t *fields;
  char* errorMessage;
  khash_t(stats_chunks) *hash;
  size_t *sequence_lengths;
  char **sequence_labels;

  if (db_on) {
    db = batch->options->db;
    hash = batch->options->hash;

    sequence_labels = counters->sequence_labels;
    sequence_lengths = counters->sequence_lengths;

    prepare_statement_bam_query_fields(db, &stmt);

    sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage);
  }

  for (int i = 0; i < num_items; i++) {
    stats = array_list_get(i, batch->bam_stats);
    if (!stats) continue;
  
    bam_progress++;

    strand = stats->strand;

    // merge tmp stats into "final" stats
    counters->num_reads++;

    if (!stats->mapped) continue;

    counters->num_unique_alignments += stats->unique_alignment;
    counters->num_mapped_reads += stats->mapped;
    counters->num_mapped_reads_1 += stats->mapped_pair_1;
    counters->num_unmapped_reads_1 += stats->unmapped_pair_1;
    counters->num_mapped_reads_2 += stats->mapped_pair_2;
    counters->num_unmapped_reads_2 += stats->unmapped_pair_2;

    counters->num_As += stats->num_As;
    counters->num_Cs += stats->num_Cs;
    counters->num_Gs += stats->num_Gs;
    counters->num_Ts += stats->num_Ts;
    counters->num_Ns += stats->num_Ns;

    gc = round(100.0f * stats->num_GCs / stats->seq_length);
    counters->GC_content[gc]++;

    counters->single_end = stats->single_end;

    counters->num_unique_alignments_strand[strand] += stats->unique_alignment;
    counters->num_mapped_reads_strand[strand] += stats->mapped;
    
    counters->num_indels += stats->num_indels;
    counters->indels_acc += stats->indels_length;

    if (stats->num_errors < NUM_ERRORS_STATS) {
      counters->num_errors[stats->num_errors]++;
    } else {
      counters->num_errors[NUM_ERRORS_STATS]++;
    }

    len = stats->seq_length;
    if (len > counters->max_alignment_length) counters->max_alignment_length = len;
    if (len < counters->min_alignment_length) counters->min_alignment_length = len;
    
    len = stats->isize;
    if (len > counters->max_insert_size) counters->max_insert_size = len;
    if (len < counters->min_insert_size) counters->min_insert_size = len;
    counters->insert_size_acc += len;
    
    quality = stats->quality;
    if (quality > counters->max_quality) counters->max_quality = quality;
    if (quality < counters->min_quality) counters->min_quality = quality;
    counters->quality_acc += quality;
    counters->quality[quality]++;
    
    // coverage addition
    bam1 = array_list_get(i, bam1s);

    seq_id = bam1->core.tid;
    start = bam1->core.pos;
    end = start + bam1->core.l_qseq;

    for (int p = start; p < end; p++) {
      counters->sequence_depths_per_nt[seq_id][p]++;
    }

    // save into db: bam query fields
    if (db_on) {
      fields = bam_query_fields_new(bam1_qname(bam1), sequence_labels[seq_id], 
				    sequence_lengths[seq_id], strand, start + 1, end + 1, 
				    (uint32_t) bam1->core.flag, stats->quality, 
				    stats->num_errors, stats->num_indels, stats->indels_length,
				    stats->isize);
      insert_statement_bam_query_fields(fields, stmt, db);

      update_chunks_hash(fields->chr, fields->chr_length, BAM_CHUNKSIZE,
			 fields->start, fields->end, hash);
      bam_query_fields_free(fields);
    } 
  } // for each item

  //  if (bam_progress % 500000 == 0) {
  //    LOG_INFO_F("%i reads processed !\n", bam_progress);
  //  }

  if (db_on) {
    sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errorMessage);

    sqlite3_finalize(stmt);
  }

  // free memory
  bam_stats_wf_batch_free(batch);
}

//--------------------------------------------------------------------
// workflow description
//--------------------------------------------------------------------

void stats_bam(stats_options_t *opts) {

  bam_file_t *bam_file = bam_fopen(opts->in_filename);
  
  size_t ref_length = 0;
  int num_targets = bam_file->bam_header_p->n_targets;

  stats_counters_t *counters = stats_counters_new();

  counters->sequence_labels = (char **) calloc(num_targets, sizeof(char *));
  counters->sequence_depths_per_nt = (uint16_t **) calloc(num_targets, sizeof(uint16_t *));
  //  counters->sequence_depths_per_nt = (int **) calloc(num_targets, sizeof(int *));
  counters->sequence_lengths = (size_t *) calloc(num_targets, sizeof(size_t));
  counters->depth_per_sequence = (double *) calloc(num_targets, sizeof(size_t));

  size_t size = 0;
  for (int i = 0; i < num_targets; i++) {
    ref_length += bam_file->bam_header_p->target_len[i];

    counters->sequence_labels[i] = strdup(bam_file->bam_header_p->target_name[i]);
    counters->sequence_depths_per_nt[i] = (uint16_t *) calloc(bam_file->bam_header_p->target_len[i], sizeof(uint16_t));
    counters->sequence_lengths[i] = bam_file->bam_header_p->target_len[i];
  }

  counters->ref_length = ref_length;
  counters->num_sequences = num_targets;

  if (opts->db_on) {
    // create db table and hash
    opts->hash = kh_init(stats_chunks);
    create_stats_db(opts->out_dbname, BAM_CHUNKSIZE, create_bam_query_fields, &opts->db);                
  }

  // bam_stats_options
  bam_stats_options_t *bam_stats_options = bam_stats_options_new(NULL, NULL);

  //------------------------------------------------------------------
  // workflow management
  //
  bam_stats_wf_input_t *wf_input = bam_stats_wf_input_new(opts,
							  bam_file,
							  bam_stats_options,
							  counters);
  
  // create and initialize workflow
  workflow_t *wf = workflow_new();

  workflow = wf;

  workflow_stage_function_t stage_functions[] = {bam_stats_worker};
  char *stage_labels[] = {"BAM stats worker"};
  workflow_set_stages(1, stage_functions, stage_labels, wf);
  
  // optional producer and consumer functions
  workflow_set_producer((workflow_producer_function_t *)bam_stats_producer, "BAM stats producer", wf);
  workflow_set_consumer((workflow_consumer_function_t *)bam_stats_consumer, "BAM stats consumer", wf);
  
  workflow_run_with(opts->num_threads, wf_input, wf);
  
  // free memory
  workflow_free(wf);
  bam_stats_wf_input_free(wf_input);
  //
  // end of workflow management
  //------------------------------------------------------------------

  // free bam stats 
  //  bam_stats_options_free(bam_stats_options);

  // compute coverage  
  size_t nt_depth, unmapped_nts = 0;
  size_t seq_len, acc_per_sequence = 0, acc = 0;
  for (int i = 0; i < num_targets; i++) {
    seq_len = counters->sequence_lengths[i];
    acc_per_sequence = 0;
    for (int j = 0; j < seq_len; j++) {
      nt_depth = ((int) counters->sequence_depths_per_nt[i][j]);
      if (nt_depth) {
	acc_per_sequence += nt_depth;
      } else {
	unmapped_nts++;
      }
    }
    acc += acc_per_sequence;
    counters->depth_per_sequence[i] = 1.0f * acc_per_sequence / counters->sequence_lengths[i];
  }
  
  counters->depth = 1.0f * acc / ref_length;
  counters->unmapped_nts = unmapped_nts;

  counters->num_nucleotides = counters->num_As + counters->num_Cs + 
    counters->num_Gs + counters->num_Ts + counters->num_Ns;

  if (opts->db_on) {
    // after all records insertions, insert touched chunks from hash,
    // and after all inserts then create index
    insert_chunk_hash(BAM_CHUNKSIZE, opts->hash, opts->db);
    create_stats_index(create_bam_index, opts->db);
  }

  // report
  report_stats(opts->in_filename, opts->out_dirname, opts->db, counters);

  printf("\n\n");
  printf("RESULTS\n");
  printf("=================================================\n");
  printf("Report files and images were stored in '%s' directory\n", opts->out_dirname);
  if (opts->region_table) {
    printf("\nFiltering: enabled (regions: %s)\n", 
	   (opts->region_list ? opts->region_list : opts->gff_region_filename));
    printf("\tSo, statistics were computed for %lu of %lu alignments.\n",
	   counters->num_passed, counters->num_passed + counters->num_failed);
  } else {
    printf("\nFiltering: disabled\n");
    printf("\tSo, statistics were computed for the whole input file.\n");
  }
  printf("=================================================\n");

  // free memory and close file
  stats_counters_free(counters);
  bam_fclose(bam_file);

  if (opts->db_on) {
    // finally, close db and free hash
    close_stats_db(opts->db, opts->hash);
  }
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
