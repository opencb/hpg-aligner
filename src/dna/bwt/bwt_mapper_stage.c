#include "bwt_mapper_stage.h"

//--------------------------------------------------------------------

void cut_adapter(char *adapter, int adapter_length, fastq_read_t *read);

//--------------------------------------------------------------------
// bwt mapper
//--------------------------------------------------------------------

int bwt_single_mapper(void *data) {
/*
  #ifdef _TIMING
  struct timeval stop, start;
  #endif

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif
  
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;

  int num_seeds = wf_batch->options->num_seeds;

  
  sa_mapping_batch_t *mapping_batch = wf_batch->mapping_batch;
  mapping_batch->options = wf_batch->options;

  sa_index3_t *sa_index = (sa_index3_t *) wf_batch->sa_index;
  
  int bam_format = mapping_batch->bam_format;

  size_t num_reads = mapping_batch->num_reads;

  float max_score;

  // smith-waterman parameters
  float match_score = wf_batch->options->match;
  float mismatch_penalty = wf_batch->options->mismatch;
  float gap_open_penalty = -1.0f * wf_batch->options->gap_open;
  float gap_extend_penalty = -1.0f * wf_batch->options->gap_extend;
  
  // CAL management


  cal_mng_t *cal_mng;
  array_list_t *cal_list;

  int sw_post_read_counter = 0;
  int sw_post_read[num_reads];

  array_list_t *cal_lists[num_reads];
  array_list_t *sw_prepare_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  fastq_read_t *read;

  cal_mng = cal_mng_new(sa_index->genome);
  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_OTHER] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif

  // for each read, create cals and prepare sw
  for (int i = 0; i < num_reads; i++) {
    read = array_list_get(i, mapping_batch->fq_reads);
    fastq_read_revcomp(read);

    if (wf_batch->options->adapter) {
      cut_adapter(wf_batch->options->adapter, wf_batch->options->adapter_length, read);
    }

    // 1) extend using mini-sw from suffix
    cal_list = create_cals(num_seeds, read, mapping_batch, sa_index, cal_mng);

    if (array_list_size(cal_list) > 0) {

    	select_best_cals(read, &cal_list);

      	if (array_list_size(cal_list) <= 0) {
			suffix_mng_search_read_cals(read, num_seeds, sa_index, cal_list, cal_mng->suffix_mng);
			if (array_list_size(cal_list) > 0) {
		  		select_best_cals(read, &cal_list);
			} else {
		  		mapping_batch->status[i] = 3; // invalid
			}
      	}

      	//fill_seed_gaps(cal_list, read, sa_index);
      	clean_cals(&cal_list, read, sa_index);
      	if (array_list_size(cal_list) <= 0) {
			mapping_batch->status[i] = 4; //clean cals
      	}

      	// 2) prepare Smith-Waterman to fill in the gaps
      	if (prepare_sw(read, sw_prepare_list, mapping_batch, sa_index, cal_list)) {
			sw_post_read[sw_post_read_counter++] = i;
      	}
    } else {
      mapping_batch->status[i] = 1; // no cals
    }

    cal_lists[i] = cal_list;
  }

  // 3) run SW to fill
  if (array_list_size(sw_prepare_list) > 0) {
    execute_sw(sw_prepare_list, mapping_batch);
    post_process_sw(sw_post_read_counter, sw_post_read, cal_lists, mapping_batch);
  }
  array_list_free(sw_prepare_list, (void *) NULL);

  // 4) prepare mappings for the writer
  for (int i = 0; i < num_reads; i++) {
    cal_list = cal_lists[i];
    read = array_list_get(i, mapping_batch->fq_reads);
    
    if (array_list_size(cal_list) > 0) {
      // filter by score
      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      max_score = get_max_score(cal_list, match_score, mismatch_penalty,
				gap_open_penalty, gap_extend_penalty);
      #ifdef _VERBOSE
      printf("\n******* after SW> max. score = %0.2f (read %s)\n", max_score, read->id);
      #endif
      filter_cals_by_max_score(max_score, &cal_list);
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_FILTER_BY_NUM_MISMATCHES] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif
    }
    
    // if BAM format, create alignments structures
    if (bam_format) {
      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      create_alignments(cal_list, read, bam_format, mapping_batch->mapping_lists[i]);
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_CREATE_ALIGNMENTS] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif
      
      // free cal list and clear seed manager for next read
      array_list_free(cal_list, (void *) NULL);
    } else {
      if (mapping_batch->mapping_lists[i]) {
	array_list_free(mapping_batch->mapping_lists[i], (void *) NULL);
	mapping_batch->mapping_lists[i] = cal_list;
      }
    }
  } // end of for reads
  
  // free memory
  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif
  cal_mng_free(cal_mng);
  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_OTHER] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif
  */
  return -1;
}

//--------------------------------------------------------------------

int bwt_pair_mapper(void *data) {
/*
  #ifdef _TIMING
  struct timeval stop, start;
  #endif

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif
  
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;

  int infer_insert;
  int pair_min_distance = wf_batch->options->pair_min_distance;
  int pair_max_distance = wf_batch->options->pair_max_distance;

  int num_seeds = wf_batch->options->num_seeds;
  
  sa_mapping_batch_t *mapping_batch = wf_batch->mapping_batch;
  mapping_batch->options = wf_batch->options;
  if (pair_min_distance > 0 && pair_max_distance > 0) {
    infer_insert = 0;
    mapping_batch->pair_min_distance = pair_min_distance; 
    mapping_batch->pair_max_distance = pair_max_distance; 
  } else {
    infer_insert = 1;
  }

  sa_index3_t *sa_index = (sa_index3_t *) wf_batch->sa_index;
  
  int bam_format = mapping_batch->bam_format;

  size_t num_reads = mapping_batch->num_reads;

  // CAL management
  cal_mng_t *cal_mng;
  array_list_t *cal_list;

  int sw_post_read_counter = 0;
  int sw_post_read[num_reads];

  array_list_t *cal_lists[num_reads];
  array_list_t *sw_prepare_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  fastq_read_t *read;

  cal_mng = cal_mng_new(sa_index->genome);
  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_OTHER] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif

  // for each read, create cals and prepare sw
  for (int i = 0; i < num_reads; i++) {
  	read = array_list_get(i, mapping_batch->fq_reads);
    fastq_read_revcomp(read);

    if (wf_batch->options->adapter) {
      cut_adapter(wf_batch->options->adapter, wf_batch->options->adapter_length, read);
    }

    // 1) extend using mini-sw from suffix
    cal_list = create_cals(num_seeds, read, mapping_batch, sa_index, cal_mng);

    if (array_list_size(cal_list) > 0) {

      //printf("------------------> %s:%i: before select_best_cals (%i):\n", __FILE__, __LINE__, array_list_size(cal_list));
      //for (int i = 0; i < array_list_size(cal_list); i++) { seed_cal_print(array_list_get(i, cal_list)); }

      select_best_cals(read, &cal_list);

      //printf("-----------------> %s:%i: after select_best_cals (%i):\n", __FILE__, __LINE__, array_list_size(cal_list));
      //for (int i = 0; i < array_list_size(cal_list); i++) { seed_cal_print(array_list_get(i, cal_list)); }

      if (array_list_size(cal_list) <= 0) {
		suffix_mng_search_read_cals(read, num_seeds, sa_index, cal_list, cal_mng->suffix_mng);
		if (array_list_size(cal_list) > 0) {
	  		select_best_cals(read, &cal_list);
		} else {
	  		mapping_batch->status[i] = 3; // invalid
		}
      }

      clean_cals(&cal_list, read, sa_index);
      if (array_list_size(cal_list) <= 0) {
		if (mapping_batch->status[i] == 0) {
	  		mapping_batch->status[i] = 4; //clean cals
		}
      } else {
		if (((seed_cal_t *)array_list_get(0, cal_list))->mapq == 0) {
	  		mapping_batch->status[i] = 2; // mapq = 0
		}
      }
    } else {
      mapping_batch->status[i] = 1; // no cals
    }
    cal_lists[i] = cal_list;
  }


  // 2) filter cals by pairs
  if (infer_insert) {
    infer_insert_size(&pair_min_distance, &pair_max_distance, 
		      num_reads, cal_lists);
    mapping_batch->pair_min_distance = pair_min_distance; 
    mapping_batch->pair_max_distance = pair_max_distance; 
  }
  filter_cals_by_pair_mode(pair_min_distance, pair_max_distance, 
    			   num_reads, cal_lists);
  
  check_pairs(cal_lists, sa_index, mapping_batch, cal_mng);

  // 3) prepare Smith-Waterman to fill in the gaps
  for (int i = 0; i < num_reads; i++) {

    read = array_list_get(i, mapping_batch->fq_reads);
    cal_list = cal_lists[i];

    if (array_list_size(cal_list) > 0) {
      if (prepare_sw(read, sw_prepare_list, mapping_batch, sa_index, cal_list)) {
	sw_post_read[sw_post_read_counter++] = i;
      }
    }
  }

  // 4) run SW to fill
  if (array_list_size(sw_prepare_list) > 0) {
    execute_sw(sw_prepare_list, mapping_batch);
    post_process_sw(sw_post_read_counter, sw_post_read, cal_lists, mapping_batch);
  }
  array_list_free(sw_prepare_list, (void *) NULL);

  // 5) prepare mappings for the writer
  for (int i = 0; i < num_reads; i++) {
    cal_list = cal_lists[i];
    read = array_list_get(i, mapping_batch->fq_reads);

    //printf("-----------------> %s:%i: before create_alignments (%lu):\n", __FILE__, __LINE__, array_list_size(cal_list));
    //for (int i = 0; i < array_list_size(cal_list); i++) { seed_cal_print(array_list_get(i, cal_list)); }

    // create alignments structures
    #ifdef _TIMING
    gettimeofday(&start, NULL);
    #endif
    create_alignments(cal_list, read, bam_format, mapping_batch->mapping_lists[i]);
    #ifdef _TIMING
    gettimeofday(&stop, NULL);
    mapping_batch->func_times[FUNC_CREATE_ALIGNMENTS] += 
      ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
    #endif
      
    // free cal list and clear seed manager for next read
    array_list_free(cal_list, (void *) NULL);
  } // end of for reads
  
  complete_pairs(mapping_batch);

  // free memory
  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif
  cal_mng_free(cal_mng);
  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_OTHER] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif
  */
  return -1;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
