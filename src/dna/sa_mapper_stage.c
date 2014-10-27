#include "sa_mapper_stage.h"

//--------------------------------------------------------------------

#define SW_FLANK       10
#define SW_LEFT_FLANK  10
#define SW_RIGHT_FLANK 10

#define SW_LEFT_FLANK_EX  20
#define SW_RIGHT_FLANK_EX 20

#define CIGAR_FROM_GAP  1
#define CIGAR_FROM_SEED 2

#define MAX_GAP_LENGTH 2000

//--------------------------------------------------------------------
// process_right_side & append_seed_linked_list
//--------------------------------------------------------------------

void cut_adapter(char *adapter, int adapter_length, fastq_read_t *read);

int generate_cals_from_suffixes(int strand, fastq_read_t *read,
				int read_pos, int suffix_len, size_t low, size_t high, 
				sa_index3_t *sa_index, cal_mng_t *cal_mng
                                #ifdef _TIMING
				, sa_mapping_batch_t *mapping_batch
                                #endif
				);


void process_right_side(seed_t *new_item, linked_list_t *seed_list) {
  int overlap;
  linked_list_iterator_t* itr = linked_list_iterator_new(seed_list);
  seed_t *aux_item;
  seed_t *item = (seed_t *)linked_list_iterator_curr(itr);

  while (item != new_item) {
    //continue loop...
    linked_list_iterator_next(itr);
    item = linked_list_iterator_curr(itr);
  } // end while

  linked_list_iterator_next(itr);
  item = linked_list_iterator_curr(itr);
  while (item != NULL) {
    if (item->read_start > new_item->read_end) {
      break;
    } else if (item->read_end <= new_item->read_end) {
      aux_item = linked_list_iterator_remove(itr);
      if (aux_item) seed_free(aux_item);
    } else {
      overlap = new_item->read_end - item->read_start + 1;
      if (item->read_end - item->read_start - overlap > 18) {
	seed_ltrim_read(overlap + 2, item);
	linked_list_iterator_next(itr);
      } else {
	aux_item = linked_list_iterator_remove(itr);
	if (aux_item) seed_free(aux_item);
      }
    }

    //continue loop...
    item = linked_list_iterator_curr(itr);
  } // end while

  
  linked_list_iterator_free(itr);
}

//--------------------------------------------------------------------

void append_seed_linked_list(seed_cal_t *cal,
			     seed_t *new_item) {

  int process = 0;
  linked_list_t* seed_list = cal->seed_list;
  size_t read_start = new_item->read_start;
  size_t read_end = new_item->read_end; 
  size_t genome_start = new_item->genome_start;
  size_t genome_end = new_item->genome_end;

  int i, j, overlap;
  int seeds_to_delete, inserted = 0;
  seed_t *item, *aux_item;
  linked_list_iterator_t* itr = linked_list_iterator_new(seed_list);
  
  if (linked_list_size(seed_list) <= 0) {
    linked_list_insert(new_item, seed_list);
  } else {
    size_t num_seeds = linked_list_size(seed_list);

    for (i = 0; i < num_seeds; i++) {
      item = linked_list_get(i, seed_list);

      if (read_start < item->read_start) {
	if (read_end < item->read_start) {
	  /**
	        new |------|
	                        |-------|
	  */
	  linked_list_insert_at(i, new_item, seed_list);
	  process = 1;
	  break;
	} else if (read_end < item->read_end) {
	  /**
	        new   |-----------|
	                        |------------|         
	  */
	  overlap = read_end - item->read_start + 1;
	  if (read_end - read_start - overlap > 18) {
	    seed_rtrim_read(overlap + 2, new_item);
	    linked_list_insert_at(i, new_item, seed_list);
	  } else {
	    seed_free(new_item);
	  }
	  process = 1;
	  break;
	} else if (read_end == item->read_end) {
	  if (genome_end != item->genome_end) {
	    cal->invalid = 1;
	  }
	  /**
	        new   |--------------------|
	                      |------------|         
	  */
	  linked_list_insert_at(i, new_item, seed_list);
	  if (i == 0) cal->start = genome_start;
	  if (linked_list_size(seed_list) <= i + 1) {
	    linked_list_remove_last(seed_list);
	  } else {
	    linked_list_remove_at(i + 1, seed_list);
	  }
	  seed_free(item);
	  process = 1;
	  break;
	} else {
	  /**
	        new   |------------------------|
	                      |--------| |--|        
	  */
	  linked_list_insert_at(i, new_item, seed_list);
	  if (linked_list_size(seed_list) <= i + 1) {
	    linked_list_remove_last(seed_list);
	  } else {
	    linked_list_remove_at(i + 1, seed_list);
	  }
	  seed_free(item);
	  process_right_side(new_item, seed_list);
	  process = 1;
	  break;
	}
      } else if (read_start == item->read_start) {
	if (genome_start != item->genome_start) {
	  cal->invalid = 1;
	}
	if (read_end < item->read_end) {
	  /**
	        new   |-----------|
	              |----------------|         
	  */
	  seed_free(new_item);
	  process = 1;
	  break;
	} else if (read_end == item->read_end) {
	  if (genome_end != item->genome_end) {
	    cal->invalid = 1;
	  }
	  /**
	        new   |------------|
	              |------------|         
	  */
	  seed_free(new_item);
	  process = 1;
	  break;
	} else {
	  /**
	        new   |--------------------|
	              |--------|   |--|   |----|         
	  */
	  linked_list_insert_at(i, new_item, seed_list);
	  if (linked_list_size(seed_list) <= i + 1) {
	    linked_list_remove_last(seed_list);
	  } else {
	    linked_list_remove_at(i + 1, seed_list);
	  }
	  seed_free(item);
	  process_right_side(new_item, seed_list);
	  process = 1;
	  break;
	}
      } else if (read_start < item->read_end) {
	if (read_end < item->read_end) {
	  /**
	        new      |--------|
	              |----------------|         
	  */
	  seed_free(new_item);
	  process = 1;
	  break;
	} else if (read_end == item->read_end) {
	  if (genome_end != item->genome_end) {
	    cal->invalid = 1;
	  }
	  /**
	        new      |---------|
	              |------------|         
	  */
	  seed_free(new_item);
	  process = 1;
	  break;
	} else {
	  /**
	        new       |-----------------|
	              |--------|    |--|        |------|        
	  */
	  overlap = item->read_end - read_start + 1;
	  if (read_end - read_start - overlap > 18) {
	    seed_ltrim_read(overlap + 2, new_item);
	    if (linked_list_size(seed_list) <= i + 1) {
	      linked_list_insert_last(new_item, seed_list);
	    } else {
	      linked_list_insert_at(i + 1, new_item, seed_list);
	    }
	    process_right_side(new_item, seed_list);
	  } else {
	    seed_free(new_item);
	  }
	  process = 1;
	  break;
	}
      }
    }
    if (!process) {
      linked_list_insert_last(new_item, seed_list);
    }

    cal->start = ((seed_t*) linked_list_get_first(seed_list))->genome_start;
    cal->end = ((seed_t *) linked_list_get_last(seed_list))->genome_end;
  }

  linked_list_iterator_free(itr);
}

//--------------------------------------------------------------------
// cal_mng_t functions
//--------------------------------------------------------------------

cal_mng_t * cal_mng_new(sa_genome3_t *genome) {

  int num_chroms = genome->num_chroms;

  linked_list_t **cals_lists = (linked_list_t **) malloc (sizeof(linked_list_t *) * num_chroms);
  for (unsigned int i = 0; i < num_chroms; i++) {
    cals_lists[i] = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
  }

  cal_mng_t *p = (cal_mng_t *) calloc(1, sizeof(cal_mng_t));
  p->read_length = 10;
  p->min_read_area = 100;
  p->max_read_area = 0;
  p->num_chroms = num_chroms;
  p->cals_lists = cals_lists;

  p->suffix_mng = suffix_mng_new(genome);

  return p;
}

//--------------------------------------------------------------------

void cal_mng_free(cal_mng_t *p) {
  if (p) {
    if (p->cals_lists) {
      for (unsigned int i = 0; i < p->num_chroms; i++) {
	if (p->cals_lists[i]) {
	  linked_list_free(p->cals_lists[i], (void *)seed_cal_free);
	}
      }
      free(p->cals_lists);
    }
    if (p->suffix_mng) suffix_mng_free(p->suffix_mng);

    free(p);
  }
}

//--------------------------------------------------------------------

void cal_mng_simple_free(cal_mng_t *p) {
  if (p) {
    if (p->cals_lists) {
      for (unsigned int i = 0; i < p->num_chroms; i++) {
	if (p->cals_lists[i]) {
	  linked_list_free(p->cals_lists[i], (void *)NULL);
	}
      }
      free(p->cals_lists);
    }
    free(p);
  }
}

//--------------------------------------------------------------------

void cal_mng_simple_clear(cal_mng_t *p) {
  int list_size;
  linked_list_item_t *item;
  linked_list_t *list;
  cal_t *cal;
  if (p) {
    if (p->cals_lists) {
      for (unsigned int i = 0; i < p->num_chroms; i++) {
	list = p->cals_lists[i];
	if (list) {
	  item = list->first;
	  while (item) {
	    cal = item->item;
	    linked_list_clear(cal->sr_list, (void *)seed_region_simple_free);
	    cal_simple_free(cal);
	    item = item->next;
	  }
	  linked_list_clear(p->cals_lists[i], (void *)NULL);
	}
      }
    }
  }
}

//--------------------------------------------------------------------

void cal_mng_clear(cal_mng_t *p) {
  if (p) {
    if (p->cals_lists) {
      for (unsigned int i = 0; i < p->num_chroms; i++) {
	if (p->cals_lists[i]) {
	  linked_list_clear(p->cals_lists[i], (void *)seed_cal_free);
	}
      }
    }
  }
}

//--------------------------------------------------------------------


void cal_mng_update(seed_t *seed, fastq_read_t *read, cal_mng_t *p) {
  if (p->cals_lists) {
    seed_cal_t *cal;
    seed_t *s_last;
    linked_list_t *seed_list;
    linked_list_t *cal_list = p->cals_lists[seed->chromosome_id];
    if (cal_list) {
      int r_gap, g_gap;
      #ifdef _VERBOSE
      printf("\t\t\tinsert this seed to the CAL manager:\n");
      print_seed("\t\t\t", seed);
      #endif
      if (linked_list_size(cal_list) <= 0) {
	// create CAL and insert it into the CAL manager
	seed_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	linked_list_insert(seed, seed_list);
	cal = seed_cal_new(seed->chromosome_id, seed->strand, 
			   seed->genome_start, seed->genome_end, seed_list);
	cal->read_area = seed->read_end - seed->read_start + 1;
	cal->num_mismatches = seed->num_mismatches + seed->num_open_gaps + seed->num_extend_gaps;
	cal->read = read;
	linked_list_insert(cal, cal_list);
      } else {
	// insert (by order)
	linked_list_iterator_t* itr = linked_list_iterator_new(cal_list);
	seed_cal_t *item = (seed_cal_t *) linked_list_iterator_curr(itr);
	while (item != NULL) {
	  #ifdef _VERBOSE
	  printf("---> merging with this CAL?\n");
	  seed_cal_print(item);
	  #endif
	  s_last = linked_list_get_last(item->seed_list);
	  if (s_last->suf_read_start != seed->suf_read_start &&
	      s_last->suf_read_end != seed->suf_read_end) {
	    r_gap = abs(seed->read_start - s_last->read_end);
	    g_gap = abs(seed->genome_start - s_last->genome_end);
	    if (g_gap <= read->length && r_gap <= read->length) {
	      if (abs(r_gap - g_gap) < 200) {
		append_seed_linked_list(item, seed);
		break;
	      }
	    }
	  }
	  if (seed->genome_start < item->start) {
	    seed_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	    linked_list_insert(seed, seed_list);
	    cal = seed_cal_new(seed->chromosome_id, seed->strand, 
			       seed->genome_start, seed->genome_end, seed_list);
	    cal->read_area = seed->read_end - seed->read_start + 1;
	    cal->num_mismatches = seed->num_mismatches + seed->num_open_gaps + seed->num_extend_gaps;
	    cal->read = read;
	    linked_list_iterator_insert(cal, itr);
	    linked_list_iterator_prev(itr);
	    break;
	  }

	  //continue loop...
	  linked_list_iterator_next(itr);
	  item = linked_list_iterator_curr(itr);
	}
	if (item == NULL) {
	  // create CAL and insert it into the CAL manager
	  seed_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	  linked_list_insert(seed, seed_list);
	  cal = seed_cal_new(seed->chromosome_id, seed->strand, 
			     seed->genome_start, seed->genome_end, seed_list);
	  cal->read_area = seed->read_end - seed->read_start + 1;
	  cal->num_mismatches = seed->num_mismatches + seed->num_open_gaps + seed->num_extend_gaps;
	  cal->read = read;
	  linked_list_insert_last(cal, cal_list);
	}
	linked_list_iterator_free(itr);
      }
    }
  }
}

//--------------------------------------------------------------------

int cal_mng_find(int strand, unsigned int chrom, size_t start, size_t end, cal_mng_t *p) {
  #ifdef _VERBOSE1
  printf("\t\t***** searching CAL: chrom %u: %lu-%lu\n", chrom, start, end);
  #endif
  int found_cal = 0;
  if (p->cals_lists) {
    linked_list_t *cal_list = p->cals_lists[chrom];
    if (cal_list) {
      seed_cal_t *cal;
      for (linked_list_item_t *item = cal_list->first; 
	   item != NULL; 
	   item = item->next) {
	cal = item->item;
	#ifdef _VERBOSE1
	printf("\t\t\t***** searching CAL: suf. seed %lu-%lu is included in cal %c:%i:%lu-%lu\n", 
	       start, end, (cal->strand == 0 ? '+' : '-'), cal->chromosome_id, cal->start, cal->end);
        #endif
	if (cal->strand == strand && cal->start <= start && cal->end >= end) {
	  found_cal = 1;
	  break;
	} 
	if (cal->start > end) {
	  break;
	}
      }
    }
  }
  #ifdef _VERBOSE1
  printf("\t\t\t\t***** searching CAL: found_cal = %i\n", found_cal);
  #endif
  return found_cal;
}

//--------------------------------------------------------------------

void cal_mng_to_array_list(int min_read_area, array_list_t *out_list, cal_mng_t *p) {
  seed_t *first, *last;
  seed_cal_t *cal;
  linked_list_iterator_t itr;

  #ifdef _VERBOSE
  printf("-----> cal_mng_to_array_list\n");
  #endif

  if (p->cals_lists) {
    linked_list_t *cal_list;
    for (unsigned int i = 0; i < p->num_chroms; i++) {
      cal_list = p->cals_lists[i];
      while (cal = (seed_cal_t *) linked_list_remove_last(cal_list)) {
	#ifdef _VERBOSE
	seed_cal_print(cal);
	#endif
	first = linked_list_get_first(cal->seed_list);
	last = linked_list_get_last(cal->seed_list);
	cal->start = first->genome_start;
	cal->end = last->genome_end;
	seed_cal_update_info(cal);
	if (cal->read_area >= min_read_area &&
	    cal->num_open_gaps < (0.05f * cal->read->length) &&
	    cal->num_mismatches < (0.09f * cal->read->length) ) {
	  array_list_insert(cal, out_list);
	} else {
	  // free CAL
	  seed_cal_free(cal);
	}
      }
    }
  }
}

//--------------------------------------------------------------------

void cal_mng_select_best(int read_area, array_list_t *valid_list, array_list_t *invalid_list, 
			 cal_mng_t *p) {
  seed_cal_t *cal;
  linked_list_iterator_t itr;

  if (p->cals_lists) {
    linked_list_t *cal_list;
    for (unsigned int i = 0; i < p->num_chroms; i++) {
      cal_list = p->cals_lists[i];
      while (cal = (seed_cal_t *) linked_list_remove_last(cal_list)) {
	if (p->min_read_area <= read_area && cal->read_area <= read_area) {
	  array_list_insert(cal, valid_list);
	} else {
	  array_list_insert(cal, invalid_list);
	}
      }
    }
  }
}


//--------------------------------------------------------------------
//
//--------------------------------------------------------------------

array_list_t *search_mate_cal_by_prefixes(seed_cal_t *cal, fastq_read_t *read,
					  sa_index3_t *sa_index, sa_mapping_batch_t *batch,
					  cal_mng_t *cal_mng) {
  array_list_t *cal_list = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  char *seq;

  unsigned int chromosome;
  int start, end;
  int read_pos, read_inc, read_end_pos;

  int mapq = cal->mapq;

  int num_seeds = batch->options->num_seeds;

  int min_distance = batch->options->pair_min_distance;
  int max_distance = batch->options->pair_max_distance;
  
  read_end_pos = read->length - sa_index->k_value;
  read_inc = read->length / num_seeds;
  if (read_inc < sa_index->k_value / 2) {
    read_inc = sa_index->k_value / 2;
  }
  
  chromosome = cal->chromosome_id;
  start = cal->start - max_distance - read->length;
  end = cal->end + max_distance + read->length;
  seq = (cal->strand ? read->sequence : read->revcomp);
  
  size_t num_prefixes, low, high;
  size_t g_start_suf, g_end_suf;
  unsigned int chrom;

  for (read_pos = 0; read_pos < read_end_pos; read_pos += read_inc)  {	

    num_prefixes = search_prefix(&seq[read_pos], &low, &high, sa_index, 0);
    if (num_prefixes <= 0) continue;

    for (size_t i = low; i <= high; i++) {
      chrom = (unsigned int) sa_index->CHROM[i];
      if (chrom == chromosome) {
	g_start_suf = sa_index->SA[i] - sa_index->genome->chrom_offsets[chrom];
	g_end_suf = g_start_suf + sa_index->k_value - 1;
	
	if (start <= g_start_suf && end >= g_end_suf) {
	  generate_cals_from_suffixes(1 - cal->strand, read, read_pos,
				      sa_index->k_value, i, i,
				      sa_index, cal_mng
                                      #ifdef _TIMING
				      , batch
                                      #endif
				      );
	}
      }
    }
  }

  cal_mng_to_array_list(read->length / 3, cal_list, cal_mng);

  if (array_list_size(cal_list) > 0) {
    select_best_cals(read, &cal_list);

    if (array_list_size(cal_list) <= 0) {
      suffix_mng_search_read_cals_by_region(read, num_seeds, sa_index, 1 - cal->strand,
					  chromosome, start, end, cal_list, cal_mng->suffix_mng);
      if (array_list_size(cal_list) > 0) {
      	select_best_cals(read, &cal_list);
      }
    }

    for (int kk = 0; kk < array_list_size(cal_list); kk++) { 
      cal = array_list_get(kk, cal_list);
      if (cal->mapq > 0) {
	cal->mapq = mapq;
      }
    }
  }
  
  return cal_list;
}

//--------------------------------------------------------------------
int is_valid_cal_pair(seed_cal_t *cal1, seed_cal_t *cal2,
		      int min_distance, int max_distance, 
		      int *out_distance) {
  int distance, valid_pair = 0;
  if (cal1->chromosome_id == cal2->chromosome_id) {
    if (cal2->start > cal1->end) {
      distance = cal2->end - cal1->start + 1;
    } else {
      distance = cal1->end - cal2->start + 1;
    }
    if (distance >= min_distance && distance <= max_distance) {
      valid_pair = 1;
      *out_distance = distance;
    }
  }
  return valid_pair;
}

//--------------------------------------------------------------------

void check_pairs(array_list_t **cal_lists, sa_index3_t *sa_index,
		 sa_mapping_batch_t *batch, cal_mng_t *cal_mng) {
  int found, inserted, max_read_area, read_area;

  int score, first_score, second_score;
  int distance, valid_pair, list_size, list1_size, list2_size;
  seed_cal_t *cal, *mate_cal, *cal1, *cal2;
  array_list_t *list, *mate_list, *mate1_list, *mate2_list, *list1, *list2;
  fastq_read_t *read, *read1, *read2;
  size_t mate_list_size, num_cals, num_reads = array_list_size(batch->fq_reads);

  int min_distance = batch->options->pair_min_distance;
  int max_distance = batch->options->pair_max_distance;

  for (int i = 0; i < num_reads; i += 2) {
    read1 = array_list_get(i, batch->fq_reads);
    read2 = array_list_get(i+1, batch->fq_reads);

    list1 = cal_lists[i];
    list2 = cal_lists[i+1];

    list1_size = array_list_size(list1);
    list2_size = array_list_size(list2);

    if (list1_size == 0 && list2_size == 0) continue;

    // no mate #1
    if (list1_size == 0) {
      // search mate by prefixes
      array_list_t *new_list = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
      for (int i2 = 0; i2 < list2_size; i2++) {
	cal2 = array_list_get(i2, list2);
	if (i2 == 0 && cal2->mapq < 10) break;

	list = search_mate_cal_by_prefixes(cal2, read1, sa_index, batch, cal_mng);
	if (array_list_size(list) > 0) {
	  for (int kk = 0; kk < array_list_size(list); kk++) { 
	    cal = array_list_get(kk, list); 
	    if (is_valid_cal_pair(cal1, cal, min_distance, max_distance, &distance)) {
	      array_list_insert(cal, new_list);
	      array_list_set(kk, NULL, list);
	    }
	  }
	}
	array_list_free(list, (void *) seed_cal_free);
      }
      if (array_list_size(new_list) > 0) {
	batch->status[i] = 5; // found mate #2 from mate #1
	cal_lists[i] = new_list;
	array_list_free(list1, (void *) NULL);
      } else {
	array_list_free(new_list, (void *) NULL);
      }
      continue;
    }

    // no mate #2
    if (list2_size == 0) {
      // search mate by prefixes
      array_list_t *new_list = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
      for (int i1 = 0; i1 < list1_size; i1++) {
	cal1 = array_list_get(i1, list1);
	if (i1 == 0 && cal1->mapq < 10) break;
	list = search_mate_cal_by_prefixes(cal1, read2, sa_index, batch, cal_mng);
	if (array_list_size(list) > 0) {
	  for (int kk = 0; kk < array_list_size(list); kk++) { 
	    cal = array_list_get(kk, list);
	    if (is_valid_cal_pair(cal1, cal, min_distance, max_distance, &distance)) {
	      array_list_insert(cal, new_list);
	      array_list_set(kk, NULL, list);
	    }
	  }
	}
	array_list_free(list, (void *) seed_cal_free);
      }
      if (array_list_size(new_list) > 0) {
	batch->status[i+1] = 6; // found mate #1 from mate #2
	cal_lists[i+1] = new_list;
	array_list_free(list2, (void *) NULL);
      } else {
	array_list_free(new_list, (void *) NULL);
      }
      continue;
    }

    // check if we have valid pairs
    score = 0;
    first_score = 0;
    second_score = 0;
    valid_pair = 0;
    for (int i1 = 0; i1 < list1_size; i1++) {
      cal1 = array_list_get(i1, list1);
      for (int i2 = 0; i2 < list2_size; i2++) {
	cal2 = array_list_get(i2, list2);
	if (is_valid_cal_pair(cal1, cal2, min_distance, max_distance, &distance)) {
	  // score management
	  score = cal1->score + cal2->score;
	  if (score < first_score) {
	    if (score > second_score) {
	      second_score = score;
	    }
	  } else if (score > first_score) {
	    if (first_score > 0) {
	      second_score = first_score;
	    }
	    first_score = score;
	  } else {
	    second_score = score;
	  }
	  
	  valid_pair++;
	}
      }
    }
    if (valid_pair) {
      continue;
    }


    if (list1_size > 1 && list2_size > 1) {
      array_list_clear(list1, (void *) seed_cal_free);
      array_list_clear(list2, (void *) seed_cal_free);
      continue;
    }

    // no valid pairs, then search mate by region
    max_read_area = 0;
    score = 0;
    first_score = 0;
    second_score = 0;

    mate1_list = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    mate2_list = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

    for (int k = 0; k < 2; k++) {
      if  (k == 0) {
	list_size = list1_size;
	list = list1;
	read = read2;
      } else {
	list_size = list2_size;
	list = list2;
	read = read1;
      }
      for (int ii = 0; ii < list_size; ii++) {
	inserted = 0;
	cal = array_list_get(ii, list);
	mate_list = search_mate_cal_by_prefixes(cal, read, sa_index, batch, cal_mng);
	mate_list_size = array_list_size(mate_list);
	int num_found_mates = 0;
	int found_mates[mate_list_size];
	for (int jj = 0; jj < mate_list_size; jj++) { 
	  mate_cal = array_list_get(jj, mate_list);
	  
	  // score management
	  score = cal->score + mate_cal->score;
	  if (score < first_score) {
	    if (score > second_score) {
	      second_score = score;
	    }
	  } else if (score > first_score) {
	    if (first_score > 0) {
	      second_score = first_score;
	    }
	    first_score = score;
	  } else {
	    second_score = score;
	  }
	  
	  // read area management
	  read_area = cal->read_area + mate_cal->read_area;
	  
	  if (read_area > max_read_area) {
	    // remove the current mates from the mate lists, and insert this pair

	    // in mate1_list, we insert the new mate found
	    array_list_clear(mate1_list, (void *) seed_cal_free);
	    array_list_set(jj, NULL, mate_list);
	    array_list_insert(mate_cal, mate1_list);

	    // now we update mate2_list
	    array_list_clear(mate2_list, (void *) NULL);
	    array_list_insert(cal, mate2_list);
	    inserted = 1;

	    // update max read area
	    max_read_area = read_area;	      
	  } else if (read_area == max_read_area) {

	    // in mate1_list, we insert the new mate found
	    array_list_set(jj, NULL, mate_list);
	    array_list_insert(mate_cal, mate1_list);

	    // now, we update mate2_list
	    if (!inserted) {
	      array_list_insert(cal, mate2_list);
	      inserted = 1;
	    }
	  }
	}
	// free memory
	array_list_free(mate_list, (void *) seed_cal_free);
      }
    }
   
    if (max_read_area == 0) {
      // free memory
      array_list_free(mate1_list, (void *) seed_cal_free);
      array_list_free(mate2_list, (void *) NULL);
      continue;
    }

    // update cal lists with the found mate lists
    list_size = array_list_size(mate2_list);
    for (int ii = 0; ii < list1_size; ii++) { 
      cal1 = array_list_get(ii, list1); 
      found = 0;
      for (int jj = 0; jj < list_size; jj++) { 
	cal2 = array_list_get(jj, mate2_list);
	if (cal1 == cal2) {
	  found = 1;
	  break;
	}
      }
      if (found) {
	array_list_set(ii, NULL, list1);	
      }
    }
    array_list_free(list1, (void *) seed_cal_free);

    list_size = array_list_size(mate2_list);
    for (int ii = 0; ii < list2_size; ii++) { 
      cal1 = array_list_get(ii, list2); 
      found = 0;
      for (int jj = 0; jj < list_size; jj++) { 
	cal2 = array_list_get(jj, mate2_list);
	if (cal1 == cal2) {
	  found = 1;
	  break;
	}
      }
      if (found) {
	array_list_set(ii, NULL, list2);	
      }
    }
    array_list_free(list2, (void *) seed_cal_free);

    cal_lists[i] = mate1_list;
    cal_lists[i+1] = mate2_list;

    batch->status[i] = 7; 
    batch->status[i+1] = 7;
  }
}


//--------------------------------------------------------------------
// generate cal from an exact read
//--------------------------------------------------------------------

void generate_cals_from_exact_read(int strand, fastq_read_t *read,
				   size_t low, size_t high, sa_index3_t *sa_index, 
				   cal_mng_t *cal_mng) {
  size_t g_start, g_end;
  unsigned int chrom;

  seed_cal_t *cal;
  cigar_t *cigar;
  seed_t *seed;
  
  for (size_t suff = low; suff <= high; suff++) {
    chrom = (unsigned int) sa_index->CHROM[suff];
    g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];
    g_end = g_start + read->length - 1;

    //    seed_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    seed = seed_new(0, read->length - 1, g_start, g_end);
    seed->chromosome_id = chrom;
    seed->strand = strand;
    cigar_append_op(read->length, '=', &seed->cigar);

    cal_mng_update(seed, read, cal_mng);
  }
}

//--------------------------------------------------------------------
// generate cals extending suffixes to left and right side 
// bt using the mini-sw
//--------------------------------------------------------------------

int generate_cals_from_suffixes(int strand, fastq_read_t *read,
				int read_pos, int suffix_len, size_t low, size_t high, 
				sa_index3_t *sa_index, cal_mng_t *cal_mng
                                #ifdef _TIMING
				, sa_mapping_batch_t *mapping_batch
                                #endif
				) {

  #ifdef _TIMING
  struct timeval stop, start;
  #endif

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif

  size_t r_start_suf, r_end_suf, g_start_suf, g_end_suf;
  size_t r_start, r_end, r_len, g_start, g_end, g_len;
  int found_cal, diff, max_map_len = 0;
  unsigned int chrom;

  float score;
  alig_out_t alig_out;
  cigar_init(&alig_out.cigar);

  cigar_t cigar;
  seed_t *seed;

  char *g_seq, *r_seq;
  r_seq = (strand ? read->revcomp : read->sequence);
  
  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_INIT_CALS_FROM_SUFFIXES] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif

  for (size_t suff = low; suff <= high; suff++) {
    #ifdef _TIMING
    gettimeofday(&start, NULL);
    #endif
    chrom = (unsigned int) sa_index->CHROM[suff];
    if (chrom < 0) {
      printf("chrom %c %i %u is < 0\n", chrom, chrom, chrom);
      exit(-1);
    }

    // extend suffix to right side
    r_start_suf = read_pos;
    r_end_suf = r_start_suf + suffix_len - 1;
    
    g_start_suf = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];
    g_end_suf = g_start_suf + suffix_len - 1;

    #ifdef _TIMING
    gettimeofday(&stop, NULL);
    mapping_batch->func_times[FUNC_SET_POSITIONS] += 
      ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
    #endif

    // skip suffixes, 
    // if found cal for this suffix, then next suffix
    #ifdef _TIMING
    gettimeofday(&start, NULL);
    #endif
    found_cal = cal_mng_find(strand, chrom, g_start_suf, g_end_suf, cal_mng);
    #ifdef _TIMING
    gettimeofday(&stop, NULL);
    mapping_batch->func_times[FUNC_SKIP_SUFFIXES] += 
      ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
    #endif
    if (found_cal) continue;


    //#ifdef _VERBOSE
    //printf("%s:%i:\t\tsuffix at [%lu|%lu-%lu|%lu] %c chrom %s\n", __FILE__, __LINE__,
    //	   g_start_suf, r_start_suf, r_end_suf, g_end_suf, (strand == 0 ? '+' : '-'), 
    //	   sa_index->genome->chrom_names[chrom]);
    //#endif

    seed = seed_new(r_start_suf, r_end_suf, g_start_suf, g_end_suf);

    // extend suffix to left side, if necessary
    if (r_start_suf > 0) {
      r_start = 0;
      r_end = r_start_suf - 1;
      r_len = r_start_suf;

      g_len = r_len + 5;
      g_end = g_start_suf - 1;
      g_start = g_end - g_len;

      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      g_seq = &sa_index->genome->S[g_start + sa_index->genome->chrom_offsets[chrom] + 1];
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_SET_REF_SEQUENCE] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif

      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      score = doscadfun_inv(r_seq, r_len, g_seq, g_len, MISMATCH_PERC,
			    &alig_out);
      //printf("%s:%i********** doscadfun_inv (score = %0.2f) cigar: %s\n", 
      //	     __FILE__, __LINE__, score, cigar_to_string(&alig_out.cigar));
			    
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_MINI_SW_LEFT_SIDE] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif
      if (score > 0.0f) {

	// update seed
	seed->num_mismatches += alig_out.mismatch;
	seed->num_open_gaps += alig_out.gap_open;
	seed->num_extend_gaps += alig_out.gap_extend;

	seed->read_start -= alig_out.map_len1;
	seed->genome_start -= alig_out.map_len2;

	// if there's a mini-gap then try to fill the mini-gap
	if (seed->read_start > 0 && seed->read_start < 5) {
          #ifdef _TIMING
	  gettimeofday(&start, NULL);
          #endif
	  g_seq = &sa_index->genome->S[seed->genome_start + sa_index->genome->chrom_offsets[chrom] - seed->read_start];
	  for (size_t k1 = 0, k2 = 0; k1 < seed->read_start; k1++, k2++) {
	    if (r_seq[k1] != g_seq[k2]) {
	      seed->num_mismatches++;
	      cigar_append_op(1, 'X', &seed->cigar);
	    } else {
	      cigar_append_op(1, '=', &seed->cigar);
	    }
	  }
	  
	  seed->genome_start -= seed->read_start;
	  seed->read_start = 0;
          #ifdef _TIMING
	  gettimeofday(&stop, NULL);
	  mapping_batch->func_times[FUNC_SEED_NEW] += 
	    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
          #endif
	}


	// update cigar with the sw output
	if (alig_out.cigar.num_ops > 0) {
	  cigar_concat(&alig_out.cigar, &seed->cigar);
	}
      }
    }
    // update cigar with the suffix length
    cigar_append_op(suffix_len, '=', &seed->cigar);

    // extend suffix to left side, if necessary
    if (r_end_suf < read->length - 1) {
      r_start = r_end_suf + 1;
      r_end = read->length - 1;
      r_len = r_end - r_start + 1;

      g_len = r_len + 5;
      g_start = g_end_suf + 1;
      g_end = g_start + g_len;

      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      g_seq = &sa_index->genome->S[g_start + sa_index->genome->chrom_offsets[chrom]];
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_SET_REF_SEQUENCE] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif

      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      score = doscadfun(&r_seq[r_start], r_len, g_seq, g_len, MISMATCH_PERC,
			&alig_out);
      //printf("%s:%i********** doscadfun (score = %0.2f) cigar: %s\n", 
      //	     __FILE__, __LINE__, score, cigar_to_string(&alig_out.cigar));
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_MINI_SW_RIGHT_SIDE] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif
      if (score > 0.0f) {

	// update seed
	seed->num_mismatches += alig_out.mismatch;
	seed->num_open_gaps += alig_out.gap_open;
	seed->num_extend_gaps += alig_out.gap_extend;

	seed->read_end += alig_out.map_len1;
	seed->genome_end += alig_out.map_len2;

	// update cigar with the sw output
	if (alig_out.cigar.num_ops > 0) {
	  cigar_concat(&alig_out.cigar, &seed->cigar);
	}
	// if there's a mini-gap then try to fill the mini-gap
	diff = read->length - seed->read_end - 1;
	if (diff > 0 && diff < 5) {
          #ifdef _TIMING
	  gettimeofday(&start, NULL);
          #endif
	  g_seq = &sa_index->genome->S[seed->genome_end + sa_index->genome->chrom_offsets[chrom] + 1];
	  for (size_t k1 = seed->read_end + 1, k2 = 0; k1 < read->length; k1++, k2++) {
	    if (r_seq[k1] != g_seq[k2]) {
	      seed->num_mismatches++;
	      cigar_append_op(1, 'X', &seed->cigar);
	    } else {
	      cigar_append_op(1, '=', &seed->cigar);
	    }
	  }

	  seed->read_end += diff;
	  seed->genome_end += diff;
          #ifdef _TIMING
	  gettimeofday(&stop, NULL);
	  mapping_batch->func_times[FUNC_SEED_NEW] += 
	    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
          #endif
	}
      }
    }

    //printf("%s:%c:\t seed final cigar: %s\n", __FILE__, __LINE__, cigar_to_string(&seed->cigar));

    // update CAL manager with this seed
    if (seed->read_end - seed->read_start + 1 > 20) {
      seed->strand = strand;
      seed->chromosome_id = chrom;
      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif

      cal_mng_update(seed, read, cal_mng);

      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_CAL_MNG_INSERT] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif
    } else {
      // free seed
      seed_free(seed);
    }
  }

  cigar_clean(&alig_out.cigar);

  return 0;
}

//--------------------------------------------------------------------
// create_cals function:
//    search prefix -> search longer suffix -> extend suffix
//--------------------------------------------------------------------

array_list_t *create_cals(int num_seeds, fastq_read_t *read,
			  sa_mapping_batch_t *mapping_batch, 
			  sa_index3_t *sa_index, cal_mng_t *cal_mng) {

  #ifdef _TIMING
  struct timeval stop, start;
  #endif

  double prefix_time, suffix_time;
  size_t suffix_len, num_suffixes;
  char *r_seq = read->sequence;

  size_t low, high;

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif

  array_list_t *cal_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  cal_mng->min_read_area = read->length;
  cal_mng->read_length = read->length;

  int max_read_area;
  int read_pos, read_inc;
  
  read_inc = read->length / num_seeds;
  if (read_inc < sa_index->k_value / 2) {
    read_inc = sa_index->k_value / 2;
  }

  // fill in the CAL manager structure
  int read_end_pos = read->length - sa_index->k_value;
  int extra_seed = (read->length - sa_index->k_value) % read_inc;

  #ifdef _VERBOSE	  
  printf("\n\n====>>>> STEP ONE <<<<====\n");
  display_cmp_sequences(read, sa_index);
  #endif

  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_OTHER] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif

  // first step, searching mappings in both strands
  // distance between seeds >= prefix value (sa_index->k_value)
  for (int strand = 0; strand < 2; strand++) {
    #ifdef _VERBOSE
    printf("=======> STRAND %c (read_end_pos = %i, %s %s)\n", 
	   (strand == 0 ? '+' : '-'), read_end_pos, read->id, read->sequence);
    #endif

    for (read_pos = 0; read_pos < read_end_pos; )  {	
      #ifdef _VERBOSE	  
      printf("\tread pos. = %lu\n", read_pos);
      #endif

      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      num_suffixes = search_suffix(&r_seq[read_pos], sa_index->k_value, 
				   MAX_NUM_SUFFIXES, sa_index, 
				   &low, &high, &suffix_len
                                   #ifdef _TIMING
				   , &prefix_time, &suffix_time
                                   #endif
				   );
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_SEARCH_SUFFIX] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  

      mapping_batch->func_times[FUNC_SEARCH_PREFIX] += prefix_time;
      mapping_batch->func_times[FUNC_SEARCH_SA] +=  suffix_time;
      #endif
      
      #ifdef _VERBOSE	  
      printf("\t\tnum. suffixes = %lu (suffix length = %lu)\n", num_suffixes, suffix_len);
      #endif

      if (num_suffixes < MAX_NUM_SUFFIXES && suffix_len) {
        #ifdef _VERBOSE	  
	display_suffix_mappings(strand, read_pos, suffix_len, low, high, sa_index);
        #endif 
	
	// exact search
	if (suffix_len == read->length) {
          #ifdef _TIMING
	  gettimeofday(&start, NULL);
          #endif
	  generate_cals_from_exact_read(strand, read, low, high, 
					sa_index, cal_mng);
          #ifdef _TIMING
	  gettimeofday(&stop, NULL);
	  mapping_batch->func_times[FUNC_CALS_FROM_EXACT_READ] += 
	    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
          #endif
	  break;
	} else {
          #ifdef _TIMING
	  gettimeofday(&start, NULL);
          #endif
	  generate_cals_from_suffixes(strand, read,
				      read_pos, suffix_len, low, high, sa_index, cal_mng
                                      #ifdef _TIMING
				      , mapping_batch
                                      #endif
				      );
          #ifdef _TIMING
	  gettimeofday(&stop, NULL);
	  mapping_batch->func_times[FUNC_CALS_FROM_SUFFIXES] += 
	    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
          #endif
	}
      }
      read_pos += read_inc;
    } // end of for read_pos
    
    if (suffix_len != read->length && extra_seed) {
      read_pos = read->length - sa_index->k_value;
      #ifdef _VERBOSE	  
      printf("\tread pos. = %lu\n", read_pos);
      #endif

      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      num_suffixes = search_suffix(&r_seq[read_pos], sa_index->k_value, 
				   MAX_NUM_SUFFIXES, sa_index, 
				   &low, &high, &suffix_len
                                   #ifdef _TIMING
				   , &prefix_time, &suffix_time
                                   #endif
				   );
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_SEARCH_SUFFIX] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  

      mapping_batch->func_times[FUNC_SEARCH_PREFIX] += prefix_time;
      mapping_batch->func_times[FUNC_SEARCH_SA] +=  suffix_time;
      #endif
      
      #ifdef _VERBOSE	  
      printf("\t\tnum. suffixes = %lu (suffix length = %lu)\n", num_suffixes, suffix_len);
      #endif
      if (num_suffixes < MAX_NUM_SUFFIXES && suffix_len) {
        #ifdef _VERBOSE	  
	display_suffix_mappings(strand, read_pos, suffix_len, low, high, sa_index);
        #endif 
	
        #ifdef _TIMING
	gettimeofday(&start, NULL);
        #endif
	read_pos += generate_cals_from_suffixes(strand, read,
						read_pos, suffix_len, low, high, sa_index, cal_mng
                                                #ifdef _TIMING
						, mapping_batch
                                                #endif
						);
        #ifdef _TIMING
	gettimeofday(&stop, NULL);
	mapping_batch->func_times[FUNC_CALS_FROM_SUFFIXES] += 
	  ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
        #endif
      }
    }

    // update cal list from cal manager
    #ifdef _TIMING
    gettimeofday(&start, NULL);
    #endif

    cal_mng_to_array_list(read->length / 3, cal_list, cal_mng);
    #ifdef _TIMING
    gettimeofday(&stop, NULL);
    mapping_batch->func_times[FUNC_CAL_MNG_TO_LIST] += 
      ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
    #endif
    
    // next, - strand
    r_seq = read->revcomp;
  } // end of for strand
  

  #ifdef _VERBOSE
  printf("\t***> max_read_area = %i\n", max_read_area);
  #endif
  
  return cal_list;
}

//--------------------------------------------------------------------
//
//--------------------------------------------------------------------

void fill_seed_gaps(array_list_t *cal_list, fastq_read_t *read, sa_index3_t *sa_index) {

  seed_t *prev_seed, *seed;
  linked_list_item_t *prev_item, *item;

  seed_cal_t *cal;
  size_t num_seeds, num_cals = array_list_size(cal_list);

  float score;
  char *g_seq, *r_seq;
  alig_out_t alig_out;
  size_t r_start, r_end, r_len, g_start, g_end, g_len;
  cigar_init(&alig_out.cigar);


  for (int i = 0; i < num_cals; i++) {
    cal = array_list_get(i, cal_list);
    num_seeds = linked_list_size(cal->seed_list);

    // first seed
    seed = linked_list_get_first(cal->seed_list);

    // middle seeds
    prev_seed = seed;
    for (int i = 1; i < num_seeds; i++) {
      seed = linked_list_get(i, cal->seed_list);

      // update previous seed
      prev_seed = seed;
    }

    // last seed
  }

  cigar_clean(&alig_out.cigar);
}

//--------------------------------------------------------------------

inline static int is_invalid_cigar(cigar_t *cigar) {
  int name, value, count = 0;
  int min = SW_FLANK * 2;

  for (int i = 0; i < cigar->num_ops; i++) {
    cigar_get_op(i, &value, &name, cigar);
    if (value > SW_FLANK && name == '=') {
      count += value;
      if (count > min) return 0;
    }
  }  
  return (count > min ? 0 : 1);
}


int check_gap_lengths(int max_gap_length, array_list_t *cal_list) {
  int invalid = 0;
  int gap_read_len, gap_genome_len;
  size_t gap_read_start, gap_read_end;
  size_t gap_genome_start, gap_genome_end;
  size_t seed_count, num_seeds, num_cals = array_list_size(cal_list);

  seed_cal_t *cal;
  seed_t *prev_seed, *seed;
  linked_list_item_t *item;

  for (int i = 0; i < num_cals; i++) {
    cal = array_list_get(i, cal_list);
    if (cal->invalid) {
      invalid = 1;
      continue;
    }
    num_seeds = cal->seed_list->size;
    if (num_seeds <= 0) {
      invalid = 1;
      cal->invalid = 1;
      continue;
    }

    // first seed
    seed = linked_list_get_first(cal->seed_list);
    if (seed->read_start > max_gap_length) {
      invalid = 1;
      cal->invalid = 1;
      continue;
    }
    // seeds at the middle positions
    prev_seed = seed;
    for (item = cal->seed_list->first->next; item != NULL; item = item->next) {
      seed = item->item;

      // gap in read
      gap_read_start = prev_seed->read_end + 1;
      gap_read_end = seed->read_start - 1;
      gap_read_len = gap_read_end - gap_read_start + 1;    
      if (abs(gap_read_len) > max_gap_length) {
	invalid = 1;
	cal->invalid = 1;
	continue;      
      }

      // gap in genome
      gap_genome_start = prev_seed->genome_end + 1;
      gap_genome_end = seed->genome_start - 1;
      gap_genome_len = gap_genome_end - gap_genome_start + 1;
      if (abs(gap_genome_len) > max_gap_length) {
	invalid = 1;
	cal->invalid = 1;
	continue;      
      }
      
      prev_seed = seed;
    }
    // last seed
    seed = linked_list_get_last(cal->seed_list);
    if (cal->read->length - seed->read_end > max_gap_length) {
      invalid = 1;
      cal->invalid = 1;
      continue;      
    }
  }

  return invalid;
}


void clean_cals(array_list_t **list, fastq_read_t *read, sa_index3_t *sa_index) {

  seed_t *prev_seed, *seed;
  linked_list_item_t *prev_item, *item;
  linked_list_iterator_t* itr;

  int trim, invalid, overlap;
  seed_cal_t *cal;
  array_list_t *cal_list = *list;
  size_t num_cals = array_list_size(cal_list);

  invalid = 0;
  for (int i = 0; i < num_cals; i++) {
    cal = array_list_get(i, cal_list);
    
    trim = 0;
    while (1) {
      item = cal->seed_list->first;
      if (item) {
	seed = item->item;
	if (is_invalid_cigar(&seed->cigar)) {
	  linked_list_remove_item(item, cal->seed_list);
	  seed_free(seed);
	  if (cal->seed_list->size <= 0) {
	    invalid = 1;
	    cal->invalid = 1;
	  }
	} else {
	  break;
	}
      } else {
	break;
      }
    }
    if (!item) continue;
    
    prev_item = item;
    while ((item = prev_item->next) != NULL) {
      prev_seed = prev_item->item;
      seed = item->item;

      if (is_invalid_cigar(&seed->cigar)) {
	linked_list_remove_item(item, cal->seed_list);
	seed_free(seed);
	if (cal->seed_list->size <= 0) {
	  invalid = 1;
	  cal->invalid = 1;
	}
      } else {
	if ( (prev_seed->read_end >= seed->read_start &&
	      prev_seed->genome_end >= seed->genome_start) ) {
	  
	  if (prev_seed->genome_end >= seed->genome_start) {
	    overlap = prev_seed->genome_end - seed->genome_start + 1;
	    
	    trim = 1;
	    seed_rtrim_read(overlap, prev_seed);
	    seed_ltrim_read(overlap, seed);
	    prev_item = item;
	  } else {
	    linked_list_remove_item(item, cal->seed_list);
	    seed_free(seed);
	    if (cal->seed_list->size <= 0) {
	      invalid = 1;
	      cal->invalid = 1;
	    }
	  }
	} else {
	  prev_item = item;
	}
      }
    }

    if (trim) {
      for (item = cal->seed_list->first; item != NULL; ) {
	seed = item->item;
	prev_item = item->next;
	if (is_invalid_cigar(&seed->cigar)) {
	  linked_list_remove_item(item, cal->seed_list);
	  seed_free(seed);
	  if (cal->seed_list->size <= 0) {
	    invalid = 1;
	    cal->invalid = 1;
	  }
	}
	item = prev_item;
      }
    }
  }
  
  if (invalid || check_gap_lengths(MAX_GAP_LENGTH, cal_list)) {
    array_list_t *new_cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    num_cals = array_list_size(cal_list);
    for (int i = 0; i < num_cals; i++) {
      cal = array_list_get(i, cal_list);
      if (!cal->invalid) {
	array_list_insert(cal, new_cal_list);
	array_list_set(i, NULL, cal_list);
      }
    }
    array_list_free(cal_list, (void *) seed_cal_free);
    *list = new_cal_list;
  }
}

//--------------------------------------------------------------------
// step three:
//    fill gaps using sw
//--------------------------------------------------------------------

void update_left_side_seed(int min_flank, seed_t *seed) {
  int q_flank = 0, r_flank = 0;
  // look at left-side cigar
  cigar_t *cigar = &seed->cigar;
  int found, op_value, op_name, trim = 0, num_ops = cigar->num_ops;

  found = 0;
  for (int i = 0; i < num_ops; i++) {
    cigar_get_op(i, &op_value, &op_name, cigar);
    if (op_name == '=' && op_value > min_flank) {
      found = 1;
      break;
    } else {
      trim++;
      if (op_name == '=' || op_name == 'X') {
	q_flank += op_value;
	r_flank += op_value;
      } else if (op_name == 'I') {
	q_flank += op_value;
      } else if (op_name == 'D') {
	r_flank += op_value;
      }
    }
  }  
  if (trim) {
    seed->read_start += q_flank;
    seed->genome_start += r_flank;
    cigar_ltrim_ops(trim, cigar);
  }
}

void update_right_side_seed(int min_flank, seed_t *seed) {
  int q_flank = 0, r_flank = 0;
  // look at right-side cigar
  cigar_t *cigar = &seed->cigar;
  int op_value, op_name, trim = 0, num_ops = cigar->num_ops;

  for (int i = num_ops - 1; i >= 0; i--) {
    cigar_get_op(i, &op_value, &op_name, cigar);
    if (op_name == '=' && op_value > min_flank) {
      break;
    } else {
      trim++;
      if (op_name == '=' || op_name == 'X') {
	q_flank += op_value;
	r_flank += op_value;
      } else if (op_name == 'I') {
	q_flank += op_value;
      } else if (op_name == 'D') {
	r_flank += op_value;
      }
    }
  }  
  if (trim) {
    seed->read_end -= q_flank;
    seed->genome_end -= r_flank;
    cigar_rtrim_ops(trim, cigar);
  }
}

int prepare_sw(fastq_read_t *read,   array_list_t *sw_prepare_list,
	       sa_mapping_batch_t *mapping_batch, sa_index3_t *sa_index, 
	       array_list_t *cal_list) {
  size_t seed_count, num_seeds, num_cals = array_list_size(cal_list);

  #ifdef _VERBOSE	  
  printf("\n\n====>>>> STEP THREE : prepare_sw <<<<====\n");
  #endif

  #ifdef _TIMING
  struct timeval stop, start;
  #endif

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif

  char *seq, *ref;
  seed_t *prev_seed, *seed;
  int gap_len;
  int gap_read_len, gap_genome_len;
  size_t start, end;
  size_t gap_read_start, gap_read_end;
  size_t gap_genome_start, gap_genome_end;

  linked_list_item_t *item;
  sw_prepare_t *sw_prepare;

  seed_cal_t *cal;
  cigar_t *cigar;
  cigarset_t *cigarset;

  int query_flank, ref_flank;
  int num_sw, num_total_sw = 0;

  for (int i = 0; i < num_cals; i++) {
    num_sw = 0;
    cal = array_list_get(i, cal_list);
    cal->num_mismatches = 0;

    #ifdef _VERBOSE
    seed_cal_print(cal);
    #endif

    // if not seeds, then next cal
    num_seeds = cal->seed_list->size;
    if (num_seeds <= 0) continue;

    seed_cal_merge_seeds(cal);
    if (cigar_get_length(&cal->cigar) == read->length) {
      continue;
    }
    
    // first seed
    seed = linked_list_get_first(cal->seed_list);

    // cal cigar
    num_seeds = cal->seed_list->size;
    cigarset = cigarset_new(num_seeds * 2 + 1);
    cal->cigarset = cigarset;

    if (seed->read_start > 0) {
      #ifdef _VERBOSE
      print_seed("-----> before updating first left-side seed: ", seed);
      #endif
      update_left_side_seed(SW_RIGHT_FLANK, seed);
      #ifdef _VERBOSE
      print_seed("-------> after updating first left-side seed: ", seed);
      #endif

      gap_genome_start = seed->genome_start - seed->read_start - 1 - SW_LEFT_FLANK_EX;
      gap_genome_end = seed->genome_start + SW_RIGHT_FLANK;
      ref = sa_genome_get_sequence(cal->chromosome_id, gap_genome_start, gap_genome_end, sa_index->genome);
      
      seq = get_subsequence((cal->strand ? read->revcomp : read->sequence), 
			    0, seed->read_start + SW_RIGHT_FLANK);
      
      sw_prepare = sw_prepare_new(seq, ref, 0, SW_RIGHT_FLANK, FIRST_SW);
      sw_prepare->seed_region = (seed_region_t *)seed;
      sw_prepare->cal = (cal_t *)cal;
      sw_prepare->read = read;
      array_list_insert(sw_prepare, sw_prepare_list);
      num_sw++;

      cigarset_info_set(CIGAR_FROM_GAP, 0, NULL, NULL, &cigarset->info[0]);
    } else {
      cigarset_info_set(0, 0, NULL, NULL, &cigarset->info[0]);
    }
    // seeds at the middle positions
    cal->num_mismatches += seed->num_mismatches;
    cal->num_open_gaps += seed->num_open_gaps;
    cal->num_extend_gaps += seed->num_extend_gaps;
    prev_seed = seed;
    cigarset_info_set(CIGAR_FROM_SEED, 0, &seed->cigar, seed, &cigarset->info[1]);
    seed_count = 0;
    for (item = cal->seed_list->first->next; item != NULL; item = item->next) {
      seed_count++;
      seed = item->item;
      cal->num_mismatches += seed->num_mismatches;
      cal->num_open_gaps += seed->num_open_gaps;
      cal->num_extend_gaps += seed->num_extend_gaps;

      #ifdef _VERBOSE
      print_seed("-----> before updating middle right-side seed: ", prev_seed);
      #endif
      update_right_side_seed(SW_LEFT_FLANK, prev_seed);
      #ifdef _VERBOSE
      print_seed("-------> after updating middle right-side seed: ", prev_seed);
      #endif

      if (prev_seed->read_end <= prev_seed->read_start) {
	cal->invalid = 1;
	break;
      }

      #ifdef _VERBOSE
      print_seed("-----> before updating middle left-side seed: ", seed);
      #endif
      update_left_side_seed(SW_RIGHT_FLANK, seed);
      #ifdef _VERBOSE
      print_seed("-------> after updating middle left-side seed: ", seed);
      #endif

      if (seed->read_end <= seed->read_start) {
	cal->invalid = 1;
	break;
      }

      gap_read_start = prev_seed->read_end + 1;
      gap_read_end = seed->read_start - 1;
      gap_read_len = gap_read_end - gap_read_start + 1;

      gap_genome_start = prev_seed->genome_end + 1;
      gap_genome_end = seed->genome_start - 1;
      gap_genome_len = gap_genome_end - gap_genome_start + 1;

      if (gap_read_len < 0) {
	// deletion
	start = gap_genome_start - SW_LEFT_FLANK - abs(gap_read_len);
	end = gap_genome_end + SW_RIGHT_FLANK + abs(gap_read_len);
	if (end <= start) {
	  cal->invalid = 1;
	  break;
	}

	seq = get_subsequence((cal->strand ? read->revcomp : read->sequence), 
			      gap_read_start - SW_LEFT_FLANK - abs(gap_read_len), 
			      gap_read_len + SW_LEFT_FLANK + SW_RIGHT_FLANK + (2*abs(gap_read_len)));

	ref = sa_genome_get_sequence(cal->chromosome_id, start, end, sa_index->genome);

	cigarset_info_set(CIGAR_FROM_GAP, abs(gap_read_len), NULL, NULL, &cigarset->info[seed_count * 2]);
      } else if (gap_genome_len < 0) {
	// insertion
	start = gap_genome_start - SW_LEFT_FLANK - abs(gap_genome_len);
	end = gap_genome_end + SW_RIGHT_FLANK + abs(gap_genome_len);
	if (end <= start) {
	  cal->invalid = 1;
	  break;
	}

	seq = get_subsequence((cal->strand ? read->revcomp : read->sequence), 
			      gap_read_start - SW_LEFT_FLANK, gap_read_len + SW_LEFT_FLANK + SW_RIGHT_FLANK);

	ref = sa_genome_get_sequence(cal->chromosome_id, start, end, sa_index->genome);
            
	cigarset_info_set(CIGAR_FROM_GAP, 0, NULL, NULL, &cigarset->info[seed_count * 2]);
      } else {

        #ifdef _VERBOSE1
	print_seed("", prev_seed);
	print_seed("", seed);
	printf("read id = %s\n", read->id);
	printf("genome gap (start, end) = (%lu, %lu), len = %i\n", gap_genome_start, gap_genome_end, gap_genome_len);
	printf("read gap (start, end) = (%lu, %lu), len = %i\n", gap_read_start, gap_read_end, gap_read_len);
	exit(-1);
        #endif

	seq = get_subsequence((cal->strand ? read->revcomp : read->sequence), 
			      gap_read_start - SW_LEFT_FLANK, gap_read_len + SW_LEFT_FLANK + SW_RIGHT_FLANK);

	ref = sa_genome_get_sequence(cal->chromosome_id, gap_genome_start - SW_LEFT_FLANK, 
				     gap_genome_end + SW_RIGHT_FLANK, sa_index->genome);

	cigarset_info_set(CIGAR_FROM_GAP, 0, NULL, NULL, &cigarset->info[seed_count * 2]);
      }
      // prepare MIDDLE_SW
      sw_prepare = sw_prepare_new(seq, ref, 0, 0, MIDDLE_SW);

      sw_prepare->seed_region = (seed_region_t *)(seed_count * 2);
      sw_prepare->cal = (cal_t *)cal;
      sw_prepare->read = read;
      array_list_insert(sw_prepare, sw_prepare_list);
      num_sw++;

      prev_seed = seed;
      cigarset_info_set(CIGAR_FROM_SEED, 0, &seed->cigar, seed, &cigarset->info[seed_count * 2 + 1]);
    }

    if (cal->invalid) {
      cigarset_free(cigarset);
      cal->cigarset = NULL;
      continue;
    }

    // last seed
    seed = linked_list_get_last(cal->seed_list);
    if (seed->read_end < read->length - 1) {
      #ifdef _VERBOSE
      print_seed("-----> before updating last right-side seed: ", seed);
      #endif
      update_right_side_seed(SW_LEFT_FLANK, seed);
      #ifdef _VERBOSE
      print_seed("-------> after updating last right-side seed: ", seed);
      #endif
	
      gap_genome_start = seed->genome_end - SW_LEFT_FLANK + 1;
      gap_genome_end = gap_genome_start + (read->length - seed->read_end) + SW_LEFT_FLANK_EX;

      if (gap_genome_end <= gap_genome_start) {
	printf("in last seed gap_genome_end <= gap_genome_start\n");
	exit(-1);
      }

      ref = sa_genome_get_sequence(cal->chromosome_id, gap_genome_start, gap_genome_end, sa_index->genome);
      
      seq = get_subsequence((cal->strand ? read->revcomp : read->sequence), 
			    seed->read_end - SW_LEFT_FLANK + 1, 
			    read->length + SW_LEFT_FLANK - seed->read_end);
      
      sw_prepare = sw_prepare_new(seq, ref, 0, 0, LAST_SW);
      sw_prepare->seed_region = (seed_region_t *)(num_seeds * 2);
      sw_prepare->cal = (cal_t *)cal;
      sw_prepare->read = read;
      array_list_insert(sw_prepare, sw_prepare_list);
      num_sw++;
      
      cigarset_info_set(CIGAR_FROM_GAP, 0, NULL, NULL, &cigarset->info[num_seeds * 2]);
    } else {
      cigarset_info_set(0, 0, NULL, NULL, &cigarset->info[num_seeds * 2]);
    }

    num_total_sw += num_sw;
  }
  
  return num_total_sw;
}

//--------------------------------------------------------------------

void execute_sw(array_list_t *sw_prepare_list, sa_mapping_batch_t *mapping_batch) {

  #ifdef _TIMING
  struct timeval stop, start;
  #endif

  sw_prepare_t *sw_prepare;

  seed_cal_t *cal;
  cigar_t *cigar;
  cigarset_t *cigarset;

  // apply smith-waterman
  sw_optarg_t sw_optarg;
  sw_optarg_init(10, 0.5, 5, -4, &sw_optarg);
  
  size_t sw_count = array_list_size(sw_prepare_list); 
  char *q[sw_count], *r[sw_count];
  for (int i = 0; i < sw_count; i++) {
    sw_prepare = array_list_get(i, sw_prepare_list);
    q[i] = sw_prepare->query;
    r[i] = sw_prepare->ref;

    #ifdef _VERBOSE
    printf("\t\t%s:%i:to SW:\n", __FILE__, __LINE__);
    printf("\t\t%i: query: %s\n", i, q[i]);
    printf("\t\t%i: ref. : %s\n", i, r[i]);
    printf("\n");
    #endif
  }
  
  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_PRE_SW] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif

  sw_multi_output_t *sw_output = sw_multi_output_new(sw_count);
  smith_waterman_mqmr(q, r, sw_count, &sw_optarg, 1, sw_output);
  #ifdef _VERBOSE
  printf("\t\t%s:%i:SW output:\n", __FILE__, __LINE__);
  sw_multi_output_save(sw_count, sw_output, stdout);
  #endif

  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_SW] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif

  // process Smith-Waterman output
  seed_t *seed;
  char *seq, *ref, *query_map, *ref_map;
  int op_name, op_value, diff, len, r_nt_mapped;
  int left_flank, right_flank, query_start, ref_start;
  for (int i = 0; i < sw_count; i++) {
    sw_prepare = array_list_get(i, sw_prepare_list);

    cal = (seed_cal_t *) sw_prepare->cal;

    if (cal->invalid) {
      // free memory
      free(sw_prepare->query);
      free(sw_prepare->ref);
      sw_prepare_free(sw_prepare);
      continue;
    }
    
    cigarset = cal->cigarset;
    cigar = cigar_new_empty();

    query_map = sw_output->query_map_p[i];
    ref_map = sw_output->ref_map_p[i];

    // nt mapped in reference
    r_nt_mapped = 0;
    len = strlen(ref_map);

    left_flank = sw_prepare->left_flank;
    right_flank = sw_prepare->right_flank;
    query_start = sw_output->query_start_p[i];
    ref_start = sw_output->ref_start_p[i];
    diff = query_start - ref_start;

    // check initial positions
    if (sw_prepare->ref_type == FIRST_SW) {
      seed = (seed_t *)sw_prepare->seed_region;
      if (query_start > 0) {
	cigar_append_op(query_start, 'S', cigar);
      }
      for(int j = 0; j < len; j++) {
	if (ref_map[j] != '-') {
	  r_nt_mapped++;
	}
      }
      cal->start = seed->genome_start - (r_nt_mapped - right_flank) + query_start;
    } else {
      if (query_start > 0) {
	if (ref_start > 0) {
	  if (diff == 0) {
	    cigar_append_op(query_start, '=', cigar);      
	  } else if (diff > 0) {
	    cigar_append_op(ref_start, '=', cigar);      
	    cigar_append_op(diff, 'I', cigar);      
	  } else {
	    cigar_append_op(query_start, '=', cigar);      
	    cigar_append_op(abs(diff), 'D', cigar);      
	  }
	} else {
	  cigar_append_op(query_start, 'I', cigar);      
	}
      } else if (ref_start > 0) {
	//cigar_append_op(ref_start, '=', cigar);      
      }
    }

    // scan map to complete cigar
    op_value = 0;
    op_name = '=';
    for(int i = 0; i < len; i++) {
      if (query_map[i] == '-') {
	// deletion (in the query)
	if (op_name != 'D' && op_value > 0) {
	  cigar_append_op(op_value, op_name, cigar);
	  op_value = 0;
	}
	op_value++;
	op_name = 'D';
      } else if (ref_map[i] == '-') {
	// insertion (in the query)
	if (op_name != 'I' && op_value > 0) {
	  cigar_append_op(op_value, op_name, cigar);
	  op_value = 0;
	}
	op_value++;
	op_name = 'I';
      } else if (ref_map[i] == query_map[i]) {
	if (op_name != '=' && op_value > 0) {
	  cigar_append_op(op_value, op_name, cigar);
	  op_value = 0;
	}
	op_value++;
	op_name = '=';
      } else {
	if (op_name != 'X' && op_value > 0) {
	  cigar_append_op(op_value, op_name, cigar);
	  op_value = 0;
	}
	op_value++;
	op_name = 'X';
      }
    }
    cigar_append_op(op_value, op_name, cigar);


    size_t gap_count = 0;
    if (sw_prepare->ref_type == FIRST_SW) {
      gap_count = 0;
    } else {
      gap_count = (size_t)sw_prepare->seed_region;
    }
    cigarset->info[gap_count].active = CIGAR_FROM_GAP;
    cigarset->info[gap_count].cigar = cigar;

    #ifdef _VERBOSE
    printf("************** sw_count: %i, for gap %i cigar %s\n", 
	   i, gap_count, cigar_to_string(cigar));
    #endif

    // free memory
    free(sw_prepare->query);
    free(sw_prepare->ref);
    sw_prepare_free(sw_prepare);
  }

  // free memory
  sw_multi_output_free(sw_output);

  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_POST_SW] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif
}

//--------------------------------------------------------------------

void post_process_sw(int sw_post_read_counter, int *sw_post_read,   
		     array_list_t **cal_lists, sa_mapping_batch_t *mapping_batch) {

  #ifdef _TIMING
  struct timeval stop, start;
  #endif

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif

  int num_cals, cigar_type;
  array_list_t *cal_list;

  seed_t *seed;
  seed_cal_t *cal;
  cigar_t *cigar, *aux_cigar;
  cigarset_t *cigarset;
  int cigar_len, op_value, op_name;

  for (int j = 0; j < sw_post_read_counter; j++) {

    cal_list = cal_lists[sw_post_read[j]];
    num_cals = array_list_size(cal_list);
    
    // re-construct CIGAR for each CAL (from seed and gap CIGARS)
    for (int i = 0; i < num_cals; i++) {
      cal = array_list_get(i, cal_list);

      if (cal->seed_list->size <= 0 || cal->invalid) continue;
      
      #ifdef _VERBOSE
      printf("\t\t\t\tcal score = %0.2f\n", cal->score);
      seed_cal_print(cal);
      #endif
      cigarset = cal->cigarset;
      cigar = &cal->cigar;
      //printf("cigarset size = %i\n", cigarset->size);
      for (int j = 0; j < cigarset->size; j++) {

	cigar_type = cigarset->info[j].active;
	if (cigar_type > 0) {
	  //printf("************** cigar type %i : %s\n",
	  //	 cigar_type, cigar_to_string(cigarset->info[j].cigar));

	  if (cigar_type == CIGAR_FROM_SEED) {
	    seed = cigarset->info[j].seed;
	    aux_cigar = cigarset->info[j].cigar;
	    // CIGAR_FROM_SEED
	    if (seed->read_start > 0) {
	      cigar_get_op(0, &op_value, &op_name, aux_cigar);
	      cigar_set_op(0, op_value - SW_LEFT_FLANK, op_name, aux_cigar);
	    }
	    if (seed->read_end < cal->read->length - 1) {
	      cigar_get_op(aux_cigar->num_ops - 1, &op_value, &op_name, aux_cigar);
	      cigar_set_op(aux_cigar->num_ops - 1, op_value - SW_RIGHT_FLANK, op_name, aux_cigar);
	    }
            #ifdef _VERBOSE
	    printf("************** gap %i -> concat SEED cigar %s into %s\n",
		   j, cigar_to_string(aux_cigar), cigar_to_string(cigar));
            #endif
	    cigar_concat(aux_cigar, cigar);
	  } else {
	    // CIGAR_FROM_GAP
            #ifdef _VERBOSE
	    printf("************** gap %i -> concat GAP cigar %s into %s\n",
		   j, cigar_to_string(cigarset->info[j].cigar), cigar_to_string(cigar));
            #endif
	    aux_cigar = cigarset->info[j].cigar;
	    if (cigarset->info[j].overlap) {
	      cigar_get_op(aux_cigar->num_ops - 1, &op_value, &op_name, aux_cigar);
	      cigar_set_op(aux_cigar->num_ops - 1, op_value - cigarset->info[j].overlap, op_name, aux_cigar);
	    }
	    cigar_concat(aux_cigar, cigar);
	    cigar_free(aux_cigar);
	  }
	}
      }
      cigar_len = cigar_get_length(cigar);
      if (cigar_len < cal->read->length) {
	if (seed && seed->read_end == cal->read->length - 1) {
	  cigar_append_op(cal->read->length - cigar_len, '=', cigar);
	} else {
	  cigar_append_op(cal->read->length - cigar_len, 'S', cigar);
	}
      }
      //printf("************** CAL cigar: %s\n", cigar_to_string(cigar));
    }
  }

  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_POST_SW] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif
}

//--------------------------------------------------------------------
// sa mapper
//--------------------------------------------------------------------

int sa_single_mapper(void *data) {

  #ifdef _TIMING
  struct timeval stop, start;
  #endif

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif
  
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;

  int num_seeds = wf_batch->options->num_seeds;
  int min_cal_size = wf_batch->options->min_cal_size;
  
  sa_mapping_batch_t *mapping_batch = wf_batch->mapping_batch;
  mapping_batch->options = wf_batch->options;

  sa_index3_t *sa_index = (sa_index3_t *) wf_batch->sa_index;
  
  int bam_format = mapping_batch->bam_format;

  size_t num_reads = mapping_batch->num_reads;
  int max_read_area, min_num_mismatches;
  float max_score;

  // smith-waterman parameters
  float match_score = wf_batch->options->match;
  float mismatch_penalty = wf_batch->options->mismatch;
  float gap_open_penalty = -1.0f * wf_batch->options->gap_open;
  float gap_extend_penalty = -1.0f * wf_batch->options->gap_extend;
  
  // CAL management
  size_t num_cals;
  seed_cal_t *cal;
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
  
  return -1;
}

//--------------------------------------------------------------------

int sa_pair_mapper(void *data) {

  #ifdef _TIMING
  struct timeval stop, start;
  #endif

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif
  
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;

  int pair_mode = wf_batch->options->pair_mode;
  int pair_min_distance = wf_batch->options->pair_min_distance;
  int pair_max_distance = wf_batch->options->pair_max_distance;

  int num_seeds = wf_batch->options->num_seeds;
  int min_cal_size = wf_batch->options->min_cal_size;
  
  sa_mapping_batch_t *mapping_batch = wf_batch->mapping_batch;
  mapping_batch->options = wf_batch->options;

  sa_index3_t *sa_index = (sa_index3_t *) wf_batch->sa_index;
  
  int bam_format = mapping_batch->bam_format;

  size_t num_reads = mapping_batch->num_reads;
  int max_read_area, min_num_mismatches;
  float max_score;

  // smith-waterman parameters
  float match_score = wf_batch->options->match;
  float mismatch_penalty = wf_batch->options->mismatch;
  float gap_open_penalty = -1.0f * wf_batch->options->gap_open;
  float gap_extend_penalty = -1.0f * wf_batch->options->gap_extend;

  // CAL management
  size_t num_cals;
  seed_cal_t *cal;
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

    //printf("\t\t%s:%i: after create_cals:\n", __FILE__, __LINE__);
    //for(int i = 0; i < array_list_size(cal_list); i++) { seed_cal_print(array_list_get(i, cal_list)); }

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
  filter_cals_by_pair_mode(pair_mode, pair_min_distance, pair_max_distance, 
    			   num_reads, cal_lists);
  
  check_pairs(cal_lists, sa_index, mapping_batch, cal_mng);

  // 3) prepare Smith-Waterman to fill in the gaps
  for (int i = 0; i < num_reads; i++) {

    read = array_list_get(i, mapping_batch->fq_reads);
    cal_list = cal_lists[i];

    //printf("---> %s:%s:%i: before prepare_sw:\n", __FILE__, __func__, __LINE__);
    //for(int i = 0; i < array_list_size(cal_list); i++) { seed_cal_print(array_list_get(i, cal_list)); }

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
  
  return -1;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
