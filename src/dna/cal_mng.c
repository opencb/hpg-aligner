#include "cal_mng.h"

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
  
  int i, overlap;

  seed_t *item;
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
  for (unsigned short int i = 0; i < num_chroms; i++) {
    cals_lists[i] = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
  }

  cal_mng_t *p = (cal_mng_t *) calloc(1, sizeof(cal_mng_t));
  p->read_length = 10;
  p->min_read_area = 100;
  p->max_read_area = 0;
  p->num_chroms = num_chroms;
  p->cals_lists = cals_lists;

  memset(p->active_mask, 0, sizeof(p->active_mask));
  p->num_active = 0;

  p->suffix_mng = suffix_mng_new(genome);

  return p;
}

//--------------------------------------------------------------------

void cal_mng_free(cal_mng_t *p) {
  if (p) {
    if (p->cals_lists) {
      for (unsigned short int i = 0; i < p->num_chroms; i++) {
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
      for (unsigned short int i = 0; i < p->num_chroms; i++) {
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

  unsigned short int chrom;
  linked_list_item_t *item;
  linked_list_t *list;
  cal_t *cal;
  if (p) {
    if (p->cals_lists) {     
      for (int i = 0; i < p->num_active; i++) {
	chrom = p->active[i];
	list = p->cals_lists[chrom];
	if (list->size > 0) {
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
      p->num_active = 0;
      memset(p->active_mask, 0, sizeof(p->active_mask));
    }
  }
}

//--------------------------------------------------------------------

void cal_mng_clear(cal_mng_t *p) {
  if (p) {
    if (p->cals_lists) {
      linked_list_t *list;
      unsigned short int chrom;
      for (int i = 0; i < p->num_active; i++) {
	chrom = p->active[i];
	list = p->cals_lists[chrom];
	if (list->size > 0) {
	  linked_list_clear(list, (void *)seed_cal_free);
	}
      }
      p->num_active = 0;
      memset(p->active_mask, 0, sizeof(p->active_mask));
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

      if (!p->active_mask[seed->chromosome_id]) {
	p->active[p->num_active++] = seed->chromosome_id;	
	p->active_mask[seed->chromosome_id] = 1;
      }
      
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

int cal_mng_find(int strand, unsigned short int chrom, size_t start, size_t end, cal_mng_t *p) {
  #ifdef _VERBOSE1
  printf("\t\t***** searching CAL: chrom %u: %lu-%lu\n", chrom, start, end);
  #endif
  int found_cal = 0;
  if (p->cals_lists) {
    linked_list_t *cal_list = p->cals_lists[chrom];
    if (cal_list->size > 0) {
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


  #ifdef _VERBOSE
  printf("-----> cal_mng_to_array_list\n");
  #endif

  if (p->cals_lists) {
    linked_list_t *cal_list;
    unsigned short int chrom;
    for (int i = 0; i < p->num_active; i++) {
      chrom = p->active[i];
      cal_list = p->cals_lists[chrom];
      if (cal_list->size > 0) {
	while ((cal = (seed_cal_t *) linked_list_remove_last(cal_list))) {
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
}

//--------------------------------------------------------------------

void cal_mng_select_best(int read_area, array_list_t *valid_list, array_list_t *invalid_list, 
			 cal_mng_t *p) {
  seed_cal_t *cal;
  
  if (p->cals_lists) {
    linked_list_t *cal_list;
    unsigned short int chrom;
    for (int i = 0; i < p->num_active; i++) {
      chrom = p->active[i];
      cal_list = p->cals_lists[chrom];
      if (cal_list->size > 0) {
	while ((cal = (seed_cal_t *) linked_list_remove_last(cal_list))) {
	  if (p->min_read_area <= read_area && cal->read_area <= read_area) {
	    array_list_insert(cal, valid_list);
	  } else {
	    array_list_insert(cal, invalid_list);
	  }
	}
      }
    }
  }
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
