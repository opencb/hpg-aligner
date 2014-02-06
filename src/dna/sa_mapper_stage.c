#include "sa_mapper_stage.h"

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
  return p;
}

//--------------------------------------------------------------------

void cal_mng_free(cal_mng_t *p) {
  if (p) {
    if (p->cals_lists) {
      for (unsigned int i = 0; i < p->num_chroms; i++) {
	if (p->cals_lists[i]) {
	  linked_list_free(p->cals_lists[i], seed_cal_free);
	}
      }
      free(p->cals_lists);
    }
    free(p);
  }
}

//--------------------------------------------------------------------

void cal_mng_clear(cal_mng_t *p) {
  if (p) {
    if (p->cals_lists) {
      for (unsigned int i = 0; i < p->num_chroms; i++) {
	if (p->cals_lists[i]) {
	  linked_list_clear(p->cals_lists[i], seed_cal_free);
	}
      }
    }
  }
}

//--------------------------------------------------------------------

void cal_mng_update(seed_t *seed, fastq_read_t *read, cal_mng_t *p) {
  int is_used = 0;
  if (p->cals_lists) {
    seed_cal_t *cal;
    seed_t *s_last;
    linked_list_t *seed_list;
    linked_list_t *cal_list = p->cals_lists[seed->chromosome_id];
    if (cal_list) {
      #ifdef _VERBOSE
      printf("\t\t\tinsert this seed to the CAL manager:\n");
      print_seed("\t\t\t", seed);
      #endif

      //      if (cal->read_area < p->min_read_area) {
	//	printf("****************** set read min. area from %i to %i\n", p->min_read_area, cal->read_area);
      //	p->min_read_area = cal->read_area;
      //      }
      if (linked_list_size(cal_list) <= 0) {
	// create CAL and insert it into the CAL manager
	is_used = 1;
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
	seed_cal_t *item = (cal_t *) linked_list_iterator_curr(itr);
	while (item != NULL) {
	  #ifdef _VERBOSE
	  printf("---> merging with this CAL?\n");
	  seed_cal_print(item);
	  #endif
	  //	  assert(cal->end > item->start);
	  s_last = linked_list_get_last(item->seed_list);
	  if (seed->read_start - s_last->read_end < read->length && 
	      seed->genome_start - s_last->genome_end < read->length) {
	    is_used = 1;
	    linked_list_insert_last(seed, item->seed_list);
	    item->end = seed->genome_end;
	    item->read_area += (seed->read_end - seed->read_start + 1);
	    item->num_mismatches += seed->num_mismatches + seed->num_open_gaps + seed->num_extend_gaps;

	    //	    cal_free(cal);
            #ifdef _VERBOSE
	    printf("---> yes, merging CAL, result:\n");
	    seed_cal_print(item);
	    #endif
	    break;
	  } else {
	    if (seed->genome_end < item->start) {
	      // create CAL and insert it into the CAL manager
	      is_used = 1;
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
	  }
	  //continue loop...
	  linked_list_iterator_next(itr);
	  item = linked_list_iterator_curr(itr);
	}
	if (item == NULL) {
	  // create CAL and insert it into the CAL manager
	  is_used = 1;
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
  if (!is_used) {
    //    seed_region_free(seed);
  }
}

//--------------------------------------------------------------------

int cal_mng_find(int chrom, size_t start, size_t end, cal_mng_t *p) {
  #ifdef _VERBOSE1
  printf("\t\t***** searching CAL: chrom %i: %lu-%lu\n", chrom, start, end);
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
	printf("\t\t\t***** searching CAL: comparing %lu-%lu to cal item %lu-%lu\n", 
	       start, end, cal->start, cal->end);
	#endif
	if (cal->start <= start && cal->end >= end) {
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

void cal_mng_to_array_list(int read_area, array_list_t *out_list, cal_mng_t *p) {
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
	if ((cal->end - cal->start) >= read_area) {
	  array_list_insert(cal, out_list);
	} else {
	  // free CAL
	  seed_cal_free(cal);
	}
	/*
	if (p->min_read_area <= read_area && cal->read_area <= read_area) {
	  array_list_insert(cal, out_list);
	} else {
	  cal_free(cal);
	}
	*/
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
// generate cal from an exact read
//--------------------------------------------------------------------

void generate_cals_from_exact_read(int strand, fastq_read_t *read,
				   size_t low, size_t high, sa_index3_t *sa_index, 
				   array_list_t *cal_list) {  
  size_t g_start, g_end;
  int chrom;

  seed_cal_t *cal;
  cigar_t *cigar;
  seed_t *seed;
  linked_list_t *seed_list;
  
  for (size_t suff = low; suff <= high; suff++) {
    chrom = sa_index->CHROM[suff];
    g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];
    g_end = g_start + read->length - 1;
    seed_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    seed = seed_new(0, read->length - 1, g_start, g_end);
    cigar_append_op(read->length, 'M', &seed->cigar);
    linked_list_insert(seed, seed_list);
    cal = seed_cal_new(chrom, strand, g_start, g_end, seed_list);
    cal->read_area = read->length;
    cal->num_mismatches = 0;
    cal->read = read;
    array_list_insert(cal, cal_list);
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
  int found_cal, chrom, diff, max_map_len = 0;

  float score;
  alig_out_t alig_out;

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
    chrom = sa_index->CHROM[suff];

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
    found_cal = cal_mng_find(chrom, g_start_suf, g_end_suf, cal_mng);
    #ifdef _TIMING
    gettimeofday(&stop, NULL);
    mapping_batch->func_times[FUNC_SKIP_SUFFIXES] += 
      ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
    #endif
    if (found_cal) continue;


    #ifdef _VERBOSE
    printf("\t\tsuffix at [%lu|%lu-%lu|%lu] %c chrom %s\n",
	   g_start_suf, r_start_suf, r_end_suf, g_end_suf, (strand == 0 ? '+' : '-'), 
	   sa_index->genome->chrom_names[chrom]);
    #endif

    seed = seed_new(r_start_suf, r_end_suf, g_start_suf, g_end_suf);

    // extend suffix to left side, if necessary
    if (r_start_suf > 0) {
      r_start = 0;
      r_end = r_start_suf - 1;
      r_len = r_start_suf;

      g_len = r_len + 5;
      g_end = g_start_suf - 1;
      g_start = g_end - g_len;
      //      printf("(g_start, g_end) = (%lu, %lu)\n", g_start, g_end);
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
			    
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_MINI_SW_LEFT_SIDE] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif
      if (1) {
      //      if (score > 0 || (alig_out.map_len1 + suffix_len) > 20) {
	// update seed
	seed->num_mismatches += alig_out.mismatch;
	seed->num_open_gaps += alig_out.gap_open;
	seed->num_extend_gaps += alig_out.gap_extend;

	seed->read_start -= alig_out.map_len1;
	seed->genome_start -= alig_out.map_len2;

	// if there's a mini-gap then try to fill the mini-gap
	if (seed->read_start > 0 && seed->read_start < 10) {
          #ifdef _TIMING
	  gettimeofday(&start, NULL);
          #endif
	  g_seq = &sa_index->genome->S[seed->genome_start + sa_index->genome->chrom_offsets[chrom] - seed->read_start];
	  for (size_t k1 = 0, k2 = 0; k1 < seed->read_start; k1++, k2++) {
	    if (r_seq[k1] != g_seq[k2]) {
	      seed->num_mismatches++;
	    }
	  }
	  cigar_append_op(seed->read_start, 'M', &seed->cigar);
	  
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
    cigar_append_op(suffix_len, 'M', &seed->cigar);

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
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_MINI_SW_RIGHT_SIDE] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif
      if (1) {
      //      if (score > 0 || (alig_out.map_len1 + suffix_len) > 20) {
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
	if (diff > 0 && diff < 10) {
          #ifdef _TIMING
	  gettimeofday(&start, NULL);
          #endif
	  g_seq = &sa_index->genome->S[seed->genome_end + sa_index->genome->chrom_offsets[chrom] + 1];
	  for (size_t k1 = seed->read_end + 1, k2 = 0; k1 < read->length; k1++, k2++) {
	    if (r_seq[k1] != g_seq[k2]) {
	      seed->num_mismatches++;
	    }
	  }
	  cigar_append_op(diff, 'M', &seed->cigar);

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

  //  return max_map_len + suffix_len + 1;
  //  return (sa_index->k_value * 2); // / 2);
  //  return (sa_index->k_value * 4);
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
  int suffix_len, num_suffixes;
  char *r_seq = read->sequence;

  size_t low, high;

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif

  array_list_t *cal_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  cal_mng->min_read_area = read->length;
  cal_mng->read_length = read->length;

  int max_read_area;// = read->length * MISMATCH_PERC;
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

  //    memset(saved_pos, 0, sizeof(saved_pos));

  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_OTHER] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif

  // first step, searching mappings in both strands
  // distance between seeds >= prefix value (sa_index->k_value)
  for (int strand = 0; strand < 2; strand++) {
    #ifdef _VERBOSE	  
    printf("=======> STRAND %c (read_end_pos = %i)\n", (strand == 0 ? '+' : '-'), read_end_pos);
    #endif

    for (read_pos = 0; read_pos < read_end_pos; )  {	
      // save this position and search suffixes from this read position
      //	saved_pos[strand][read_pos] = 1;
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
	//display_suffix_mappings(strand, read_pos, suffix_len, low, high, sa_index);
        #endif 
	
	// exact search
	if (suffix_len == read->length) {
          #ifdef _TIMING
	  gettimeofday(&start, NULL);
          #endif
	  generate_cals_from_exact_read(strand, read,
					low, high, sa_index, cal_list);
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
    
    if (extra_seed) {
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
	//display_suffix_mappings(strand, read_pos, suffix_len, low, high, sa_index);
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
    cal_mng_to_array_list((read->length / 3), cal_list, cal_mng);
    #ifdef _TIMING
    gettimeofday(&stop, NULL);
    mapping_batch->func_times[FUNC_CAL_MNG_TO_LIST] += 
      ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
    #endif
    
    // next, - strand
    r_seq = read->revcomp;
  } // end of for strand
  
  //  printf("**************** filter min_read_area: = %i, num_cals = %i\n", 
  //	 cal_mng->min_read_area, array_list_size(cal_list));
  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif
  max_read_area = get_max_read_area(cal_list);
  if (max_read_area > read->length) max_read_area = read->length;
  filter_cals_by_max_read_area(max_read_area, &cal_list);
  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_FILTER_BY_READ_AREA] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif

  #ifdef _VERBOSE
  printf("\t***> max_read_area = %i\n", max_read_area);
  #endif
  
  return cal_list;
}

//--------------------------------------------------------------------
// step two:
//    search prefix + extend using min-sw
//--------------------------------------------------------------------
/*
void step_two(fastq_read_t *read, char *revcomp_seq,
	      sa_mapping_batch_t *mapping_batch, 
	      sa_index3_t *sa_index, cal_mng_t *cal_mng,
	      array_list_t *valid_cal_list, array_list_t *invalid_cal_list) {
  #ifdef _TIMING
  struct timeval stop, start;
  #endif

  int suffix_len, num_suffixes;
  char *r_seq = read->sequence;

  size_t low, high;

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif

  cal_mng->min_read_area = 100;

  int max_read_area = read->length * MISMATCH_PERC;
  int read_pos, read_end_pos, read_inc = sa_index->k_value / 2;
 
  #ifdef _VERBOSE	  
  printf("\n\n====>>>> STEP TWO <<<<====\n");
  display_cmp_sequences(read, revcomp_seq, sa_index);
  #endif

  // fill in the CAL manager structure
  read_end_pos = read->length - sa_index->k_value;
  //    memset(saved_pos, 0, sizeof(saved_pos));

  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_OTHER] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif

  // first step, searching mappings in both strands
  // distnce between seeds >= prefix value (sa_index->k_value)
  for (int strand = 0; strand < 2; strand++) {
    #ifdef _VERBOSE	  
    printf("=======> STRAND %c\n", (strand == 0 ? '+' : '-'));
    #endif

    for (read_pos = 0; read_pos < read_end_pos; )  {	
      // save this position and search suffixes from this read position
      //	saved_pos[strand][read_pos] = 1;
      #ifdef _VERBOSE	  
      printf("\tread pos. = %lu\n", read_pos);
      #endif

      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      num_suffixes = search_prefix(&r_seq[read_pos], &low, &high, sa_index, 0);
      suffix_len = sa_index->k_value;
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_SEARCH_PREFIX] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif
      
      #ifdef _VERBOSE	  
      printf("\t\tnum. suffixes = %lu (suffix length = %lu)\n", num_suffixes, suffix_len);
      #endif
      if (num_suffixes && num_suffixes < MAX_NUM_SUFFIXES && suffix_len) {
        #ifdef _VERBOSE	  
	//display_suffix_mappings(strand, read_pos, suffix_len, low, high, sa_index);
        #endif 
	
	assert(suffix_len != read->length);

        #ifdef _TIMING
	gettimeofday(&start, NULL);
        #endif
	read_pos += generate_cals_from_suffixes(strand, read, revcomp_seq,
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
      } else {
	read_pos += read_inc;
      }
    } // end of for read_pos
    
      // update cal list from cal manager
    #ifdef _TIMING
    gettimeofday(&start, NULL);
    #endif
    cal_mng_select_best(max_read_area, valid_cal_list, invalid_cal_list, cal_mng);
    #ifdef _TIMING
    gettimeofday(&stop, NULL);
    mapping_batch->func_times[FUNC_CAL_MNG_TO_LIST] += 
      ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
    #endif
    
    // next, - strand
    r_seq = revcomp_seq;
  } // end of for strand

  //#ifdef _TIMING
  //gettimeofday(&start, NULL);
  //#endif
  //filter_cals_by_max_read_area(cal_mng->min_read_area, &cal_list);
  //#ifdef _TIMING
  //gettimeofday(&stop, NULL);
  //mapping_batch->func_times[FUNC_FILTER_BY_READ_AREA] += 
  //  ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  //#endif

}
*/

//--------------------------------------------------------------------
//
//--------------------------------------------------------------------

void clean_cals(array_list_t *cal_list, fastq_read_t *read, sa_index3_t *sa_index) {

  seed_t *prev_seed, *seed;
  linked_list_item_t *prev_item, *item;

  seed_cal_t *cal;
  size_t num_cals = array_list_size(cal_list);

  for (int i = 0; i < num_cals; i++) {
    cal = array_list_get(i, cal_list);
    prev_item = cal->seed_list->first;
    while ((item = prev_item->next) != NULL) {
      prev_seed = prev_item->item;
      seed = item->item;

      if (prev_seed->read_end >= seed->read_start &&
      	  prev_seed->genome_end >= seed->genome_start) {

      //      if (prev_seed->read_start == seed->read_start ||
      //	  prev_seed->read_end == seed->read_end) {

	printf("----> TO CLEAN (read %s):\n", read->id);
	display_cmp_sequences(read, sa_index);
	seed_cal_print(cal);

	linked_list_remove_item(item, cal->seed_list);
	seed_free(seed);

	//	printf("----> after remove:\n");
	//	seed_cal_print(cal);
      } else {
	prev_item = item;
	//	printf("----> keep it\n");
	//	display_cmp_sequences(read, sa_index);
	//	seed_cal_print(cal);
      }
    }
  }
}

//--------------------------------------------------------------------
// step three:
//    fill gaps using sw
//--------------------------------------------------------------------

void step_three0(fastq_read_t *read, sa_mapping_batch_t *mapping_batch, 
		sa_index3_t *sa_index, array_list_t *cal_list) {

  size_t seed_count, num_seeds, num_cals = array_list_size(cal_list);
  seed_cal_t *cal;
  seed_t *prev_seed, *seed;

  for (int i = 0; i < num_cals; i++) {
    cal = array_list_get(i, cal_list);
    cal->num_mismatches = 0;

    // if not seeds, then next cal
    num_seeds = cal->seed_list->size;
    if (num_seeds <= 0) continue;

    if (num_seeds == 1) mapping_batch->counters[0]++;
    if (num_seeds > 1) mapping_batch->counters[1]++;

    // first seed
    seed = linked_list_get_first(cal->seed_list);
    cal->num_open_gaps = seed->num_open_gaps;
    if (seed->read_start > 0) {
      cal->start -= seed->read_start;
      // SW here?
    }

    // middle seeds
    linked_list_item_t *item;
    prev_seed = seed;
    for (item = cal->seed_list->first->next; item != NULL; item = item->next) {
      seed = item->item;
      
      cal->num_open_gaps += seed->num_open_gaps;

      if (prev_seed->read_end >= seed->read_start) {
	mapping_batch->counters[2]++;
      }

      prev_seed = seed;
    }
    // last seed
    seed = linked_list_get_last(cal->seed_list);
    if (seed->read_end < read->length - 1) {
    }
    if (num_seeds > 1) {
      cal->num_open_gaps += seed->num_open_gaps;
    }
    if (cal->num_open_gaps == 0) mapping_batch->counters[3]++;
    if (cal->num_open_gaps == 1) mapping_batch->counters[4]++;
    if (cal->num_open_gaps > 1) mapping_batch->counters[5]++;


    cal->AS = 254;
    cigar_append_op(read->length, 'M', &cal->cigar);
  }
}

//--------------------------------------------------------------------

#define SW_LEFT_FLANK 5
#define SW_RIGHT_FLANK 5

#define SW_LEFT_FLANK_EX 20
#define SW_RIGHT_FLANK_EX 20

#define CIGAR_FROM_GAP  1
#define CIGAR_FROM_SEED 2

void update_left_side_seed(int min_flank, seed_t *seed) {
  int q_flank = 0, r_flank = 0;
  // look at left-side cigar
  cigar_t *cigar = &seed->cigar;
  int op_value, op_name, trim = 0, num_ops = cigar->num_ops;

  for (int i = 0; i < num_ops; i++) {
    cigar_get_op(i, &op_value, &op_name, cigar);
    if (op_name == 'M' && op_value > min_flank) {
      break;
    } else {
      trim++;
      if (op_name == 'M') {
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
    if (op_name == 'M' && op_value > min_flank) {
      break;
    } else {
      trim++;
      if (op_name == 'M') {
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
  int gap_read_len, gap_genome_len;
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

    // first seed
    seed = linked_list_get_first(cal->seed_list);
    //    if (seed->read_start > 0) {
    //      cal->start -= seed->read_start;
    //    }

    // cal cigar
    cigarset = cigarset_new(num_seeds * 2 + 1);
    cal->info = (void *) cigarset;

    if (seed->read_start > 0) {
      #ifdef _VERBOSE
      print_seed("-----> before updating first left-side seed: ", seed);
      #endif
      update_left_side_seed(SW_RIGHT_FLANK, seed);
      #ifdef _VERBOSE
      print_seed("-----> after updating first left-side seed: ", seed);
      #endif

      gap_genome_start = seed->genome_start - seed->read_start - 1 - SW_LEFT_FLANK_EX;
      gap_genome_end = seed->genome_start + SW_RIGHT_FLANK; //ref_flank; //SW_RIGHT_FLANK;
      ref = sa_genome_get_sequence(cal->chromosome_id, gap_genome_start, gap_genome_end, sa_index->genome);
      
      seq = get_subsequence((cal->strand ? read->revcomp : read->sequence), 
			    0, seed->read_start + SW_RIGHT_FLANK); //query_flank); //SW_RIGHT_FLANK);
      
      sw_prepare = sw_prepare_new(seq, ref, 0, SW_RIGHT_FLANK, FIRST_SW);
      sw_prepare->seed_region = seed;
      sw_prepare->cal = cal;
      sw_prepare->read = read;
      array_list_insert(sw_prepare, sw_prepare_list);
      num_sw++;

      //printf("case 4: %s\n", read->id);
      //exit(-1);

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
      print_seed("-----> before updating middle right-side seed: ", prev_seed);
      #endif
      #ifdef _VERBOSE
      print_seed("-----> before updating middle lefth-side seed: ", seed);
      #endif
      update_left_side_seed(SW_RIGHT_FLANK, seed);
      #ifdef _VERBOSE
      print_seed("-----> before updating middle left-side seed: ", seed);
      #endif


      gap_read_start = prev_seed->read_end + 1;
      gap_read_end = seed->read_start - 1;
      gap_read_len = gap_read_end - gap_read_start + 1;

      gap_genome_start = prev_seed->genome_end + 1;
      gap_genome_end = seed->genome_start - 1;
      gap_genome_len = gap_genome_end - gap_genome_start + 1;

      if (gap_read_len <= 0) {
	// deletion
	//	printf("gap_read_len = %i\n", gap_read_len);
	//	exit(-1);
	//	cigarset->active[seed_count * 2] = CIGAR_FROM_GAP; //1;
	//	cigarset->cigars[seed_count * 2] = cigar_new(gap_genome_len, 'D');	
	//	cal->num_open_gaps += 1;
	//	cal->num_extend_gaps += (gap_genome_len - 1);
	//	prev_seed = seed;

	seq = get_subsequence((cal->strand ? read->revcomp : read->sequence), 
			      gap_read_start - SW_LEFT_FLANK - abs(gap_read_len), 
			      gap_read_len + SW_LEFT_FLANK + SW_RIGHT_FLANK + (2*abs(gap_read_len)));

	ref = sa_genome_get_sequence(cal->chromosome_id, gap_genome_start - SW_LEFT_FLANK - abs(gap_read_len), 
				     gap_genome_end + SW_RIGHT_FLANK + abs(gap_read_len), sa_index->genome);

	
	cigarset_info_set(CIGAR_FROM_GAP, abs(gap_read_len), NULL, NULL, &cigarset->info[seed_count * 2]);

	//	printf("case 1: %s\n", read->id);
	//	exit(-1);
      } else if (gap_genome_len <= 0) {
	// insertion
	//	cigarset->active[seed_count * 2] = CIGAR_FROM_GAP; //1;
	//	cigarset->cigars[seed_count * 2] = cigar_new(gap_genome_len, 'I');	
	//	cal->num_open_gaps += 1;
	//	cal->num_extend_gaps += (gap_genome_len - 1);

	seq = get_subsequence((cal->strand ? read->revcomp : read->sequence), 
			      gap_read_start - SW_LEFT_FLANK, gap_read_len + SW_LEFT_FLANK + SW_RIGHT_FLANK);

	ref = sa_genome_get_sequence(cal->chromosome_id, gap_genome_start - SW_LEFT_FLANK - abs(gap_genome_len), 
				     gap_genome_end + SW_RIGHT_FLANK + abs(gap_genome_len), sa_index->genome);
            
	cigarset_info_set(CIGAR_FROM_GAP, 0, NULL, NULL, &cigarset->info[seed_count * 2]);

	//	printf("case 2: %s\n", read->id);
	//	exit(-1);
      } else {

        #ifdef _VERBOSE1
	print_seed("", prev_seed);
	print_seed("", seed);
	printf("read id = %s\n", read->id);
	printf("genome gap (start, end) = (%lu, %lu), len = %i\n", gap_genome_start, gap_genome_end, gap_genome_len);
	printf("read gap (start, end) = (%lu, %lu), len = %i\n", gap_read_start, gap_read_end, gap_read_len);
	exit(-1);
        #endif

	//assert(gap_read_len > 0);
	seq = get_subsequence((cal->strand ? read->revcomp : read->sequence), 
			      gap_read_start - SW_LEFT_FLANK, gap_read_len + SW_LEFT_FLANK + SW_RIGHT_FLANK);

	//assert(gap_genome_len > 0);
	ref = sa_genome_get_sequence(cal->chromosome_id, gap_genome_start - SW_LEFT_FLANK, 
				     gap_genome_end + SW_RIGHT_FLANK, sa_index->genome);

	cigarset_info_set(CIGAR_FROM_GAP, 0, NULL, NULL, &cigarset->info[seed_count * 2]);

	//      printf("case 3: %s\n", read->id);
	//      exit(-1);
      }
      // prepare MIDDLE_SW
      sw_prepare = sw_prepare_new(seq, ref, 0, 0, MIDDLE_SW);
      sw_prepare->seed_region = seed_count * 2;
      sw_prepare->cal = cal;
      sw_prepare->read = read;
      array_list_insert(sw_prepare, sw_prepare_list);
      num_sw++;

      prev_seed = seed;
      cigarset_info_set(CIGAR_FROM_SEED, 0, &seed->cigar, seed, &cigarset->info[seed_count * 2 + 1]);
    }
    // last seed
    seed = linked_list_get_last(cal->seed_list);
    if (seed->read_end < read->length - 1) {
      #ifdef _VERBOSE
      print_seed("-----> before updating last right-side seed: ", seed);
      #endif
      update_right_side_seed(SW_LEFT_FLANK, seed);
      #ifdef _VERBOSE
      print_seed("-----> after updating last right-side seed: ", seed);
      #endif

      gap_genome_start = seed->genome_end - SW_LEFT_FLANK + 1;
      gap_genome_end = gap_genome_start + (read->length - seed->read_end) + SW_LEFT_FLANK_EX;
      ref = sa_genome_get_sequence(cal->chromosome_id, gap_genome_start, gap_genome_end, sa_index->genome);

      seq = get_subsequence((cal->strand ? read->revcomp : read->sequence), 
			    seed->read_end - SW_LEFT_FLANK + 1, 
			    read->length + SW_LEFT_FLANK - seed->read_end);

      sw_prepare = sw_prepare_new(seq, ref, 0, 0, LAST_SW);
      sw_prepare->seed_region = num_seeds * 2;
      sw_prepare->cal = cal;
      sw_prepare->read = read;
      array_list_insert(sw_prepare, sw_prepare_list);
      num_sw++;

      cigarset_info_set(CIGAR_FROM_GAP, 0, NULL, NULL, &cigarset->info[num_seeds * 2]);

      //printf("case 5: %s\n", read->id);
      //exit(-1);
    } else {
      cigarset_info_set(0, 0, NULL, NULL, &cigarset->info[num_seeds * 2]);
    }

    if (num_sw == 0) {
      seed_cal_set_cigar_by_seed(seed, cal);
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
    cal = sw_prepare->cal;
    
    cigarset = cal->info;
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
      seed = sw_prepare->seed_region;
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
	printf("case query_start > 0: read %s\n", cal->read->id);
	//	exit(-1);
	if (ref_start > 0) {
	  printf("case query_start and ref_start > 0: read %s\n", cal->read->id);
	  //	  exit(-1);
	  if (diff == 0) {
	    cigar_append_op(query_start, 'M', cigar);      
	    //assert(1 == 0);
	  } else if (diff > 0) {
	    cigar_append_op(ref_start, 'M', cigar);      
	    cigar_append_op(diff, 'I', cigar);      
	    //	  assert(1 == 0);
	  } else {
	    cigar_append_op(query_start, 'M', cigar);      
	    cigar_append_op(abs(diff), 'D', cigar);      
	  }
	} else {
	  cigar_append_op(query_start, 'I', cigar);      
	}
      } else if (ref_start > 0) {
	//	printf("case ref_start > 0\n");
	//	exit(-1);
	//cigar_append_op(ref_start, 'D', cigar);      
	//      assert(1 == 0);
      }
    }

    // scan map to complete cigar
    op_value = 0;
    op_name = 'M';
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
      } else {
	if (op_name != 'M' && op_value > 0) {
	  cigar_append_op(op_value, op_name, cigar);
          #ifdef _VERBOSE
	  printf("****************** (value, name) = (%i, %c) : op. change: current cigar: %s\n", 
		 op_value, op_name, cigar_to_string(cigar));
          #endif
	  op_value = 0;
	}
	op_value++;
	op_name = 'M';
      }
    }
    cigar_append_op(op_value, op_name, cigar);

    /*
    if (sw_output->query_start_p[i] > 0 && sw_output->ref_start_p[i] > 0) {
      diff = sw_output->query_start_p[i] - sw_output->ref_start_p[i];
      if (diff < 0) {
	//	cigar_append_op(abs(diff), 'D', cigar);
	//	op_value = sw_output->query_start_p[i];
	//	op_name = 'M';
      } else if (diff > 0) {
	cigar_append_op(abs(diff), 'I', cigar);
	op_value = sw_output->query_start_p[i];
	op_name = 'M';
	cal->num_open_gaps++;
	cal->num_extend_gaps += (diff - 1);
      } else {
	op_value = sw_output->query_start_p[i];
	op_name = 'M';
      }
      if (op_value) {
	seq = &q[i][abs(diff)];
	ref = &r[i][abs(diff)];
	for(int j = 0; j < op_value; j++) {
	  if (seq[j] != ref[j]) {
	    cal->num_mismatches++;	
	  }
	}
      }
    } else if (sw_output->query_start_p[i] > 0) {
      op_value = sw_output->query_start_p[i];
      op_name = 'I';
    } else {
      //      op_value = sw_output->ref_start_p[i];
      //      op_name = 'D';
    }

    len = strlen(sw_output->query_map_p[i]);
    for(int j = 0; j < len; j++) {
      if (sw_output->query_map_p[i][j] == '-') {
	// deletion (in the query)
	if (op_name != 'D' && op_value > 0) {
	  cigar_append_op(op_value, op_name, cigar);
	  if (op_name == 'I') {
	    cal->num_open_gaps++;
	    cal->num_extend_gaps += (op_value - 1);
	  }
	  op_value = 0;
	  op_name = 'D';
	}
	op_value++;
	//cal->num_mismatches++;	
      } else if (sw_output->ref_map_p[i][j] == '-') {
	// insertion (in the query)
	if (op_name != 'I' && op_value > 0) {
	  cigar_append_op(op_value, op_name, cigar);
	  if (op_name == 'D') {
	    cal->num_open_gaps++;
	    cal->num_extend_gaps += (op_value - 1);
	  }
	  op_value = 0;
	  op_name = 'I';
	}
	op_value++;
	//	cal->num_mismatches++;	
      } else {
	if (op_name != 'M' && op_value > 0) {
	  cigar_append_op(op_value, op_name, cigar);
	  cal->num_open_gaps++;
	  cal->num_extend_gaps += (op_value - 1);
          #ifdef _VERBOSE
	  printf("****************** (value, name) = (%i, %c) : op. change: current cigar: %s\n", 
		 op_value, op_name, cigar_to_string(cigar));
          #endif
	  op_value = 0;
	  op_name = 'M';
	}
	op_value++;
	if (sw_output->query_map_p[i][j] != sw_output->ref_map_p[i][j]) {
	  cal->num_mismatches++;	
	}
      }
    }
    if (op_value > 0) {
      if (op_name != 'M') {
	cal->num_open_gaps++;
	cal->num_extend_gaps += (op_value - 1);
      }
      cigar_append_op(op_value, op_name, cigar);
    }
*/

    int gap_count = 0;
    if (sw_prepare->ref_type == FIRST_SW) {
      gap_count = 0;
    } else {
      gap_count = (int)sw_prepare->seed_region;
    }
    cigarset->info[gap_count].active = CIGAR_FROM_GAP;
    cigarset->info[gap_count].cigar = cigar;
    //    cigarset->info[gap_count].overlap = 0;
    #ifdef _VERBOSE
    printf("************** sw_count: %i, for gap %i cigar %s\n", 
	   i, gap_count, cigar_to_string(cigar));
    #endif

    // free
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
		     array_list_t **cal_lists) {

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
      cal->score = ((-4.0f) * cal->num_mismatches) 
	+ ((-10.0f) * cal->num_open_gaps) 
	+ ((-0.5f) * cal->num_extend_gaps);
      
      cal->AS = 254;
      
      #ifdef _VERBOSE
      printf("\t\t\t\tcal score = %0.2f\n", cal->score);
      #endif
      cigarset = cal->info;
      cigar = &cal->cigar;
      for (int j = 0; j < cigarset->size; j++) {
	cigar_type = cigarset->info[j].active;
	if (cigar_type > 0) {
	  if (cigar_type == CIGAR_FROM_SEED) {
	    seed = cigarset->info[j].seed;
	    aux_cigar = cigarset->info[j].cigar;
	    // CIGAR_FROM_SEED
	    if (seed->read_start > 0) {
	      cigar_get_op(0, &op_value, &op_name, aux_cigar);
	      if (op_value < SW_LEFT_FLANK || op_name != 'M') {
		printf("read = %s: attention: op (%i, %c) sw_left_flank!!!\n", cal->read->id, op_value, op_name);
		//		exit(-1);
	      }
	      cigar_set_op(0, op_value - SW_LEFT_FLANK, op_name, aux_cigar);
	    }
	    if (seed->read_end < cal->read->length - 1) {
	      cigar_get_op(aux_cigar->num_ops - 1, &op_value, &op_name, aux_cigar);
	      if (op_value < SW_RIGHT_FLANK || op_name != 'M') {
		printf("read = %s: attention: op (%i, %c) sw_right_flank!!!\n", cal->read->id, op_value, op_name);
		//		exit(-1);
	      }
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
	      if (op_value < cigarset->info[j].overlap || op_name != 'M') {
		printf("read = %s: attention OVERLAP: op (%i, %c) overlap (%i)!!!\n", cal->read->id, op_value, op_name, cigarset->info[j].overlap);
		exit(-1);
	      }
	      cigar_set_op(aux_cigar->num_ops - 1, op_value - cigarset->info[j].overlap, op_name, aux_cigar);
	    }
	    cigar_concat(aux_cigar, cigar);
	    cigar_free(aux_cigar);
	  }
	}
      }
      cigar_len = cigar_get_length(cigar);
      if (cigar_len > cal->read->length) {
	printf("cigar_len (%i : %s) > cal->read->length (%i): read %s\n", 
	       cigar_len, cigar_to_string(cigar), cal->read->length, cal->read->id);
	//	exit(-1);
      }
      if (cigar_len < cal->read->length) {
	if (seed && seed->read_end == cal->read->length - 1) {
	  cigar_append_op(cal->read->length - cigar_len, 'M', cigar);
	} else {
	  cigar_append_op(cal->read->length - cigar_len, 'S', cigar);
	}
      }
      cigarset_free(cigarset);
      cal->info = NULL;
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

int sa_mapper(void *data) {

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
  sa_index3_t *sa_index = (sa_index3_t *) wf_batch->sa_index;
 
  size_t num_reads = mapping_batch->num_reads;
  int max_read_area, min_num_mismatches;
  float max_score;

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

    // 1) extend using mini-sw from suffix
    cal_list = create_cals(num_seeds, read, mapping_batch, sa_index, cal_mng);

    if (array_list_size(cal_list) > 0) {
      min_num_mismatches = get_min_num_mismatches(cal_list);
      #ifdef _VERBOSE
      printf("\t*** before SW> min_num_mismatches = %i\n", min_num_mismatches);
      #endif
      
      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      filter_cals_by_max_num_mismatches(min_num_mismatches, &cal_list);
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_FILTER_BY_NUM_MISMATCHES] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif

      clean_cals(cal_list, read, sa_index);

      // 2) prepare Smith-Waterman to fill in the gaps
      if (prepare_sw(read, sw_prepare_list, mapping_batch, sa_index, cal_list)) {
	sw_post_read[sw_post_read_counter++] = i;
      }
    }
    cal_lists[i] = cal_list;
  }

  // 3) run SW to fill
  if (array_list_size(sw_prepare_list) > 0) {
    execute_sw(sw_prepare_list, mapping_batch);
    post_process_sw(sw_post_read_counter, sw_post_read, cal_lists);
  }
  array_list_free(sw_prepare_list, (void *) NULL);

  // 4) for each read, create alignments
  for (int i = 0; i < num_reads; i++) {
    cal_list = cal_lists[i];
    read = array_list_get(i, mapping_batch->fq_reads);

    if (array_list_size(cal_list) > 0) {
      // filter by score
      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      max_score = get_max_score(cal_list);
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

    // create alignments structures
    #ifdef _TIMING
    gettimeofday(&start, NULL);
    #endif
    create_alignments(cal_list, read, mapping_batch->mapping_lists[i]);
    #ifdef _TIMING
    gettimeofday(&stop, NULL);
    mapping_batch->func_times[FUNC_CREATE_ALIGNMENTS] += 
      ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
    #endif
    
    // free cal list and clear seed manager for next read
    array_list_free(cal_list, (void *) NULL);
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
