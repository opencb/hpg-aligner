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
	  linked_list_free(p->cals_lists[i], cal_free_ex);
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
	  linked_list_clear(p->cals_lists[i], cal_free_ex);
	}
      }
    }
  }
}

//--------------------------------------------------------------------

void cal_mng_update(seed_t *seed, cal_mng_t *p) {
  int is_used = 0;
  if (p->cals_lists) {
    cal_t *cal;
    seed_t *s_last;
    linked_list_t *sr_list;
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
	sr_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	linked_list_insert(seed, sr_list);
	cal = cal_new(seed->chromosome_id, seed->strand, 
		      seed->genome_start, seed->genome_end, 1, sr_list, NULL);
	cal->read_area = seed->read_end - seed->read_start + 1;
	cal->num_mismatches = seed->num_mismatches;
	linked_list_insert(cal, cal_list);
      } else {
	// insert (by order)
	linked_list_iterator_t* itr = linked_list_iterator_new(cal_list);
	cal_t *item = (cal_t *) linked_list_iterator_curr(itr);
	while (item != NULL) {
	  #ifdef _VERBOSE1
	  printf("---> merging with this CAL?\n");
	  cal_print(item);
	  #endif
	  //	  assert(cal->end > item->start);
	  s_last = linked_list_get_last(item->sr_list);
	  if (seed->read_start - s_last->read_end < 100 && 
	      seed->genome_start - s_last->genome_end < 100) {
	    is_used = 1;
	    linked_list_insert_last(seed, item->sr_list);
	    item->end = seed->genome_end;
	    item->read_area += (seed->read_end - seed->read_start + 1);
	    item->num_mismatches += seed->num_mismatches;

	    //	    cal_free(cal);
            #ifdef _VERBOSE
	    printf("---> yes, merging CAL, result:\n");
	    cal_print(item);
	    #endif
	    break;
	  } else {
	    if (seed->genome_end < item->start) {
	      // create CAL and insert it into the CAL manager
	      is_used = 1;
	      sr_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	      linked_list_insert(seed, sr_list);
	      cal = cal_new(seed->chromosome_id, seed->strand, seed->genome_start, seed->genome_end, 
			    1, sr_list, NULL);
	      cal->read_area = seed->read_end - seed->read_start + 1;
	      cal->num_mismatches = seed->num_mismatches;
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
	  sr_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	  linked_list_insert(seed, sr_list);
	  cal = cal_new(seed->chromosome_id, seed->strand, seed->genome_start, seed->genome_end, 1, sr_list, NULL);
	  cal->read_area = seed->read_end - seed->read_start + 1;
	  cal->num_mismatches = seed->num_mismatches;
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
      cal_t *cal;
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
  cal_t *cal;
  linked_list_iterator_t itr;

  if (p->cals_lists) {
    linked_list_t *cal_list;
    for (unsigned int i = 0; i < p->num_chroms; i++) {
      cal_list = p->cals_lists[i];
      while (cal = (cal_t *) linked_list_remove_last(cal_list)) {
	if ((cal->end - cal->start) >= read_area) {
	  array_list_insert(cal, out_list);
	} else {
	  // free CAL
	  cal_free_ex(cal);
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
  cal_t *cal;
  linked_list_iterator_t itr;

  if (p->cals_lists) {
    linked_list_t *cal_list;
    for (unsigned int i = 0; i < p->num_chroms; i++) {
      cal_list = p->cals_lists[i];
      while (cal = (cal_t *) linked_list_remove_last(cal_list)) {
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

void generate_cals_from_exact_read(int strand, fastq_read_t *read, char *revcomp,
				   size_t low, size_t high, sa_index3_t *sa_index, 
				   array_list_t *cal_list) {  
  size_t g_start, g_end;
  int chrom;

  cal_t *cal;
  cigar_t *cigar;
  seed_t *seed;
  linked_list_t *sr_list;
  
  for (size_t suff = low; suff <= high; suff++) {
    chrom = sa_index->CHROM[suff];
    g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];
    g_end = g_start + read->length - 1;
    sr_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    seed = seed_new(0, read->length - 1, g_start, g_end);
    cigar_append_op(read->length, 'M', &seed->cigar);
    linked_list_insert(seed, sr_list);
    cal = cal_new(chrom, strand, g_start, g_end, 1, sr_list, NULL);
    cal->read_area = read->length;
    cal->num_mismatches = 0;
    array_list_insert(cal, cal_list);
  }
}

//--------------------------------------------------------------------
// generate cals extending suffixes to left and right side 
// bt using the mini-sw
//--------------------------------------------------------------------

int generate_cals_from_suffixes(int strand, fastq_read_t *read, char *revcomp,
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
  r_seq = (strand ? revcomp : read->sequence);
  
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
      cal_mng_update(seed, cal_mng);
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
  return (sa_index->k_value / 2);
  //  return (sa_index->k_value * 4);
}

//--------------------------------------------------------------------
// step one:
//    search prefix + search suffix + extend using min-sw
//--------------------------------------------------------------------

array_list_t *step_one(fastq_read_t *read, char *revcomp_seq,
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
  int read_pos, read_end_pos, read_inc = sa_index->k_value / 2;
 
  #ifdef _VERBOSE	  
  printf("\n\n====>>>> STEP ONE <<<<====\n");
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
	  generate_cals_from_exact_read(strand, read, revcomp_seq,
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
	}
      } else {
	read_pos += read_inc;
      }
    } // end of for read_pos
    
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
    r_seq = revcomp_seq;
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
// step three:
//    fill gaps using sw
//--------------------------------------------------------------------

void step_three(fastq_read_t *read, char *revcomp_seq, 
		sa_mapping_batch_t *mapping_batch, sa_index3_t *sa_index, 
		array_list_t *cal_list) {
  size_t seed_count, num_seeds, num_cals = array_list_size(cal_list);

  #ifdef _VERBOSE	  
  printf("\n\n====>>>> STEP THREE <<<<====\n");
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
  array_list_t *sw_prepare_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  cal_t *cal;
  cigar_t *cigar;
  cigarset_t *cigarset;

  for (int i = 0; i < num_cals; i++) {
    cal = array_list_get(i, cal_list);
    #ifdef _VERBOSE
    cal_print(cal);
    #endif

    // if not seeds, then next cal
    num_seeds = cal->sr_list->size;
    if (num_seeds <= 0) continue;

    // cal cigar
    cigarset = cigarset_new(num_seeds * 2 + 1);
    cal->info = (void *) cigarset;

    // first seed
    seed = linked_list_get_first(cal->sr_list);

    if (seed->read_start > 0) {
      gap_genome_start = seed->genome_start - seed->read_start - 1;// - 10;
      gap_genome_end = seed->genome_start - 1;
      ref = sa_genome_get_sequence(cal->chromosome_id, gap_genome_start, gap_genome_end, sa_index->genome);
      
      seq = get_subsequence((cal->strand ? revcomp_seq : read->sequence), 
			    0, seed->read_start + 1);
      
      sw_prepare = sw_prepare_new(seq, ref, 0, 0, FIRST_SW);
      sw_prepare->seed_region = 0;
      sw_prepare->cal = cal;
      sw_prepare->read = read;
      array_list_insert(sw_prepare, sw_prepare_list);

      cigarset->active[0] = 1;
    }
    // seeds at the middle positions
    prev_seed = seed;
    cigarset->active[1] = 2;
    cigarset->cigars[1] = &seed->cigar;
    seed_count = 0;
    for (item = cal->sr_list->first->next; item != NULL; item = item->next) {
      seed_count++;
      seed = item->item;

      gap_read_start = prev_seed->read_end + 1;
      gap_read_end = seed->read_start - 1;
      gap_read_len = gap_read_end - gap_read_start + 1;

      gap_genome_start = prev_seed->genome_end + 1;
      gap_genome_end = seed->genome_start - 1;
      gap_genome_len = gap_genome_end - gap_genome_start + 1;

      if (gap_read_len <= 0) {
	// deletion
	cigarset->active[seed_count * 2] = 1;
	cigarset->cigars[seed_count * 2] = cigar_new(gap_genome_len, 'D');	
	prev_seed = seed;
	continue;
      }

      if (gap_genome_len <= 0) {
	// insertion
	cigarset->active[seed_count * 2] = 1;
	cigarset->cigars[seed_count * 2] = cigar_new(gap_genome_len, 'I');	
	prev_seed = seed;
	continue;
      }

      #ifdef _VERBOSE1
      print_seed("", prev_seed);
      print_seed("", seed);
      printf("read id = %s\n", read->id);
      printf("genome gap (start, end) = (%lu, %lu), len = %i\n", gap_genome_start, gap_genome_end, gap_genome_len);
      printf("read gap (start, end) = (%lu, %lu), len = %i\n", gap_read_start, gap_read_end, gap_read_len);
      exit(-1);
      #endif

      assert(gap_read_len > 0);
      seq = get_subsequence((cal->strand ? revcomp_seq : read->sequence), 
			    gap_read_start, gap_read_len);

      assert(gap_genome_len > 0);
      ref = sa_genome_get_sequence(cal->chromosome_id, gap_genome_start, gap_genome_end, sa_index->genome);
      
      
      sw_prepare = sw_prepare_new(seq, ref, 0, 0, MIDDLE_SW);
      sw_prepare->seed_region = seed_count * 2;
      sw_prepare->cal = cal;
      sw_prepare->read = read;
      array_list_insert(sw_prepare, sw_prepare_list);
      cigarset->active[seed_count * 2] = 1;

      prev_seed = seed;
      cigarset->active[seed_count * 2 + 1] = 2;
      cigarset->cigars[seed_count * 2 + 1] = &seed->cigar;
    }
    // last seed
    seed = linked_list_get_last(cal->sr_list);
    if (seed->read_end < read->length - 1) {
      gap_genome_start = seed->genome_end + 1;
      gap_genome_end = gap_genome_start + (read->length - seed->read_end);// + 10;
      ref = sa_genome_get_sequence(cal->chromosome_id, gap_genome_start, gap_genome_end, sa_index->genome);

      seq = get_subsequence((cal->strand ? revcomp_seq : read->sequence), 
			    seed->read_end + 1, read->length - seed->read_end - 1);

      sw_prepare = sw_prepare_new(seq, ref, 0, 0, LAST_SW);
      sw_prepare->seed_region = num_seeds * 2;
      sw_prepare->cal = cal;
      sw_prepare->read = read;
      array_list_insert(sw_prepare, sw_prepare_list);

      cigarset->active[num_seeds * 2] = 1;
    }
  }

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
  int op_name, op_value, diff, len;
  for (int i = 0; i < sw_count; i++) {
    sw_prepare = array_list_get(i, sw_prepare_list);
    cal = sw_prepare->cal;
    
    cigarset = cal->info;
    cigar = cigar_new_empty();
    op_value = 0;
    op_name = 'M';
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
	  op_value = 0;
	  op_name = 'D';
	}
	op_value++;
	cal->num_mismatches++;	
      } else if (sw_output->ref_map_p[i][j] == '-') {
	// insertion (in the query)
	if (op_name != 'I' && op_value > 0) {
	  cigar_append_op(op_value, op_name, cigar);
	  op_value = 0;
	  op_name = 'I';
	}
	op_value++;
	cal->num_mismatches++;	
      } else {
	if (op_name != 'M' && op_value > 0) {
	  cigar_append_op(op_value, op_name, cigar);
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
      cigar_append_op(op_value, op_name, cigar);
    }
    cigarset->active[(int)sw_prepare->seed_region] = 1;
    cigarset->cigars[(int)sw_prepare->seed_region] = cigar;
    #ifdef _VERBOSE
    printf("************** sw_count: %i, cigar for gap %i: %s\n", 
	   i, sw_prepare->seed_region, cigar_to_string(cigar));
    #endif

    // free
    sw_prepare_free(sw_prepare);
  }

  for (int i = 0; i < num_cals; i++) {
    cal = array_list_get(i, cal_list);
    cigar = cigar_new_empty();
    cigarset = cal->info;
    for (int j = 0; j < cigarset->num_cigars; j++) {
      if (cigarset->active[j] > 0) {
	#ifdef _VERBOSE
	printf("************** gap %i -> concat cigar %s into %s\n",
	       j, cigar_to_string(cigarset->cigars[j]), cigar_to_string(cigar));
	#endif
	cigar_concat(cigarset->cigars[j], cigar);
	if (cigarset->active[j] == 1) {
	  cigar_free(cigarset->cigars[j]);
	}
      }
    }
    cigarset_free(cigarset);
    cal->info = cigar;
  }
  // free memory
  sw_multi_output_free(sw_output);
  array_list_free(sw_prepare_list, (void *) NULL);

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
  
  sa_mapping_batch_t *mapping_batch = wf_batch->mapping_batch;
  sa_index3_t *sa_index = (sa_index3_t *) wf_batch->sa_index;
 
  size_t num_reads = mapping_batch->num_reads;
  int max_read_area, min_num_mismatches;
  
  // CAL management
  size_t num_cals;
  cal_t *cal;
  cal_mng_t *cal_mng;
  array_list_t *cal_list;

  fastq_read_t *read;

  //  int saved_pos[2][1024];
  // TODO !!! 20 = min. cal size
  uint min_cal_size = 20;
  cal_mng = cal_mng_new(sa_index->genome);
  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_OTHER] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif

  // for each read
  for (int i = 0; i < num_reads; i++) {
    read = array_list_get(i, mapping_batch->fq_reads);
    //printf("read id = %s\n", read->id);
    //max_num_mismatches = read->length * MISMATCH_PERC;

    // step one:: extend using mini-sw from suffix
    cal_list = step_one(read, mapping_batch->revcomp_seqs[i], mapping_batch, sa_index, cal_mng);

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

      // step three (Smith-Waterman to fill in the gaps)
      step_three(read, mapping_batch->revcomp_seqs[i], mapping_batch, sa_index, cal_list);

      min_num_mismatches = get_min_num_mismatches(cal_list);
      #ifdef _VERBOSE
      printf("\t*** after SW> min_num_mismatches = %i\n", min_num_mismatches);
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
    }

    //    printf("******************* cal list size = %i\n", array_list_size(cal_list));
    //    exit(-1);
    /*
    if ((num_cals = array_list_size(cal_list)) <= 0) {
      // step two: extend using mini-sw from prefix
      array_list_t *invalid_cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
      step_two(read, mapping_batch->revcomp_seqs[i], mapping_batch, sa_index, cal_mng, 
	       cal_list, invalid_cal_list);

      //      printf("********************** after step_two num. valid = %i, invalid = %i\n", 
      //	     array_list_size(cal_list), array_list_size(invalid_cal_list));

      if (array_list_size(cal_list) > 0) {
	array_list_free(invalid_cal_list, (void *) cal_free);
      } else {
	// step three: apply smith-waterman
	array_list_free(cal_list, (void *) NULL);	
	//step_three(read, mapping_batch->revcomp_seqs[i], mapping_batch, sa_index, invalid_cal_list);
	cal_list = invalid_cal_list;

        #ifdef _TIMING
	gettimeofday(&start, NULL);
        #endif
	filter_cals_by_max_read_area(max_read_area, &cal_list);
        #ifdef _TIMING
	gettimeofday(&stop, NULL);
	mapping_batch->func_times[FUNC_FILTER_BY_READ_AREA] += 
	  ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
        #endif
      }
    */

    /*
    if (array_list_size(cal_list) > 0) {
      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      filter_min_read_area_cals(&cal_list);
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_FILTER_BY_READ_AREA] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif
    }
    */
    // fill in gaps on CALs and create alignments structures
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
