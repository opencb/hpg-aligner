#include "cal_seeker.h"

//------------------------------------------------------------------------------------
// functions to handle gaps between mapped seeds
//    - fill_gaps
//    - merge_seed_regions
//------------------------------------------------------------------------------------

void display_sr_lists(char *msg, mapping_batch_t *mapping_batch) {

  fastq_read_t *read;
  array_list_t *fq_batch = mapping_batch->fq_batch;

  size_t read_index;

  cal_t *cal;
  array_list_t *cal_list = NULL;
  size_t num_cals, num_targets = mapping_batch->num_targets;

  seed_region_t *s, *prev_s, *new_s;
  linked_list_iterator_t* itr;

  LOG_DEBUG_F("%s\n", msg);

  // debugging....
  for (size_t i = 0; i < num_targets; i++) {
    read_index = mapping_batch->targets[i];
    read = (fastq_read_t *) array_list_get(read_index, fq_batch);
    
    LOG_DEBUG_F("Read %s\n", read->id);
    
    cal_list = mapping_batch->mapping_lists[read_index];
    num_cals = array_list_size(cal_list);
    
    if (num_cals <= 0) continue;
    
    // processing each CAL from this read
    for(size_t j = 0; j < num_cals; j++) {
      
      // get cal and read index
      cal = array_list_get(j, cal_list);
      LOG_DEBUG_F("\tCAL #%i of %i (strand %i, chr %i:%lu-%lu), sr_list size = %i\n", 
		  j, num_cals, cal->strand, cal->chromosome_id - 1, cal->start, cal->end, cal->sr_list->size);
      itr = linked_list_iterator_new(cal->sr_list);
      s = (seed_region_t *) linked_list_iterator_curr(itr);
      while (s != NULL) {
	LOG_DEBUG_F("\t\t%s (dist. %i)\t[%i|%i - %i|%i]\n", 
		    (s->info ? new_cigar_code_string((cigar_code_t *) s->info) : ">>>>>> gap"),
		    (s->info ? ((cigar_code_t *) s->info)->distance : -1),
		    s->genome_start, s->read_start, s->read_end, s->genome_end);
	linked_list_iterator_next(itr);
	s = linked_list_iterator_curr(itr);
      }
      linked_list_iterator_free(itr);
    }
  }
}
//------------------------------------------------------------------------------------

void fill_gaps(mapping_batch_t *mapping_batch, sw_optarg_t *sw_optarg, 
	       genome_t *genome, int min_gap, int min_distance) {

  int sw_count = 0;

  fastq_read_t *read;
  array_list_t *fq_batch = mapping_batch->fq_batch;

  size_t read_index, read_len;

  cal_t *cal;
  array_list_t *cal_list = NULL;
  size_t num_cals, num_targets = mapping_batch->num_targets;

  char *revcomp_seq = NULL;

  seed_region_t *s, *prev_s, *new_s;
  linked_list_iterator_t* itr;

  cigar_code_t *cigar_code;

  size_t start, end;
  size_t gap_read_start, gap_read_end, gap_read_len;
  size_t gap_genome_start, gap_genome_end, gap_genome_len;

  int left_flank, right_flank;
  sw_prepare_t *sw_prepare;
  array_list_t *sw_prepare_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  char *query,  *ref;
  int distance, first, last;

  //  LOG_DEBUG("\n\n P R E   -   P R O C E S S\n");

  // initialize query and reference sequences to Smith-Waterman
  for (size_t i = 0; i < num_targets; i++) {

    read_index = mapping_batch->targets[i];
    read = (fastq_read_t *) array_list_get(read_index, fq_batch);
    
    cal_list = mapping_batch->mapping_lists[read_index];
    num_cals = array_list_size(cal_list);
    
    if (num_cals <= 0) continue;

    read_len = read->length;

    min_distance = read_len*0.2;

    LOG_DEBUG_F(">>>>> read %s\n", read->id);
    //    printf(">>>>> read %s\n", read->id);

    // processing each CAL from this read
    for(size_t j = 0; j < num_cals; j++) {

      // get cal and read index
      cal = array_list_get(j, cal_list);
      LOG_DEBUG_F("CAL #%i of %i (strand %i), sr_list size = %i, sr_duplicate_list size = %i\n", 
		  j, num_cals, cal->strand, cal->sr_list->size, cal->sr_duplicate_list->size);

      prev_s = NULL;
      itr = linked_list_iterator_new(cal->sr_list);
      s = (seed_region_t *) linked_list_iterator_curr(itr);
      while (s != NULL) {
	{
	  // for debugging
	  size_t start = s->genome_start;// + 1;
	  size_t end = s->genome_end;// + 1;
	  size_t len = end - start + 1;
	  //	  printf(":::::::::: %lu - %lu = %i ::::::::::::\n", end, start, len );
	  char *ref = (char *) malloc((len + 1) * sizeof(char));
	  genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
					    &start, &end, genome);
	  ref[len] = '\0';
	  //
	  LOG_DEBUG_F("\tseed: [%i|%i - %i|%i] %s (len = %i)\n", 
		      s->genome_start, s->read_start, s->read_end, s->genome_end, ref, len);
	  free(ref);
	}

	// set the cigar for the current region
	gap_read_len = s->read_end - s->read_start + 1;
	cigar_code = cigar_code_new();
	cigar_code_append_op(cigar_op_new(gap_read_len, 'M'), cigar_code);
	s->info = (void *) cigar_code;

	cigar_code = NULL;
	sw_prepare = NULL;

	if ((prev_s == NULL && s->read_start != 0) || (prev_s != NULL)) {
	  distance = 0;
	  mapping_batch->num_gaps++;
	  if (prev_s == NULL) {
	    // gap at the first position
	    gap_read_start = 0;
	    gap_read_end = s->read_start - 1;

	    gap_genome_start = s->genome_start - s->read_start;
	    gap_genome_end = s->genome_start - 1;

	    gap_read_len = gap_read_end - gap_read_start + 1;
	    gap_genome_len = gap_genome_end - gap_genome_start + 1;

	    cal->start = gap_genome_start;

	    assert(gap_read_len != 0);
	    assert(gap_genome_len != 0);

	    if (gap_read_len > min_gap) {
	      // the gap is too big, may be there's another CAL to cover it
	      cigar_code = cigar_code_new();
	      cigar_code_append_op(cigar_op_new(gap_read_len, 'H'), cigar_code);	      
	    } else {
	      left_flank = 0;
	      right_flank = DOUBLE_FLANK;
	    }
	  } else {
	    assert(prev_s->read_end < s->read_start);

	    // gap in a middle position
	    gap_read_start = prev_s->read_end + 1;
	    gap_read_end = s->read_start - 1;

	    gap_genome_start = prev_s->genome_end + 1;
	    gap_genome_end = s->genome_start - 1;

	    gap_read_len = gap_read_end - gap_read_start + 1;
	    gap_genome_len = gap_genome_end - gap_genome_start + 1;

	    LOG_DEBUG_F("gap (read, genome) = (%i, %i)\n", gap_read_len, gap_genome_len);

	    if (gap_genome_len == 0) { printf("#@#: %s\n", read->id); }
	    assert(gap_genome_len != 0);

	    if (gap_read_len == 0) {
	      // there's a deletion just between two consecutives seeds
	      cigar_code = (cigar_code_t *)prev_s->info;

	      cigar_code_append_op(cigar_op_new(gap_genome_len, 'D'), cigar_code);
	      cigar_code->distance += gap_genome_len;

	      cigar_code_append_op(cigar_op_new(s->read_end - s->read_start + 1, 'M'), cigar_code);
	      cigar_code->distance += ((cigar_code_t *)s->info)->distance;

	      prev_s->read_end = s->read_end;
	      prev_s->genome_end = s->genome_end;

	      LOG_DEBUG_F("prev cigar = %s\n", new_cigar_code_string((cigar_code_t *)prev_s->info));

	      // continue loop...
	      linked_list_iterator_remove(itr);
	      s = linked_list_iterator_curr(itr);
	      continue;
	    }
	      
	    left_flank = SINGLE_FLANK;
	    right_flank = SINGLE_FLANK;
	  }

	  if (!cigar_code) {
	    // we have to try to fill this gap and get a cigar
	    if (gap_read_len == gap_genome_len) {
	      //    1) first, for from  begin -> end, and begin <- end
	      start = gap_genome_start;// + 1;
	      end = gap_genome_end;// + 1;
	      first = -1;
	      last = -1;
	      ref = (char *) malloc((gap_genome_len + 5) * sizeof(char));
	      genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
						&start, &end, genome);
	      // handle strand -
	      if (cal->strand) {
		if (revcomp_seq == NULL) {
		  revcomp_seq = strdup(read->sequence);
		  seq_reverse_complementary(revcomp_seq, read_len);
		}
		query = &revcomp_seq[gap_read_start];
	      } else {
		query = &read->sequence[gap_read_start];
	      }
	      
	      for (int k = 0; k < gap_read_len; k++) {
		if (query[k] != ref[k]) {
		  distance++;
		  if (first == -1) first = k;
		  last = k;
		}
	      }

	      if (distance < min_distance) {
		cigar_code = cigar_code_new();
		cigar_code_append_op(cigar_op_new(gap_read_len, 'M'), cigar_code);
		cigar_code_inc_distance(distance, cigar_code);
	      }
	    }
	    if (!cigar_code) {
	      //    2) second, prepare SW to run

	      // get query sequence, revcomp if necessary
	      size_t read_start = gap_read_start - left_flank;
	      size_t read_end = gap_read_end + right_flank;
	      int gap_read_len_ex = read_end - read_start + 1;
	      query = (char *) malloc((gap_read_len_ex + 1) * sizeof(char));
	      // handle strand -
	      if (cal->strand) {
		if (revcomp_seq == NULL) {
		  revcomp_seq = strdup(read->sequence);
		  seq_reverse_complementary(revcomp_seq, read_len);
		}
		memcpy(query, &revcomp_seq[read_start], gap_read_len_ex);
	      } else {
		memcpy(query, &read->sequence[read_start], gap_read_len_ex);
	      }
	      query[gap_read_len_ex] = '\0';
	      
	      // get ref. sequence
	      size_t genome_start = gap_genome_start - left_flank;// + 1;
	      size_t genome_end = gap_genome_end + right_flank;// + 1;
	      int gap_genome_len_ex = genome_end - genome_start + 1;
	      ref = (char *) malloc((gap_genome_len_ex + 1) * sizeof(char));;
	      genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
						&genome_start, &genome_end, genome);	      
	      ref[gap_genome_len_ex] = '\0';

	      if (prev_s == NULL) {
		sw_prepare = sw_prepare_new(query, ref, left_flank, right_flank, FIRST_SW);
	      } else {
		sw_prepare = sw_prepare_new(query, ref, left_flank, right_flank, MIDDLE_SW);
	      }

	      array_list_insert(sw_prepare, sw_prepare_list);
	      
	      // increase counter
	      sw_count++;	  

	      LOG_DEBUG_F("query: %s\n", query);
	      LOG_DEBUG_F("ref  : %s\n", ref);
	      LOG_DEBUG_F("dist.: %i (min. %i) of %i (first = %i, last = %i)\n", 
			  distance, min_distance, gap_read_len, first, last);
	      LOG_DEBUG_F("\tto SW (read %lu-%lu, genome %lu-%lu) = (%i, %i): read %s\n", 
			  gap_read_start, gap_read_end, gap_genome_start, gap_genome_end,
			  gap_read_end - gap_read_start + 1, gap_genome_end - gap_genome_start + 1, 
			  read->id);

	    }
	  }
	  
	  // insert gap in the list
	  new_s = seed_region_new(gap_read_start, gap_read_end, gap_genome_start, gap_genome_end, 0, 0, 0);
	  new_s->info = (void *) cigar_code;
	  linked_list_iterator_insert(new_s, itr);

	  if (sw_prepare) {
	    sw_prepare->seed_region = new_s;
	    sw_prepare->cal = cal;
	    sw_prepare->read = read;
	  }
	}

	// continue loop...
	prev_s = s;
	linked_list_iterator_next(itr);
	s = linked_list_iterator_curr(itr);
      }

      // check for a gap at the last position
      sw_prepare = NULL;
      if (prev_s != NULL && prev_s->read_end < read_len - 1) { 
	cigar_code = NULL;
	mapping_batch->num_gaps++;
	//	mapping_batch->num_sws++;
	//	mapping_batch->num_ext_sws++;

	// gap at the last position
	gap_read_start = prev_s->read_end + 1;
	gap_read_end = read_len - 1;
	gap_read_len = gap_read_end - gap_read_start + 1;

	assert(gap_read_len != 0);

	gap_genome_len = gap_read_len;
	gap_genome_start = prev_s->genome_end + 1;
	gap_genome_end = gap_genome_start + gap_genome_len - 1;

	cal->end = gap_genome_end;

	assert(gap_genome_len != 0);

	//	LOG_DEBUG_F("\t\tgap_read_len = %i, gap_genome_len = %i\n", gap_read_len, gap_genome_len);
	//	LOG_DEBUG_F("\t\t%i : [%lu|%lu - %lu|%lu]\n", 
	//		    sw_count, gap_genome_start, gap_read_start, gap_read_end, gap_genome_end);

	if (gap_read_len > min_gap) {
	  // the gap is too big, may be there's another CAL to cover it
	  cigar_code = cigar_code_new();
	  cigar_code_append_op(cigar_op_new(gap_read_len, 'H'), cigar_code);	      
	} else {
	  // we have to try to fill this gap and get a cigar
	  
	  //    1) first, for from  begin -> end, and begin <- end
	  start = gap_genome_start;// + 1;
	  end = gap_genome_end;// + 1;
	  first = -1;
	  last = -1;
	  ref = (char *) malloc((gap_genome_len + 1) * sizeof(char));;
	  genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
					    &start, &end, genome);
	  // handle strand -
	  if (cal->strand) {
	    if (revcomp_seq == NULL) {
	      revcomp_seq = strdup(read->sequence);
	      seq_reverse_complementary(revcomp_seq, read_len);
	    }
	    query = &revcomp_seq[gap_read_start];
	  } else {
	    query = &read->sequence[gap_read_start];
	  }
	  
	  distance = 0;
	  for (int k = 0; k < gap_read_len; k++) {
	    if (query[k] != ref[k]) {
	      distance++;
	      if (first == -1) first = k;
	      last = k;
	    }
	  }
	  if (distance < min_distance) {
	    cigar_code = cigar_code_new();
	    cigar_code_append_op(cigar_op_new(gap_read_len, 'M'), cigar_code);
	    cigar_code_inc_distance(distance, cigar_code);
	  } else {
	    //    2) second, prepare SW to run

	    left_flank = DOUBLE_FLANK;
	    right_flank = 0;
	    
	    // get query sequence, revcomp if necessary
	    size_t read_start = gap_read_start - left_flank;
	    size_t read_end = gap_read_end + right_flank;
	    int gap_read_len_ex = read_end - read_start + 1;
	    query = (char *) malloc((gap_read_len_ex + 1) * sizeof(char));
	    // handle strand -
	    if (cal->strand) {
	      if (revcomp_seq == NULL) {
		revcomp_seq = strdup(read->sequence);
		seq_reverse_complementary(revcomp_seq, read_len);
	      }
	      memcpy(query, &revcomp_seq[read_start], gap_read_len_ex);
	    } else {
	      memcpy(query, &read->sequence[read_start], gap_read_len_ex);
	    }
	    query[gap_read_len_ex] = '\0';
	    
	    // get ref. sequence
	    size_t genome_start = gap_genome_start - left_flank;// + 1;
	    size_t genome_end = gap_genome_end + right_flank;// + 1;
	    int gap_genome_len_ex = genome_end - genome_start + 1;
	    ref = (char *) malloc((gap_genome_len_ex + 1) * sizeof(char));;
	    genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
					      &genome_start, &genome_end, genome);
	    query[gap_genome_len_ex] = '\0';

	    sw_prepare = sw_prepare_new(query, ref, left_flank, right_flank, LAST_SW);
	    array_list_insert(sw_prepare, sw_prepare_list);
	    
	    // increase counter
	    sw_count++;	  

	    LOG_DEBUG_F("query: %s\n", query);
	    LOG_DEBUG_F("ref  : %s\n", ref);
	    LOG_DEBUG_F("dist.: %i (min. %i) of %i (first = %i, last = %i)\n", 
			distance, min_distance, gap_read_len, first, last);
	    LOG_DEBUG_F("\tto SW (read %lu-%lu, genome %lu-%lu) = (%i, %i): read %s\n", 
			gap_read_start, gap_read_end, gap_genome_start, gap_genome_end,
			gap_read_end - gap_read_start + 1, gap_genome_end - gap_genome_start + 1, 
			read->id);
	  }
	}
	
	// insert gap in the list
	new_s = seed_region_new(gap_read_start, gap_read_end, gap_genome_start, gap_genome_end, 0, 0, 0);
	new_s->info = (void *) cigar_code;
	linked_list_insert_last(new_s, cal->sr_list);

	if (sw_prepare) {
	  sw_prepare->seed_region = new_s;
	  sw_prepare->cal = cal;
	  sw_prepare->read = read;
	}
      }
      linked_list_iterator_free(itr);      
    }

    // free memory
    if (revcomp_seq) {
      free(revcomp_seq);
      revcomp_seq = NULL;
    }
  }

  //  display_sr_lists("ATER pre-process in fill_gaps", mapping_batch);

  LOG_DEBUG_F("\nR U N   S W (sw_count = %i, sw_prepare_list size = %i)\n", sw_count, array_list_size(sw_prepare_list));
  assert(sw_count == array_list_size(sw_prepare_list));

  char *q[sw_count], *r[sw_count];
  for (int i = 0; i < sw_count; i++) {
    sw_prepare = array_list_get(i, sw_prepare_list);
    q[i] = sw_prepare->query;
    r[i] = sw_prepare->ref;
  }
  sw_multi_output_t *output = sw_multi_output_new(sw_count);

  // run Smith-Waterman
  smith_waterman_mqmr(q, r, sw_count, sw_optarg, 1, output);
  
  LOG_DEBUG("P O S T   -   P R O C E S S\n");
  cigar_op_t* cigar_op;
  cigar_code_t *cigar_c;
  for (int i = 0; i < sw_count; i++) {
    sw_prepare = array_list_get(i, sw_prepare_list);
    s = sw_prepare->seed_region;

    int read_gap_len = s->read_end - s->read_start + 1;
    int genome_gap_len = s->genome_end - s->genome_start + 1;

    int read_gap_len_ex = read_gap_len_ex + sw_prepare->left_flank + sw_prepare->right_flank;
    int genome_gap_len_ex = genome_gap_len_ex + sw_prepare->left_flank + sw_prepare->right_flank;

    LOG_DEBUG_F("\tgap (read %lu-%lu, genome %lu-%lu) = (%i, %i): read %s\n", 
		s->read_start, s->read_end, s->genome_start, s->genome_end,
		read_gap_len, genome_gap_len, sw_prepare->read->id);
    LOG_DEBUG_F("\tflanks (left, right) = (%i, %i)\n", sw_prepare->left_flank, sw_prepare->right_flank);
    LOG_DEBUG_F("\tquery : %s\n", sw_prepare->query);
    LOG_DEBUG_F("\tref   : %s\n", sw_prepare->ref);
    LOG_DEBUG_F("\tmquery: %s (start %i)\n", output->query_map_p[i], output->query_start_p[i]);
    LOG_DEBUG_F("\tmref  : %s (start %i)\n", output->ref_map_p[i], output->ref_start_p[i]);

    cigar_code_t *cigar_c = generate_cigar_code(output->query_map_p[i], output->ref_map_p[i],
						strlen(output->query_map_p[i]), output->query_start_p[i],
						output->ref_start_p[i], read_gap_len, genome_gap_len,
						&distance, sw_prepare->ref_type);
    LOG_DEBUG_F("\tscore : %0.2f, cigar: %s (distance = %i)\n", 
		output->score_p[i], new_cigar_code_string(cigar_c), distance);

    /*
    if (output->query_start_p[i] > 0 && output->ref_start_p[i] > 0 && 
	output->query_start_p[i] != output->ref_start_p[i]) { 
      LOG_DEBUG("both map start points > 0 and are different lengths");
      exit(-1);
    }
    */
    //    assert(output->query_start_p[i] == 0);
    //    assert(output->ref_start_p[i] == 0);

    cigar_op = cigar_code_get_op(0, cigar_c);
    if (cigar_op) {
      if (cigar_op->name == 'H') {
	if (output->ref_start_p[i] == 0) { 
	  cigar_op->name = 'I';
	} else {
	  cigar_op->name = 'M';
	}
      } else if (cigar_op->name == '=') cigar_op->name = 'M';
    }

    cigar_op = cigar_code_get_last_op(cigar_c);
    if (cigar_op && cigar_op->name == 'H') cigar_op->name = 'I';

    LOG_DEBUG_F("gap_read_len = %i, cigar_code_length (%s) = %i\n", 
		read_gap_len, new_cigar_code_string(cigar_c), cigar_code_nt_length(cigar_c));
    assert(read_gap_len == cigar_code_nt_length(cigar_c));

    /*
    if (cigar_code_get_num_ops(cigar_c) > 2) {
      if (sw_prepare->left_flank > 0) {
	cigar_op = cigar_code_get_op(0, cigar_c);
	assert(cigar_op->number >= sw_prepare->left_flank && cigar_op->name == 'M');
	cigar_op->number -= sw_prepare->left_flank;
      }
      if (sw_prepare->right_flank > 0) {
	cigar_op = cigar_code_get_last_op(cigar_c);
	assert(cigar_op->number >= sw_prepare->right_flank && cigar_op->name == 'M');
	cigar_op->number -= sw_prepare->right_flank;
      }
      init_cigar_string(cigar_c);
      LOG_DEBUG_F("\tnew cigar: %s\n", new_cigar_code_string(cigar_c));
    } else {
      assert(cigar_code_get_num_ops(cigar_c) == 1);
      if (sw_prepare->right_flank > 0) {
	cigar_op = cigar_code_get_last_op(cigar_c);
	assert(cigar_op->number >= sw_prepare->right_flank && cigar_op->name == 'M');
	cigar_op->number -= (sw_prepare->left_flank + sw_prepare->right_flank);
	if (cigar_op->number > read_gap_len) {
	  cigar_code_append_op(cigar_op_new(cigar_op->number - read_gap_len, 'D'), cigar_c);
	} else if (cigar_op->number < read_gap_len) {
	  cigar_code_append_op(cigar_op_new(read_gap_len - cigar_op->number, 'I'), cigar_c);
	} else{
	  init_cigar_string(cigar_c);
	}
	//	LOG_DEBUG_F("\tnew cigar: %s\n", new_cigar_code_string(cigar_c));
      }
    }
    */
    // and now set the cigar for this gap
    s->info = (void *) cigar_c;

    // free
    sw_prepare_free(sw_prepare);
  }

  display_sr_lists("END of fill_gaps", mapping_batch);
    
  // free memory
  sw_multi_output_free(output);
  array_list_free(sw_prepare_list, (void *) NULL);
}

//------------------------------------------------------------------------------------

void fill_end_gaps(mapping_batch_t *mapping_batch, sw_optarg_t *sw_optarg, 
		   genome_t *genome, int min_H, int min_distance) {

  int sw_count = 0;

  fastq_read_t *fq_read;
  array_list_t *fq_batch = mapping_batch->fq_batch;

  size_t read_index, read_len;

  cal_t *cal;
  array_list_t *cal_list = NULL;
  size_t num_cals, num_targets = mapping_batch->num_targets;

  char *seq, *revcomp_seq = NULL;

  seed_region_t *s;

  cigar_op_t *cigar_op;
  cigar_code_t *cigar_code;

  size_t start, end;
  size_t gap_read_start, gap_read_end, gap_read_len;
  size_t gap_genome_start, gap_genome_end, gap_genome_len;

  int first, last, mode, distance, flank = 5;
  sw_prepare_t *sw_prepare;
  array_list_t *sw_prepare_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  
  char *ref, *query;

  // initialize query and reference sequences to Smith-Waterman
  for (size_t i = 0; i < num_targets; i++) {

    read_index = mapping_batch->targets[i];
    fq_read = (fastq_read_t *) array_list_get(read_index, fq_batch);

    cal_list = mapping_batch->mapping_lists[read_index];
    num_cals = array_list_size(cal_list);
    
    if (num_cals <= 0) continue;

    read_len = fq_read->length;
    revcomp_seq = NULL;

    // processing each CAL from this read
    for(size_t j = 0; j < num_cals; j++) {

      // get cal and read index
      cal = array_list_get(j, cal_list);
      if (cal->sr_list->size == 0) continue;

      sw_prepare = NULL;
      s = (seed_region_t *) linked_list_get_first(cal->sr_list);
      cigar_code = (cigar_code_t *) s->info;
      LOG_DEBUG_F("CAL #%i of %i (strand %i), sr_list size = %i, cigar = %s (distance = %i)\n", 
		  j, num_cals, cal->strand, cal->sr_list->size, new_cigar_code_string(cigar_code), cigar_code->distance);
      
      for (int k = 0; k < 2; k++) {
	mode = NONE_POS;
	if (k == 0) {
	  if ((cigar_op = cigar_code_get_op(0, cigar_code)) &&
	      cigar_op->name == 'H' && cigar_op->number > min_H) {
	    LOG_DEBUG_F("%i%c\n", cigar_op->number, cigar_op->name);

	    mode = BEGIN_POS;
	    gap_read_start = 0;
	    gap_read_end = cigar_op->number - 1;
	    gap_genome_start = s->genome_start;
	    gap_genome_end = gap_genome_start + cigar_op->number - 1;
	  }
	} else {
	  if ((cigar_op = cigar_code_get_last_op(cigar_code)) &&
	      cigar_op->name == 'H' && cigar_op->number > min_H) {
	    LOG_DEBUG_F("%i%c\n", cigar_op->number, cigar_op->name);

	    mode = END_POS;
	    gap_read_start = read_len - cigar_op->number;
	    gap_read_end = read_len - 1;
	    gap_genome_end = s->genome_end;
	    gap_genome_start = gap_genome_end - cigar_op->number + 1;
	  }
	}
	    
	if (mode == NONE_POS) continue;

	// get query sequence, revcomp if necessary
	if (cal->strand) {
	  if (revcomp_seq == NULL) {
	    revcomp_seq = strdup(fq_read->sequence);
	    seq_reverse_complementary(revcomp_seq, read_len);
	  }
	  seq = revcomp_seq;
	} else {
	  seq = fq_read->sequence;
	}

	gap_read_len = gap_read_end - gap_read_start + 1;
	/*	
	char *query = (char *) malloc((gap_read_len + 1) * sizeof(char));
	memcpy(query, seq, gap_len);
	query[gap_read_len] = '\0';
	*/

	// get ref. sequence
	start = gap_genome_start;// + 1;
	end = gap_genome_end;// + 1;
	gap_genome_len = end - start + 1;
	ref = (char *) malloc((gap_genome_len + 1) * sizeof(char));
	genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
					  &start, &end, genome);
	ref[gap_genome_len] = '\0';
	
	first = -1; 
	last = -1;
	distance = 0;
	for (int k = 0, k1 = gap_read_start; k < gap_read_len; k++, k1++) {
	  if (seq[k1] != ref[k]) {
	    distance++;
	    if (first == -1) first = k;
	    last = k;
	  }
	  //	  LOG_DEBUG_F("k = %i, k.read = %i: %c - %c : distance = %i, (first, last) = (%i, %i)\n", 
	  //		      k, k1, seq[k1], ref[k], distance, first, last);
	}

	if (distance < min_distance) {
	  cigar_op->name = 'M';
	  cigar_code->distance += distance;
	  free(ref);
	  continue;
	} else {
	  //	  LOG_DEBUG_F("query: %s\n", &seq[gap_read_start]);
	  //	  LOG_DEBUG_F("ref. : %s\n", ref);
	  LOG_FATAL_F("here we must run SW: distance = %i: first = %i, last = %i, gaps (read, genome) = (%i, %i)\n", 
		      distance, first, last, gap_read_len, gap_genome_len);
	}

	// we must run the SW algorithm
	

	//	sw_prepare = sw_prepare_new(0, 0, 0, 0);
	//	sw_prepare_sequences( cal, genome, sw_prepare);
	//	array_list_insert(sw_prepare, sw_prepare_list);
	//	sw_count++;
      }
    }
  }
  LOG_DEBUG_F("sw_count = %i\n", sw_count);


  // debugging....
  for (size_t i = 0; i < num_targets; i++) {
    read_index = mapping_batch->targets[i];
    fq_read = (fastq_read_t *) array_list_get(read_index, fq_batch);

    LOG_DEBUG_F("Read %s\n", fq_read->id);
    
    cal_list = mapping_batch->mapping_lists[read_index];
    num_cals = array_list_size(cal_list);
    
    if (num_cals <= 0) continue;

    for(size_t j = 0; j < num_cals; j++) {

      // get cal and read index
      cal = array_list_get(j, cal_list);
      if (cal->sr_list->size == 0) continue;

      sw_prepare = NULL;
      s = (seed_region_t *) linked_list_get_first(cal->sr_list);
      cigar_code = (cigar_code_t *) s->info;
      LOG_DEBUG_F("\tCAL #%i of %i (strand %i), sr_list size = %i, cigar = %s (distance = %i)\n", 
		  j, num_cals, cal->strand, cal->sr_list->size, new_cigar_code_string(cigar_code), cigar_code->distance);
    }
  }
}

//------------------------------------------------------------------------------------

void merge_seed_regions(mapping_batch_t *mapping_batch) {
  linked_list_item_t *list_item;
  cal_t *cal;
  seed_region_t *s, *s_first;
  cigar_code_t *cigar_code, *cigar_code_prev;
  cigar_op_t *cigar_op, *cigar_op_prev;
  int num_ops;
  int op;
  array_list_t *cals_list;
  fastq_read_t *fq_read;
  linked_list_iterator_t itr;

  register size_t num_targets = mapping_batch->num_targets;
  register size_t num_cals;
  register int i;
  register size_t t;

  for (t = 0; t < num_targets; t++) {
    cals_list = mapping_batch->mapping_lists[mapping_batch->targets[t]];
    num_cals = array_list_size(cals_list);
    fq_read = array_list_get(mapping_batch->targets[t], mapping_batch->fq_batch);

    //LOG_DEBUG_F("Read %s\n", fq_read->id);

    for (i = 0; i < num_cals; i++) {
      cal = array_list_get(i, cals_list);
      //LOG_DEBUG_F("\tCAL %i:\n", i);
      linked_list_iterator_init(cal->sr_list, &itr);

      s_first = linked_list_iterator_curr(&itr);      

      if (s_first) {
	cigar_code_prev = (cigar_code_t *)s_first->info;
	s = linked_list_iterator_next(&itr);
	while (s) {
	  //LOG_DEBUG_F("\t\tItem [%lu|%i - %i|%lu]: \n", s->genome_start, s->read_start, s->read_end, s->genome_end);
	  cigar_code = (cigar_code_t *)s->info;
	  if (cigar_code) { //TODO: delete
	    num_ops = array_list_size(cigar_code->ops);
	    for (op = 0, cigar_op = array_list_get(op, cigar_code->ops); 
		 op < num_ops;
		 op++, cigar_op = array_list_get(op, cigar_code->ops)) {
	      cigar_code_append_op(cigar_op, cigar_code_prev);	    
	    }
	    cigar_code_prev->distance += cigar_code->distance;
	  } 
	  
	  s_first->read_end = s->read_end;
	  s_first->genome_end = s->genome_end;
	  
	  linked_list_iterator_remove(&itr);	
	  s = linked_list_iterator_curr(&itr);
	}
	cal->info = (void *)cigar_code_prev;
      } else {
	cal->info = NULL;
	//LOG_DEBUG("\t\tLINKED LIST EMPTY");
      }
      /*LOG_DEBUG_F("\t\tItem [%lu|%i - %i|%lu]: Distance(%i) %s\n", s->genome_start, s->read_start,
	s->read_end, s->genome_end, cigar_code->distance, new_cigar_code_string(cigar_code));*/
    }
  }
}



int apply_caling_rna(cal_seeker_input_t* input, batch_t *batch) {

  LOG_DEBUG("========= APPLY CALING RNA =========\n");

  struct timeval start, end;
  double time;
  //if (time_on) { start_timer(start); }

  metaexons_t *metaexons = input->metaexons;
  bwt_optarg_t *bwt_optarg = input->bwt_optarg;
  bwt_index_t *bwt_index = input->index;
  cal_optarg_t *cal_optarg = input->cal_optarg;
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  array_list_t *allocate_cals;
  size_t num_cals, select_cals, total_cals = 0;
  size_t num_batches = 0, num_reads_unmapped = 0, num_without_cals = 0;
  size_t total_reads = 0;
  size_t num_targets, target_pos, total_targets, extra_target_pos;
  fastq_read_t *read;
  genome_t *genome = input->genome;
  unsigned int num_chromosomes = genome->num_chromosomes;
  size_t *targets_aux;
  int min_seeds, max_seeds;
  int seed_size = input->cal_optarg->seed_size;
  array_list_t *cal_list, *list;
  cal_t *cal;
  //array_list_t *region_list;
  region_t *bwt_region_back, *bwt_region_forw;
  linked_list_t *linked_list;
  seed_region_t *seed_region_start, *seed_region_end, *seed_region;
  int gap_nt, anchor_nt;
  bwt_anchor_t *bwt_anchor_back, *bwt_anchor_forw;

  num_targets = mapping_batch->num_targets;
  total_targets = 0;
  extra_target_pos = 0;
  total_reads += num_targets;
  target_pos = 0;

  mapping_batch->extra_stage_do = 1;

  /*  int t, target;
  for (t = 0; t < num_targets; t++) {
    target = mapping_batch->targets[t];
    mapping_batch->mapping_lists[target]->size = 0;
  }
  return RNA_POST_PAIR_STAGE;
  */

  array_list_t *region_list = array_list_new(1000, 
					     1.25f, 
					     COLLECTION_MODE_ASYNCHRONIZED);
  
  extern pthread_mutex_t mutex_sp;
  //extern size_t TOTAL_READS_SEEDING, TOTAL_READS_SEEDING2;

  //pthread_mutex_lock(&mutex_sp);
  //TOTAL_READS_SEEDING += num_targets;
  //pthread_mutex_unlock(&mutex_sp);

  for (size_t i = 0; i < num_targets; i++) {
    read = array_list_get(mapping_batch->targets[i], mapping_batch->fq_batch); 
    
    //printf("From CAL Seeker %s\n", read->id);
    list = mapping_batch->mapping_lists[mapping_batch->targets[i]];
    
    //if (array_list_get_flag(region_list) == 0 || 
    //	array_list_get_flag(region_list) == 2) {
    //We have normal and extend seeds (anchors)
    max_seeds = (read->length / 15)*2 + 10;      
    //printf("%i\n", input->cal_optarg->min_cal_size);
    num_cals = bwt_generate_cals(read->sequence, 
				 seed_size, 
				 bwt_optarg,
				 cal_optarg,
				 bwt_index, 
				 list, 
				 num_chromosomes);


    // if we want to seed with 24-length seeds,
    if (num_cals == 0) {
      //printf("No Cals seeding...\n");
      
      //pthread_mutex_lock(&mutex_sp);
      //extern size_t seeds_1err;
      //seeds_1err++;
      //pthread_mutex_unlock(&mutex_sp);

      int seed_size = 24;
      //First, Delete old regions
      array_list_clear(region_list, (void *)region_bwt_free);
      
      //Second, Create new regions with seed_size 24 and 1 Mismatch
      
      bwt_map_inexact_seeds_seq(read->sequence, seed_size, seed_size/2,
				bwt_optarg, bwt_index, 
				region_list);
      
      max_seeds = (read->length / 15)*2 + 10;
      //int prev_min_cal = input->cal_optarg->min_cal_size;
      //input->cal_optarg->min_cal_size = seed_size + seed_size / 2;
      //printf("NO CALS, new seeds %lu\n", array_list_size(region_list));

      num_cals = bwt_generate_cal_list_linked_list(region_list,
						   input->cal_optarg,
						   &min_seeds, &max_seeds,
						   genome->num_chromosomes + 1,
						   list, read->length,
						   cal_optarg->min_cal_size,
						   0);

      //input->cal_optarg->min_cal_size = prev_min_cal;

      //pthread_mutex_lock(&mutex_sp);
      //TOTAL_READS_SEEDING2++;
      //pthread_mutex_unlock(&mutex_sp);

    } 

    array_list_clear(region_list, (void *)region_bwt_free);

    //filter-incoherent CALs
    int founds[num_cals], found = 0;
    for (size_t j = 0; j < num_cals; j++) {
      founds[j] = 0;
      cal = array_list_get(j, list);
      LOG_DEBUG_F("\tcal %i of %i: sr_list size = %i (cal->num_seeds = %i) %i:%lu-%lu\n", 
		  j, num_cals, cal->sr_list->size, cal->num_seeds,
		  cal->chromosome_id, cal->start, cal->end);
      if (cal->sr_list->size > 0) {
	int start = 0;
	size_t genome_start = 0;
	int first = 1;
	for (linked_list_item_t *list_item = cal->sr_list->first; list_item != NULL; list_item = list_item->next) {
	  seed_region_t *s = list_item->item;
	  
	  LOG_DEBUG_F("\t\t:: star %lu > %lu s->read_start\n", start, s->read_start);
	  LOG_DEBUG_F("\t\t:: read_star %lu > read_end %lu \n", s->read_start, s->read_end);
	  if (start > s->read_start || s->read_start >= s->read_end) {
	    LOG_DEBUG("\t\t\t:: remove\n");
	    found++;
	    founds[j] = 1;
	  }

	  if (!first && 
	      ((s->genome_start < genome_start) || 
	      (s->genome_start - genome_start) > 2*read->length)) {
	    //printf("Remove (genome_start = %i s->genome_start = %i)\n", genome_start, s->genome_start);
	    //cal_print(cal);
	    found++;
	    founds[j] = 1;
	  }

	  first = 0;
	  start = s->read_end + 1;
	  genome_start = s->genome_end + 1;
	}
      } else {
	found++;
	founds[j] = 1;
      }
    }

    if (found) {
      min_seeds = 100000;
      max_seeds = 0;
      cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
      for (size_t j = 0; j < num_cals; j++) {
	if (!founds[j]) {
	  cal = array_list_get(j, list);
	  cal->num_seeds = cal->sr_list->size;
	  if (cal->num_seeds > max_seeds) max_seeds = cal->num_seeds;
	  if (cal->num_seeds < min_seeds) min_seeds = cal->num_seeds;
	  array_list_insert(cal, cal_list);
	  array_list_set(j, NULL, list);
	}
      }
      array_list_free(list, (void *) cal_free);
      num_cals = array_list_size(cal_list);
      list = cal_list;
    }

    mapping_batch->mapping_lists[mapping_batch->targets[i]] = list;
    num_cals = array_list_size(list);

    int max = 100;
    if (num_cals > max) {
      select_cals = num_cals - max;
      for(int j = num_cals - 1; j >= max; j--) {
	cal_free(array_list_remove_at(j, mapping_batch->mapping_lists[mapping_batch->targets[i]]));
      }
    }

    //mapping_batch->targets[target_pos++] = mapping_batch->targets[i];
    //} //else if (num_cals > 0) {

    mapping_batch->targets[target_pos++] = mapping_batch->targets[i];

    
    /* printf("<<<<<===== CAL SERVER =====>>>>>\n"); */
    /* for (int c = 0; c < array_list_size(mapping_batch->mapping_lists[mapping_batch->targets[i]]); c++) { */
    /*   cal_t *cal_aux = array_list_get(c, mapping_batch->mapping_lists[mapping_batch->targets[i]]); */
    /*   cal_print(cal_aux); */
    /* } */
    /* printf("<<<<<===== CAL SERVER END =====>>>>>\n"); */
       
  }

  mapping_batch->num_targets = target_pos;

  array_list_free(region_list, NULL);

  //if (time_on) { stop_timer(start, end, time); timing_add(time, CAL_SEEKER, timing); }

  LOG_DEBUG("========= APPLY CALING RNA END =========\n");

  //  return RNA_STAGE;
  if (batch->mapping_mode == RNA_MODE) {
    return RNA_STAGE;
  }

  if (batch->pair_input->pair_mng->pair_mode != SINGLE_END_MODE) {
    return PRE_PAIR_STAGE;
  } else if (batch->mapping_batch->num_targets > 0) {
    return SW_STAGE;
  }
  
  return DNA_POST_PAIR_STAGE;

}


    /*} else {
      //We have double anchors with smaller distance between they
      //printf("Easy case... Two anchors and same distance between read gap and genome distance\n");
      num_cals = 0;
      for (int a = array_list_size(region_list) - 1; a >= 0; a -= 2) {
	bwt_anchor_back = array_list_remove_at(a, region_list);
	bwt_anchor_forw = array_list_remove_at(a - 1, region_list);

	linked_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);

	
	//Seed for the first anchor
	anchor_nt = bwt_anchor_forw->end - bwt_anchor_forw->start;
	//printf("\t seed0[%i-%i][%lu-%lu]\n", 0, anchor_nt - 1,
	//     bwt_anchor_forw->start, bwt_anchor_forw->end);
	seed_region_start = seed_region_new(0, anchor_nt - 1,
					    bwt_anchor_forw->start, bwt_anchor_forw->end, 0);

	//Seed for the first anchor
	gap_nt = read->length - (anchor_nt + (bwt_anchor_back->end - bwt_anchor_back->start));
	//printf("\t gap_nt = %i, anchor_nt = %i\n", gap_nt, anchor_nt);
	//printf("\t seed1[%i-%i][%lu-%lu]\n", anchor_nt + gap_nt, read->length - 1, 
	//     bwt_anchor_back->start + 1, bwt_anchor_back->end);
	seed_region_end = seed_region_new(anchor_nt + gap_nt, read->length - 1,
				      bwt_anchor_back->start + 1, bwt_anchor_back->end, 1);

	//The reference distance is 0 and the read distance not
	//The read distance is 0 and the reference distance not
	//if (seed_region_start->genome_end > seed_region_end->genome_start || 
	//  seed_region_start->read_end > seed_region_end->read_start) { 
	//array_list_clear(region_list, NULL);
	//continue;
	if (seed_region_end->genome_start - seed_region_start->genome_end < 5 || 
	    seed_region_end->read_start - seed_region_start->read_end < 5) {
	  seed_region_start->genome_end -= 5;
	  seed_region_start->read_end -= 5;
	  seed_region_end->genome_start += 5;
	  seed_region_end->read_start += 5;
	}

	linked_list_insert(seed_region_start, linked_list);
	linked_list_insert_last(seed_region_end, linked_list);

	cal = cal_new(bwt_anchor_forw->chromosome + 1,
		      bwt_anchor_forw->strand,
		      bwt_anchor_forw->start,
		      bwt_anchor_back->end + 1,
		      2,
		      linked_list,
		      linked_list_new(COLLECTION_MODE_ASYNCHRONIZED));
	array_list_insert(cal, list);
	num_cals++;
      }
      }*/


    //array_list_free(mapping_batch->mapping_lists[mapping_batch->targets[i]], region_bwt_free);
    //mapping_batch->mapping_lists[mapping_batch->targets[i]] = allocate_cals;

    /*if (num_cals > MAX_RNA_CALS) {
      if (!mapping_batch->extra_stage_do) {
	mapping_batch->extra_targets[extra_target_pos] = mapping_batch->targets[i];
	mapping_batch->extra_stage_id[extra_target_pos++] = MAX_CALS;
	array_list_clear(mapping_batch->mapping_lists[mapping_batch->targets[i]], cal_free);	
      } else {
	select_cals = num_cals - MAX_RNA_CALS;
	for(size_t j = num_cals - 1; j >= MAX_RNA_CALS; j--) {
	  cal_free(array_list_remove_at(j, mapping_batch->mapping_lists[mapping_batch->targets[i]]));
	}
	mapping_batch->targets[target_pos++] = mapping_batch->targets[i];
      }
    }else if (num_cals > 0) {
      mapping_batch->targets[target_pos++] = mapping_batch->targets[i];
    }else {
      if (!mapping_batch->extra_stage_do) {
	mapping_batch->extra_targets[extra_target_pos] = mapping_batch->targets[i];
	mapping_batch->extra_stage_id[extra_target_pos++] = NO_CALS;
      }
      }*/
  //printf("APPLY CAL SEEKER DONE!\n");
/*
  // this code does not offer great improvements !!!
  // should we comment it ?
  if (mapping_batch->extra_stage_do) {
    //printf("Go to original targets & Fusion...\n");
    targets_aux = mapping_batch->targets;
    mapping_batch->targets = mapping_batch->extra_targets;
    mapping_batch->extra_targets = targets_aux;
    mapping_batch->num_targets = mapping_batch->num_extra_targets;
   
    for (size_t i = 0; i < target_pos; i++) {
      mapping_batch->targets[mapping_batch->num_targets++] = mapping_batch->extra_targets[i];
    }
    //mapping_batch->num_targets += target_pos;
    mapping_batch->num_extra_targets = 0;
    mapping_batch->extra_stage_do = 0;
    /*printf("Total Targets %i\n", mapping_batch->num_targets);
    printf("\t--->TARGETS: ", mapping_batch->extra_stage);
    for (int i = 0; i < mapping_batch->num_targets; i++){
      printf("%i,", mapping_batch->targets[i]);
    }
    printf("\n");
    
  }else if (extra_target_pos) {
    targets_aux = mapping_batch->targets;
    mapping_batch->targets = mapping_batch->extra_targets;
    mapping_batch->extra_targets = targets_aux;
    
    mapping_batch->num_targets = extra_target_pos;
    mapping_batch->num_extra_targets = target_pos;
    //printf("Go back to stage Seeding. Extra Targets = %i, Targets = %i\n", mapping_batch->num_extra_targets, mapping_batch->num_targets );
    mapping_batch->extra_stage_do = 1;

    return SEEDING_STAGE;
  } 

  //return RNA_PREPROCESS_STAGE;

}
*/

//====================================================================================
// apply_caling
//====================================================================================
int apply_caling(cal_seeker_input_t* input, batch_t *batch) {
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  array_list_t *list = NULL;
  size_t read_index, num_cals;
  int min_seeds, max_seeds;
  int min_limit;

  cal_t *cal;
  array_list_t *cal_list;

  fastq_read_t *read;

  bwt_optarg_t *bwt_optarg = input->bwt_optarg;
  bwt_index_t *bwt_index = input->index;
  size_t num_chromosomes = input->genome->num_chromosomes + 1;
  size_t num_targets = mapping_batch->num_targets;
  size_t *targets = mapping_batch->targets;
  size_t new_num_targets = 0;
  array_list_t *region_list;
  bwt_anchor_t *bwt_anchor_back, *bwt_anchor_forw;
  linked_list_t *linked_list;
  int anchor_nt, gap_nt;
  seed_region_t *seed_region_start, *seed_region_end;
  //max_seeds = input->cal_optarg->num_seeds;
  
  //  size_t *new_targets = (size_t *) calloc(num_targets, sizeof(size_t));
  
  // set to zero
  mapping_batch->num_to_do = 0;

  for (size_t i = 0; i < num_targets; i++) {

    read_index = targets[i];
    read = array_list_get(read_index, mapping_batch->fq_batch); 
    region_list = mapping_batch->mapping_lists[read_index];
    // for debugging
    //    LOG_DEBUG_F("%s\n", ((fastq_read_t *) array_list_get(read_index, mapping_batch->fq_batch))->id);
    
    if (!list) {
      list = array_list_new(1000, 
			    1.25f, 
			    COLLECTION_MODE_ASYNCHRONIZED);
    }


    if (array_list_get_flag(region_list) == 0 || 
	array_list_get_flag(region_list) == 2) {
      //We have normal and extend seeds (anchors)
      max_seeds = (read->length / 15)*2 + 10;
      num_cals = bwt_generate_cal_list_linked_list(region_list,
						   input->cal_optarg,
						   &min_seeds, &max_seeds,
						   num_chromosomes,
						   list, read->length,
						   input->cal_optarg->min_cal_size, 0);
    } else {
      //We have double anchors with smaller distance between they
      //printf("Easy case... Two anchors and same distance between read gap and genome distance\n");
      num_cals = 0;
      for (int a = array_list_size(region_list) - 1; a >= 0; a -= 2) {
	max_seeds = 2;
	min_seeds = 2;
	bwt_anchor_back = array_list_remove_at(a, region_list);
	bwt_anchor_forw = array_list_remove_at(a - 1, region_list);

	linked_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);

	
	//Seed for the first anchor
	anchor_nt = bwt_anchor_forw->end - bwt_anchor_forw->start;
	//printf("\t seed0[%i-%i][%lu-%lu]\n", 0, anchor_nt - 1,
	//     bwt_anchor_forw->start, bwt_anchor_forw->end);
	seed_region_start = seed_region_new(0, anchor_nt - 1,
					    bwt_anchor_forw->start, bwt_anchor_forw->end, 0, 0, 0);

	//Seed for the first anchor
	gap_nt = read->length - (anchor_nt + (bwt_anchor_back->end - bwt_anchor_back->start));
	//printf("\t gap_nt = %i, anchor_nt = %i\n", gap_nt, anchor_nt);
	//printf("\t seed1[%i-%i][%lu-%lu]\n", anchor_nt + gap_nt, read->length - 1, 
	//     bwt_anchor_back->start + 1, bwt_anchor_back->end);
	seed_region_end = seed_region_new(anchor_nt + gap_nt, read->length - 1,
					  bwt_anchor_back->start + 1, bwt_anchor_back->end, 1, 0, 0);

	//The reference distance is 0 and the read distance not
	//The read distance is 0 and the reference distance not
	//if (seed_region_start->genome_end > seed_region_end->genome_start || 
	//  seed_region_start->read_end > seed_region_end->read_start) { 
	//array_list_clear(region_list, NULL);
	//continue;
	if (seed_region_end->genome_start - seed_region_start->genome_end < 5 || 
	    seed_region_end->read_start - seed_region_start->read_end < 5) {
	  seed_region_start->genome_end -= 5;
	  seed_region_start->read_end -= 5;
	  seed_region_end->genome_start += 5;
	  seed_region_end->read_start += 5;
	}

	linked_list_insert(seed_region_start, linked_list);
	linked_list_insert_last(seed_region_end, linked_list);

	cal = cal_new(bwt_anchor_forw->chromosome + 1,
		      bwt_anchor_forw->strand,
		      bwt_anchor_forw->start,
		      bwt_anchor_back->end + 1,
		      2,
		      linked_list,
		      linked_list_new(COLLECTION_MODE_ASYNCHRONIZED));
	array_list_insert(cal, list);
	num_cals++;
      }
    }

    // for debugging
    LOG_DEBUG_F("read %s : num. cals = %i, min. seeds = %i, max. seeds = %i\n", 
		read->id, num_cals, min_seeds, max_seeds);


    /*    if (num_cals == 0) {
      int seed_size = 24;
      //First, Delete old regions
      array_list_clear(mapping_batch->mapping_lists[read_index], region_bwt_free);
      //Second, Create new regions with seed_size 24 and 1 Mismatch
      bwt_map_inexact_seeds_seq(read->sequence, seed_size, seed_size/2,
				bwt_optarg, bwt_index, 
				mapping_batch->mapping_lists[read_index]);

      num_cals = bwt_generate_cal_list_linked_list(mapping_batch->mapping_lists[mapping_batch->targets[i]], 
						   input->cal_optarg,
						   &min_seeds, &max_seeds,
						   num_chromosomes,
						   list, read->length);
						   }*/

    /*
    for (size_t j = 0; j < num_cals; j++) {
      cal = array_list_get(j, list);
      LOG_DEBUG_F("\tchr: %i, strand: %i, start: %lu, end: %lu, num_seeds = %i, num. regions = %lu\n", 
		  cal->chromosome_id, cal->strand, cal->start, cal->end, cal->num_seeds, cal->sr_list->size);
    }
    */
    //    printf("min_seeds = %i, max_seeds = %i, min_limit = %i, num_cals = %i\n", 
    //	   min_seeds, max_seeds, min_limit, array_list_size(list));

    // filter incoherent CALs
    int founds[num_cals], found = 0;
    for (size_t j = 0; j < num_cals; j++) {
      founds[j] = 0;
      cal = array_list_get(j, list);
      LOG_DEBUG_F("\tcal %i of %i: sr_list size = %i (cal->num_seeds = %i) %i:%lu-%lu\n", 
		  j, num_cals, cal->sr_list->size, cal->num_seeds,
		  cal->chromosome_id, cal->start, cal->end);
      if (cal->sr_list->size > 0) {
	int start = 0;
	for (linked_list_item_t *list_item = cal->sr_list->first; list_item != NULL; list_item = list_item->next) {
	  seed_region_t *s = list_item->item;
	  
	  LOG_DEBUG_F("\t\t:: star %lu > %lu s->read_start\n", start, s->read_start);
	  if (start > s->read_start) {
	    LOG_DEBUG("\t\t\t:: remove\n");
	    found++;
	    founds[j] = 1;
	  }
	  start = s->read_end + 1;
	}
      } else {
	found++;
	founds[j] = 1;
      }
    }
    if (found) {
      min_seeds = 100000;
      max_seeds = 0;
      cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
      for (size_t j = 0; j < num_cals; j++) {
	if (!founds[j]) {
	  cal = array_list_get(j, list);
	  cal->num_seeds = cal->sr_list->size;
	  if (cal->num_seeds > max_seeds) max_seeds = cal->num_seeds;
	  if (cal->num_seeds < min_seeds) min_seeds = cal->num_seeds;
	  array_list_insert(cal, cal_list);
	  array_list_set(j, NULL, list);
	}
      }
      array_list_free(list, (void *) cal_free);
      num_cals = array_list_size(cal_list);
      list = cal_list;
    }
  
    //    LOG_FATAL_F("num. cals = %i, min. seeds = %i, max. seeds = %i\n", num_cals, min_seeds, max_seeds);
    // filter CALs by the number of seeds

    cal_list = list;
    list = NULL;
    /*
    int min_limit = input->cal_optarg->min_num_seeds_in_cal;

    if (min_limit < 0) min_limit = max_seeds;
    //    min_limit -= 3;
    
    if (min_seeds == max_seeds || min_limit <= min_seeds) {
      cal_list = list;
      list = NULL;
    } else {
      cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
      for (size_t j = 0; j < num_cals; j++) {
	cal = array_list_get(j, list);
	if (cal->num_seeds >= min_limit) {
	  array_list_insert(cal, cal_list);
	  array_list_set(j, NULL, list);
	}
      }
      array_list_clear(list, (void *) cal_free);
      num_cals = array_list_size(cal_list);
    }
    */
    if (num_cals > MAX_CALS) {
      for (size_t j = num_cals - 1; j >= MAX_CALS; j--) {
	cal = (cal_t *) array_list_remove_at(j, cal_list);
	cal_free(cal);
      }
      num_cals = array_list_size(cal_list);
    }
    
    //    LOG_DEBUG_F("num. cals = %i, MAX_CALS = %i\n", num_cals, MAX_CALS);

    if (num_cals > 0 && num_cals <= MAX_CALS) {
      array_list_set_flag(2, cal_list);
      targets[new_num_targets++] = read_index;

      /*
      int count1 = 0, count2 = 0;
      // count number of sw to do

      // method #1
      //      printf("method #1\n");
      seed_region_t *s, *prev_s;
      linked_list_iterator_t* itr;
      for (size_t j = 0; j < num_cals; j++) {
	prev_s = NULL;
	cal = array_list_get(j, cal_list);
	itr = linked_list_iterator_new(cal->sr_list);
	s = (seed_region_t *) linked_list_iterator_curr(itr);
	while (s != NULL) {
	  if ((prev_s == NULL && s->read_start != 0) || (prev_s != NULL)) {
	    //	    printf("\t\t\tcase 1\n");
	    count1++;
	  }
	  prev_s = s;
	  linked_list_iterator_next(itr);
	  s = linked_list_iterator_curr(itr);
	}
	if (prev_s != NULL && prev_s->read_end < read->length - 1) { 
	  count1++;
	  //	  printf("\t\t\tcase 2 (%i < %i)\n", prev_s->read_end, read->length - 1);
	}
	linked_list_iterator_free(itr);
      }

      // method #2
      printf("method #2\n");
      for (size_t j = 0; j < num_cals; j++) {
	cal = array_list_get(j, cal_list);
	printf("\t: %i\n", j);
	if (cal->sr_list->size > 0) {
	  int start = 0;
	  for (linked_list_item_t *list_item = cal->sr_list->first; list_item != NULL; list_item = list_item->next) {
	    seed_region_t *s = list_item->item;
	    printf("\t\t[%i|%i - %i|%i]\n", s->genome_start, s->read_start, s->read_end, s->genome_end);
	    if (s->read_start != start) {
	      count2++;
	    }
	    start = s->read_end + 1;
	  }
	  if (start < read->length) { 
	    count2++;
	  }
	}
      }
      printf("count #1 = %i, count #2 = %i\n", count1, count2);
      assert(count1 == count2);

      mapping_batch->num_to_do += count1;
*/

      // we have to free the region list
      array_list_free(mapping_batch->mapping_lists[read_index], (void *) region_bwt_free);
      mapping_batch->mapping_lists[read_index] = cal_list;
    } else {
      array_list_set_flag(0, mapping_batch->mapping_lists[read_index]);
      // we have to free the region list
      array_list_clear(mapping_batch->mapping_lists[read_index], (void *) region_bwt_free);
      if (cal_list) array_list_free(cal_list, (void *) cal_free);
      if (list) array_list_clear(list, (void *) cal_free);
    }

    /*    
    cal_list = list;
    list = NULL;
    array_list_set_flag(2, cal_list);
    //    mapping_batch->num_to_do += num_cals;
    targets[new_num_targets++] = read_index;
    
    // we have to free the region list
    array_list_free(mapping_batch->mapping_lists[read_index], (void *) region_bwt_free);
    mapping_batch->mapping_lists[read_index] = cal_list;
    */
    /*
    // filter CALs by the number of seeds
    int min_limit = input->cal_optarg->min_num_seeds_in_cal;
    if (min_limit < 0) min_limit = max_seeds;

    printf("min_seeds = %i, max_seeds = %i, min_limit = %i, num_cals = %i\n", 
	   min_seeds, max_seeds, min_limit, array_list_size(list));
    
    if (min_seeds == max_seeds || min_limit <= min_seeds) {
      cal_list = list;
      list = NULL;
    } else {
      cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
      for (size_t j = 0; j < num_cals; j++) {
	cal = array_list_get(j, list);
	if (cal->num_seeds >= min_limit) {
	  array_list_insert(cal, cal_list);
	  array_list_set(j, NULL, list);
	}
      }
      array_list_clear(list, (void *) cal_free);
      num_cals = array_list_size(cal_list);
      printf("************, num_cals = %i\n", num_cals);
    }

    if (num_cals > MAX_CALS) {
      for (size_t j = num_cals - 1; j >= MAX_CALS; j--) {
	cal = (cal_t *) array_list_remove_at(j, cal_list);
	cal_free(cal);
      }
      num_cals = array_list_size(cal_list);
    }

    if (num_cals > 0 && num_cals <= MAX_CALS) {
      array_list_set_flag(2, cal_list);
      mapping_batch->num_to_do += num_cals;
      targets[new_num_targets++] = read_index;
      
      // we have to free the region list
      array_list_free(mapping_batch->mapping_lists[read_index], (void *) region_bwt_free);
      mapping_batch->mapping_lists[read_index] = cal_list;
    } else {
      array_list_set_flag(0, mapping_batch->mapping_lists[read_index]);
      // we have to free the region list
      array_list_clear(mapping_batch->mapping_lists[read_index], (void *) region_bwt_free);
      if (cal_list) array_list_free(cal_list, (void *) cal_free);
      if (list) array_list_clear(list, (void *) cal_free);
    }
    */
  } // end for 0 ... num_targets

  // update batch
  mapping_batch->num_targets = new_num_targets;

  //  LOG_DEBUG_F("num. SW to do: %i\n", 	mapping_batch->num_to_do);

  //  exit(-1);

  // free memory
  if (list) array_list_free(list, NULL);

  if (batch->mapping_mode == RNA_MODE) {
    return RNA_STAGE;
  }

  if (batch->pair_input->pair_mng->pair_mode != SINGLE_END_MODE) {
    return PRE_PAIR_STAGE;
  } else if (batch->mapping_batch->num_targets > 0) {
    return SW_STAGE;
  }
  
  return DNA_POST_PAIR_STAGE;
}

//====================================================================================
// apply_caling bs
//====================================================================================

array_list_t *filter_cals(size_t num_cals, size_t read_length, array_list_t *list) {
  cal_t *cal;
  int min_seeds, max_seeds;
  array_list_t *cal_list;
  size_t select_cals;

  //filter-incoherent CALs
  int founds[num_cals], found = 0;
  for (size_t j = 0; j < num_cals; j++) {
    founds[j] = 0;
    cal = array_list_get(j, list);
    LOG_DEBUG_F("\tcal %i of %i: sr_list size = %i (cal->num_seeds = %i) %i:%lu-%lu\n", 
		j, num_cals, cal->sr_list->size, cal->num_seeds,
		cal->chromosome_id, cal->start, cal->end);
    if (cal->sr_list->size > 0) {
      int start = 0;
      size_t genome_start = 0;
      int first = 1;
      for (linked_list_item_t *list_item = cal->sr_list->first; list_item != NULL; list_item = list_item->next) {
	seed_region_t *s = list_item->item;
	
	LOG_DEBUG_F("\t\t:: star %lu > %lu s->read_start\n", start, s->read_start);
	LOG_DEBUG_F("\t\t:: read_star %lu > read_end %lu \n", s->read_start, s->read_end);
	if (start > s->read_start || s->read_start >= s->read_end) {
	  LOG_DEBUG("\t\t\t:: remove\n");
	  found++;
	  founds[j] = 1;
	}
	
	if (!first && 
	    ((s->genome_start < genome_start) || 
	     (s->genome_start - genome_start) > 2 * read_length)) {
	  //printf("Remove (genome_start = %i s->genome_start = %i)\n", genome_start, s->genome_start);
	  //cal_print(cal);
	  found++;
	  founds[j] = 1;
	}
	
	first = 0;
	start = s->read_end + 1;
	genome_start = s->genome_end + 1;
      }
    } else {
      found++;
      founds[j] = 1;
    }
  }
  
  if (found) {
    min_seeds = 100000;
    max_seeds = 0;
    cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    for (size_t j = 0; j < num_cals; j++) {
      if (!founds[j]) {
	cal = array_list_get(j, list);
	cal->num_seeds = cal->sr_list->size;
	if (cal->num_seeds > max_seeds) max_seeds = cal->num_seeds;
	if (cal->num_seeds < min_seeds) min_seeds = cal->num_seeds;
	array_list_insert(cal, cal_list);
	array_list_set(j, NULL, list);
      }
    }
    array_list_free(list, (void *) cal_free);
    num_cals = array_list_size(cal_list);
    list = cal_list;
  }
  
  num_cals = array_list_size(list);
  
  int max = 100;
  if (num_cals > max) {
    select_cals = num_cals - max;
    for(int j = num_cals - 1; j >= max; j--) {
      cal_free(array_list_remove_at(j, list));
    }
  }
 
  return list;
}

//------------------------------------------------------------------------------------

int apply_caling_bs(cal_seeker_input_t* input, batch_t *batch) {

  LOG_DEBUG("========= APPLY CALING BS =========\n");
  
  struct timeval start, end;
  double time;
  //if (time_on) { start_timer(start); }

  metaexons_t *metaexons = input->metaexons;
  bwt_optarg_t *bwt_optarg = input->bwt_optarg;
  bwt_index_t *bwt_index = input->index;
  bwt_index_t *bwt_index2 = input->index2;
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  array_list_t *allocate_cals;
  size_t num_cals, total_cals = 0;
  size_t num_batches = 0, num_reads_unmapped = 0, num_without_cals = 0;
  size_t max_seeds, total_reads = 0;
  size_t num_targets, target_pos = 0, target_pos2 = 0;
  fastq_read_t *read, *read2;
  genome_t *genome = input->genome;
  size_t *targets_aux, target_index;
  int seed_size = input->cal_optarg->seed_size;
  array_list_t *list;
  region_t *bwt_region_back, *bwt_region_forw;
  linked_list_t *linked_list;
  seed_region_t *seed_region_start, *seed_region_end, *seed_region;
  int gap_nt, anchor_nt;
  bwt_anchor_t *bwt_anchor_back, *bwt_anchor_forw;

  size_t *targets = mapping_batch->targets;
  size_t *targets2 = mapping_batch->targets2;

  num_targets = mapping_batch->num_targets;
  total_reads += num_targets;

  mapping_batch->extra_stage_do = 1;

  //extern pthread_mutex_t mutex_sp;
  //extern size_t TOTAL_READS_SEEDING, TOTAL_READS_SEEDING2;

  //pthread_mutex_lock(&mutex_sp);
  //TOTAL_READS_SEEDING += num_targets;
  //pthread_mutex_unlock(&mutex_sp);

  for (size_t i = 0; i < num_targets; i++) {
    target_index = targets[i];
    
    // GA reads
    LOG_DEBUG("searching CALs for GA reads\n");

    read = array_list_get(target_index, mapping_batch->GA_rev_fq_batch);
    read2 = array_list_get(target_index, mapping_batch->GA_fq_batch);
   
    //printf("From CAL Seeker %s\n", read->id);
    list = mapping_batch->mapping_lists[target_index];

    //    if (array_list_get_flag(list) == DOUBLE_ANCHORS) {
    //      printf("******************************* double anchors\n");
    //    } else {
    //printf("---------------------------> SEEDING 1\nindex1 %s\tindex2 %s\n", bwt_index->nucleotides, bwt_index2->nucleotides);

      max_seeds = (read->length / 15) * 2 + 10;      
      num_cals = bwt_generate_cals_bs(read->sequence, read2->sequence, seed_size, 
				      bwt_optarg, bwt_index2, bwt_index, list);
      //    }
    
    for (int c = 0; c < array_list_size(list); c++) {
      cal_print(array_list_get(c, list));
    }

    // filter CALs
    list = filter_cals(num_cals, read->length, list);

    for (int c = 0; c < array_list_size(list); c++) {
      cal_print(array_list_get(c, list));
    }


    // and update targets
    mapping_batch->mapping_lists[target_index] = list;
    if (array_list_size(list) > 0) {
      mapping_batch->targets[target_pos++] = target_index;
    }


    // CT reads
    LOG_DEBUG("searching CALs for CT reads\n");

    read = array_list_get(target_index, mapping_batch->CT_rev_fq_batch);
    read2 = array_list_get(target_index, mapping_batch->CT_fq_batch);
   
    //printf("From CAL Seeker %s\n", read->id);
    list = mapping_batch->mapping_lists2[target_index];

    printf("---------------------------> SEEDING 2\n");
    
    //    if (array_list_get_flag(list) == DOUBLE_ANCHORS) {
    //    } else {
      max_seeds = (read->length / 15) * 2 + 10;      
      num_cals = bwt_generate_cals_bs(read->sequence, read2->sequence, seed_size, 
				      bwt_optarg, bwt_index, bwt_index2, list);
      //    }
    
    for (int c = 0; c < array_list_size(list); c++) {
      cal_print(array_list_get(c, list));
    }

    // filter CALs
    list = filter_cals(num_cals, read->length, list);

    for (int c = 0; c < array_list_size(list); c++) {
      cal_print(array_list_get(c, list));
    }

    // and update targets
    mapping_batch->mapping_lists2[target_index] = list;
    if (array_list_size(list) > 0) {
      mapping_batch->targets2[target_pos2++] = target_index;
    }
  } // end of main loop (targets)

  // updating number of targets for the next stage
  mapping_batch->num_targets = target_pos;
  mapping_batch->num_targets2 = target_pos2;

  //if (time_on) { stop_timer(start, end, time); timing_add(time, CAL_SEEKER, timing); }

  LOG_DEBUG("========= END OF APPLYING CALING BS =========\n");

  if (batch->pair_input->pair_mng->pair_mode != SINGLE_END_MODE) {
    return BS_PRE_PAIR_STAGE;
  } else if ( (batch->mapping_batch->num_targets  > 0) ||
	      (batch->mapping_batch->num_targets2 > 0)   ) {
    return BS_SW_STAGE;
  }
  
  return BS_POST_PAIR_STAGE;
}

//------------------------------------------------------------------------------------


/* int apply_caling_bs_OOOOLLDDDDD(cal_seeker_input_t* input, batch_t *batch) { */

/*   //printf("APPLY CALLING BS...\n"); */

/*   mapping_batch_t *mapping_batch = batch->mapping_batch; */
/*   array_list_t *list = NULL; */
/*   size_t read_index, num_cals, min_seeds, max_seeds; */
/*   int min_limit; */

/*   cal_t *cal; */
/*   array_list_t *cal_list; */

/*   size_t num_chromosomes = input->genome->num_chromosomes + 1; */

/*   size_t num_targets = mapping_batch->num_targets; */
/*   size_t *targets = mapping_batch->targets; */
/*   size_t new_num_targets = 0; */
  
/*   // new variables for bisulfite */
/*   size_t num_targets2 = mapping_batch->num_targets2; */
/*   size_t *targets2 = mapping_batch->targets2; */
/*   size_t new_num_targets2 = 0; */
/*   array_list_t *list2 = NULL; */
/*   size_t num_cals2; */
/*   array_list_t *cal_list2; */
/*   size_t read_index2; */

/*   // set to zero */
/*   mapping_batch->num_to_do = 0; */
/*   mapping_batch->num_to_do2 = 0; */


/*   //////////////////////////////// */
/*   /\* */
/*   size_t reads_mapp  = 0; */
/*   size_t reads_mapp2 = 0; */
/*   size_t reads_no_mapp  = 0; */
/*   size_t reads_no_mapp2 = 0; */
/*   size_t reads_discard  = 0; */
/*   size_t reads_discard2 = 0; */

/*   size_t reads_cals = 0; */
/*   size_t reads_cals_acum = 0; */
/*   *\/ */
/*   //////////////////////////////// */


/*   //printf("targets 1 %lu\ntargets 2 %lu\n", num_targets, num_targets2); */

/*   //printf("----->primera lista\n"); */
/*   for (size_t i = 0; i < num_targets; i++) { */

/*     read_index = targets[i]; */

/*     // for debugging */
/*     //    LOG_DEBUG_F("%s\n", ((fastq_read_t *) array_list_get(read_index, mapping_batch->fq_batch))->id); */
    
/*     if (!list) { */
/*       list = array_list_new(1000,  */
/* 			    1.25f,  */
/* 			    COLLECTION_MODE_ASYNCHRONIZED); */
/*     } */

/*     // optimized version */
/*     num_cals = bwt_generate_cal_list_linkedlist(mapping_batch->mapping_lists[read_index],  */
/* 						input->cal_optarg, */
/* 						&min_seeds, &max_seeds, */
/* 						num_chromosomes, */
/* 						list); */
/*     //printf("read %lu\tcals1 = %lu\n", read_index, num_cals); */

/*     /\* */
/*     // for debugging */
/*     LOG_DEBUG_F("num. cals = %i, min. seeds = %i, max. seeds = %i\n", num_cals, min_seeds, max_seeds); */
/*     for (size_t j = 0; j < num_cals; j++) { */
/*       cal = array_list_get(j, list); */
/*       LOG_DEBUG_F("\tchr: %i, strand: %i, start: %lu, end: %lu, num. seeds = %lu, flank: (start, end) = (%lu, %lu)\n",  */
/* 		  cal->chromosome_id, cal->strand, cal->start, cal->end, cal->num_seeds, cal->flank_start, cal->flank_end); */
/*     } */
/*     *\/ */

/*     // filter CALs by the number of seeds */
/*     int min_limit = input->cal_optarg->min_num_seeds_in_cal; */

/*     if (min_limit < 0) min_limit = max_seeds; */

/*     if (min_seeds == max_seeds || min_limit <= min_seeds) { */
/*       cal_list = list; */
/*       list = NULL; */

/*     } else { */
/*       cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED); */

/*       ///////////////// */
/*       /\* */
/*       if (num_cals > 0) */
/* 	reads_cals++; */
/*       *\/ */
/*       ///////////////// */

/*       for (size_t j = 0; j < num_cals; j++) { */
/* 	cal = array_list_get(j, list); */
/* 	//filter cals with few seeds */
/* 	if (cal->num_seeds >= min_limit) { */

/* 	  ///////////////// */
/* 	  //reads_cals_acum++; */
/* 	  ///////////////// */


/* 	  //	  LOG_DEBUG_F("\t\tchr: %i, strand: %i, start: %lu, end: %lu, num. seeds = %lu (min. limit %i seeds), flank: (start, end) = (%lu, %lu)\n",  */
/* 	  //		      cal->chromosome_id, cal->strand, cal->start, cal->end, cal->num_seeds, min_limit, cal->flank_start, cal->flank_end); */
	  
/* 	  array_list_insert(cal, cal_list); */
/* 	  array_list_set(j, NULL, list); */
/* 	} */
/*       } */
/*       array_list_clear(list, (void *) cal_free); */
/*       num_cals = array_list_size(cal_list); */
/*     } */

/*     if (num_cals > MAX_CALS) { */
/*       for (size_t j = num_cals - 1; j >= MAX_CALS; j--) { */
/* 	cal = (cal_t *) array_list_remove_at(j, cal_list); */
/* 	cal_free(cal); */
/*       } */
/*       num_cals = array_list_size(cal_list); */
/*     } */

/*     if (num_cals > 0 && num_cals <= MAX_CALS) {   */
/*       ///////////// */
/*       //reads_mapp++; */
/*       ///////////// */

/*       array_list_set_flag(2, cal_list); */
/*       mapping_batch->num_to_do += num_cals; */
/*       targets[new_num_targets++] = read_index; */
      
/*       // we have to free the region list */
/*       array_list_free(mapping_batch->mapping_lists[read_index], (void *) region_bwt_free); */
/*       mapping_batch->mapping_lists[read_index] = cal_list; */
/*     } else { */
/*       ///////////// */
/*       //reads_discard++; */
/*       ///////////// */
/*       array_list_set_flag(0, mapping_batch->mapping_lists[read_index]); */
/*       // we have to free the region list */
/*       array_list_clear(mapping_batch->mapping_lists[read_index], (void *) region_bwt_free); */
/*       if (cal_list) array_list_free(cal_list, (void *) cal_free); */
/*       if (list) array_list_clear(list, (void *) cal_free); */
/*     } */

/*     //printf("read %lu\tcals  = %lu\n", read_index, num_cals); */
/*   } // end for 0 ... num_seqs */


/*   //printf("----->segunda lista\n"); */
/*   // add for bisulfite */
/*   for (size_t i = 0; i < num_targets2; i++) { */

/*     read_index2 = targets2[i]; */
/*     //printf("\nread %lu\n", read_index2); */

/*     // for debugging */
/*     //    LOG_DEBUG_F("%s\n", ((fastq_read_t *) array_list_get(read_index, mapping_batch->fq_batch))->id); */
    
/*     if (!list2) { */
/*       list2 = array_list_new(1000,  */
/* 			     1.25f,  */
/* 			     COLLECTION_MODE_ASYNCHRONIZED); */
/*     } */

/*     // optimized version */
/*     num_cals2 = bwt_generate_cal_list_linkedlist(mapping_batch->mapping_lists2[read_index2],  */
/* 						 input->cal_optarg, */
/* 						 &min_seeds, &max_seeds, */
/* 						 num_chromosomes, */
/* 						 list2); */
/*     //printf("read %lu\tcals2 = %lu\n", read_index2, num_cals2); */

/*     /\* */
/*     // for debugging */
/*     LOG_DEBUG_F("num. cals = %i, min. seeds = %i, max. seeds = %i\n", num_cals2, min_seeds, max_seeds); */
/*     for (size_t j = 0; j < num_cals2; j++) { */
/*       cal = array_list_get(j, list); */
/*       LOG_DEBUG_F("\tchr: %i, strand: %i, start: %lu, end: %lu, num. seeds = %lu, flank: (start, end) = (%lu, %lu)\n",  */
/* 		  cal->chromosome_id, cal->strand, cal->start, cal->end, cal->num_seeds, cal->flank_start, cal->flank_end); */
/*     } */
/*     *\/ */

/*     // filter CALs by the number of seeds */
/*     int min_limit = input->cal_optarg->min_num_seeds_in_cal; */

/*     if (min_limit < 0) min_limit = max_seeds; */

/*     if (min_seeds == max_seeds || min_limit <= min_seeds) { */
/*       //printf("read %lu\tborrar listas\n", read_index2); */
/*       cal_list2 = list2; */
/*       list2 = NULL; */
/*       //printf("read %lu\tborrar listas\n", read_index2); */

/*     } else { */
/*       //printf("read %lu\t1recortar lista\n", read_index2); */
/*       cal_list2 = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED); */

/*       ///////////////// */
/*       /\* */
/*       if (num_cals2 > 0) */
/* 	reads_cals++; */
/*       *\/ */
/*       ///////////////// */

/*       for (size_t j = 0; j < num_cals2; j++) { */
/* 	cal = array_list_get(j, list2); */
/* 	//printf("read %lu\tnum_seeds %lu\tmin_limit %i\n", read_index2, cal->num_seeds, min_limit); */
/* 	//filter cals with few seeds */
/* 	if (cal->num_seeds >= min_limit) { */

/* 	  ///////////////// */
/* 	  //reads_cals_acum++; */
/* 	  ///////////////// */

/* 	  //printf("keep cal %lu, with %lu seed (of %lu needed)\n", j, cal->num_seeds, min_limit); */
/* 	  //	  LOG_DEBUG_F("\t\tchr: %i, strand: %i, start: %lu, end: %lu, num. seeds = %lu (min. limit %i seeds), flank: (start, end) = (%lu, %lu)\n",  */
/* 	  //		      cal->chromosome_id, cal->strand, cal->start, cal->end, cal->num_seeds, min_limit, cal->flank_start, cal->flank_end); */
	  
/* 	  //printf("pre  callist2 size %lu\n", array_list_size(cal_list2)); */
/* 	  array_list_insert(cal, cal_list2); */
/* 	  //printf("post callist2 size %lu\n", array_list_size(cal_list2)); */
/* 	  array_list_set(j, NULL, list2); */
/* 	  //printf("list2 size %lu\n", array_list_size(list2)); */
/* 	} */
/*       } */
/*       array_list_clear(list2, (void *) cal_free); */
/*       num_cals2 = array_list_size(cal_list2); */
/*       //printf("read %lu\t2recortar lista\n", read_index2); */
/*     } */

/*     if (num_cals2 > MAX_CALS) { */
/*       //printf("read %lu\tacortar lista\n", read_index2); */
/*       for (size_t j = num_cals2 - 1; j >= MAX_CALS; j--) { */
/* 	cal = (cal_t *) array_list_remove_at(j, cal_list2); */
/* 	cal_free(cal); */
/*       } */
/*       num_cals2 = array_list_size(cal_list2); */
/*       //printf("read %lu\tacortar lista\n", read_index2); */
/*     } */

/*     //printf("read %lu\tcals2 = %lu\n", read_index2, num_cals2); */

/*     if (num_cals2 > 0 && num_cals2 <= MAX_CALS) { */
/*       ///////////// */
/*       //reads_mapp2++; */
/*       ///////////// */

/*       //printf("read %lu\tactualizar lista\n", read_index2); */
/*       array_list_set_flag(2, cal_list2); */
/*       mapping_batch->num_to_do2 += num_cals2; */
/*       targets2[new_num_targets2++] = read_index2; */
      
/*       // we have to free the region list */
/*       array_list_free(mapping_batch->mapping_lists2[read_index2], (void *) region_bwt_free); */
/*       mapping_batch->mapping_lists2[read_index2] = cal_list2; */
/*       //printf("read %lu\tactualizar lista\n", read_index2); */
/*     } else { */
/*       ///////////// */
/*       //reads_discard2++; */
/*       ///////////// */

/*       //printf("read %lu\tdescartar lista\n", read_index2); */
/*       array_list_set_flag(0, mapping_batch->mapping_lists2[read_index2]); */
/*       // we have to free the region list */
/*       array_list_clear(mapping_batch->mapping_lists2[read_index2], (void *) region_bwt_free); */
/*       if (cal_list2) array_list_free(cal_list2, (void *) cal_free); */
/*       if (list2) array_list_clear(list2, (void *) cal_free); */
/*       //printf("read %lu\tdescartar lista\n", read_index2); */
/*     } */

/*     //printf("read %lu\tcals2 = %lu\n", read_index2, num_cals2); */

/*     //printf("read %lu\tcals1 = %lu\tsize mapps  %i\n", read_index, num_cals, array_list_size(mapping_batch->mapping_lists[read_index])); */
/*     //printf("read %lu\tcals2 = %lu\tsize mapps2 %lu\n", read_index2, num_cals2, array_list_size(mapping_batch->mapping_lists2[read_index2])); */

    
/*   } // end for 0 ... num_seqs */
/*   // end add for bisulfite */

/*   /\* */
/*   printf("CAL_seek1   \t%3lu\thave CAL (to SW)\t%3lu\thave no CALs     \t%3lu\n",  */
/* 	 num_targets, reads_mapp, reads_discard); */
/*   printf("CAL_seek2   \t%3lu\thave CAL (to SW)\t%3lu\thave no CALs     \t%3lu\n",  */
/* 	 num_targets2, reads_mapp2, reads_discard2); */
/*   *\/ */
/*   /\* */
/*   printf("2 reads CAL   \t%3lu\ttotal CALs      \t%3lu\tCALs promedio    \t%6.2f\n",  */
/* 	 reads_cals, reads_cals_acum, 1.0 * reads_cals_acum / reads_cals); */
/*   *\/ */

/*   //printf("new targets1 = %lu, new targets2 = %lu\n", new_num_targets, new_num_targets2); */

/*   // update batch */
/*   mapping_batch->num_targets = new_num_targets; */

/*   // free memory */
/*   if (list) array_list_free(list, NULL); */

/*   // added for bisulfite */
/*   // update batch */
/*   mapping_batch->num_targets2 = new_num_targets2; */

/*   // free memory */
/*   if (list2) array_list_free(list2, NULL); */
/*   // end added for bisulfite */

/*   //printf("End cal stage\nSW_STAGE = %lu\n", SW_STAGE); */

/*   //return SW_STAGE; */

/*   if (batch->pair_input->pair_mng->pair_mode != SINGLE_END_MODE) { */
/*     //printf("return PRE_PAIR_STAGE\n"); */
/*     return BS_PRE_PAIR_STAGE; */
/*   } else if (batch->mapping_batch->num_targets > 0 || batch->mapping_batch->num_targets2 > 0) { */
/*     //printf("return SW_STAGE\n"); */
/*     return BS_SW_STAGE; */
/*   } */

/*   //printf("return POST_PAIR_STAGE\n"); */
/*   return BS_POST_PAIR_STAGE; */
/* } */

//------------------------------------------------------------------------------------

/* int apply_caling_bs_OOOOLDDDDDDD(cal_seeker_input_t* input, batch_t *batch) { */

/*   //printf("APPLY CALLING BS...\n"); */

/*   mapping_batch_t *mapping_batch = batch->mapping_batch; */
/*   array_list_t *list = NULL; */
/*   size_t read_index, num_cals, min_seeds, max_seeds; */
/*   int min_limit; */

/*   cal_t *cal; */
/*   array_list_t *cal_list; */

/*   size_t num_chromosomes = input->genome->num_chromosomes + 1; */

/*   size_t num_targets = mapping_batch->num_targets; */
/*   size_t *targets = mapping_batch->targets; */
/*   size_t new_num_targets = 0; */
  
/*   // new variables for bisulfite */
/*   size_t num_targets2 = mapping_batch->num_targets2; */
/*   size_t *targets2 = mapping_batch->targets2; */
/*   size_t new_num_targets2 = 0; */
/*   array_list_t *list2 = NULL; */
/*   size_t num_cals2; */
/*   array_list_t *cal_list2; */
/*   size_t read_index2; */

/*   // set to zero */
/*   mapping_batch->num_to_do = 0; */
/*   mapping_batch->num_to_do2 = 0; */


/*   //////////////////////////////// */
/*   /\* */
/*   size_t reads_mapp  = 0; */
/*   size_t reads_mapp2 = 0; */
/*   size_t reads_no_mapp  = 0; */
/*   size_t reads_no_mapp2 = 0; */
/*   size_t reads_discard  = 0; */
/*   size_t reads_discard2 = 0; */

/*   size_t reads_cals = 0; */
/*   size_t reads_cals_acum = 0; */
/*   *\/ */
/*   //////////////////////////////// */


/*   //printf("targets 1 %lu\ntargets 2 %lu\n", num_targets, num_targets2); */

/*   //printf("----->primera lista\n"); */
/*   for (size_t i = 0; i < num_targets; i++) { */

/*     read_index = targets[i]; */

/*     // for debugging */
/*     //    LOG_DEBUG_F("%s\n", ((fastq_read_t *) array_list_get(read_index, mapping_batch->fq_batch))->id); */
    
/*     if (!list) { */
/*       list = array_list_new(1000,  */
/* 			    1.25f,  */
/* 			    COLLECTION_MODE_ASYNCHRONIZED); */
/*     } */

/*     // optimized version */
/*     num_cals = bwt_generate_cal_list_linkedlist(mapping_batch->mapping_lists[read_index],  */
/* 						input->cal_optarg, */
/* 						&min_seeds, &max_seeds, */
/* 						num_chromosomes, */
/* 						list); */
/*     //printf("read %lu\tcals1 = %lu\n", read_index, num_cals); */

/*     /\* */
/*     // for debugging */
/*     LOG_DEBUG_F("num. cals = %i, min. seeds = %i, max. seeds = %i\n", num_cals, min_seeds, max_seeds); */
/*     for (size_t j = 0; j < num_cals; j++) { */
/*       cal = array_list_get(j, list); */
/*       LOG_DEBUG_F("\tchr: %i, strand: %i, start: %lu, end: %lu, num. seeds = %lu, flank: (start, end) = (%lu, %lu)\n",  */
/* 		  cal->chromosome_id, cal->strand, cal->start, cal->end, cal->num_seeds, cal->flank_start, cal->flank_end); */
/*     } */
/*     *\/ */

/*     // filter CALs by the number of seeds */
/*     int min_limit = input->cal_optarg->min_num_seeds_in_cal; */

/*     if (min_limit < 0) min_limit = max_seeds; */

/*     if (min_seeds == max_seeds || min_limit <= min_seeds) { */
/*       cal_list = list; */
/*       list = NULL; */

/*     } else { */
/*       cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED); */

/*       ///////////////// */
/*       /\* */
/*       if (num_cals > 0) */
/* 	reads_cals++; */
/*       *\/ */
/*       ///////////////// */

/*       for (size_t j = 0; j < num_cals; j++) { */
/* 	cal = array_list_get(j, list); */
/* 	//filter cals with few seeds */
/* 	if (cal->num_seeds >= min_limit) { */

/* 	  ///////////////// */
/* 	  //reads_cals_acum++; */
/* 	  ///////////////// */


/* 	  //	  LOG_DEBUG_F("\t\tchr: %i, strand: %i, start: %lu, end: %lu, num. seeds = %lu (min. limit %i seeds), flank: (start, end) = (%lu, %lu)\n",  */
/* 	  //		      cal->chromosome_id, cal->strand, cal->start, cal->end, cal->num_seeds, min_limit, cal->flank_start, cal->flank_end); */
	  
/* 	  array_list_insert(cal, cal_list); */
/* 	  array_list_set(j, NULL, list); */
/* 	} */
/*       } */
/*       array_list_clear(list, (void *) cal_free); */
/*       num_cals = array_list_size(cal_list); */
/*     } */

/*     if (num_cals > MAX_CALS) { */
/*       for (size_t j = num_cals - 1; j >= MAX_CALS; j--) { */
/* 	cal = (cal_t *) array_list_remove_at(j, cal_list); */
/* 	cal_free(cal); */
/*       } */
/*       num_cals = array_list_size(cal_list); */
/*     } */

/*     if (num_cals > 0 && num_cals <= MAX_CALS) {   */
/*       ///////////// */
/*       //reads_mapp++; */
/*       ///////////// */

/*       array_list_set_flag(2, cal_list); */
/*       mapping_batch->num_to_do += num_cals; */
/*       targets[new_num_targets++] = read_index; */
      
/*       // we have to free the region list */
/*       array_list_free(mapping_batch->mapping_lists[read_index], (void *) region_bwt_free); */
/*       mapping_batch->mapping_lists[read_index] = cal_list; */
/*     } else { */
/*       ///////////// */
/*       //reads_discard++; */
/*       ///////////// */
/*       array_list_set_flag(0, mapping_batch->mapping_lists[read_index]); */
/*       // we have to free the region list */
/*       array_list_clear(mapping_batch->mapping_lists[read_index], (void *) region_bwt_free); */
/*       if (cal_list) array_list_free(cal_list, (void *) cal_free); */
/*       if (list) array_list_clear(list, (void *) cal_free); */
/*     } */

/*     //printf("read %lu\tcals  = %lu\n", read_index, num_cals); */
/*   } // end for 0 ... num_seqs */


/*   //printf("----->segunda lista\n"); */
/*   // add for bisulfite */
/*   for (size_t i = 0; i < num_targets2; i++) { */

/*     read_index2 = targets2[i]; */
/*     //printf("\nread %lu\n", read_index2); */

/*     // for debugging */
/*     //    LOG_DEBUG_F("%s\n", ((fastq_read_t *) array_list_get(read_index, mapping_batch->fq_batch))->id); */
    
/*     if (!list2) { */
/*       list2 = array_list_new(1000,  */
/* 			     1.25f,  */
/* 			     COLLECTION_MODE_ASYNCHRONIZED); */
/*     } */

/*     // optimized version */
/*     num_cals2 = bwt_generate_cal_list_linkedlist(mapping_batch->mapping_lists2[read_index2],  */
/* 						 input->cal_optarg, */
/* 						 &min_seeds, &max_seeds, */
/* 						 num_chromosomes, */
/* 						 list2); */
/*     //printf("read %lu\tcals2 = %lu\n", read_index2, num_cals2); */

/*     /\* */
/*     // for debugging */
/*     LOG_DEBUG_F("num. cals = %i, min. seeds = %i, max. seeds = %i\n", num_cals2, min_seeds, max_seeds); */
/*     for (size_t j = 0; j < num_cals2; j++) { */
/*       cal = array_list_get(j, list); */
/*       LOG_DEBUG_F("\tchr: %i, strand: %i, start: %lu, end: %lu, num. seeds = %lu, flank: (start, end) = (%lu, %lu)\n",  */
/* 		  cal->chromosome_id, cal->strand, cal->start, cal->end, cal->num_seeds, cal->flank_start, cal->flank_end); */
/*     } */
/*     *\/ */

/*     // filter CALs by the number of seeds */
/*     int min_limit = input->cal_optarg->min_num_seeds_in_cal; */

/*     if (min_limit < 0) min_limit = max_seeds; */

/*     if (min_seeds == max_seeds || min_limit <= min_seeds) { */
/*       //printf("read %lu\tborrar listas\n", read_index2); */
/*       cal_list2 = list2; */
/*       list2 = NULL; */
/*       //printf("read %lu\tborrar listas\n", read_index2); */

/*     } else { */
/*       //printf("read %lu\t1recortar lista\n", read_index2); */
/*       cal_list2 = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED); */

/*       ///////////////// */
/*       /\* */
/*       if (num_cals2 > 0) */
/* 	reads_cals++; */
/*       *\/ */
/*       ///////////////// */

/*       for (size_t j = 0; j < num_cals2; j++) { */
/* 	cal = array_list_get(j, list2); */
/* 	//printf("read %lu\tnum_seeds %lu\tmin_limit %i\n", read_index2, cal->num_seeds, min_limit); */
/* 	//filter cals with few seeds */
/* 	if (cal->num_seeds >= min_limit) { */

/* 	  ///////////////// */
/* 	  //reads_cals_acum++; */
/* 	  ///////////////// */

/* 	  //printf("keep cal %lu, with %lu seed (of %lu needed)\n", j, cal->num_seeds, min_limit); */
/* 	  //	  LOG_DEBUG_F("\t\tchr: %i, strand: %i, start: %lu, end: %lu, num. seeds = %lu (min. limit %i seeds), flank: (start, end) = (%lu, %lu)\n",  */
/* 	  //		      cal->chromosome_id, cal->strand, cal->start, cal->end, cal->num_seeds, min_limit, cal->flank_start, cal->flank_end); */
	  
/* 	  //printf("pre  callist2 size %lu\n", array_list_size(cal_list2)); */
/* 	  array_list_insert(cal, cal_list2); */
/* 	  //printf("post callist2 size %lu\n", array_list_size(cal_list2)); */
/* 	  array_list_set(j, NULL, list2); */
/* 	  //printf("list2 size %lu\n", array_list_size(list2)); */
/* 	} */
/*       } */
/*       array_list_clear(list2, (void *) cal_free); */
/*       num_cals2 = array_list_size(cal_list2); */
/*       //printf("read %lu\t2recortar lista\n", read_index2); */
/*     } */

/*     if (num_cals2 > MAX_CALS) { */
/*       //printf("read %lu\tacortar lista\n", read_index2); */
/*       for (size_t j = num_cals2 - 1; j >= MAX_CALS; j--) { */
/* 	cal = (cal_t *) array_list_remove_at(j, cal_list2); */
/* 	cal_free(cal); */
/*       } */
/*       num_cals2 = array_list_size(cal_list2); */
/*       //printf("read %lu\tacortar lista\n", read_index2); */
/*     } */

/*     //printf("read %lu\tcals2 = %lu\n", read_index2, num_cals2); */

/*     if (num_cals2 > 0 && num_cals2 <= MAX_CALS) { */
/*       ///////////// */
/*       //reads_mapp2++; */
/*       ///////////// */

/*       //printf("read %lu\tactualizar lista\n", read_index2); */
/*       array_list_set_flag(2, cal_list2); */
/*       mapping_batch->num_to_do2 += num_cals2; */
/*       targets2[new_num_targets2++] = read_index2; */
      
/*       // we have to free the region list */
/*       array_list_free(mapping_batch->mapping_lists2[read_index2], (void *) region_bwt_free); */
/*       mapping_batch->mapping_lists2[read_index2] = cal_list2; */
/*       //printf("read %lu\tactualizar lista\n", read_index2); */
/*     } else { */
/*       ///////////// */
/*       //reads_discard2++; */
/*       ///////////// */

/*       //printf("read %lu\tdescartar lista\n", read_index2); */
/*       array_list_set_flag(0, mapping_batch->mapping_lists2[read_index2]); */
/*       // we have to free the region list */
/*       array_list_clear(mapping_batch->mapping_lists2[read_index2], (void *) region_bwt_free); */
/*       if (cal_list2) array_list_free(cal_list2, (void *) cal_free); */
/*       if (list2) array_list_clear(list2, (void *) cal_free); */
/*       //printf("read %lu\tdescartar lista\n", read_index2); */
/*     } */

/*     //printf("read %lu\tcals2 = %lu\n", read_index2, num_cals2); */

/*     //printf("read %lu\tcals1 = %lu\tsize mapps  %i\n", read_index, num_cals, array_list_size(mapping_batch->mapping_lists[read_index])); */
/*     //printf("read %lu\tcals2 = %lu\tsize mapps2 %lu\n", read_index2, num_cals2, array_list_size(mapping_batch->mapping_lists2[read_index2])); */

    
/*   } // end for 0 ... num_seqs */
/*   // end add for bisulfite */

/*   /\* */
/*   printf("CAL_seek1   \t%3lu\thave CAL (to SW)\t%3lu\thave no CALs     \t%3lu\n",  */
/* 	 num_targets, reads_mapp, reads_discard); */
/*   printf("CAL_seek2   \t%3lu\thave CAL (to SW)\t%3lu\thave no CALs     \t%3lu\n",  */
/* 	 num_targets2, reads_mapp2, reads_discard2); */
/*   *\/ */
/*   /\* */
/*   printf("2 reads CAL   \t%3lu\ttotal CALs      \t%3lu\tCALs promedio    \t%6.2f\n",  */
/* 	 reads_cals, reads_cals_acum, 1.0 * reads_cals_acum / reads_cals); */
/*   *\/ */

/*   //printf("new targets1 = %lu, new targets2 = %lu\n", new_num_targets, new_num_targets2); */

/*   // update batch */
/*   mapping_batch->num_targets = new_num_targets; */

/*   // free memory */
/*   if (list) array_list_free(list, NULL); */

/*   // added for bisulfite */
/*   // update batch */
/*   mapping_batch->num_targets2 = new_num_targets2; */

/*   // free memory */
/*   if (list2) array_list_free(list2, NULL); */
/*   // end added for bisulfite */

/*   //printf("End cal stage\nSW_STAGE = %lu\n", SW_STAGE); */

/*   //return SW_STAGE; */

/*   if (batch->pair_input->pair_mng->pair_mode != SINGLE_END_MODE) { */
/*     //printf("return PRE_PAIR_STAGE\n"); */
/*     return BS_PRE_PAIR_STAGE; */
/*   } else if (batch->mapping_batch->num_targets > 0 || batch->mapping_batch->num_targets2 > 0) { */
/*     //printf("return SW_STAGE\n"); */
/*     return BS_SW_STAGE; */
/*   } */

/*   //printf("return POST_PAIR_STAGE\n"); */
/*   return BS_POST_PAIR_STAGE; */
/* } */

//------------------------------------------------------------------------------------
// cal_seeker_input functions: init
//------------------------------------------------------------------------------------

void cal_seeker_input_init(list_t *regions_list, cal_optarg_t *cal_optarg, 
			   list_t* write_list, unsigned int write_size, 
			   list_t *sw_list, list_t *pair_list, 
			   genome_t *genome, bwt_optarg_t *bwt_optarg, 
			   bwt_index_t *index, metaexons_t *metaexons, 
			   cal_seeker_input_t *input) {
  input->regions_list = regions_list;
  input->cal_optarg = cal_optarg;
  input->batch_size = write_size;
  input->sw_list = sw_list;
  input->pair_list = pair_list;
  input->write_list = write_list;
  input->genome = genome;
  input->bwt_optarg = bwt_optarg;
  input->index = index;
  input->metaexons = metaexons;
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
