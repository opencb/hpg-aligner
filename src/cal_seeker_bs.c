#include "cal_seeker.h"

//------------------------------------------------------------------------------------
// functions to handle gaps between mapped seeds
//    - fill_gaps
//    - merge_seed_regions
//------------------------------------------------------------------------------------

void display_sr_lists_bs(char *msg, mapping_batch_t *mapping_batch, int bs_id) {

  fastq_read_t *read1;
  fastq_read_t *read2;
  array_list_t *fq_batch1;
  array_list_t *fq_batch2;

  size_t read_index;

  cal_t *cal;
  array_list_t *cal_list = NULL;
  array_list_t **mapping_lists;
  size_t num_cals, num_targets;
  size_t *targets;

  if (bs_id == 0) {
    mapping_lists = mapping_batch->mapping_lists;
    targets = mapping_batch->targets;
    num_targets = mapping_batch->num_targets;
    fq_batch1 = mapping_batch->GA_fq_batch;
    fq_batch2 = mapping_batch->GA_rev_fq_batch;
  } else {
    mapping_lists = mapping_batch->mapping_lists2;
    targets = mapping_batch->targets2;
    num_targets = mapping_batch->num_targets2;
    fq_batch1 = mapping_batch->CT_fq_batch;
    fq_batch2 = mapping_batch->CT_rev_fq_batch;
  }

  seed_region_t *s, *prev_s, *new_s;
  linked_list_iterator_t* itr;

  LOG_DEBUG_F("%s\n", msg);

  // debugging....
  for (size_t i = 0; i < num_targets; i++) {
    read_index = targets[i];
    read1 = (fastq_read_t *) array_list_get(read_index, fq_batch1);
    read2 = (fastq_read_t *) array_list_get(read_index, fq_batch2);
    
    LOG_DEBUG_F("Read %s\n", read1->id);
    
    cal_list = mapping_lists[read_index];
    num_cals = array_list_size(cal_list);
    
    if (num_cals <= 0) continue;
    
    // processing each CAL from this read
    for(size_t j = 0; j < num_cals; j++) {
      /*
      if (cal->strand) {
	LOG_DEBUG_F("Read %s\n", read1->sequence);
      } else {
	LOG_DEBUG_F("Read %s\n", read2->sequence);
      }
      */
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

void fill_gaps_bs(mapping_batch_t *mapping_batch, sw_optarg_t *sw_optarg, 
		  genome_t *genome1, genome_t *genome2, int min_gap, int min_distance,
		  int bs_id) {

  array_list_t **mapping_lists;
  size_t num_targets;
  size_t *targets;
  array_list_t *fq_batch;
  array_list_t *fq_batch1;
  array_list_t *fq_batch2;

  if (bs_id == 0) {
    mapping_lists = mapping_batch->mapping_lists;
    num_targets = mapping_batch->num_targets;
    targets = mapping_batch->targets;
    fq_batch1 = mapping_batch->GA_fq_batch;
    fq_batch2 = mapping_batch->GA_rev_fq_batch;
  } else {
    mapping_lists = mapping_batch->mapping_lists2;
    num_targets = mapping_batch->num_targets2;
    targets = mapping_batch->targets2;
    fq_batch1 = mapping_batch->CT_fq_batch;
    fq_batch2 = mapping_batch->CT_rev_fq_batch;
  }

  int sw_count = 0;

  fastq_read_t *read;
  fastq_read_t *read1;
  fastq_read_t *read2;

  size_t read_index, read_len;

  cal_t *cal;
  array_list_t *cal_list = NULL;
  size_t num_cals;

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

    read_index = targets[i];
    //    read = (fastq_read_t *) array_list_get(read_index, fq_batch);
    read1 = (fastq_read_t *) array_list_get(read_index, fq_batch1);
    read2 = (fastq_read_t *) array_list_get(read_index, fq_batch2);
    
    cal_list = mapping_lists[read_index];
    num_cals = array_list_size(cal_list);
    
    if (num_cals <= 0) continue;

    read_len = read1->length;

    min_distance = read_len*0.2;

    LOG_DEBUG_F(">>>>> read %s\n", read1->id);
    //    printf(">>>>> read %s\n", read->id);

    // processing each CAL from this read
    for(size_t j = 0; j < num_cals; j++) {

      // get cal and read index
      cal = array_list_get(j, cal_list);
      cal_print(cal);
      //      LOG_DEBUG_F("CAL #%i of %i (strand %i), sr_list size = %i, sr_duplicate_list size = %i\n", 
      //		  j, num_cals, cal->strand, cal->sr_list->size, cal->sr_duplicate_list->size);

      prev_s = NULL;
      itr = linked_list_iterator_new(cal->sr_list);
      s = (seed_region_t *) linked_list_iterator_curr(itr);
      while (s != NULL) {
	/*
	{
	  // for debugging
	  size_t start = s->genome_start;// + 1;
	  size_t end = s->genome_end;// + 1;
	  size_t len = end - start + 1;
	  //	  printf(":::::::::: %lu - %lu = %i ::::::::::::\n", end, start, len );
	  char *ref = (char *) malloc((len + 1) * sizeof(char));
	  genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
					    &start, &end, genome1);
	  ref[len] = '\0';
	  //
	  LOG_DEBUG_F("\tseed: [%i|%i - %i|%i] %s (len = %i)\n", 
		      s->genome_start, s->read_start, s->read_end, s->genome_end, ref, len);
	  free(ref);
	}
	*/
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
	      printf("**************************\n");
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
	      // handle strand -
	      if (cal->strand) {
		genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
						  &start, &end, genome2);
		query = &read2->sequence[gap_read_start];
		//query = &read2->sequence[read2->length - gap_read_end];
	      } else {
		genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
						  &start, &end, genome1);
		query = &read1->sequence[gap_read_start];
	      }
	      //query[gap_read_len] = '\0';

	      //query = &read->sequence[gap_read_start];
	      
	      for (int k = 0; k < gap_read_len; k++) {
		if (query[k] != ref[k]) {
		  distance++;
		  if (first == -1) first = k;
		  last = k;
		}
	      }

	      LOG_DEBUG_F("query    : %s (%lu-%lu %c)\n", query, gap_read_start, gap_read_end, (cal->strand) ? '-':'+');
	      LOG_DEBUG_F("ref      : %s (%i:%lu-%lu)\n", ref, cal->chromosome_id, start, end);
	      LOG_DEBUG_F("distance : %i\n", distance);

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
	      memcpy(query, &read->sequence[read_start], gap_read_len_ex);
	      query[gap_read_len_ex] = '\0';

	      // get ref. sequence
	      size_t genome_start = gap_genome_start - left_flank;// + 1;
	      size_t genome_end = gap_genome_end + right_flank;// + 1;
	      int gap_genome_len_ex = genome_end - genome_start + 1;
	      ref = (char *) malloc((gap_genome_len_ex + 1) * sizeof(char));;
	      // handle strand -
	      if (cal->strand) {
		genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
						  &genome_start, &genome_end, genome2);	      
	      } else {
		genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
						  &genome_start, &genome_end, genome1);	      
	      }
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
	  // handle strand -
	  if (cal->strand) {
	    genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
					      &start, &end, genome2);
	    query = &read2->sequence[gap_read_start];
	    //query = &read2->sequence[read2->length - gap_read_end];
	  } else {
	    genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
					      &start, &end, genome1);
	    query = &read1->sequence[gap_read_start];
	  }
	  //query[gap_read_len] = '\0';
	  
	  //query = &read->sequence[gap_read_start];
	  
	  distance = 0;
	  for (int k = 0; k < gap_read_len; k++) {
	    if (query[k] != ref[k]) {
	      distance++;
	      if (first == -1) first = k;
	      last = k;
	    }
	  }

	  LOG_DEBUG_F("query    : %s (%lu-%lu %c)\n", query, gap_read_start, gap_read_end, (cal->strand) ? '-':'+');
	  LOG_DEBUG_F("ref      : %s (%i:%lu-%lu)\n", ref, cal->chromosome_id, start, end);
	  LOG_DEBUG_F("distance : %i\n", distance);

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
	    //memcpy(query, &read->sequence[read_start], gap_read_len_ex);
	    
	    // get ref. sequence
	    size_t genome_start = gap_genome_start - left_flank;// + 1;
	    size_t genome_end = gap_genome_end + right_flank;// + 1;
	    int gap_genome_len_ex = genome_end - genome_start + 1;
	    ref = (char *) malloc((gap_genome_len_ex + 1) * sizeof(char));;
	    // handle strand -
	    if (cal->strand) {
	      genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
						&genome_start, &genome_end, genome2);
	      memcpy(query, &read2->sequence[read_start], gap_read_len_ex);
	    } else {
	      genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
						&genome_start, &genome_end, genome1);
	      memcpy(query, &read1->sequence[read_start], gap_read_len_ex);
	    }
	    ref[gap_genome_len_ex] = '\0';
	    query[gap_read_len_ex] = '\0';
	    
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
  
  LOG_DEBUG_F("R U N   S W (sw_count = %i, sw_prepare_list size = %i)\n", sw_count, array_list_size(sw_prepare_list));
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

  //display_sr_lists("END of fill_gaps", mapping_batch);
  display_sr_lists_bs("END of fill_gaps BS", mapping_batch, bs_id);
    
  // free memory
  sw_multi_output_free(output);
  array_list_free(sw_prepare_list, (void *) NULL);
}

//------------------------------------------------------------------------------------

void fill_end_gaps_bs(mapping_batch_t *mapping_batch, sw_optarg_t *sw_optarg, 
		      genome_t *genome1, genome_t *genome2, int min_H, int min_distance,
		      int bs_id) {
  array_list_t **mapping_lists;
  size_t num_targets;
  size_t *targets;
  array_list_t *fq_batch;
  array_list_t *fq_batch1;
  array_list_t *fq_batch2;

  fastq_read_t *read1, *read2;

  if (bs_id == 0) {
    mapping_lists = mapping_batch->mapping_lists;
    num_targets = mapping_batch->num_targets;
    targets = mapping_batch->targets;
    fq_batch1 = mapping_batch->GA_fq_batch;
    fq_batch2 = mapping_batch->GA_rev_fq_batch;
  } else {
    mapping_lists = mapping_batch->mapping_lists2;
    num_targets = mapping_batch->num_targets2;
    targets = mapping_batch->targets2;
    fq_batch1 = mapping_batch->CT_fq_batch;
    fq_batch2 = mapping_batch->CT_rev_fq_batch;
  }

  int sw_count = 0;

  size_t read_index, read_len;

  cal_t *cal;
  array_list_t *cal_list = NULL;
  size_t num_cals;

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

    read_index = targets[i];
    //    read = (fastq_read_t *) array_list_get(read_index, fq_batch);
    read1 = (fastq_read_t *) array_list_get(read_index, fq_batch1);
    read2 = (fastq_read_t *) array_list_get(read_index, fq_batch2);

    cal_list = mapping_lists[read_index];
    num_cals = array_list_size(cal_list);
    
    if (num_cals <= 0) continue;

    read_len = read1->length;
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
	// get ref. sequence
	start = gap_genome_start;// + 1;
	end = gap_genome_end;// + 1;
	gap_genome_len = end - start + 1;
	ref = (char *) malloc((gap_genome_len + 1) * sizeof(char));
	if (cal->strand) {
	  seq = read2->sequence;
	  genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
					    &start, &end, genome2);
	} else {
	  seq = read1->sequence;
	  genome_read_sequence_by_chr_index(ref, 0, cal->chromosome_id - 1, 
					    &start, &end, genome1);
	}

	gap_read_len = gap_read_end - gap_read_start + 1;
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
    read_index = targets[i];
    read1 = (fastq_read_t *) array_list_get(read_index, fq_batch1);
    
    cal_list = mapping_lists[read_index];
    num_cals = array_list_size(cal_list);

    LOG_DEBUG_F("Read %s: num_cals = %lu\n", read1->id, num_cals);
    
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

void merge_seed_regions_bs(mapping_batch_t *mapping_batch, int bs_id) {

  array_list_t **mapping_lists;
  register size_t num_targets;
  size_t *targets;

  if (bs_id == 0) {
    mapping_lists = mapping_batch->mapping_lists;
    num_targets = mapping_batch->num_targets;
    targets = mapping_batch->targets;
  } else {
    mapping_lists = mapping_batch->mapping_lists2;
    num_targets = mapping_batch->num_targets2;
    targets = mapping_batch->targets2;
  }

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

  register size_t num_cals;
  register int i;
  register size_t t;

  for (t = 0; t < num_targets; t++) {
    cals_list = mapping_lists[targets[t]];
    num_cals = array_list_size(cals_list);
    fq_read = array_list_get(targets[t], mapping_batch->fq_batch);

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
    }
  }
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
