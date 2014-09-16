#include "sa_dna_commons.h"


//--------------------------------------------------------------------
// commons
//--------------------------------------------------------------------

void seed_free(seed_t *p) {
  if (p) {
    cigar_clean(&p->cigar);
    free(p);
  }
}

//--------------------------------------------------------------------

void seed_cal_free(seed_cal_t *p) {
  if (p) {

    //    printf("freeing cal of read %s\n", p->read->id);
    //    seed_cal_print(p);

    cigar_clean(&p->cigar);
    if (p->seed_list) linked_list_free(p->seed_list, (void *) seed_free);
    if (p->cigarset) cigarset_free(p->cigarset);
    free(p);
  }
}

//--------------------------------------------------------------------
// utils
//--------------------------------------------------------------------

#ifdef _TIMING
void init_func_names() {
  strcpy(func_names[0], "search_suffix");
  strcpy(func_names[1], "search_prefix");
  strcpy(func_names[2], "search_sa");
  strcpy(func_names[3], "generate_cals_from_exact_read");
  strcpy(func_names[4], "generate_cals_from_suffixes");
  strcpy(func_names[5], "init_cals_from_suffixes");
  strcpy(func_names[6], "set_positions");
  strcpy(func_names[7], "set_reference_sequence");
  strcpy(func_names[8], "skip_suffixes");
  strcpy(func_names[9], "mini_sw_right_side");
  strcpy(func_names[10], "mini_sw_left_side");
  strcpy(func_names[11], "seed_new");
  strcpy(func_names[12], "seed_list_insert");
  strcpy(func_names[13], "cal_new");
  strcpy(func_names[14], "cal_mng_insert");
  strcpy(func_names[15], "cal_mng_to_array_list");
  strcpy(func_names[16], "filter_by_read_area");
  strcpy(func_names[17], "filter_by_num_mismatches");
  strcpy(func_names[18], "sw_pre_processing");
  strcpy(func_names[19], "sw_execution");
  strcpy(func_names[20], "sw_post_processing");
  strcpy(func_names[21], "other functions");
  strcpy(func_names[22], "create_alignments");
}
#endif

//--------------------------------------------------------------------

float get_max_score(array_list_t *cal_list,
		    float match_score, float mismatch_penalty,
		    float gap_open_penalty, float gap_extend_penalty) {
  seed_cal_t *cal;
  int num_matches, num_mismatches, num_open_gaps, num_extend_gaps;
  int num_cals = array_list_size(cal_list);
  float max_score = -1000000.0f;

  for (int j = 0; j < num_cals; j++) {
    cal = array_list_get(j, cal_list);

    cal->score = cigar_compute_score(match_score, mismatch_penalty, 
				     gap_open_penalty, gap_extend_penalty, &cal->cigar);
    //cal->score = cigar_compute_score(5.0f, -4.0f, -10.0f, -0.05f, &cal->cigar);
    cal->cigar_len = cigar_get_length(&cal->cigar);

    //    printf("(invalid, read length, cigar_len, score) = (%i, %i, %i, %0.2f)\n",
    //	     cal->invalid, cal->read->length, cal->cigar_len, cal->score);

    //seed_cal_print(cal);

    if (!cal->invalid && 
	(cal->cigar_len == cal->read->length) && 
	(cal->num_open_gaps < (0.05f * cal->read->length)) &&
	(cal->num_mismatches < (0.09f * cal->read->length)) &&
	(cal->score > max_score)) {
      max_score = cal->score;
    }
  }
  //printf("----> max. score = %0.2f\n", max_score);
  return max_score;
}

//--------------------------------------------------------------------

int get_min_num_mismatches(array_list_t *cal_list) {
  seed_cal_t *cal;
  size_t num_cals = array_list_size(cal_list);
  int num_mismatches, min_num_mismatches = 100000;
  for (size_t j = 0; j < num_cals; j++) {
    cal = array_list_get(j, cal_list);

    //    printf("previous mapq = %i -> updated (num-cals = %i) = %i\n", cal->mapq, num_cals, cal->mapq + (61 - num_cals));

    /*
    if (cal->mapq > 1000) {
      cal->mapq = (61 - num_cals) - log(cal->mapq);
    } else {
      cal->mapq = (61 - num_cals);
    }
    */
    cal->mapq = (61 - num_cals);
    if (cal->mapq > 60) cal->mapq = 60; 
    if (cal->mapq <= 0) cal->mapq = 1; 

    //    num_mismatches = cal->num_mismatches;
    num_mismatches = cal->num_mismatches + cal->num_open_gaps + cal->num_extend_gaps;
    if (num_mismatches < min_num_mismatches) {
      min_num_mismatches = num_mismatches;
    }
  }
  return min_num_mismatches;
}

//--------------------------------------------------------------------

int get_max_read_area(array_list_t *cal_list) {
  seed_cal_t *cal;
  size_t num_cals = array_list_size(cal_list);
  int max_read_area = 0;
  //printf("------> begin get max read area (num_cals = %i)\n", num_cals);
  //printf("num-cals %04i\n", num_cals);
  for (size_t j = 0; j < num_cals; j++) {
    cal = array_list_get(j, cal_list);
    cal->mapq = num_cals;
    //cal->mapq = 61 - num_cals;
    //if (cal->mapq <= 0) cal->mapq = 1; 
    //    cal->mapq = (num_cals == 1 ? 60 : 1);
    //seed_cal_update_info(cal);
    //seed_cal_print(cal);
    if (cal->read_area > max_read_area) {
      max_read_area = cal->read_area;
    }
  }

  //printf("<------ begin get max read area = %i\n", max_read_area);
  return max_read_area;
}

//--------------------------------------------------------------------

void filter_cals_by_max_read_area(int max_read_area, array_list_t **list) {
  seed_cal_t *cal;
  array_list_t *cal_list = *list;
  size_t num_cals = array_list_size(cal_list);
  array_list_t *new_cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  for (size_t j = 0; j < num_cals; j++) {
    cal = array_list_get(j, cal_list);
    //    seed_cal_print(cal);
    if (cal->read_area >= max_read_area) {
      array_list_insert(cal, new_cal_list);
      array_list_set(j, NULL, cal_list);
      //      break;
    }
  }
  array_list_free(cal_list, (void *) seed_cal_free);
  *list = new_cal_list;
}

//--------------------------------------------------------------------

void filter_cals_by_min_read_area(int read_area, array_list_t **list) {
  seed_cal_t *cal;
  array_list_t *cal_list = *list;
  size_t num_cals = array_list_size(cal_list);
  array_list_t *new_cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  for (size_t j = 0; j < num_cals; j++) {
    cal = array_list_get(j, cal_list);
    if (cal->read_area - cal->num_mismatches <= read_area) {
      array_list_insert(cal, new_cal_list);
      array_list_set(j, NULL, cal_list);
    }
  }
  array_list_free(cal_list, (void *) seed_cal_free);
  *list = new_cal_list;
}

//--------------------------------------------------------------------

void filter_cals_by_max_score(float score, array_list_t **list) {
  seed_cal_t *cal;
  array_list_t *cal_list = *list;
  size_t num_cals = array_list_size(cal_list);
  array_list_t *new_cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  for (size_t j = 0; j < num_cals; j++) {
    cal = array_list_get(j, cal_list);
    if (!cal->invalid && 
	(cal->cigar_len == cal->read->length) && 
	(cal->num_open_gaps < (0.05f * cal->read->length)) &&
	(cal->num_mismatches < (0.09f * cal->read->length)) &&
	(cal->score >= score)) {
      array_list_insert(cal, new_cal_list);
      array_list_set(j, NULL, cal_list);
    }
  }
  array_list_free(cal_list, (void *) seed_cal_free);
  *list = new_cal_list;
}

//--------------------------------------------------------------------

void filter_cals_by_max_num_mismatches(int min_num_mismatches, array_list_t **list) {
  seed_cal_t *cal;
  array_list_t *cal_list = *list;
  size_t num_cals = array_list_size(cal_list);
  array_list_t *new_cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  int num_mismatches;

  for (size_t j = 0; j < num_cals; j++) {
    cal = array_list_get(j, cal_list);
    //    num_mismatches = cal->num_mismatches;
    num_mismatches = cal->num_mismatches + cal->num_open_gaps + cal->num_extend_gaps;
    if (num_mismatches <= min_num_mismatches) {
      array_list_insert(cal, new_cal_list);
      array_list_set(j, NULL, cal_list);
    }
  }
  array_list_free(cal_list, (void *) seed_cal_free);
  *list = new_cal_list;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

typedef struct sort_cal {
  int index;
  seed_cal_t *cal;
} sort_cal_t;

//--------------------------------------------------------------------

int sort_cal_cmp(const void *a, const void *b) {
  seed_cal_t *cal_a = ((sort_cal_t*) a)->cal;
  seed_cal_t *cal_b = ((sort_cal_t*) b)->cal;

  return (cal_b->score - cal_a->score);
}

//--------------------------------------------------------------------

int sort_cal_cmp2(const void *a, const void *b) {
  seed_cal_t *cal_a = ((sort_cal_t*) a)->cal;
  seed_cal_t *cal_b = ((sort_cal_t*) b)->cal;

  int res = cal_b->read_area - cal_a->read_area;
  if (res == 0) {
    res = cal_a->num_mismatches - cal_b->num_mismatches;
  }
  return res;
}

//--------------------------------------------------------------------

int compute_num_hits(int threshold, int num_cals, sort_cal_t *sort_cals) {
  int num_hits = 1;
  seed_cal_t *best_cal = sort_cals[0].cal;
  for (int i = 1; i < num_cals; i++) {
    /*
    printf("abs(best_cal->score %i - sort_cals[i].cal->score %i) = %i\n",
	   (int)best_cal->score, (int)sort_cals[i].cal->score, abs((int)best_cal->score - (int)sort_cals[i].cal->score));
    */
    if (abs((int)best_cal->score - (int)sort_cals[i].cal->score) < threshold) {
      num_hits++;
    }
  }
  return num_hits;
}

//--------------------------------------------------------------------

int compute_mapq(seed_cal_t *best_cal, seed_cal_t *second_cal, int num_hits) {
  int mapq = 0;
  int second_score = (second_cal == NULL) ? 18 : second_cal->score;

  if (second_score >= best_cal->score) return 0;

  double identity, tmp;
  identity = 1. - (double) (best_cal->read_area - best_cal->score) / 5.0 / best_cal->read_area;
  tmp = best_cal->read_area < 50 ? 1. : 3.9120 / log(best_cal->read_area);
  tmp *= identity * identity;
  mapq = (int) (6.02 * (best_cal->score - second_score) * tmp * tmp + .499);
  //mapq = (int) (4.02 * (best_cal->score - second_score) * tmp * tmp + .499);

  if (num_hits > 1) mapq -= (int)(4.343 * log(num_hits) + .499);
  
  if (mapq > 60) mapq = 60;
  if (mapq < 0) mapq = 0;

  return mapq;
}

//--------------------------------------------------------------------

void select_best_cals(fastq_read_t *read, array_list_t **cal_list) {

  int mapq, num_hits = 1;
  int num_cals = array_list_size(*cal_list);

  if (num_cals == 0) return;

  seed_cal_t *cal, *first_cal, *second_cal = NULL;
  sort_cal_t *sort_cal;

  sort_cal_t sort_cals[num_cals];
  for (int i = 0; i < num_cals; i++) {
    cal = array_list_get(i, *cal_list);
    sort_cals[i].index = i;
    sort_cals[i].cal = cal;
  }

  if (num_cals > 1) {
    /*
    printf("before sorting...\n");
    for (int i = 0; i < num_cals; i++) {
      seed_cal_print(sort_cals[i].cal);
    }
    */
    qsort(sort_cals, num_cals, sizeof(sort_cal_t), sort_cal_cmp);
    second_cal = sort_cals[1].cal;

    num_hits = compute_num_hits(7, num_cals, sort_cals);
  }
  first_cal = sort_cals[0].cal;

  mapq = compute_mapq(first_cal, second_cal, num_hits);

  /*
  printf("--------------------------------------------------------------------\n");
  printf("after sorting and trimming (mapq = %i, first score = %0.2f, second score = %0.2f, num_hits = %i) %s ...\n", 
	 mapq, (first_cal == NULL ? -99999.9f : first_cal->score), (second_cal == NULL ? -99999.9f : second_cal->score), num_hits,
	 first_cal->read->sequence);
  printf("--------------------------------------------------------------------\n");
  printf("first CAL:\n");
  seed_cal_print(first_cal);
  if (second_cal) {
    printf("second CAL:\n");
    seed_cal_print(second_cal);
    for (int i = 2; i < num_hits; i++) {
      seed_cal_print(sort_cals[i].cal);
      if (i == 9) {
	if (num_hits - (i + 1) > 0) {
	  printf("\t...\n\t...and %i more\n", num_hits - (i + 1));
	}
	break;
      }
    }
  }
  */

  if (first_cal->invalid) {
    array_list_clear(*cal_list, (void *) seed_cal_free);
  } else {
    for (int i = 0; i < num_hits; i++) {
      cal = sort_cals[i].cal;
      if ( cal->invalid &&
	   ((abs(cal->start - first_cal->end) <= (2 * cal->read->length)) ||
	     (abs(cal->end - first_cal->start) <= (2 * cal->read->length))) ) {
	array_list_clear(*cal_list, (void *) seed_cal_free);
	return;
      }
    }

    if (num_cals > 1) {
      if (mapq > 0) num_hits = 1;
      array_list_t *new_cal_list = array_list_new(num_hits, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
      for (int i = 0; i < num_hits; i++) {
	sort_cals[i].cal->mapq = mapq;
	sort_cals[i].cal->num_hits = num_hits;

	array_list_insert(sort_cals[i].cal, new_cal_list);	
	array_list_set(sort_cals[i].index, NULL, *cal_list);
      }
      array_list_free(*cal_list, (void *) seed_cal_free);
      *cal_list = new_cal_list;      
    } else {
      sort_cals[0].cal->mapq = mapq;
      sort_cals[0].cal->num_hits = num_hits;
    }
  }
  /*
  num_cals = array_list_size(*cal_list);
  for (int i = 0; i < num_cals; i++) {
    seed_cal_print(array_list_get(i, *cal_list));
  }
  */
}
  
//--------------------------------------------------------------------

void select_best_cals2(fastq_read_t *read, array_list_t **cal_list) {

  int max_read_area, min_num_mismatches;

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif
  max_read_area = get_max_read_area(*cal_list);
  if (max_read_area > read->length) max_read_area = read->length;
  filter_cals_by_max_read_area(max_read_area, cal_list);
  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_FILTER_BY_READ_AREA] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif

  #ifdef _VERBOSE
  printf("\t***> max_read_area = %i\n", max_read_area);
  #endif

  min_num_mismatches = get_min_num_mismatches(*cal_list);
  #ifdef _VERBOSE
  printf("\t*** before SW> min_num_mismatches = %i\n", min_num_mismatches);
  #endif
      
  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif
  filter_cals_by_max_num_mismatches(min_num_mismatches, cal_list);
  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_FILTER_BY_NUM_MISMATCHES] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif
}

//--------------------------------------------------------------------

void create_bam_alignments(array_list_t *cal_list, fastq_read_t *read, 
			   array_list_t *mapping_list) {

  // CAL
  seed_cal_t *cal;
  uint num_cals = array_list_size(cal_list);

  // alignments
  alignment_t *alignment;

  if (num_cals <= 0) {
    // no CALs -> no alignment
    return;
  }

  int i, cigar_len;
  linked_list_item_t *list_item; 


  for (i = 0; i < num_cals; i++) {
    cal = array_list_get(i, cal_list);
    
    #ifdef _VERBOSE	  
    printf("--> CAL #%i (cigar %s):\n", i, cigar_to_string(&cal->cigar));
    seed_cal_print(cal);
    #endif
    
    if (!cal->invalid && ((cigar_len = cigar_get_length(&cal->cigar) <= read->length))) {
      cal->AS = (cal->score * 253 / (read->length * 5));
      
      // create the aligments
      alignment = alignment_new();	       
      alignment_init_single_end(strdup(read->id), strdup(read->sequence), strdup(read->quality), 
				cal->strand, cal->chromosome_id, cal->start,
				cigar_to_string(&cal->cigar), cal->cigar.num_ops, cal->AS, 1, (num_cals > 1),
				0, 0, alignment);  
      
      array_list_insert(convert_to_bam(alignment, 33), mapping_list);
      alignment_free(alignment);	 
    }
    
    // free memory
    seed_cal_free(cal);
  }
}

//--------------------------------------------------------------------

void create_alignments(array_list_t *cal_list, fastq_read_t *read, 
		       int bam_format, array_list_t *mapping_list) {

  // CAL
  seed_cal_t *cal;
  uint num_cals = array_list_size(cal_list);

  // alignments
  alignment_t *alignment;

  if (num_cals <= 0) {
    // no CALs -> no alignment
    return;
  }

  cigar_t *cigar;
  char *sequence, *revcomp, *quality;

  char *seq, *p, *optional_fields, *cigar_string, *cigar_M_string;
  int AS, optional_fields_length, num_mismatches, num_cigar_ops, len;
  linked_list_item_t *list_item; 

  for (int i = 0; i < num_cals; i++) {
    cal = array_list_get(i, cal_list);
    
    #ifdef _VERBOSE	  
    printf("--> CAL #%i (cigar %s):\n", i, cigar_to_string(&cal->cigar));
    seed_cal_print(cal);
    #endif

    if (read->adapter) {
      // sequences and cigar
      len = read->length + abs(read->adapter_length);
      sequence = (char *) malloc(len + 1);
      revcomp = (char *) malloc(len + 1);
      quality = (char *) malloc(len + 1);
      cigar = cigar_new_empty();

      if (read->adapter_length < 0) {
	strcpy(quality, read->adapter_quality);
	strcat(quality, read->quality);
      } else {
	strcpy(quality, read->quality);
	strcat(quality, read->adapter_quality);
      }

      if ( (cal->strand == 1 && 
	    ((read->adapter_strand == 0 && read->adapter_length > 0) || 
	     (read->adapter_strand == 1 && read->adapter_length < 0)))
	   ||
	   (cal->strand == 0 && 
	    ((read->adapter_strand == 0 && read->adapter_length < 0) ||
	     (read->adapter_strand == 1 && read->adapter_length > 0))) ) {
	strcpy(sequence, read->adapter);
	strcat(sequence, read->sequence);
	strcpy(revcomp, read->adapter_revcomp);
	strcat(revcomp, read->revcomp);
	
	cigar_append_op(abs(read->adapter_length), 'S', cigar);
	cigar_concat(&cal->cigar, cigar);
      } else {
	strcpy(sequence, read->sequence);
	strcat(sequence, read->adapter);
	strcpy(revcomp, read->revcomp);
	strcat(revcomp, read->adapter_revcomp);
	
	cigar_concat(&cal->cigar, cigar);
	cigar_append_op(read->adapter_length, 'S', cigar);
      }
      sequence[len] = 0; 
      revcomp[len] = 0; 
      quality[len] = 0; 
    } else {
      // sequences and cigar
      sequence = read->sequence;
      revcomp = read->revcomp;
      quality = read->quality;
      cigar = &cal->cigar;
    }

    cigar_string = cigar_to_string(cigar);
    cigar_M_string = cigar_to_M_string(&num_mismatches, &num_cigar_ops, cigar);
    len = strlen(cigar_string);

    optional_fields_length = 100 + len;

    cal->AS = (cal->score * 253 / (read->length * 5));
    AS = (int) cal->score;


    optional_fields = (char *) calloc(optional_fields_length, sizeof(char));
    if (bam_format) {
      p = optional_fields;
      /*
      sprintf(p, "ASi");
      p += 3;
      memcpy(p, &AS, sizeof(int));
      p += sizeof(int);
      */

      sprintf(p, "NHi");
      p += 3;
      memcpy(p, &num_cals, sizeof(int));
      p += sizeof(int);

      sprintf(p, "NMi");
      p += 3;
      memcpy(p, &cal->num_mismatches, sizeof(int));
      p += sizeof(int);

      /*
      sprintf(p, "XCZ");
      p += 3;
      memcpy(p, cigar_string, len);
      p += len;
      */

      optional_fields_length = p - optional_fields;

    } else {
      sprintf(optional_fields, "NM:i:%i", num_mismatches);
      //sprintf(optional_fields, "NM:i:%i\tXC:Z:%s", num_mismatches, cigar_string);
    }
    
    free(cigar_string);

    if (cal->strand) {
      seq = revcomp;
    } else {
      seq = sequence;
    }

    // create the aligments
    alignment = alignment_new();	       
    alignment_init_single_end(strdup(read->id), strdup(seq), strdup(quality), 
			      cal->strand, cal->chromosome_id, cal->start,
			      cigar_M_string, num_cigar_ops, cal->AS, 1, (num_cals > 1),
			      optional_fields_length, optional_fields, alignment);  
    alignment->mapq = cal->mapq;
    alignment->mate_chromosome = 0;
    alignment->mate_position = 0;

    array_list_insert(alignment, mapping_list);
    
    // free memory
    seed_cal_free(cal);
    if (read->adapter) {
      // sequences and cigar
      free(sequence);
      free(revcomp);
      free(quality);
      cigar_free(cigar);
    }
  }
}

//--------------------------------------------------------------------

void display_suffix_mappings(int strand, size_t r_start, size_t suffix_len, 
			     size_t low, size_t high, sa_index3_t *sa_index) {
  int chrom;
  size_t r_end, g_start, g_end;
  for (size_t suff = low; suff <= high; suff++) {
    r_end = r_start + suffix_len - 1;
    chrom = sa_index->CHROM[suff];
    g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];
    g_end = g_start + suffix_len - 1;
    printf("\t\t[%lu|%lu-%lu|%lu] %c chrom %s\n",
	   g_start, r_start, r_end, g_end, (strand == 0 ? '+' : '-'), 
	   sa_index->genome->chrom_names[chrom]);
  }
}

//--------------------------------------------------------------------

void print_seed(char *msg, seed_t *s) {
  printf("%s%c:%i[%lu|%lu - %lu|%lu] suf[%lu|%lu - %lu|%lu] (cigar (len = %i): %s, num. mismatches = %i)\n",  msg, (s->strand == 0 ? '+' : '-'),
	 s->chromosome_id, s->genome_start, s->read_start, s->read_end, s->genome_end,
	 s->suf_genome_start, s->suf_read_start, s->suf_read_end, s->suf_genome_end,
	 cigar_get_length(&s->cigar), /*""*/cigar_to_string(&s->cigar), s->num_mismatches);
}

//--------------------------------------------------------------------

void display_sequence(uint j, sa_index3_t *index, uint len) {
  char *p = &index->genome->S[index->SA[j]];
  char chrom = index->CHROM[j];
  for (int i = 0; i < len; i++) {
    printf("%c", *p);
    p++;
  }
  printf("\t%u\t%s:%lu\n", index->SA[j], 
	 index->genome->chrom_names[chrom], index->SA[j] - index->genome->chrom_offsets[chrom]);
}

//--------------------------------------------------------------------

char *get_subsequence(char *seq, size_t start, size_t len) {
  char *subseq = (char *) malloc((len + 1) * sizeof(char));
  memcpy(subseq, seq + start, len);
  subseq[len] = 0;
  return subseq;
}
//--------------------------------------------------------------------

void display_cmp_sequences(fastq_read_t *read, sa_index3_t *sa_index) {
  size_t pos, chrom, strand;
  char *ref, *seq, *chrom_str, *aux, *p1, *p2;
  
  aux = strdup(read->id);
  p1 = strstr(aux, "_");
  *p1 = 0;
  chrom_str = strdup(aux);
  for (chrom = 0; chrom < sa_index->genome->num_chroms; chrom++) {
    if (strcmp(chrom_str, sa_index->genome->chrom_names[chrom]) == 0) {
      break;
    }
  }
  p2 = strstr(p1 + 1, "_");
  *p2 = 0;
  pos = atol(p1 + 1);
  
  p1 = strstr(p2 + 1, "_");
  p2 = strstr(p1 + 1, "_");
  *p2 = 0;
  strand = atoi(p1 + 1);
  
  free(aux);
  free(chrom_str);
  
  printf("\n======> %s\n", read->id);
  for (int i = 0; i < read->length; i++) {
    if (i % 10 == 0) {
      printf("%i", i / 10);
    } else {
      printf(" ");
    }
  }
  printf("\n");
  for (int i = 0; i < read->length; i++) {
    printf("%i", i % 10);
  }
  printf("\n");
  if (strand) {
    seq = read->revcomp;
  } else {
   seq = read->sequence;
  }
  printf("%s\n", seq);
  ref = &sa_index->genome->S[pos + sa_index->genome->chrom_offsets[chrom] - 1];
  for (int i = 0; i < read->length; i++) {
    if (seq[i] == ref[i]) {
      printf("|");
    } else {
      printf("x");
    }
  }
  printf("\n");
  for (int i = 0; i < read->length; i++) {
    printf("%c", ref[i]);
  }
  printf("\n");
  //  printf("chrom = %i\n", chrom);
}

//--------------------------------------------------------------------
// Support for paired mode
//--------------------------------------------------------------------

void filter_cals_by_pair_mode(int pair_mode, int pair_min_distance, int pair_max_distance, 
			      int num_lists, array_list_t **cal_lists) {

  int diff, pairs, size1, size2;
  seed_cal_t *cal1, *cal2;
  array_list_t *cal_list1, *cal_list2;
  array_list_t *new_cal_list1 = NULL, *new_cal_list2 = NULL;
  for (int i = 0; i < num_lists; i += 2) {
    cal_list1 = cal_lists[i];
    cal_list2 = cal_lists[i+1];

    // prepare new_cal_list1
    if ((size1 = array_list_size(cal_list1)) > 0) {
      if (!new_cal_list1) {
	new_cal_list1 = array_list_new(3 * size1, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
      }
    }

    // prepare new_cal_list2
    if ((size2 = array_list_size(cal_list2)) > 0) {
      if (!new_cal_list2) {
	new_cal_list2 = array_list_new(3 * size2, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
      }
    }

    /*
    printf("-----> PAIR 1:\n");
    for (int i1 = 0; i1 < size1; i1++) {
      cal1 = array_list_get(i1, cal_list1);
      seed_cal_print(cal1);
    }
    printf("-----> PAIR 2:\n");
    for (int i2 = 0; i2 < size2; i2++) {
      cal2 = array_list_get(i2, cal_list2);
      seed_cal_print(cal2);
    }
    */

    if (size1 > 0 && size2 > 0) {
      int pair1[size1], pair2[size2];
      for (int j = 0; j < size1; j++) pair1[j] = 0;
      for (int j = 0; j < size2; j++) pair2[j] = 0;

      pairs = 0;
      for (int i1 = 0; i1 < size1; i1++) {
	cal1 = array_list_get(i1, cal_list1);
	for (int i2 = 0; i2 < size2; i2++) {
	  cal2 = array_list_get(i2, cal_list2);

	  if ((cal1->chromosome_id == cal2->chromosome_id)) {
	    //&&
	    //	      ((cal1->strand != cal2->strand && pair_mode == PAIRED_END_MODE) ||
	    //	       (cal1->strand == cal2->strand && pair_mode == MATE_PAIR_MODE))) {

	    diff = abs(cal2->start - cal1->start) + 1;

	    if (diff >= pair_min_distance && diff <= pair_max_distance) {
	      pair1[i1] = 1;
	      pair2[i2] = 1;
	      pairs++;

	      //printf(">>>>>>>>>> pair %i : distance = %i\n", pairs, diff);
	      //seed_cal_print(cal1);
	      //seed_cal_print(cal2);
	    }
	  }
	}
      }

      //printf("*********** pairs found: %i *****************\n", pairs);

      // update lists
      if (pairs) {
	for (int i1 = 0; i1 < size1; i1++) {
	  if (pair1[i1]) {
	    cal1 = array_list_get(i1, cal_list1);
	    array_list_insert(cal1, new_cal_list1);
	    array_list_set(i1, NULL, cal_list1);
	  }
	}
	array_list_free(cal_list1, (void *) seed_cal_free);
	cal_lists[i] = new_cal_list1;
	new_cal_list1 = NULL;

	for (int i2 = 0; i2 < size2; i2++) {
	  if (pair2[i2]) {
	    cal2 = array_list_get(i2, cal_list2);
	    array_list_insert(cal2, new_cal_list2);
	    array_list_set(i2, NULL, cal_list2);
	  }
	}
	array_list_free(cal_list2, (void *) seed_cal_free);
	cal_lists[i+1] = new_cal_list2;
	new_cal_list2 = NULL;

	if (pairs == 1) {
	  //	  if (cal2->num_hits < 100 && cal1->num_hits < 100 && 
	  if (cal2->mapq == 0 && cal1->mapq == 0) {
	    cal1->mapq = 0; //0; //61;
	    cal2->mapq = 0; //0; //61;
	  }
	}
      }

      /*
      printf("*******-----> MATES #1:\n");
      size1 = array_list_size(cal_lists[i]);
      for (int i1 = 0; i1 < size1; i1++) {
	cal1 = array_list_get(i1, cal_lists[i]);
	seed_cal_print(cal1);
      }
      printf("*******-----> MATES #2:\n");
      size2 = array_list_size(cal_lists[i+1]);
      for (int i2 = 0; i2 < size2; i2++) {
	cal2 = array_list_get(i2, cal_lists[i+1]);
	seed_cal_print(cal2);
      }
      */
    }
  }

  if (new_cal_list1) array_list_free(new_cal_list1, (void *) NULL);
  if (new_cal_list2) array_list_free(new_cal_list2, (void *) NULL);
}

//--------------------------------------------------------------------

array_list_t *create_list(size_t *valid_items, size_t num_valids, array_list_t *list) {
  void *item;
  int num = 0;

  size_t num_items = array_list_size(list);
  //int flag = array_list_get_flag(list);

  array_list_t *new_list = array_list_new(num_valids, 
					  1.25f, 
					  COLLECTION_MODE_ASYNCHRONIZED);
  //array_list_set_flag(flag, new_list);

  for (int k = 0; k < num_items; k++) {
    if (valid_items[k] == 1) {
      array_list_insert(array_list_get(k, list), new_list);
      array_list_set(k, NULL, list);
    }
  }

  //  if (flag == 1) {
  array_list_free(list, (void *) alignment_free);
    //  } else {
    //    array_list_free(list, (void *) seed_cal_free);
    //  }

  return new_list;
}

//--------------------------------------------------------------------

void complete_pairs(sa_mapping_batch_t *batch) {

  size_t num_items1, num_items2, num_pairs, num_reads = array_list_size(batch->fq_reads);

  int distance;
  int min_distance = batch->options->pair_min_distance;
  int max_distance = batch->options->pair_max_distance;
  int pair_mode = batch->options->pair_mode;
  int report_only_paired = batch->options->report_only_paired;

  array_list_t *list1, *list2;

  alignment_t *alig1, *alig2;
  size_t mapped1_counter = 0, mapped2_counter = 0;
  size_t allocated_mapped1 = 100, allocated_mapped2 = 100;
  size_t *mapped1 = (size_t *) malloc(allocated_mapped1 * sizeof(size_t));
  size_t *mapped2 = (size_t *) malloc(allocated_mapped2 * sizeof(size_t));

  short int chr1, chr2, strand1, strand2;
  size_t end1, start2;

  int pair_found;

  int num_best1, num_best2;
  float best_score, best_score1, best_score2, score;

  pair_t *pair, *new_pair;
  linked_list_t *pair_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
  linked_list_iterator_t *pair_list_itr = linked_list_iterator_new(pair_list);
  int num_hits, counter_hits;
  int n_best = batch->options->report_n_best;
  int n_hits = batch->options->report_n_hits;
  int all = batch->options->report_all;
  int best = batch->options->report_best;

  if (n_best) {
    num_hits = n_best;
  } else  if (n_hits) {
    num_hits = n_hits;
  } else {
    num_hits = 10000;
  }

  for (int i = 0; i < num_reads; i += 2) {
    num_pairs = 0;
    best_score = 0;

    //    num_best1 = 0;
    //    num_best2 = 0;
    //    best_score1 = 0;
    //    best_score2 = 0;

    //fastq_read_t *read = array_list_get(i, batch->fq_reads);
    //printf("read id %s\n", read->id);

    list1 = batch->mapping_lists[i];
    list2 = batch->mapping_lists[i+1];

    num_items1 = 0;
    if (list1 != NULL)  num_items1 = array_list_size(list1);
    num_items2 = 0;
    if (list2 != NULL)  num_items2 = array_list_size(list2);

    //printf("%i - %i\n", num_items1, num_items2);

    if (num_items1 > 0 && num_items2 > 0) {
      // initalizes memory and counters
      mapped1_counter = 0;
      if (allocated_mapped1 < num_items1) {
	free(mapped1);
	mapped1 = (size_t *) malloc(num_items1 * sizeof(size_t));
	allocated_mapped1 = num_items1;
      }
      memset(mapped1, 0, num_items1 * sizeof(size_t));

      mapped2_counter = 0;
      if (allocated_mapped2 < num_items2) {
	free(mapped2);
	mapped2 = (size_t *) malloc(num_items2 * sizeof(size_t));
	allocated_mapped2 = num_items2;
      }
      memset(mapped2, 0, num_items2 * sizeof(size_t));

      // search for pairs properly aligned
      for (size_t j1 = 0; j1 < num_items1; j1++) {
	alig1 = (alignment_t *) array_list_get(j1, list1);
	chr1 = alig1->chromosome;
	strand1 = alig1->seq_strand;
	end1 = alig1->position;
	//printf("Item %i Pair1 [chr %i - start %i] cigar = %s, score = %0.2f\n", j1, chr1, end1, alig1->cigar, alig1->map_quality);	
	for (size_t j2 = 0; j2 < num_items2; j2++) {
	  //if (mapped2[j2] == 1) continue;
	  alig2 = (alignment_t *) array_list_get(j2, list2);
	  chr2 = alig2->chromosome;
	  strand2 = alig2->seq_strand;
	  start2 = alig2->position;
	  // computes distance between alignments,
	  // is a valid distance ?
	  if (start2 > end1) {
	    //  _____________   _____________  //
	    // [____ALIG1____] [____ALIG2____] //
	    if (alig2->map_len == 0) {
	      generate_alignment_len(alig2);	      
	    }
	    distance = (start2 + alig2->map_len) - end1;
	  } else {
	    //  _____________   _____________  //
	    // [____ALIG2____] [____ALIG1____] //
	    if (alig1->map_len == 0) {
	      generate_alignment_len(alig1);
	    }
	    distance = (end1 + alig1->map_len) - start2;
	  }
	  //distance = (start2 > end1 ? start2 - end1 : end1 - start2); // abs //
	  //printf("\t\tItem %i Pair2 [chr %i - start %i], cigar = %s, score = %0.2f\n", 
	  //	 j2, chr2, start2, alig1->cigar, alig2->map_quality);
	  //	  printf("\t\t\t*** chr1: %i == chr2: %i; str1: %i == str2: %i; distance = %lu, min_distance = %i, max_distance = %i, strand1 = %i, strand2 = %i, pair_mode = %i\n",  chr1, chr2, strand1, strand2, distance, min_distance, max_distance, strand1, strand2, pair_mode);

	  
	  if ( (chr1 == chr2) &&
	       (distance >= min_distance) && (distance <= max_distance)) { //&&
	    //((strand1 != strand2 && pair_mode == PAIRED_END_MODE) ||
	       //(strand1 == strand2 && pair_mode == MATE_PAIR_MODE )   ) ) {
	    //printf("pair!\n");
	    // order proper pairs by best score
	    // create the new pair
	    score = 0.5f * (alig1->map_quality + alig2->map_quality);
	    if (score > best_score) best_score = score;
	    //printf("\t\t\t\t-----> score = %0.2f\n", score);
	    new_pair = pair_new(j1, j2, score);
	    num_pairs++;
	    // insert the new pair in the correct position
	    // acording to its score
	    linked_list_iterator_first(pair_list_itr);
	    pair = (pair_t *) linked_list_iterator_curr(pair_list_itr);
	    while (pair != NULL) {
	      if (score > pair->score) {
		linked_list_iterator_insert(new_pair, pair_list_itr);
		linked_list_iterator_prev(pair_list_itr);
		break;
	      }
	      // continue loop...
	      linked_list_iterator_next(pair_list_itr);
	      pair = linked_list_iterator_curr(pair_list_itr);
	    }
	    if (pair == NULL) {
	      linked_list_insert_last(new_pair, pair_list);
	    }
	  }
	} // end for j2
      } // end for j1

      //      printf("***** num pairs = %i\n", num_pairs);
      /*
      if (num_pairs == 1) {
	linked_list_iterator_first(pair_list_itr);
	pair = (pair_t *) linked_list_iterator_curr(pair_list_itr);
	alig1 = (alignment_t *) array_list_get(pair->index1, list1);
	if (alig1->mapq == 0) {
	  alig1->mapq = 14;
	}
	alig2 = (alignment_t *) array_list_get(pair->index2, list2);
	if (alig2->mapq == 0) {
	  alig2->mapq = 14;
	}
      }
      */
      // compute number of mappings with the best score
      num_pairs = 0;
      linked_list_iterator_first(pair_list_itr);
      pair = (pair_t *) linked_list_iterator_curr(pair_list_itr);
      while (pair != NULL) {
	if (pair->score == best_score) {
	  num_pairs++;
	} else {
	  break;
	}
	// continue loop...
	linked_list_iterator_next(pair_list_itr);
	pair = linked_list_iterator_curr(pair_list_itr);
      }

      //      printf("***** num pairs with best-score = %i (best score = %0.2f)\n", num_pairs, best_score);
      
      //batch->counters[(num_pairs > 9) ? 9 : num_pairs]++;

      // filter pairs
      counter_hits = 0;
      linked_list_iterator_first(pair_list_itr);
      pair = (pair_t *) linked_list_iterator_curr(pair_list_itr);
      while (pair != NULL) {
	if (mapped1[pair->index1] == 0 && mapped2[pair->index2] == 0) {
	  mapped1[pair->index1] = 1;
	  mapped2[pair->index2] = 1;

	  mapped1_counter++;
	  mapped2_counter++;

	  alig1 = (alignment_t *) array_list_get(pair->index1, list1);
	  alig2 = (alignment_t *) array_list_get(pair->index2, list2);

	  // set pair1 fields
	  alig1->mate_position = alig2->position;
	  alig1->mate_chromosome = alig2->chromosome;
     
	  alig1->is_paired_end = 1;
	  alig1->is_paired_end_mapped = 1;
	  alig1->is_mate_mapped = 1;
	  alig1->mate_strand = alig2->seq_strand;
	  alig1->pair_num = 1;
	  
	  // set pair2 fields
	  alig2->mate_position = alig1->position;
	  alig2->mate_chromosome = alig1->chromosome;
	    
	  alig2->is_paired_end = 1;
	  alig2->is_paired_end_mapped = 1;
	  alig2->is_mate_mapped = 1;
	  alig2->mate_strand = alig1->seq_strand;
	  alig2->pair_num = 2;

	  if (alig1->position > alig2->position) {
	    //  _____________   _____________  //
	    // [____ALIG2____] [____ALIG1____] //
	    //Alig1->template_length -
	    //Alig2->template_length +
	    alig1->template_length = alig2->position - (alig1->position + alig1->map_len);
	    alig2->template_length = (alig1->position + alig1->map_len) - alig2->position;
	  } else {
	    //  _____________   _____________  //
	    // [____ALIG1____] [____ALIG2____] //
	    //Alig1->template_length +
	    //Alig2->template_length 
	    alig1->template_length = (alig2->position + alig2->map_len) - alig1->position;
	    alig2->template_length = alig1->position - (alig2->position + alig2->map_len);
	  }

	  if (num_pairs == 1) {
	    if (alig1->mapq == 0) {
	      alig1->mapq = alig2->mapq;
	    } else if (alig2->mapq == 0) {
	      alig2->mapq = alig1->mapq;
	    }
	    /*
	  } else if (num_pairs > 1) {
	    if (best_score != second_score) {
	      printf("----> %s\t%s num pairs = %i (best score, second)  = (%0.2f, %0.2f)\n\n", 
		     alig1->query_name, num_pairs, best_score, second_score);
	    }
	    */
	  }

	  if ( (++counter_hits) >= num_hits) {
	    /*
	    if (num_pairs == 1) {
	      alig1->map_quality = 3;
	      alig2->map_quality = 3;
	    } else if (num_pairs == 2) {
	      alig1->map_quality = 2;
	      alig2->map_quality = 2;
	    } else if (num_pairs > 2 && num_pairs < 9) {
	      alig1->map_quality = 1;
	      alig2->map_quality = 1;
	    } else {
	      alig1->map_quality = 0;
	      alig2->map_quality = 0;
	    }
	    */
	    break;
	  }
	}

	// continue loop...
	linked_list_iterator_next(pair_list_itr);
	pair = linked_list_iterator_curr(pair_list_itr);
      }

      linked_list_clear(pair_list, (void *) pair_free);

      //      printf("\n------------> counter_hits = %i\n", counter_hits);
      //      printf("(all = %i, n_best = %i, n_hits = %i, best = %i)\n", all, n_best, n_hits, best);

      // check if there are unproperly aligned pairs
      if (counter_hits) {
	// remove unpaired alignments and
	// report only pair alignments found
	if (mapped1_counter != num_items1) {
	  batch->mapping_lists[i] = create_list(mapped1, mapped1_counter, list1);
	}
	if (mapped2_counter != num_items2) {
	  batch->mapping_lists[i + 1] = create_list(mapped2, mapped2_counter, list2);
	}
      } else {
	// all aligments are unpaired
	//printf("Alignments are unpaired\n");
	if (!report_only_paired && (all || n_best || n_hits || best)) {
	  size_t num_items = num_items1 + num_items2;
	  // report all mappings 
	  if (all || num_items <= n_best || num_items <= n_hits) {
	    // report all mappings 
	    update_mispaired_pairs(num_items1, num_items2, list1, list2);
	  } else if (best) {
	    filter_alignments(0, 0, 0, 1, batch->mapping_lists[i]);
	  } else if (n_hits) {
	    //select n hits from first pair
	    filter_alignments(0, 0, n_hits, 0, batch->mapping_lists[i]);
	    update_mispaired_pair(1, array_list_size(batch->mapping_lists[i]), batch->mapping_lists[i]);
	    if (num_items1 < n_hits) {
	      size_t new_n_hits = n_hits - num_items1;
	      filter_alignments(0, 0, new_n_hits, 0, batch->mapping_lists[i + 1]);
	      update_mispaired_pair(2, array_list_size(batch->mapping_lists[i + 1]), batch->mapping_lists[i + 1]);
	    } else {
	      array_list_clear(batch->mapping_lists[i + 1], (void *) alignment_free);
	    }
	  } else if (n_best) {
	    update_mispaired_pairs(num_items1, num_items2, list1, list2);
	    for (int n = num_items2 - 1; n >= 0; n--) {
	      alig2 = array_list_remove_at(n, batch->mapping_lists[i + 1]);
	      array_list_insert(alig2, batch->mapping_lists[i]);
	    }
	    filter_alignments(0, n_best, 0, 0, batch->mapping_lists[i]);
	    size_t num_items = array_list_size(batch->mapping_lists[i]);
	    for (int n = num_items - 1; n >= 0; n--) {
	      alig1 = array_list_get(n, batch->mapping_lists[i]);
	      if (alig1->pair_num == 2) {
		alig1 = array_list_remove_at(n, batch->mapping_lists[i]);
		array_list_insert(alig1, batch->mapping_lists[i + 1]);
	      }
	    }
	  }

	  //------------------------------------------
	  // temp
	  {
	    alignment_t *alig;
	    array_list_t *list = batch->mapping_lists[i];
	    int num_items = array_list_size(list);
	    for (size_t ii = 0; ii < num_items; ii++) {
	      alig = (alignment_t *) array_list_get(ii, list);
	      if (num_items == 1) {
		alig->map_quality = 3;
	      } else if (num_items == 2) {
		alig->map_quality = 2;
	      } else if (num_items > 2 && num_items < 9) {
		alig->map_quality = 1;
	      } else {
		alig->map_quality = 0;
	      }
	    }
	    list = batch->mapping_lists[i + 1];
	    num_items = array_list_size(list);
	    for (size_t ii = 0; ii < num_items; ii++) {
	      alig = (alignment_t *) array_list_get(ii, list);
	      if (num_items == 1) {
		alig->map_quality = 3;
	      } else if (num_items == 2) {
		alig->map_quality = 2;
	      } else if (num_items > 2 && num_items < 9) {
		alig->map_quality = 1;
	      } else {
		alig->map_quality = 0;
	      }
	    }
	  }
	  // end of temp
	  //------------------------------------------

	} else {
	  array_list_clear(batch->mapping_lists[i], (void *) alignment_free);
	  array_list_clear(batch->mapping_lists[i + 1], (void *) alignment_free);
	}
      }
    } else {
      //batch->counters[0]++;

      //printf("This section\n");
      // pairs are not properly aligned, only one is mapped
      if (!report_only_paired) {
	// report all, n-best or n-hits
	array_list_t *list;
	int num_pair;
	//printf("num items1 = %i, num_items2 = %i\n", num_items1, num_items2);
	if (num_items1) {
	  list = batch->mapping_lists[i];
	  num_pair = 1;
	} else {
	  list = batch->mapping_lists[i + 1];
	  num_pair = 2;
	} 

	//filter_alignments(all, n_best, n_hits, best, list);
	update_mispaired_pair(num_pair, array_list_size(list), list);   
      } else {
	// no report_unpaired option set, delete all mappings found
	array_list_clear(batch->mapping_lists[i], (void *) alignment_free);
	array_list_clear(batch->mapping_lists[i + 1], (void *) alignment_free);
      }
    }
  } // end for num_reads

  // free memory
  free(mapped1);
  free(mapped2);
  linked_list_free(pair_list, (void *) pair_free);
  linked_list_iterator_free(pair_list_itr);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
