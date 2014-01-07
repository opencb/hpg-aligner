#include "sa_dna_commons.h"


//--------------------------------------------------------------------
// commons
//--------------------------------------------------------------------

void seed_free(seed_t *p) {
  if (p) free(p);
}

//--------------------------------------------------------------------

void seed_cal_free(seed_cal_t *p) {
  if (p) {
    if (p->seed_list) linked_list_free(p->seed_list, (void *) seed_free);
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

float get_max_score(array_list_t *cal_list) {
  seed_cal_t *cal;
  size_t num_cals = array_list_size(cal_list);
  float max_score = -1000000.0f;
  for (size_t j = 0; j < num_cals; j++) {
    cal = array_list_get(j, cal_list);
    if (cal->score > max_score) {
      max_score = cal->score;
    }
  }
  return max_score;
}

//--------------------------------------------------------------------

int get_min_num_mismatches(array_list_t *cal_list) {
  seed_cal_t *cal;
  size_t num_cals = array_list_size(cal_list);
  int min_num_mismatches = 100000;
  for (size_t j = 0; j < num_cals; j++) {
    cal = array_list_get(j, cal_list);
    if (cal->num_mismatches < min_num_mismatches) {
      min_num_mismatches = cal->num_mismatches;
    }
  }
  return min_num_mismatches;
}

//--------------------------------------------------------------------

int get_max_read_area(array_list_t *cal_list) {
  seed_cal_t *cal;
  size_t num_cals = array_list_size(cal_list);
  int max_read_area = 0;
  for (size_t j = 0; j < num_cals; j++) {
    cal = array_list_get(j, cal_list);
    if (cal->read_area > max_read_area) {
      max_read_area = cal->read_area;
    }
  }

  return max_read_area;
}

//--------------------------------------------------------------------

void filter_cals_by_max_score(float score, array_list_t **list) {
  seed_cal_t *cal;
  array_list_t *cal_list = *list;
  size_t num_cals = array_list_size(cal_list);
  array_list_t *new_cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  for (size_t j = 0; j < num_cals; j++) {
    cal = array_list_get(j, cal_list);
    if (cal->score >= score) {
      array_list_insert(cal, new_cal_list);
      array_list_set(j, NULL, cal_list);
    }
  }
  array_list_free(cal_list, (void *) seed_cal_free);
  *list = new_cal_list;
}

//--------------------------------------------------------------------

void filter_cals_by_max_read_area(int read_area, array_list_t **list) {
  seed_cal_t *cal;
  array_list_t *cal_list = *list;
  size_t num_cals = array_list_size(cal_list);
  array_list_t *new_cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  for (size_t j = 0; j < num_cals; j++) {
    cal = array_list_get(j, cal_list);
    if (cal->read_area >= read_area) {
      array_list_insert(cal, new_cal_list);
      array_list_set(j, NULL, cal_list);
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
    if (cal->read_area <= read_area) {
      array_list_insert(cal, new_cal_list);
      array_list_set(j, NULL, cal_list);
    }
  }
  array_list_free(cal_list, (void *) seed_cal_free);
  *list = new_cal_list;
}

//--------------------------------------------------------------------

void filter_cals_by_max_num_mismatches(int num_mismatches, array_list_t **list) {
  seed_cal_t *cal;
  array_list_t *cal_list = *list;
  size_t num_cals = array_list_size(cal_list);
  array_list_t *new_cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  for (size_t j = 0; j < num_cals; j++) {
    cal = array_list_get(j, cal_list);
    if (cal->num_mismatches <= num_mismatches) {
      array_list_insert(cal, new_cal_list);
      array_list_set(j, NULL, cal_list);
    }
  }
  array_list_free(cal_list, (void *) seed_cal_free);
  *list = new_cal_list;
}

//--------------------------------------------------------------------

void create_alignments(array_list_t *cal_list, fastq_read_t *read, 
		       array_list_t *mapping_list) {

  // CAL
  seed_cal_t *cal;
  uint num_cals = array_list_size(cal_list);

  // alignments
  char cigar[1000];
  int num_cigar_ops;
  alignment_t *alignment;

  if (num_cals <= 0) {
    // no CALs -> no alignment
    return;
  }

  int AS;
  size_t i, pos;
  linked_list_item_t *list_item; 
  seed_t *s_first, *s_last;

  for (i = 0; i < num_cals; i++) {
    cigar[0] = 0;
    cal = array_list_get(i, cal_list);

    #ifdef _VERBOSE	  
    printf("--> CAL #%i (cigar %s):\n", i, cigar_to_string(&cal->cigar));
    cal_print(cal);
    #endif

    // computing aligments
    AS = read->length * 5;
    s_first = linked_list_get_first(cal->seed_list);
    s_last = linked_list_get_last(cal->seed_list);

    // sanity checking
    //    assert(s_first != NULL && s_last != NULL);

    if (s_first != NULL && s_last != NULL) {
      //      printf("-----> cigar_code = %x\n", s_first->info);
      //      printf("-----> cigar_code = %x\n", s_last->info);
      /*
      if (num_cals > 1) {
	printf("\t\t%s\t%i of %i\t[%i|%i - %i|%i]\n", 
	       read->id, i, num_cals, 
	       s_first->genome_start, s_first->read_start, s_first->read_end, s_first->genome_end);
	/*
	cigar_code_t *cigar = (cigar_code_t *) cal->info;
	printf("\t\t%s\t%i of %i\t[%i|%i - %i|%i] (cigar %s, distance = %i)\n", 
	       read->id, i, num_cals, 
	       s_first->genome_start, s_first->read_start, s_first->read_end, s_first->genome_end,
	       (cigar != NULL ? new_cigar_code_string(cigar) : "none"), 
	       (cigar != NULL ? cigar->distance : 0));
	*/
      //      }

      num_cigar_ops = 0;
      pos = cal->start;
      if (s_first->read_start > 0) {
	pos -= s_first->read_start;
	//sprintf(cigar, "%s%iS", cigar, s_first->read_start);
	//      num_cigar_ops++;
      }
      sprintf(cigar, "%s%iM", cigar, read->length);
      num_cigar_ops++;
      /*
      s_first->read_start = 0; // be carefull !!! to remove...now !!!
      sprintf(cigar, "%s%iM", cigar, s_last->read_end - s_first->read_start + 1);
      num_cigar_ops++;
      if (s_last->read_end < read->length - 1) {
	sprintf(cigar, "%s%iS", cigar, read->length - s_last->read_end - 1);
	num_cigar_ops++;
      }
      */
      alignment = alignment_new();	       
      alignment_init_single_end(strdup(read->id), strdup(read->sequence), strdup(read->quality), 
				cal->strand, cal->chromosome_id, 
				pos,
				strdup(cigar), num_cigar_ops, (AS * 254) / (read->length * 5), 1, (num_cals > 1),
				0, 0, alignment);  
      
      array_list_insert(alignment, mapping_list);
      //      array_list_insert(convert_to_bam(alignment, 33), mapping_list);
      
      // free memory
      //      alignment_free(alignment);
    }
    
    seed_cal_free(cal);
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
  printf("%s%c:%i[%lu|%lu - %lu|%lu] (cigar: %s, num. mismatches = %i)\n",  msg, (s->strand == 0 ? '+' : '-'),
	 s->chromosome_id, s->genome_start, s->read_start, s->read_end, s->genome_end,
	 cigar_to_string(&s->cigar), s->num_mismatches);
}

//--------------------------------------------------------------------

void display_sequence(uint j, sa_index3_t *index, uint len) {
  char *p = &index->genome->S[index->SA[j]];
  char chrom = index->CHROM[j];
  for (int i = 0; i < len; i++) {
    printf("%c", *p);
    p++;
  }
  printf("\t%u\t%s:%u\n", index->SA[j], 
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

void display_cmp_sequences(fastq_read_t *read, char *revcomp_seq, sa_index3_t *sa_index) {
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
  
  printf("\n\n======> %s\n", read->id);
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
    seq = revcomp_seq;
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
//--------------------------------------------------------------------
