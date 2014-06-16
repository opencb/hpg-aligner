#include "pair_server.h"

//------------------------------------------------------------------------------------
// this functions writes reads to disk, reads came from
// a list
//------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------

void prepare_paired_alignments(pair_server_input_t *input, mapping_batch_t *batch);

//------------------------------------------------------------------------------------

void filter_alignments(char report_all,
			      size_t report_n_best, 
			      size_t report_n_hits,
			      int report_best,
			      array_list_t *mapping_list);

//unsigned long alignment_hash_code(void *p);
//int alignment_compare(void *p1, void *p2);

inline int remove_alignments(size_t *valid_items, array_list_t *list) {
  alignment_t *alig;
  int num = 0;

  size_t num_items = array_list_size(list);

  for (int k = num_items - 1; k >= 0; k--) {
    if (valid_items[k] == 0) {
      alig = (alignment_t *) array_list_remove_at(k, list);
      alignment_free(alig);
      ++num;
    }
  }
  return num;
}

//------------------------------------------------------------------------------------

void prepare_pair_server(pair_server_input_t* input) {
  list_item_t *pair_item;
  mapping_batch_t *batch;
  LOG_DEBUG_F("pair_server (%i): START\n", omp_get_thread_num());  
  //printfinput->pair_mng->pair_mode
  while ( (pair_item = list_remove_item(input->pair_list)) != NULL ) {
    //printf("------->Extract pair server item ...\n");
    batch = (mapping_batch_t *)pair_item->data_p;

    if (input->pair_mng->pair_mode != SINGLE_END_MODE) {
      prepare_paired_alignments(input, batch);
    }
    list_insert_item(pair_item, input->write_list);
  }
  list_decr_writers(input->write_list);
  LOG_DEBUG("pair_server : END\n");
}

//------------------------------------------------------------------------------------

array_list_t *create_new_list(size_t *valid_items, size_t num_valids, array_list_t *list) {
  void *item;
  int num = 0;

  size_t num_items = array_list_size(list);
  int flag = array_list_get_flag(list);

  array_list_t *new_list = array_list_new(num_valids, 
					  1.25f, 
					  COLLECTION_MODE_ASYNCHRONIZED);
  array_list_set_flag(flag, new_list);

  for (int k = 0; k < num_items; k++) {
    if (valid_items[k] == 1) {
      array_list_insert(array_list_get(k, list), new_list);
      array_list_set(k, NULL, list);
    }
  }

  if (flag == 1) {
    array_list_free(list, (void *) alignment_free);
  } else {
    array_list_free(list, (void *) cal_free);
  }

  return new_list;
}

static void prepare_single_alignments(pair_server_input_t *input, mapping_batch_t *batch) {

  static char aux[4096];

  char *p, *optional_fields;
  int optional_fields_length;
      
  array_list_t *fq_batch = batch->fq_batch;
  fastq_read_t *fq_read;

  int num_items, distance, AS;
  size_t read_index, num_targets = batch->num_targets;

  size_t read_len, header_len, mapped_len, len, pos;
  size_t deletion_n, num_cigar_ops, shift;
  char *cigar, *header_match, *read_match, *quality_match;

  size_t mref_start, mquery_start;
  alignment_t *alignment;

  sw_output_t *sw_output;
  array_list_t *sw_list, *alignment_list;

  int c;
  char *quality;

  // convert the SW output to alignments
  for (size_t i = 0; i < num_targets; i++) {
    //    printf("pair_server.c, prepare_single_alignments: target = #%i of %i\n", i, num_targets);
    read_index = batch->targets[i];
    fq_read = (fastq_read_t *) array_list_get(read_index, fq_batch);
    
    sw_list = batch->mapping_lists[read_index];
    
    read_len = fq_read->length;

    header_len = strlen(fq_read->id);
    num_items = array_list_size(sw_list);

    alignment_list = array_list_new(1000, 
				    1.25f, 
				    COLLECTION_MODE_ASYNCHRONIZED);
    array_list_set_flag(1, alignment_list);

    //    printf("pair_server.c, prepare_single_alignments: process read #%i with %i mappings\n", 
    //	   index, num_items);

    for (size_t j = 0; j < num_items; j++) {
      sw_output = (sw_output_t *) array_list_get(j, sw_list);
      
      mapped_len = sw_output->mref_len;

      header_len = get_to_first_blank(fq_read->id, header_len, aux);
      
      header_match = (char *) malloc(sizeof(char) * (header_len + 1));
      memcpy(header_match, aux, header_len);
      header_match[header_len] = '\0';
      
      mquery_start = sw_output->mquery_start;
      mref_start = sw_output->mref_start;

      //printf("%s\n", header_match);
      //printf("\tstrand = %i, chromosome = %i\n", sw_output->strand, sw_output->chromosome - 1);
      //printf("\tmquery_start = %i, mref_start = %i\n", mquery_start, mref_start);

      pos = sw_output->ref_start + mref_start - 1;

      read_match = (char *) malloc(sizeof(char) * (mapped_len + 1));
      quality_match = (char *) malloc(sizeof(char) * (mapped_len + 1));
      c = 0;
      quality = fq_read->quality + mquery_start;
      for (int ii = 0; ii < mapped_len; ii++){
	   if (sw_output->mquery[ii] != '-') {
		read_match[c] = sw_output->mquery[ii]; 
		quality_match[c] = quality[c];
		c++;
	   }
      }
      read_match[c] = '\0'; 
      quality_match[c] = '\0'; 

      cigar =  generate_cigar_str(sw_output->mquery, 
				  sw_output->mref, 
				  mquery_start, 
				  read_len, 
				  sw_output->mref_len, 
				  &distance,
				  (int *) &num_cigar_ops);

      /*
      printf(">>>>> %s : %s, pos = %lu\n", header_match, cigar, pos);
      printf("\tref_start = %lu, ref_len = %lu, mref_start = %lu ,mref_len = %lu\n",
	     sw_output->ref_start, sw_output->ref_len, mref_start, sw_output->mref_len);
      printf("\tmquery = %s (start = %i)\n", sw_output->mquery, mquery_start);
      printf("\tmref   = %s (start = %i)\n", sw_output->mref, mref_start);
      */
      // set optional fields
      //      optional_fields_length = 0;
      //      optional_fields = NULL;
      optional_fields_length = 100;
      optional_fields = (char *) calloc(optional_fields_length, sizeof(char));

      p = optional_fields;
      AS = (int) sw_output->score;

      sprintf(p, "ASi");
      p += 3;
      memcpy(p, &AS, sizeof(int));
      p += sizeof(int);
      
      sprintf(p, "NHi");
      p += 3;
      memcpy(p, &num_items, sizeof(int));
      p += sizeof(int);

      sprintf(p, "NMi");
      p += 3;
      memcpy(p, &distance, sizeof(int));
      p += sizeof(int);
      //      *p = '\0';


      optional_fields_length = p - optional_fields;
      // create the alignment and insert into the list
      alignment = alignment_new();
      alignment_init_single_end(header_match, read_match, quality_match, 
				sw_output->strand, 
				sw_output->chromosome - 1, 
				pos,
				cigar, num_cigar_ops, sw_output->norm_score * 254, 1, (num_items > 1),
				optional_fields_length, optional_fields, alignment);

      array_list_insert(alignment, alignment_list);

      // free memory (sw output)
      sw_output_free(sw_output);
    } // end for sw items

    // free the sw list, and update the mapping list with the alignment list
    array_list_free(sw_list, NULL);
    batch->mapping_lists[read_index] = alignment_list;
  } // end for targets
  //  printf("pair_server.c, prepare_single_alignments: Done\n");
}

//------------------------------------------------------------------------------------

static void prepare_single_alignments_bs(pair_server_input_t *input, mapping_batch_t *batch) {

  static char aux[4096];

  char *p, *optional_fields;
  int optional_fields_length;
      
  array_list_t *fq_batch = batch->fq_batch;
  fastq_read_t *fq_read;

  int num_items, distance, AS;
  size_t read_index, num_targets = batch->num_targets, num_targets2 = batch->num_targets2;

  size_t read_len, header_len, mapped_len, len, pos;
  size_t deletion_n, num_cigar_ops, shift;
  char *cigar, *header_match, *read_match, *quality_match;

  size_t mref_start, mquery_start;
  alignment_t *alignment;

  sw_output_t *sw_output;
  array_list_t *sw_list, *alignment_list;

  int c;
  char *quality;

  // convert the SW output to alignments
  for (size_t i = 0; i < num_targets; i++) {
    //    printf("pair_server.c, prepare_single_alignments: target = #%i of %i\n", i, num_targets);
    LOG_DEBUG_F("------------------- list 1 item %i (%i)...\n", i, num_targets2);
    read_index = batch->targets[i];
    fq_read = (fastq_read_t *) array_list_get(read_index, fq_batch);
    
    sw_list = batch->mapping_lists[read_index];
    
    read_len = fq_read->length;

    header_len = strlen(fq_read->id);
    num_items = array_list_size(sw_list);

    alignment_list = array_list_new(1000, 
				    1.25f, 
				    COLLECTION_MODE_ASYNCHRONIZED);
    array_list_set_flag(1, alignment_list);

    //    printf("pair_server.c, prepare_single_alignments: process read #%i with %i mappings\n", 
    //	   index, num_items);

    for (size_t j = 0; j < num_items; j++) {
      sw_output = (sw_output_t *) array_list_get(j, sw_list);
      
      mapped_len = sw_output->mref_len;

      header_len = get_to_first_blank(fq_read->id, header_len, aux);
      
      header_match = (char *) malloc(sizeof(char) * (header_len + 1));
      memcpy(header_match, aux, header_len);
      header_match[header_len] = '\0';
      
      mquery_start = sw_output->mquery_start;
      mref_start = sw_output->mref_start;

      //printf("%s\n", header_match);
      //printf("\tstrand = %i, chromosome = %i\n", sw_output->strand, sw_output->chromosome - 1);
      //printf("\tmquery_start = %i, mref_start = %i\n", mquery_start, mref_start);

      pos = sw_output->ref_start + mref_start - 1;

      read_match = (char *) malloc(sizeof(char) * (mapped_len + 1));
      quality_match = (char *) malloc(sizeof(char) * (mapped_len + 1));
      c = 0;
      quality = fq_read->quality + mquery_start;
      for (int ii = 0; ii < mapped_len; ii++){
	   if (sw_output->mquery[ii] != '-') {
		read_match[c] = sw_output->mquery[ii]; 
		quality_match[c] = quality[c];
		c++;
	   }
      }
      read_match[c] = '\0'; 
      quality_match[c] = '\0'; 

      cigar =  generate_cigar_str(sw_output->mquery, 
				  sw_output->mref, 
				  mquery_start, 
				  read_len, 
				  sw_output->mref_len, 
				  &distance,
				  (int *) &num_cigar_ops);

      /*
      printf(">>>>> %s : %s, pos = %lu\n", header_match, cigar, pos);
      printf("\tref_start = %lu, ref_len = %lu, mref_start = %lu ,mref_len = %lu\n",
	     sw_output->ref_start, sw_output->ref_len, mref_start, sw_output->mref_len);
      printf("\tmquery = %s (start = %i)\n", sw_output->mquery, mquery_start);
      printf("\tmref   = %s (start = %i)\n", sw_output->mref, mref_start);
      */

      // set optional fields
      //      optional_fields_length = 0;
      //      optional_fields = NULL;
      optional_fields_length = 100;
      optional_fields = (char *) calloc(optional_fields_length, sizeof(char));

      p = optional_fields;
      AS = (int) sw_output->score;

      sprintf(p, "ASi");
      p += 3;
      memcpy(p, &AS, sizeof(int));
      p += sizeof(int);
      
      sprintf(p, "NHi");
      p += 3;
      memcpy(p, &num_items, sizeof(int));
      p += sizeof(int);

      sprintf(p, "NMi");
      p += 3;
      memcpy(p, &distance, sizeof(int));
      p += sizeof(int);
      //      *p = '\0';


      optional_fields_length = p - optional_fields;
      // create the alignment and insert into the list
      alignment = alignment_new();
      alignment_init_single_end(header_match, read_match, quality_match, 
				sw_output->strand, 
				sw_output->chromosome - 1, 
				pos,
				cigar, num_cigar_ops, sw_output->norm_score * 254, 1, (num_items > 1),
				optional_fields_length, optional_fields, alignment);

      array_list_insert(alignment, alignment_list);

      // free memory (sw output)
      sw_output_free(sw_output);
    } // end for sw items

    // free the sw list, and update the mapping list with the alignment list
    array_list_free(sw_list, NULL);
    batch->mapping_lists[read_index] = alignment_list;
  } // end for targets

/////////////////////////////////////////////////
  // convert the SW output2 to alignments
  for (size_t i = 0; i < num_targets2; i++) {
    //    printf("pair_server.c, prepare_single_alignments: target = #%i of %i\n", i, num_targets);
    LOG_DEBUG_F("------------------- list 2 item %i (%i)...\n", i, num_targets2);
    read_index = batch->targets2[i];
    fq_read = (fastq_read_t *) array_list_get(read_index, fq_batch);
    
    sw_list = batch->mapping_lists2[read_index];
    
    read_len = fq_read->length;

    header_len = strlen(fq_read->id);
    num_items = array_list_size(sw_list);

    alignment_list = array_list_new(1000, 
				    1.25f, 
				    COLLECTION_MODE_ASYNCHRONIZED);
    array_list_set_flag(1, alignment_list);

    //    printf("pair_server.c, prepare_single_alignments: process read #%i with %i mappings\n", 
    //	   index, num_items);

    for (size_t j = 0; j < num_items; j++) {
      sw_output = (sw_output_t *) array_list_get(j, sw_list);
      
      mapped_len = sw_output->mref_len;

      header_len = get_to_first_blank(fq_read->id, header_len, aux);
      
      header_match = (char *) malloc(sizeof(char) * (header_len + 1));
      memcpy(header_match, aux, header_len);
      header_match[header_len] = '\0';
      
      mquery_start = sw_output->mquery_start;
      mref_start = sw_output->mref_start;

      //printf("%s\n", header_match);
      //printf("\tstrand = %i, chromosome = %i\n", sw_output->strand, sw_output->chromosome - 1);
      //printf("\tmquery_start = %i, mref_start = %i\n", mquery_start, mref_start);

      pos = sw_output->ref_start + mref_start - 1;

      read_match = (char *) malloc(sizeof(char) * (mapped_len + 1));
      quality_match = (char *) malloc(sizeof(char) * (mapped_len + 1));
      c = 0;
      quality = fq_read->quality + mquery_start;
      for (int ii = 0; ii < mapped_len; ii++){
	   if (sw_output->mquery[ii] != '-') {
		read_match[c] = sw_output->mquery[ii]; 
		quality_match[c] = quality[c];
		c++;
	   }
      }
      read_match[c] = '\0'; 
      quality_match[c] = '\0'; 

      cigar =  generate_cigar_str(sw_output->mquery, 
				  sw_output->mref, 
				  mquery_start, 
				  read_len, 
				  sw_output->mref_len, 
				  &distance,
				  (int *) &num_cigar_ops);
      /*
      printf(">>>>> %s : %s, pos = %lu\n", header_match, cigar, pos);
      printf("\tref_start = %lu, ref_len = %lu, mref_start = %lu ,mref_len = %lu\n",
	     sw_output->ref_start, sw_output->ref_len, mref_start, sw_output->mref_len);
      printf("\tmquery = %s (start = %i)\n", sw_output->mquery, mquery_start);
      printf("\tmref   = %s (start = %i)\n", sw_output->mref, mref_start);
      */
      // set optional fields
      //      optional_fields_length = 0;
      //      optional_fields = NULL;
      optional_fields_length = 100;
      optional_fields = (char *) calloc(optional_fields_length, sizeof(char));

      p = optional_fields;
      AS = (int) sw_output->score;

      sprintf(p, "ASi");
      p += 3;
      memcpy(p, &AS, sizeof(int));
      p += sizeof(int);
      
      sprintf(p, "NHi");
      p += 3;
      memcpy(p, &num_items, sizeof(int));
      p += sizeof(int);

      sprintf(p, "NMi");
      p += 3;
      memcpy(p, &distance, sizeof(int));
      p += sizeof(int);
      //      *p = '\0';


      optional_fields_length = p - optional_fields;
      // create the alignment and insert into the list
      alignment = alignment_new();
      alignment_init_single_end(header_match, read_match, quality_match, 
				sw_output->strand, 
				sw_output->chromosome - 1, 
				pos,
				cigar, num_cigar_ops, sw_output->norm_score * 254, 1, (num_items > 1),
				optional_fields_length, optional_fields, alignment);

      array_list_insert(alignment, alignment_list);

      // free memory (sw output)
      sw_output_free(sw_output);
    } // end for sw items

    // free the sw list, and update the mapping list with the alignment list
    array_list_free(sw_list, NULL);
    batch->mapping_lists2[read_index] = alignment_list;
  } // end for targets
  //printf("pair_server.c, prepare_single_alignments: Done\n");
}

//------------------------------------------------------------------------------------

void update_mispaired_pair(int pair_num, size_t num_items, array_list_t *list) {
  alignment_t *alig;

  for (size_t i = 0; i < num_items; i++) {
    alig = (alignment_t *) array_list_get(i, list);

    // set pair fields
    alig->mate_position = 0;
    alig->mate_chromosome = 0;
    alig->template_length = 0;
     
    alig->is_paired_end = 1;
    alig->is_paired_end_mapped = 0;
    alig->is_mate_mapped = 0;
    alig->mate_strand = 0;
    alig->pair_num = pair_num;

    // tmp
    if (num_items == 1) {
      alig->map_quality = 3;
    } else if (num_items == 2) {
      alig->map_quality = 2;
    } else if (num_items > 2 && num_items < 9) {
      alig->map_quality = 1;
    } else {
      alig->map_quality = 0;
    }


//    printf("update_mispaired_pair = %i\n", alig->template_length);
  }
}

//------------------------------------------------------------------------------------

void update_mispaired_alignment(int pair_num, alignment_t *alig) {
    // set pair fields
    alig->mate_position = 0;
    alig->mate_chromosome = 0;
    alig->template_length = 0;
     
    alig->is_paired_end = 1;
    alig->is_paired_end_mapped = 0;
    alig->is_mate_mapped = 0;
    alig->mate_strand = 0;
    alig->pair_num = pair_num;

//    printf("update_mispaired_alignment = %i\n", alig->template_length);
}

//------------------------------------------------------------------------------------

void update_mispaired_pairs(size_t num_items1, size_t num_items2,
			    array_list_t *list1, array_list_t *list2) {
  alignment_t *alig;
  alignment_t *first1 = array_list_get(0, list1);
  alignment_t *first2 = array_list_get(0, list2);

  for (size_t i = 0; i < num_items1; i++) {
    alig = (alignment_t *) array_list_get(i, list1);

    // set pair1 fields
    alig->mate_position = first2->position;
    //alig->mate_position = 0;

    alig->mate_chromosome = first2->chromosome;
    //alig->mate_chromosome = 0;
    
    
    alig->template_length = 0;
    //alig->template_length = first2->position - alig->position;;
     
    alig->is_paired_end = 1;
    alig->is_paired_end_mapped = 0;
    alig->is_mate_mapped = 1;
    alig->mate_strand = first2->seq_strand;
    alig->pair_num = 1;

//    printf("update_mispaired_pairs, set pair1 fields, distance  = %i\n", alig->template_length);

  }

  for (size_t i = 0; i < num_items2; i++) {
    alig = (alignment_t *) array_list_get(i, list2);

    // set pair2 fields
    alig->mate_position = first1->position;
    //alig->mate_position = 0;

    alig->mate_chromosome = first1->chromosome;
    //alig->mate_chromosome = 0;

    //alig->template_length = first1->position - alig->position;
    alig->template_length = 0;
     
    alig->is_paired_end = 1;
    alig->is_paired_end_mapped = 0;
    alig->is_mate_mapped = 1;
    alig->mate_strand = first1->seq_strand;
    alig->pair_num = 2;

//    printf("update_mispaired_pairs, set pair2 fields, distance  = %i\n", alig->template_length);

  }
}

//------------------------------------------------------------------------------------

void prepare_paired_alignments(pair_server_input_t *input, mapping_batch_t *batch) {

  size_t num_items1, num_items2, num_reads = array_list_size(batch->fq_batch);

  int distance;
  int min_distance = input->pair_mng->min_distance;
  int max_distance = input->pair_mng->max_distance;
  int pair_mode = input->pair_mng->pair_mode;
  int report_only_paired = input->pair_mng->report_only_paired;

  array_list_t *list1, *list2;

  alignment_t *alig1, *alig2;
  size_t mapped1_counter = 0, mapped2_counter = 0;
  size_t allocated_mapped1 = 100, allocated_mapped2 = 100;
  size_t *mapped1 = (size_t *) malloc(allocated_mapped1 * sizeof(size_t));
  size_t *mapped2 = (size_t *) malloc(allocated_mapped2 * sizeof(size_t));

  short int chr1, chr2, strand1, strand2;
  size_t end1, start2;

  int pair_found;

  float score;
  pair_t *pair, *new_pair;
  int num_hits, counter_hits;
  linked_list_t *pair_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
  linked_list_iterator_t *pair_list_itr = linked_list_iterator_new(pair_list);
  int n_best = input->report_optarg->n_best;
  int n_hits = input->report_optarg->n_hits;
  int all = input->report_optarg->all;
  int best = input->report_optarg->best;

  if (n_best) {
    num_hits = n_best;
  } else  if (n_hits) {
    num_hits = n_hits;
  } else {
    num_hits = 10000;
  }

  for (int i = 0; i < num_reads; i += 2) {
    //fastq_read_t *read = array_list_get(i, batch->fq_batch);
    //printf("%s\n", read->sequence);
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
	//printf("Item %i Pair1 [chr %i - start %i]\n", j1, chr1, end1);	
	for (size_t j2 = 0; j2 < num_items2; j2++) {
	  //if (mapped2[j2] == 1) continue;
	  alig2 = (alignment_t *) array_list_get(j2, list2);
	  chr2 = alig2->chromosome;
	  strand2 = alig2->seq_strand;
	  start2 = alig2->position;
	  //printf("Item Pair2 %i [chr %i - start %i]\n", j2, chr2, start2);
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
	  
	  //printf("*** chr1: %i == chr2: %i; str1: %i == str2: %i; distance = %lu, min_distance = %i, max_distance = %i, strand1 = %i, strand2 = %i, pair_mode = %i\n",  chr1, chr2, strand1, strand2, distance, min_distance, max_distance, strand1, strand2, pair_mode);

	  if ( (chr1 == chr2) &&
	       (distance >= min_distance) && (distance <= max_distance)) { //&&
	    //((strand1 != strand2 && pair_mode == PAIRED_END_MODE) ||
	       //(strand1 == strand2 && pair_mode == MATE_PAIR_MODE )   ) ) {
	    //printf("pair!\n");
	    // order proper pairs by best score
	    // create the new pair
	    score = 0.5f * (alig1->map_quality + alig2->map_quality);
	    new_pair = pair_new(j1, j2, score);
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

	  if ( (++counter_hits) >= num_hits) {
	    break;
	  }
	}

	// continue loop...
	linked_list_iterator_next(pair_list_itr);
	pair = linked_list_iterator_curr(pair_list_itr);
      }

      linked_list_clear(pair_list, (void *) pair_free);

      // check if there are unproperly aligned pairs
      if (counter_hits) {
	/*if (report_unpaired) {
	  // report all alignments
	  for (size_t j = 0; j < num_items1; j++) {
	    if (!mapped1[j]) {
	      alig1 = (alignment_t *) array_list_get(j, list1);
	      update_mispaired_alignment(1, alig1);
	    }
	  } 
	  
	  for (size_t j = 0; j < num_items2; j++) {
	    if (!mapped2[j]) {
	      alig2 = (alignment_t *) array_list_get(j, list2);
	      update_mispaired_alignment(2, alig2);
	    }  
	  }
	*/
	//} else {
	// remove unpaired alignments and
	// report only pair alignments found
	if (mapped1_counter != num_items1) {
	  batch->mapping_lists[i] = create_new_list(mapped1, mapped1_counter, list1);
	}
	if (mapped2_counter != num_items2) {
	  batch->mapping_lists[i + 1] = create_new_list(mapped2, mapped2_counter, list2);
	}
      } else {
	// all aligments are unpaired
	//printf("Alignments are unpaired\n");
	if (!report_only_paired && (all || n_best || n_hits || best)) {
	  size_t num_items = num_items1 + num_items2;
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
	} else {
	  array_list_clear(batch->mapping_lists[i], (void *) alignment_free);
	  array_list_clear(batch->mapping_lists[i + 1], (void *) alignment_free);
	}
      }
    } else {
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
	filter_alignments(all, n_best, n_hits, best, list);
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

//------------------------------------------------------------------------------------

static inline size_t select_n_hits(array_list_t *mapping_list,
			    size_t report_n_hits) {

  array_list_t *mapping_list_filter;
  alignment_t *aux_alignment;
  int i;
  size_t num_mappings = array_list_size(mapping_list);
  
  mapping_list_filter = array_list_new(num_mappings + 1, 
				       1.25f, 
				       COLLECTION_MODE_ASYNCHRONIZED);
  
  for (i = num_mappings - 1; i >= report_n_hits; i--) {
    //printf("Remove... \n");
    aux_alignment = array_list_remove_at(i, mapping_list);
    alignment_free(aux_alignment);
    //printf("Remove ok!\n");
  }

  return array_list_size(mapping_list);
}

//-----------------------------------------------------------------------------
extern unsigned int alignmentcmp(alignment_t *alignment_1, alignment_t *alignment_2);

static inline size_t select_best_hits(array_list_t *mapping_list,
			       size_t report_n_best) {
  int j, i, z;
  size_t best_pos, array_size;
  alignment_t *best_alignment, *aux_alignment;
  array_list_t *mapping_list_filter;
  size_t num_mappings = array_list_size(mapping_list);
  int primary_delete = 0;

  mapping_list_filter = array_list_new(num_mappings + 1, 
				       1.25f, 
				       COLLECTION_MODE_ASYNCHRONIZED);
    
  //  LOG_DEBUG_F("Initial array size %i\n", array_list_size(mapping_list));
  for (j = 0; j < report_n_best; j++) {     
    best_pos = 0;
    best_alignment = array_list_get( 0, mapping_list);
    array_size = array_list_size(mapping_list);
    for (i = 1; i < array_size; i++) {
      aux_alignment = array_list_get(i, mapping_list);
      if (alignmentcmp(best_alignment, aux_alignment) == 2) {
	best_alignment = aux_alignment;
	best_pos = i;
      }
    }
    //    printf("Remove item %i\n", best_pos);
    array_list_insert(array_list_remove_at(best_pos, mapping_list), 
		      mapping_list_filter);
  }
  
  // free all mapings discarded
  array_size = array_list_size(mapping_list);
  //  printf("Array Size %i\n", array_size);
  for (j = array_size - 1; j >= 0; j--) {
    aux_alignment = array_list_remove_at(j, mapping_list);
    if (!is_secondary_alignment(aux_alignment)) { primary_delete = 1; } 
    alignment_free(aux_alignment);
  }

  //  LOG_DEBUG_F("Move %i elements to list with %i elements\n", 
  //	      array_list_size(mapping_list_filter), array_list_size(mapping_list));

  for (j = report_n_best - 1; j >= 0; j--) {
    aux_alignment = array_list_remove_at(j, mapping_list_filter);
    array_list_insert(aux_alignment, mapping_list);
  }
  
  array_list_free(mapping_list_filter, NULL);

  if (primary_delete) {
    aux_alignment = array_list_get(0, mapping_list);
    if (aux_alignment) {
      set_secondary_alignment(0, aux_alignment);
    }
  }

  //  LOG_DEBUG_F("Final array size %i\n", array_list_size(mapping_list));
  return array_list_size(mapping_list);
}

//-----------------------------------------------------------------------------

static inline void select_best (array_list_t *mapping_list) {
  size_t num_items = array_list_size(mapping_list);
  int max_quality = 0;
  alignment_t *alig;
  int num_best = 0;

  //printf("===============================================================================\n");
  for (int i = 0; i < num_items; i++) {
    alig = array_list_get(i, mapping_list);
    if (alig->map_quality >= max_quality) {
      num_best++;
      max_quality = alig->map_quality;
    }
    //alignment_print(alig);
  }
  //printf("===============================================================================\n");

  //int dbg = 0;

  //printf("items = %i\tmax quality = %i\n", num_items, max_quality);

  for (int i = num_items - 1; i >= 0; i--) {
    alig = array_list_get(i, mapping_list);
    //printf("elem = %i\tquality = %i\n", i, alig->map_quality);
    if (alig->map_quality != max_quality) {
      //printf("borrar item = %i\n", i);
      alig = array_list_remove_at(i, mapping_list);
      //printf("item %i borrado\n", i);
      alignment_free(alig);
      //printf("item %i liberado\n", i);
    } /*else {
      dbg++;
      }*/
  }

  /*
  if (dbg == 2) {
    //printf("===============================================================================\n");
    for (int i = 0; i < array_list_size(mapping_list); i++) {
      alignment_t *alig = array_list_get(i, mapping_list);
      //alignment_print(alig);
    }
    //printf("===============================================================================\n");
  }*/

}

//-----------------------------------------------------------------------------

void filter_alignments(char report_all,
		       size_t report_n_best, 
		       size_t report_n_hits,
		       int report_best,
		       array_list_t *mapping_list) {

  size_t num_mappings = array_list_size(mapping_list);
  //printf("filter\n");
  //printf("\treport all = '%c'\nreport n best = %lu\nreport hits = %lu\nreport best = %i\n", report_all, report_n_best, report_n_hits, report_best);

  if (!report_all && num_mappings) {
    //printf("select reports\n");
    if (report_best) {
      //printf("report best\n");
      // max-score
      select_best(mapping_list);
      //printf("end select best\n");
    } else if (report_n_best > 0) {
      //printf("report n best\n");
      // n-best
      if (num_mappings > report_n_best) { 
	//printf("select n best\n");
	select_best_hits(mapping_list, report_n_best);
      }
    } else if (report_n_hits > 0) {
      //printf("report n hits\n");
      // n-hits
      if (num_mappings > report_n_hits) { 	
	//printf("select n hits\n");
	select_n_hits(mapping_list, report_n_hits);        
      }
    }
  }
  //printf("end filter\n");
}

//-----------------------------------------------------------------------------

void filter_alignments_lists(char report_all, 
			     size_t report_n_best, 
			     size_t report_n_hits,
			     int report_best,
			     size_t num_lists,
			     array_list_t **mapping_lists) {
  array_list_t *mapping_list;
  //printf("num_list = %lu\n", num_lists);
  //printf("report all = '%c'\nreport best = %lu\nreport hits = %lu\n", report_all, report_n_best, report_n_hits);
  for (int i = 0; i < num_lists; i++) {
    //printf("%i\n", i);
    mapping_list = mapping_lists[i];
    filter_alignments(report_all, 
		      report_n_best, 
		      report_n_hits,
		      report_best,
		      mapping_list);
  }
  //printf("\n");
}

//====================================================================================
// main functions: apply pair and prperare alignments
//====================================================================================

int apply_pair(pair_server_input_t* input, batch_t *batch) {

  printf("START: apply_pair\n"); 
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  char *seq;
  list_t *list = NULL;
  array_list_t *fq_batch = mapping_batch->fq_batch;

  int pair_mode = input->pair_mng->pair_mode;
  size_t min_distance = input->pair_mng->min_distance;
  size_t max_distance = input->pair_mng->max_distance;
  int distance;

  //  size_t num_targets = batch->num_targets;
  size_t num_items1, num_items2, num_reads = array_list_size(fq_batch);

  int flag1, flag2;
  array_list_t *list1, *list2;

  size_t end1, start2;
  short int chr1, chr2, strand1, strand2;

  size_t mapped1_counter = 0, mapped2_counter = 0;
  size_t allocated_mapped1 = 100, allocated_mapped2 = 100;
  size_t *mapped1 = (size_t *) malloc(allocated_mapped1 * sizeof(size_t));
  size_t *mapped2 = (size_t *) malloc(allocated_mapped2 * sizeof(size_t));

  int pair_found;

  alignment_t *alig;
  cal_t *cal;
  
  int total_removed = 0, num_to_do = 0;

  for (size_t i = 0; i < num_reads; i += 2) {

    list1 = mapping_batch->mapping_lists[i];
    list2 = mapping_batch->mapping_lists[i + 1];

    flag1 = array_list_get_flag(list1);
    flag2 = array_list_get_flag(list2);

    num_items1 = 0;
    if (list1 != NULL)  num_items1 = array_list_size(list1);
    num_items2 = 0;
    if (list2 != NULL) num_items2 = array_list_size(list2);

    if (num_items1 > 1 && num_items2 > 1) {

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

      pair_found = 0;

      // search for pairs properly aligned
      for (size_t j1 = 0; j1 < num_items1; j1++) {

	if (flag1 == 1) {
	  alig = (alignment_t *) array_list_get(j1, list1);
	  chr1 = alig->chromosome;
	  strand1 = alig->seq_strand;
	  end1 = alig->position + strlen(alig->sequence);
	} else if (flag1 == 2) {
	  cal = (cal_t *) array_list_get(j1, list1);
	  chr1 = cal->chromosome_id - 1;
	  strand1 = cal->strand;
	  end1 = cal->end;
	} else {
	  printf("Error in pair_server.c, apply_sw function (pair1)\n");
	  abort();
	}
	
	//for (size_t j2 = num_items2 - 1; j2 > 0; j2--) {
	for (size_t j2 = 0; j2 < num_items2; j2++) {
	  if (mapped2[j2] == 1) continue;

	  if (flag2 == 1) {
	    alig = (alignment_t *) array_list_get(j2, list2);
	    chr2 = alig->chromosome;
	    strand2 = alig->seq_strand;
	    start2 = alig->position;
	  } else if (flag2 == 2) {
	    cal = (cal_t *) array_list_get(j2, list2);
	    chr2 = cal->chromosome_id - 1;
	    strand2 = cal->strand;
	    start2 = cal->start;
	  } else {
	    printf("Error in pair_server.c, apply_sw function (pair2)\n");
	    abort();
	  }

	  // computes distance between alignments,
	  // is a valid distance ?
	  printf("strand = %i, strand = %i\n", strand1, strand2);
	  distance = (start2 > end1 ? start2 - end1 : end1 - start2); // abs
	  if ( (chr1 == chr2) &&
	       (distance >= min_distance) && (distance <= max_distance)) { //&&
	       //((strand1 != strand2 && pair_mode == PAIRED_END_MODE) ||
	       //(strand1 == strand2 && pair_mode == MATE_PAIR_MODE )   ) ) {
	    printf("Pair!!\n");
	    mapped1[j1] = 1;
	    mapped2[j2] = 1;
	    mapped1_counter++;
	    mapped2_counter++;

	    pair_found = 1;
	    break;
	  }
	} // end for j2..num_items2
      } // end for j1..num_item1


      if (pair_found) {
	// removing no valid items
	if (mapped1_counter != num_items1) {
	  mapping_batch->mapping_lists[i] = create_new_list(mapped1, mapped1_counter, list1);
	}
	if (mapped2_counter != num_items2) {
	  mapping_batch->mapping_lists[i + 1] = create_new_list(mapped2, mapped2_counter, list2);
	}
      }      
    }
  }

  // free memory
  free(mapped1);
  free(mapped2);

  // go to the next stage
  if (batch->mapping_batch->num_targets > 0) {
    return SW_STAGE;
  }

  return DNA_POST_PAIR_STAGE;
}

//------------------------------------------------------------------------------------

int prepare_alignments(pair_server_input_t *input, batch_t *batch) {

  if (batch->mapping_mode == DNA_MODE) {
    // convert Smith-Waterman output objects to alignments
    //    prepare_single_alignments(input, batch->mapping_batch);
  }

  if (input->pair_mng->pair_mode == SINGLE_END_MODE) {
    // filter alignments by best-alignments, n-hits, all and unpaired reads
    filter_alignments_lists(input->report_optarg->all, 
			    input->report_optarg->n_best, 
			    input->report_optarg->n_hits,
			    input->report_optarg->best,
			    array_list_size(batch->mapping_batch->fq_batch),
			    batch->mapping_batch->mapping_lists);
  } else {
    // first, search for proper pairs and
    // then filter paired alignments by best-alignments, n-hits, all and unpaired reads
    prepare_paired_alignments(input, batch->mapping_batch);
  }

  //*********************************************************

  return CONSUMER_STAGE;

  //printf("pair_server.c: prepare_alignments done (pair mode = %i)\n", input->pair_mng->pair_mode);
  //  printf("pair_server.c: 1: after prepare_single_alignments\n");
}

//------------------------------------------------------------------------------------

int prepare_alignments_bs(pair_server_input_t *input, batch_t *batch) {

  //printf("---------->Prepare alignments\n");
  if (batch->mapping_mode == BS_MODE) {
    //printf("Convert\n");
    // convert Smith-Waterman output objects to alignments
    //prepare_single_alignments_bs(input, batch->mapping_batch);
    //printf("\tEnd Convert\n");
  }

  if (input->pair_mng->pair_mode == SINGLE_END_MODE) {
    //printf("single end\n");
    // filter alignments by best-alignments, n-hits, all and unpaired reads

    filter_alignments_lists(input->report_optarg->all, 
			    input->report_optarg->n_best, 
			    input->report_optarg->n_hits,
			    input->report_optarg->best,
			    array_list_size(batch->mapping_batch->fq_batch),
			    batch->mapping_batch->mapping_lists);

    //printf("  single end\n");

    filter_alignments_lists(input->report_optarg->all, 
			    input->report_optarg->n_best, 
			    input->report_optarg->n_hits,
			    input->report_optarg->best,
			    array_list_size(batch->mapping_batch->fq_batch),
			    batch->mapping_batch->mapping_lists2);

    //printf("    end single end\n");
  } else {
    // first, search for proper pairs and
    // then filter paired alignments by best-alignments, n-hits, all and unpaired reads
    prepare_paired_alignments(input, batch->mapping_batch);
  }

  //*********************************************************

  //printf("---------->pair_server.c: prepare_alignments done (pair mode = %i)\n", input->pair_mng->pair_mode);
  //  printf("pair_server.c: 1: after prepare_single_alignments\n");

  // return for the analysis without the postprocess
  return CONSUMER_STAGE;
  // return for the analysis without the postprocess
  //return BS_STATUS_STAGE;
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

void pair_server_input_init(pair_mng_t *pair_mng, report_optarg_t *report_optarg, 
			    list_t* pair_list, list_t *sw_list,
			    list_t *write_list, pair_server_input_t* input) {

  input->report_optarg = report_optarg;

  input->pair_mng = pair_mng;
  input->pair_list = pair_list;
  input->sw_list = sw_list;
  input->write_list = write_list;
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
