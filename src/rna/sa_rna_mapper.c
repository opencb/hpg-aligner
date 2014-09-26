#include "sa_rna_mapper.h"

#define MAX_DEPTH 4
//#define DEBUG 1

extern int min_intron, max_intron;

int min_score;
int MAX_ALIG;

int get_score (cigar_code_t *cc, fastq_read_t *read, float match, float mismatch, float gap_open, float gap_extend) {
  //Smith-Waterman score
  int tot_match = read->length - cc->distance;
  int tot_gap_open = 0;
  int tot_indel_extend = 0;
  int tot_indel = 0;
  for (int i = 0; i < cc->ops->size; i++) {
    cigar_op_t *op = array_list_get(i, cc->ops);
    if (op->name == 'I' || op->name == 'D') {
      tot_gap_open++;
      tot_indel_extend += op->number - 1;
      tot_indel += op->number;
    }
  }
  
  int tot_mismatch = cc->distance - tot_indel;
  if (tot_mismatch < 0) {
    tot_mismatch = 0;
  }

  float score = tot_match * match + tot_mismatch * mismatch - tot_indel_extend * gap_extend - tot_gap_open * gap_open;
  if (score < 0) { score = 0; }
  int score_f = score * 100 / (read->length * match);

  
  //printf("%s, score = (%i * %f + %i * %f + %i * %f + %i * %f)= %f/%i\n", read->id, tot_match, match, tot_mismatch, mismatch, tot_indel_extend, gap_extend, tot_gap_open, gap_open, score, score_f);

  return score_f;

}

//#define SCORE(X, Y) ((X * 0.5 - Y*0.4)*100) / (X * 0.5);


sa_batch_t *sa_batch_new(array_list_t *fq_reads) {

  fastq_read_t *read;
  size_t num_reads = array_list_size(fq_reads);

  sa_batch_t *p = (sa_batch_t *) malloc(sizeof(sa_batch_t));
  p->num_reads = num_reads;
  p->fq_reads = fq_reads;
  p->mapping_lists = (array_list_t **) malloc(num_reads * sizeof(array_list_t *));

  for (size_t i = 0; i < num_reads; i++) {
    p->mapping_lists[i] = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  }

  return p;

}

sa_batch_t *sa_batch_simple_new(array_list_t *fq_reads) {

  fastq_read_t *read;
  size_t num_reads = array_list_size(fq_reads);

  sa_batch_t *p = (sa_batch_t *) malloc(sizeof(sa_batch_t));
  p->num_reads = num_reads;
  p->fq_reads = fq_reads;

  return p;

}


void sa_batch_free(sa_batch_t *p) {
  if (p) {
    if (p->fq_reads) { 
      array_list_free(p->fq_reads, (void *) fastq_read_free);
    }
    if (p->mapping_lists) { 
      free(p->mapping_lists); 
    }
    free(p);
  }
}



cigar_code_t* search_splice_junction(sw_optarg_t *sw_optarg,
				     seed_region_t *s_prev, seed_region_t *s_next,
				     int chromosome_id, int strand, 
				     char *sequence, genome_t *genome, 
				     size_t *sp_start, size_t *sp_end,
				     int *sp_type,
				     int *distance);


//Sa Mappings Reader, 2nd round
void *sa_alignments_reader_rna(void *input) {
  sa_wf_input_t *wf_input = (sa_wf_input_t *) input;  
  sa_wf_batch_t *new_wf_batch = NULL;
  sa_wf_batch_t *curr_wf_batch = wf_input->wf_batch;  
  sa_rna_input_t *sa_rna = curr_wf_batch->data_input;
  FILE *fd = sa_rna->file1;
 
  const int MAX_READS = 100;
  fastq_batch_reader_input_t *fq_reader_input = wf_input->fq_reader_input;
  array_list_t *reads = array_list_new(MAX_READS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  size_t num_reads = 0;
  array_list_t **mapping_lists = (array_list_t **) malloc(MAX_READS * sizeof(array_list_t *));
  size_t num_items;


  pair_server_input_t *pair_input = sa_rna->pair_input;
  int pair_mode = pair_input->pair_mng->pair_mode;

  extern size_t reads_ph2;

  if (pair_mode == SINGLE_END_MODE) {
    while (1) { 
      int type = file_read_type_items(fd);
      fastq_read_t *fq_read = file_read_fastq_reads(&num_items, fd);    
      if (fq_read == NULL) { break; }
    
      array_list_insert(fq_read, reads);    
      mapping_lists[num_reads] = array_list_new(50,
						1.25f, 
						COLLECTION_MODE_ASYNCHRONIZED);
    
      sa_file_read_alignments(num_items, mapping_lists[num_reads],
			      fq_read, fd);
      array_list_set_flag(type, mapping_lists[num_reads]);
      num_reads++;

      if (num_reads >= MAX_READS) { break; }

    }
  } else {
    while (1) {
      int type = file_read_type_items(fd);
      if (type == -1) { break; }
      fastq_read_t *fq_read = file_read_fastq_reads(&num_items, fd);    
      if (fq_read == NULL) { printf("File ERROR 1\n"); exit(-1); }    

      array_list_insert(fq_read, reads);      
      mapping_lists[num_reads] = array_list_new(50,
						1.25f, 
						COLLECTION_MODE_ASYNCHRONIZED);     
      if (type == SA_PARTIAL_TYPE) {
	sa_file_read_alignments(num_items, mapping_lists[num_reads],
				fq_read, fd);
      } else {
	file_read_alignments(num_items, mapping_lists[num_reads],
			     fq_read, fd);
      }
      array_list_set_flag(type, mapping_lists[num_reads]);
      num_reads++;


      type = file_read_type_items(fd);      
      fq_read = file_read_fastq_reads(&num_items, fd);    
      if (fq_read == NULL) { printf("File ERROR 2\n"); exit(-1); }    

      array_list_insert(fq_read, reads); 
      mapping_lists[num_reads] = array_list_new(50,
						1.25f, 
						COLLECTION_MODE_ASYNCHRONIZED);
      if (type == SA_PARTIAL_TYPE) {
	sa_file_read_alignments(num_items, mapping_lists[num_reads],
				fq_read, fd);
      } else {
	file_read_alignments(num_items, mapping_lists[num_reads],
			     fq_read, fd);
      }
      array_list_set_flag(type, mapping_lists[num_reads]);
      num_reads++;

      if (num_reads >= MAX_READS) { break; }

    }    
  }

  if (num_reads) {
    sa_batch_t *sa_batch = sa_batch_simple_new(reads);
    sa_batch->mapping_lists = mapping_lists;
    sa_batch->num_reads = num_reads;
    new_wf_batch = sa_wf_batch_new(NULL,
				   curr_wf_batch->sa_index,
				   curr_wf_batch->writer_input, 
				   sa_batch,
				   curr_wf_batch->data_input);
  } else {
    array_list_free(reads, (void *)fastq_read_free);
    free(mapping_lists);
  }

  reads_ph2 += num_reads;

  return new_wf_batch;

}


//Fastq Reader
void *sa_fq_reader_rna(void *input) {
  sa_wf_input_t *wf_input = (sa_wf_input_t *) input;
  sa_wf_batch_t *new_wf_batch = NULL;
  sa_wf_batch_t *curr_wf_batch = wf_input->wf_batch;  
  fastq_batch_reader_input_t *fq_reader_input = wf_input->fq_reader_input;
  array_list_t *reads = array_list_new(fq_reader_input->batch_size, 1.25f, 
				       COLLECTION_MODE_ASYNCHRONIZED);
  
  extern size_t fd_read_bytes;

  if (fq_reader_input->gzip) {
    // Gzip fastq file
    if (fq_reader_input->flags == SINGLE_END_MODE) {
      fastq_gzread_bytes_se(reads, fq_reader_input->batch_size, fq_reader_input->fq_gzip_file1);
    } else {
      fastq_gzread_bytes_pe(reads, fq_reader_input->batch_size, fq_reader_input->fq_gzip_file1, fq_reader_input->fq_gzip_file2);
    }
  } else {
    // Fastq file
    if (fq_reader_input->flags == SINGLE_END_MODE) {
      fd_read_bytes += fastq_fread_bytes_se(reads, fq_reader_input->batch_size, fq_reader_input->fq_file1);
    } else {
      fd_read_bytes += fastq_fread_bytes_aligner_pe(reads, fq_reader_input->batch_size, 
				   fq_reader_input->fq_file1, fq_reader_input->fq_file2);
    }
  }
  
  size_t num_reads = array_list_size(reads);
  
  if (num_reads == 0) {
    array_list_free(reads, (void *)fastq_read_free);
  } else {
    sa_batch_t *sa_batch = sa_batch_new(reads);

    new_wf_batch = sa_wf_batch_new(NULL,
				   curr_wf_batch->sa_index,
				   curr_wf_batch->writer_input, 
				   sa_batch,
				   curr_wf_batch->data_input);
  }

  return new_wf_batch;

}

int sa_bam_writer_rna(void *data) {
  //=======================
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;  
  sa_batch_t *mapping_batch = (sa_batch_t *) wf_batch->mapping_batch;  
  batch_writer_input_t *writer_input = wf_batch->writer_input;
  bam_file_t *bam_file = writer_input->bam_file;     
  size_t num_reads = mapping_batch->num_reads;
  fastq_read_t *fq_read;
  size_t num_items;
  array_list_t *read_list = mapping_batch->fq_reads;
  //=======================
  
  struct timeval start, end;
  double time;
  
  //batch_t *batch = (batch_t *) data;
  array_list_t *array_list;

  //mapping_batch_t *mapping_batch = (mapping_batch_t *) batch->mapping_batch;
  
  extern size_t total_reads;
  extern size_t reads_no_map;

  extern pthread_mutex_t mutex_sp;
  
  pthread_mutex_lock(&mutex_sp);
  total_reads += num_reads;
  pthread_mutex_unlock(&mutex_sp);
  
  //
  // DNA/RNA mode
  //

  for (size_t i = 0; i < num_reads; i++) {
    num_items = array_list_size(mapping_batch->mapping_lists[i]);
    fq_read = (fastq_read_t *) array_list_get(i, read_list);
    // mapped or not mapped ?	 
    if (num_items == 0) {
      pthread_mutex_lock(&mutex_sp);
      reads_no_map++;
      pthread_mutex_unlock(&mutex_sp);      
      write_unmapped_read(fq_read, bam_file);
      if (mapping_batch->mapping_lists[i]) {
	array_list_free(mapping_batch->mapping_lists[i], NULL);
      }
    } else {
      write_mapped_read(mapping_batch->mapping_lists[i], bam_file);
    }
  }
  
  // free memory
  sa_batch_free(mapping_batch);
  
  if (wf_batch) sa_wf_batch_free(wf_batch);
  
}

//Fastq Writer
int sa_sam_writer_rna(void *data) {
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;
  sa_rna_input_t *sa_rna = wf_batch->data_input;
  pair_server_input_t *pair_input = sa_rna->pair_input;
  int pair_mode = pair_input->pair_mng->pair_mode;  
  sa_batch_t *mapping_batch = (sa_batch_t *) wf_batch->mapping_batch;

  int flag, pnext = 0, tlen = 0;
  char rnext[4] = "*\0";
  char optional_flags[512] = "\0";

  fastq_read_t *read;
  array_list_t *read_list = mapping_batch->fq_reads;

  alignment_t *alig;
  array_list_t *mapping_list;
  FILE *out_file = (FILE *) wf_batch->writer_input->bam_file;

  sa_genome3_t *genome = wf_batch->sa_index->genome;

  size_t num_reads, num_mappings, num_mate_mappings;

  num_reads = mapping_batch->num_reads;

  //extern size_t total_reads, unmapped_reads, correct_reads;
  extern size_t total_reads;
  extern size_t reads_no_map;

  extern pthread_mutex_t mutex_sp;

  //struct timeval time_free_s, time_free_e;
  //extern double time_free_alig, time_free_list, time_free_batch;

  pthread_mutex_lock(&mutex_sp);
  total_reads += num_reads;
  pthread_mutex_unlock(&mutex_sp);

  if (pair_mode != SINGLE_END_MODE) {
    /*
    // PAIR MODE
    int len;
    char *sequence, *quality;
    char *seq;
    alignment_t *alig;
    array_list_t *mate_list;

    for (size_t i = 0; i < num_reads; i++) {
      read = (fastq_read_t *) array_list_get(i, read_list);
      if (i % 2 == 0)  {
	mate_list = mapping_batch->mapping_lists[i+1];
	num_mate_mappings = array_list_size(mate_list);
      } else {
	mate_list = mapping_list;
	num_mate_mappings = num_mappings;
      }
      
      mapping_list = mapping_batch->mapping_lists[i];
      num_mappings = array_list_size(mapping_list);
      
      if (num_mappings > 0) {
	//num_mapped_reads++;
	for (size_t j = 0; j < num_mappings; j++) {
	  alig = (alignment_t *) array_list_get(j, mapping_list);
	  
	  // update alignment
	  alig->secondary_alignment = 0;
	  if (num_mate_mappings != 1) {
	    alig->is_mate_mapped = 0;
	    alig->is_paired_end_mapped = 0;
	    alig->mate_strand = 0;
	  }

	  flag = 0;
	  if (alig->is_paired_end)                              flag += BAM_FPAIRED;
	  if (alig->is_paired_end_mapped)                       flag += BAM_FPROPER_PAIR;
	  if (!alig->is_seq_mapped)                             flag += BAM_FUNMAP;   
	  if ((!alig->is_mate_mapped) && (alig->is_paired_end)) flag += BAM_FMUNMAP;
	  if (alig->mate_strand)                                flag += BAM_FMREVERSE;
	  if (alig->pair_num == 1)	                        flag += BAM_FREAD1;
	  if (alig->pair_num == 2)                              flag += BAM_FREAD2;
	  if (alig->secondary_alignment)                        flag += BAM_FSECONDARY;
	  if (alig->fails_quality_check)                        flag += BAM_FQCFAIL;
	  if (alig->pc_optical_duplicate)                       flag += BAM_FDUP;
	  if (alig->seq_strand) {                               
	    flag += BAM_FREVERSE;
	    //seq = read->revcomp;
	  }

	  //num_total_mappings++;
	  fprintf(out_file, "%s\t%i\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%s\t%s\t%s\n", 
		  read->id,
		  flag,
		  genome->chrom_names[alig->chromosome],
		  alig->position + 1,
		  alig->mapq, //60, //(alig->map_quality > 3 ? 0 : alig->map_quality),
		  alig->cigar,
		  genome->chrom_names[alig->mate_chromosome],
		  alig->mate_position + 1,
		  alig->template_length,
		  alig->sequence,
		  alig->quality,
		  (alig->optional_fields ? alig->optional_fields : "")
		  );

	}
	alignment_free(alig);	 
      } else {
	//num_unmapped_reads++;
	fprintf(out_file, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", 
		read->id,
		read->sequence,
		read->quality
		);
      }
      
      array_list_free(mapping_list, (void *) NULL);
    }
    */
  } else {
    //SINGLE MODE
    for (size_t i = 0; i < num_reads; i++) {
      read = (fastq_read_t *) array_list_get(i, read_list);
      mapping_list = mapping_batch->mapping_lists[i];
      num_mappings = array_list_size(mapping_list);
      //printf("%i.Read %s (num_mappings %i)\n", i, read->id, num_mappings);    
      if (num_mappings > 0) {
	for (size_t j = 0; j < num_mappings; j++) {
	  alig = (alignment_t *) array_list_get(j, mapping_list);
	  flag = (alig->seq_strand ? 16 : 0);

	  int score = ((read->length * 0.5 - alig->map_quality*0.4)*100) / (read->length * 0.5);

	  fprintf(out_file, "%s\t%i\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%s\t%s\t%s\n", 
		  alig->query_name,
		  flag,
		  genome->chrom_names[alig->chromosome - 1],
		  alig->position + 1,
		  255,
		  alig->cigar,
		  rnext,
		  pnext,
		  tlen,
		  alig->seq_strand ? read->revcomp : read->sequence,
		  alig->quality,
		  optional_flags
		  );

	  //start_timer(time_free_s);
	  alignment_free(alig);
	  //stop_timer(time_free_s, time_free_e, time_free_alig);

	}
      } else {
	pthread_mutex_lock(&mutex_sp);
	reads_no_map++;
	pthread_mutex_unlock(&mutex_sp);
      
	fprintf(out_file, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", 
		read->id,
		read->sequence,
		read->quality
		);
      
      }
      //start_timer(time_free_s);
      array_list_free(mapping_list, (void *) NULL);
      //stop_timer(time_free_s, time_free_e, time_free_list);
    }

  }

  // free memory
  //start_timer(time_free_s);
  sa_batch_free(mapping_batch);
  if (wf_batch) sa_wf_batch_free(wf_batch);
  //stop_timer(time_free_s, time_free_e, time_free_batch);

  return 0;

}


void insert_seed_to_cal(size_t genome_start, size_t genome_end, 
			size_t read_start, size_t read_end, cal_t *cal, int id) {  
  assert(cal);
  linked_list_t *list_p = cal->sr_list;
  unsigned char actualization = 0;
  seed_region_t *item, *item_aux, *new_item_p, *item_free;
  linked_list_iterator_t* itr = linked_list_iterator_new(list_p);
  int max_cal_distance = 0;

  //SA_DEBUG("\tSEED INSERT %lu|%i-%i|%lu:\n", genome_start, read_start, read_end, genome_end);
  if (linked_list_size(list_p) <= 0) {
    new_item_p = seed_region_new(read_start, read_end, genome_start, genome_end, id, 0, 0);
    linked_list_insert(new_item_p, list_p);
    //Insert Seed
  } else {
    item = (seed_region_t *)linked_list_iterator_curr(itr);
    while (item != NULL) {
      if (read_start < item->read_start) {
	if (read_end < item->read_start && 
	    genome_end < item->genome_start) {
	  /*********************************************
	   *    Case 1: New item insert before item.   *
           *                                           *
           *        new item     item                  *
           *       |-------| |--------|                *
           ********************************************/
	  new_item_p = seed_region_new(read_start, read_end, genome_start, genome_end, id, 0, 0);
	  linked_list_iterator_insert(new_item_p, itr);
	  /*
	  new_item_p = cal_simple_new(chromosome_id,
				      strand, start, end);
	  linked_list_iterator_insert(new_item_p, itr);
	  insert_seed_to_cal(start, end, seq_start, seq_end, new_item_p, seed_id);	  
	  linked_list_iterator_prev(itr);
	  */
	} else {
	  /********************************************
           *  Case 2: Actualization item start        *
           *           new item                       *
           *          |-------|   item                *
           *                   |--------|             *
           ********************************************/
	  item->read_start   = read_start;
	  item->genome_start = genome_start;
	  if (read_end > item->read_end) {
	    /**************************************************
             *  Case 3: Actualization item start and item end *
             *          new item                              *
             *         |------------|                         *
             *              item                              *    
             *           |--------|                           *    
             **************************************************/
	    item->read_end = read_end;
	    item->genome_end = genome_end;
	    actualization = 1;
	  }
	}
	break;
      } else {
	if (read_end <= item->read_end) {
	  /**************************************************
           *  Case 4: The new item don't insert in the list * 
           *              item                              * 
           *         |-------------|                        * 
           *             new item                           * 
           *            |--------|                          * 
           **************************************************/
	  //Insert Seed
	  break;
	} else if (item->read_end >= read_start && 
		   item->genome_end >= genome_start) {
	  /********************************************
           *  Case 5: Actualization item end          *
           *            item                          * 
           *          |-------| new item              * 
           *                 |--------|               * 
           ********************************************/
	  item->read_end = read_end;
	  item->genome_end = genome_end;
	  //insert_seed_to_cal(start, end, seq_start, seq_end, item, seed_id);
	  //Insert Seed
	  actualization = 1;
	  //item->num_seeds++;
	  break;
	}
      } // end else
 
      //continue loop...
      linked_list_iterator_next(itr);
      item = linked_list_iterator_curr(itr);
      
    } // end while

    if (item == NULL) {
      /******************************************************* 
       * Case 5: Insert new item at the end of the list      * 
       *                 item    new item                    * 
       *              |-------| |--------|                   *    
       *******************************************************/
      //new_item_p = cal_simple_new(chromosome_id,
      //			  strand, start, end);
      new_item_p = seed_region_new(read_start, read_end, genome_start, genome_end, id, 0, 0);
      linked_list_insert_last(new_item_p, list_p);
      //Insert Seed
    }
    
    if (actualization == 1) {
      linked_list_iterator_next(itr);
      item_aux = linked_list_iterator_curr(itr);
      while (item_aux != NULL) {
	if (item->read_end + max_cal_distance < item_aux->read_start) {
	  break;
	} else {
	  if (item->read_end < item_aux->read_end) {
	    item->read_end = item_aux->read_end;
	    item->genome_end = item_aux->genome_end;
	  }
	  item_free = linked_list_iterator_remove(itr);
	  if (item_free) { seed_region_simple_free(item_free); }
	  item_aux = linked_list_iterator_curr(itr);
	}                              
      }
    }
  }//end first else

  //printf("End insert and actualization\n");
  linked_list_iterator_free(itr);

}


void generate_cals(size_t chromosome_id, 
		   short int strand,
		   size_t start, size_t end, 
		   size_t seq_start, size_t seq_end,
		   linked_list_t* list_p, size_t max_cal_distance, int seed_id) {  
  unsigned char actualization = 0;
  cal_t *item, *item_aux, *new_item_p, *item_free;
  
  linked_list_iterator_t* itr = linked_list_iterator_new(list_p);

  //  printf("LINKED LIST INSERT (%i:%i)%lu|%i-%i|%lu:\n", chromosome_id, strand, start, seq_start, seq_end, end);
  
  if (linked_list_size(list_p) <= 0) {
    new_item_p = cal_simple_new(chromosome_id, 
				strand, start, end);
    linked_list_insert(new_item_p, list_p);
    insert_seed_to_cal(start, end, seq_start, seq_end, new_item_p, seed_id);
  } else {
    item = (cal_t *)linked_list_iterator_curr(itr);
    while (item != NULL) {
      if (start < item->start) {
	if (end + max_cal_distance < item->start) {
	  /*********************************************
	   *    Case 1: New item insert before item.   *
           *                                           *
           *        new item     item                  *
           *       |-------| |--------|                *
           ********************************************/
	  new_item_p = cal_simple_new(chromosome_id,
				      strand, start, end);
	  linked_list_iterator_insert(new_item_p, itr);
	  insert_seed_to_cal(start, end, seq_start, seq_end, new_item_p, seed_id);	  
	  linked_list_iterator_prev(itr);
	} else {
	  /********************************************
           *  Case 2: Actualization item start        *
           *           new item                       *
           *          |-------|   item                *
           *                   |--------|             *
           ********************************************/
	  item->start = start;
	  //item->num_seeds++;
	  if (end > item->end) {
	    /**************************************************
             *  Case 3: Actualization item start and item end *
             *          new item                              *
             *         |------------|                         *
             *              item                              *    
             *           |--------|                           *    
             **************************************************/
	    item->end = end;
	    actualization = 1;
	  }
	  insert_seed_to_cal(start, end, seq_start, seq_end, item, seed_id);
	}
	break;
      } else {
	if (end <= item->end) {
	  /**************************************************
           *  Case 4: The new item don't insert in the list * 
           *              item                              * 
           *         |-------------|                        * 
           *             new item                           * 
           *            |--------|                          * 
           **************************************************/
	  //item->num_seeds++;
	  //Insert Seed
	  insert_seed_to_cal(start, end, seq_start, seq_end, item, seed_id);
	  break;
	} else if (item->end + max_cal_distance >= start) {
	  /********************************************
           *  Case 5: Actualization item end          *
           *            item                          * 
           *          |-------| new item              * 
           *                 |--------|               * 
           ********************************************/
	  item->end = end;
	  insert_seed_to_cal(start, end, seq_start, seq_end, item, seed_id);
	  //Insert Seed
	  actualization = 1;
	  //item->num_seeds++;
	  break;
	}
      } // end else
 
      //continue loop...
      linked_list_iterator_next(itr);
      item = linked_list_iterator_curr(itr);
      
    } // end while

    if (item == NULL) {
      /******************************************************* 
       * Case 5: Insert new item at the end of the list      * 
       *                 item    new item                    * 
       *              |-------| |--------|                   *    
       *******************************************************/
      new_item_p = cal_simple_new(chromosome_id,
				  strand, start, end);      
      linked_list_insert_last(new_item_p, list_p);
      insert_seed_to_cal(start, end, seq_start, seq_end, new_item_p, seed_id);
      //Insert Seed
    }
    
    if (actualization == 1) {
      linked_list_iterator_next(itr);
      item_aux = linked_list_iterator_curr(itr);
      while (item_aux != NULL) {
	if (item->end + max_cal_distance < item_aux->start) {
	  //printf("\t\tSTOP Actualization\n");
	  break;
	} else {
	  //printf("\t\tCONTINUE Actualization. item->end=%d < item_aux->end=%d?\n", item->end, item_aux->end);
	  if (item->end < item_aux->end) {
	    //printf("\t\tActualization end value %d\n", item_aux->end);
	    item->end = item_aux->end;
	    //item->seq_end = item_aux->seq_end;	    
	    seed_region_t *seed_aux;
	    while (seed_aux = linked_list_remove_first(item_aux->sr_list)) {
	      //Insert Seed
	      insert_seed_to_cal(seed_aux->genome_start, seed_aux->genome_end, 
				 seed_aux->read_start, seed_aux->read_end, item, seed_aux->id);
	      seed_region_free(seed_aux);
	    }	    
	  }
          //printf("\t\tDelete item %d-%d\n", item_aux->start, item_aux->end);
	  item_free = linked_list_iterator_remove(itr);
	  if (item_free) { cal_simple_free(item_free); }
	  //printf("\t\tDelete OK!\n");
	  item_aux = linked_list_iterator_curr(itr);
	}                                                                       
      }
    }
  }//end first else

  //printf("End insert and actualization\n");
  linked_list_iterator_free(itr);
  /*
  printf(" ============== RESULT LIST: ============ \n");
  linked_list_item_t *item_list = list_p->first;
  cal_t *cal;
  while (item_list != NULL) {
    cal = item_list->item;
    //printf("[%lu-%lu]\n", cal->start, cal->end);
    cal_print(cal);
    item_list = item_list->next;
  }
  printf(" ============== ============= ============ \n");
  */
}

//--------------------------------------------------------------------

void cal_mng_print(cal_mng_t *p) {
  cal_t *cal;
  seed_region_t *seed_first, *seed_last;
  linked_list_iterator_t itr;
  size_t num_items;
  int items = 0;
  linked_list_item_t *list_item;

  if (p->cals_lists) {
    linked_list_t *cal_list;
    for (unsigned int i = 0; i < p->num_chroms; i++) {
      cal_list = p->cals_lists[i];
      num_items = linked_list_size(cal_list);
      if (num_items) {
	printf("CHROMOSOME %i:\n", i + 1);
	items = 1;
	linked_list_iterator_init(cal_list, &itr);
	cal = linked_list_iterator_curr(&itr);
	while (cal) {
	  printf(" [%lu-%lu](%lu) :\n", cal->start, cal->end, linked_list_size(cal->sr_list));
	  seed_first = linked_list_get_first(cal->sr_list);
	  seed_last = linked_list_get_last(cal->sr_list);

	  list_item = cal->sr_list->first;
	  while (list_item) { 
	    seed_first = list_item->item;
	    printf("\t[%lu|%i-%i|%lu] ", seed_first->genome_start, seed_first->read_start, seed_first->read_end, 
		   seed_first->genome_end);
	    list_item = list_item->next;
	  }
	  printf("\n");

	  cal = linked_list_iterator_next(&itr);
	}
	printf("\n");
      }
    }
  }

  if (!items) {
    printf("EMPTY!\n");
  }

}

/*
array_list_t *rna_seeding(fastq_read_t *read,
			  sa_index3_t *sa_index, 
			  genome_t *genome,
			  avls_list_t *avls_list,
			  metaexons_t *metaexons,
			  array_list_t *alignments_list) { 


}
*/


void cal_purge (cal_t *cal) {
  linked_list_item_t *list_item, *list_item_prev = NULL;
  seed_region_t *seed_prev = NULL, *seed_next;
  linked_list_iterator_t iter;
  //printf("///////////////////////////////// purge ................\n");
  //cal_print(cal);

  seed_prev = linked_list_get_first(cal->sr_list);
  seed_next = linked_list_get_last(cal->sr_list);

  if (seed_prev->genome_start != cal->start || 
      seed_next->genome_end != cal->end) {

    //Delete Seeds out of limits
    linked_list_iterator_init(cal->sr_list, &iter);
    seed_next = linked_list_iterator_curr(&iter);
    while (seed_next != NULL) {
      //printf("Seed: [%lu-%lu]\n", seed_next->genome_start, seed_next->genome_end);
      
      if (seed_next->genome_start < cal->start || seed_next->genome_end > cal->end) {
	seed_next = linked_list_iterator_remove(&iter);
	seed_region_free(seed_next);
	
	seed_next = linked_list_iterator_curr(&iter);
	continue;
      }
      seed_next = linked_list_iterator_next(&iter);    
    }    
  }

  linked_list_iterator_init(cal->sr_list, &iter);
  seed_next = linked_list_iterator_curr(&iter);
  while (seed_next != NULL) {
    if (seed_prev != NULL) {
      //printf("Results\n");
      //seed_prev = list_item_prev->item;
      //printf("Results 1\n");
      //seed_next = list_item->item;
      //printf("Results 2\n");
      //printf("[%i-%i], [%i-%i]\n", seed_prev->read_start, seed_prev->read_end, seed_next->read_start, seed_next->read_end);
      if ((seed_prev->read_end >= seed_next->read_start) || 
	  (seed_prev->genome_end >= seed_next->genome_start)) {
	//We have a problem! with read coords or genome coords!
	//Select the best seed
	int seed_len_prev = seed_prev->read_end - seed_prev->read_start;
	int seed_len_next = seed_next->read_end - seed_next->read_start;
	if (seed_len_prev > seed_len_next) {
	  //printf("Delete Next **************\n");
	  //Delete Seed next
	  seed_next = linked_list_iterator_remove(&iter);
	  seed_region_free(seed_next);

	  seed_next = linked_list_iterator_curr(&iter);
	  continue;
	} else {
	  //printf("Delete Prev **************\n");
	  //Delete seed prev
	  seed_prev = linked_list_iterator_prev(&iter);
	  seed_prev = linked_list_iterator_remove(&iter);
	  seed_region_free(seed_prev);

	  seed_next = linked_list_iterator_curr(&iter);	  
	}
      }
    }

    seed_prev = seed_next;
    seed_next = linked_list_iterator_next(&iter);

  }
}

int validate_cal(cal_t *cal) {
  //if (!cal) { return 0; }

  seed_region_t *seed_prev = linked_list_get_first(cal->sr_list);
  seed_region_t *seed_next = linked_list_get_last(cal->sr_list);
  linked_list_iterator_t iter;

  if (seed_prev->genome_start != cal->start || 
      seed_next->genome_end != cal->end) {
    return 0;
  }

  seed_prev = NULL;
  linked_list_iterator_init(cal->sr_list, &iter);
  seed_next = linked_list_iterator_curr(&iter);
  while (seed_next != NULL) {
    if (seed_prev) {
      if (seed_prev->read_end > seed_next->read_start) {
	//return 0;
	int dsp = seed_prev->read_end - seed_next->read_start + 1;
	if (dsp > 10) { return 0; }

	seed_prev->read_end -= dsp;
	seed_next->read_start += dsp;
	seed_prev->genome_end -= dsp;
	seed_next->genome_start += dsp;	
      }
      if (seed_prev->genome_end > seed_next->genome_start) {
	int dsp = seed_prev->genome_end - seed_next->genome_start + 1;
	if (dsp > 10) { return 0; }

	seed_prev->read_end -= dsp;
	seed_next->read_start += dsp;
	seed_prev->genome_end -= dsp;
	seed_next->genome_start += dsp;	
      }
    }
    seed_prev = seed_next;
    seed_next = linked_list_iterator_next(&iter);
  }  
  
  return 1;

}

float generate_cals_merge_score(array_list_t *cals_list, int read_length) {
  int len_cal = 0;
  size_t num_cals = array_list_size(cals_list);
  cal_t *cal;
  linked_list_iterator_t itr;  
  seed_region_t *s, *s_prev;
  int i;

  //for (i = 0; i < num_cals; i++) {    
  //cal = array_list_get(i, cals_list);
  //if (!validate_cal(cal)) { return 0.0; }
  //}

  for (i = 0; i < num_cals; i++) {    
    cal = array_list_get(i, cals_list);

    linked_list_iterator_init(cal->sr_list, &itr);
    s = (seed_region_t *) linked_list_iterator_curr(&itr);
    while (s != NULL) {
      //if (s->read_start > s->read_end) { assert(s->read_start); }
      len_cal += s->read_end - s->read_start;
      s = (seed_region_t *) linked_list_iterator_next(&itr);
    }
  }

  //printf("Score = %f\n", (float)(len_cal*100)/(float)read_length);

  return (float)(len_cal*100)/(float)read_length;

}

//==================================================================

typedef struct sa_sw_depth {
  int depth;
  char *q[MAX_DEPTH];
  char *r[MAX_DEPTH];
  //seed_region_t *seed_ref[MAX_DEPTH];
  void *item_ref[MAX_DEPTH];
  int type[MAX_DEPTH];
} sa_sw_depth_t;

void sw_process(sa_sw_depth_t *sw_depth, sw_optarg_t *sw_optarg, sw_multi_output_t *output) {

  if (sw_depth->depth == 0) { return; }

  float match = sw_optarg->subst_matrix['A']['A'];
  int distance;

  /*
  printf("depth: %i\n", sw_depth->depth);
  for (int i = 0; i < sw_depth->depth; i++) {
    printf("REF: (%i)%s\n", strlen(sw_depth->r[i]), sw_depth->r[i]);
    printf("SEQ: (%i)%s\n", strlen(sw_depth->q[i]), sw_depth->q[i]);
  }
  */

  smith_waterman_mqmr(sw_depth->q, sw_depth->r, sw_depth->depth, sw_optarg, 1, output);

  for (int i = 0; i < sw_depth->depth; i++) {
    if (sw_depth->type[i] != SJ_SW) {
      cigar_code_t *cigar_code = generate_cigar_code(output->query_map_p[i],
						     output->ref_map_p[i],
						     strlen(output->ref_map_p[i]),
						     output->query_start_p[i], output->ref_start_p[i],
						     strlen(sw_depth->q[i]), strlen(sw_depth->r[i]),
						     &distance, sw_depth->type[i]);
    
      //printf("==========================SW============================\n");
      //printf("REF: %s\n", output->ref_map_p[i]);
      //printf("SEQ: %s\n", output->query_map_p[i]);
      //printf("======================================================\n");

      //printf("/////////// cigar : %s\n", new_cigar_code_string(cigar_code));
      if (sw_depth->type[i] == FIRST_SW || sw_depth->type[i] == LAST_SW) {
	float norm_score = NORM_SCORE(output->score_p[i], strlen(sw_depth->q[i]), match);
	//printf("norm_score = %f\n", norm_score);
	if (norm_score > 0.4) {
	  ((seed_region_t *)sw_depth->item_ref[i])->info = cigar_code;
	} else {
	  array_list_clear(cigar_code->ops, (void *)cigar_op_free);
	  cigar_code_free(cigar_code);
	}
      } else {
	((seed_region_t *)sw_depth->item_ref[i])->info = cigar_code;
      }

      free(sw_depth->q[i]);
      free(sw_depth->r[i]);
    
      free(output->query_map_p[i]);
      free(output->ref_map_p[i]);
      output->query_map_p[i] = NULL;
      output->ref_map_p[i] = NULL;
    } else {
      //printf("REF: %s\n", output->ref_map_p[i]);
      //printf("SEQ: %s\n", output->query_map_p[i]);      
    }    
  }
  
}

void sw_insert_item(char *q, char *r, int type, void *item, 
		    sw_optarg_t *sw_optarg, sw_multi_output_t *output,
		    sa_sw_depth_t *sw_depth) {

  if (q != NULL && r != NULL) {
    sw_depth->q[sw_depth->depth] = strdup(q);
    sw_depth->r[sw_depth->depth] = strdup(r);
    sw_depth->type[sw_depth->depth] = type;

    sw_depth->item_ref[sw_depth->depth++] = item;
  }

  if (sw_depth->depth == 4 || 
      (q == NULL && r == NULL)) {
    sw_process(sw_depth, sw_optarg, output);
    sw_depth->depth = 0;
  }
  
}

//==================================================================

void cal_fill_gaps(cal_t *cal,
		   char *query_map,
		   genome_t *genome,
		   metaexons_t *metaexons,
		   avls_list_t *avls_list, 
		   sa_sw_depth_t *sw_depth, 
		   sw_optarg_t *sw_optarg, 
		   sw_multi_output_t *output) { 

  //cal_print(cal);
  //printf("======= FILL GAPS =======\n");
  //cigar_code_t *cigar_code;
  int max_size = 2048;
  char *reference = (char *)calloc(max_size, sizeof(char));
  int distance;
  int first, last;
  int min_distance;
  int closed;  
  avl_node_t *node_avl_start;
  avl_node_t *node_avl_end;
  int add_seed, sw_add, sp_found;
  //sw_item_t *sw_item;

  //cal_t *cal = array_list_get(i, meta_alignment->cals_list);
  seed_region_t *seed_prev = NULL, *seed_next;
  linked_list_item_t *list_item, *list_item_prev;
  size_t sp_start, sp_end;
  int distance_aux;
  
  const int FLANK = 8;
  const int MAX_GAPS = 1024;
  const int MIN_GAP = 5;
  //int *gaps_info = (int *)calloc(MAX_GAPS, sizeof(int));
  int total_gaps = 0;
  int set_sw;

  //Last. Generate cigars and smith-watermans
  list_item = cal->sr_list->first;
  seed_prev = NULL;
  while (list_item != NULL) {
    seed_next = (seed_region_t *)list_item->item;
    int num_match = seed_next->read_end - seed_next->read_start + 1;
    if (seed_prev != NULL) {
      //Close nt 
      if (seed_next->read_start <= seed_prev->read_end || 
	  seed_next->genome_start <= seed_prev->genome_end || 
	  seed_next->read_start - seed_prev->read_end <= MIN_GAP || 
	  seed_next->genome_start - seed_prev->genome_end <= MIN_GAP ) {
	seed_prev->genome_end   = seed_prev->genome_end - FLANK;
	seed_next->genome_start = seed_next->genome_start + FLANK;
	seed_prev->read_end   = seed_prev->read_end - FLANK;
	seed_next->read_start = seed_next->read_start + FLANK;
      }

      assert(seed_next->read_start > seed_prev->read_end);      
      assert(seed_next->genome_start > seed_prev->genome_end);

      size_t gap_read = seed_next->read_start - seed_prev->read_end  - 1;
      size_t gap_genome = seed_next->genome_start - seed_prev->genome_end  - 1; 
      //assert(gap_genome > 0);
      //printf("FILL GAPS ::: (gap_read = %lu), (gap_genome = %lu)\n", gap_read, gap_genome);      
      distance = 0;
      closed = 0;
      add_seed = 0;
      sw_add = 0;
      sp_found = 0;

      if (gap_read == gap_genome &&
	  gap_read == 0) {
	closed = 1;
      } else if (gap_read == gap_genome) {
	size_t genome_end   = seed_next->genome_start - 1;
	size_t genome_start = seed_prev->genome_end + 1;
	int read_end        = seed_next->read_start - 1;
	int read_start      = seed_prev->read_end + 1;
	char *query         = &query_map[read_start];
	first = -1;
	last = -1;
	//assert(genome_end - genome_start + 1 < 2048);
	if (genome_end - genome_start >= max_size) {
	  fprintf(stderr, "Error max reference\n");
	  exit(-1);

	  free(reference);
	  max_size = genome_end - genome_start + 1024;
	  reference = (char *)calloc(max_size, sizeof(char));
	}
	genome_read_sequence_by_chr_index(reference, 0, cal->chromosome_id - 1,
					  &genome_start, &genome_end, genome);
	
	//printf("[%lu|%i]-GAP-[%i|%lu]\n", genome_start, read_start, read_end, genome_end);
	for (int k = 0; k < gap_read; k++) {
	  //printf("[q:%c vs r:%c]\n", query[k], reference[k]);
	  if (k >= strlen(query)) {
	    fprintf(stderr, "Error read sequence overflow\n");
	    exit(-1);
	  }
	  if (query[k] != reference[k]) {
	    distance++;
	    if (first == -1) first = k;
	    last = k;
	  }
	}

	min_distance = (gap_genome / 3) + 2;
	//printf("Distance %i <= %i\n", distance, min_distance);
	if (distance <= min_distance) {
	  closed = 1;
	  add_seed = 1;
	}		
      }

      if (!closed) {
	int gap = gap_genome - gap_read; 
	if (gap > min_intron) {
	  //printf("Search splice");
	  int sp_type;
	  int nt = search_simple_splice_junction(seed_prev, seed_next,
						 cal->chromosome_id, cal->strand, 
						 query_map, genome, 
						 &sp_start, 
						 &sp_end,
						 &sp_type,
						 &distance_aux);
	  
	  if (!nt) {
	    nt = search_simple_splice_junction_semi_cannonical(seed_prev, seed_next,
							       cal->chromosome_id, cal->strand, 
							       query_map, genome, 
							       &sp_start, 
							       &sp_end,
							       &sp_type,
							       &distance_aux);	    
	  }
	  
	  if (nt >= min_intron && 
	      seed_prev->genome_start < sp_start && 
	      sp_end < seed_next->genome_end && 
	      sp_start < sp_end) {
	    closed = 1;
	    sp_found = 1;
	  }
	}
      }
      
      if (!closed) {
	//SMITH-WATERMAN
	if (gap_read <= 0 || gap_genome <= 0) {
	  seed_prev->genome_end   = seed_prev->genome_end - FLANK;
	  seed_next->genome_start = seed_next->genome_start + FLANK;
	  seed_prev->read_end   = seed_prev->read_end - FLANK;
	  seed_next->read_start = seed_next->read_start + FLANK;	    
	}       
	add_seed = 1;
	sw_add = 1;
      }
      
      cigar_code_t *cc = cigar_code_new();
      cigar_code_append_new_op(seed_prev->read_end - seed_prev->read_start + 1, 'M', cc);
      seed_prev->info = cc;
      //printf("********************(fill gaps)[%lu-%lu] cigar code: %s\n", seed_prev->read_start, seed_prev->read_end, new_cigar_code_string(cc));
      
      if (sp_found) {
	//printf("Insert!!!!!!!!!\n");
	cigar_code_append_new_op(sp_end - sp_start + 1, 'N', cc);
	cc->distance = distance_aux;
      }

      if (add_seed) {	
	seed_region_t *new_seed = seed_region_new(seed_prev->read_end + 1, seed_next->read_start - 1, 
						  seed_prev->genome_end + 1, seed_next->genome_start - 1, 
						  seed_prev->id + 1, 0, 0);	
	linked_list_item_t *new_item = linked_list_item_new(new_seed);

	list_item_prev->next = new_item;
	new_item->prev = list_item_prev;
	list_item->prev = new_item;
	new_item->next = list_item; 
	cal->sr_list->size++;
	
	//printf(" ::> %lu|%i-%i|%lu\n", new_seed->genome_start, new_seed->read_start, new_seed->read_end, new_seed->genome_end);
	
	if (sw_add) {
	  char reference[2048];
	  size_t genome_start = seed_prev->genome_end + 1;
	  size_t genome_end   = seed_next->genome_start - 1;
	  int read_start      = seed_prev->read_end + 1;
	  int read_end        = seed_next->read_start - 1;
	  char query[2048];

	  assert(read_start < strlen(query_map));
	  assert(read_end < strlen(query_map));
	  assert(read_start + (read_end - read_start + 1) < strlen(query_map));
	  //assert(genome_end - genome_start + 1 < 2048);

	  if (genome_end - genome_start + 1 >= 2048 || 
	      read_end - read_start + 1 <= 0) {
	    //fprintf(stderr, "%i - %i + 1 = %i\n", read_end, read_start, read_end - read_start + 1);
	    break;
	  }

	  genome_read_sequence_by_chr_index(reference, 0, cal->chromosome_id - 1,
					    &genome_start, &genome_end, genome);
	  
	  memcpy(query, &query_map[read_start],  read_end - read_start + 1);
	  query[read_end - read_start + 1] = '\0';

	  //printf("Read   : %i - %i\n", read_start, read_end);
	  //printf("Genome : %i - %i\n", genome_start, genome_end);
	  
	  //printf("SW-Query: %s\n", query);
	  //printf("SW-Refer: %s\n", reference);
	  
	  sw_insert_item(query, reference, MIDDLE_SW, new_seed, 
			 sw_optarg, output, sw_depth);
	  
	} else {
	  cigar_code_t *cc = cigar_code_new();
	  cigar_code_append_new_op(new_seed->read_end - new_seed->read_start + 1, 'M', cc);
	  cc->distance = distance;
	  new_seed->info = cc;
	  //printf("********************>>(fill gaps)[%lu-%lu] cigar code: %s\n", new_seed->read_start, new_seed->read_end, new_cigar_code_string(cc));
	}
      }
    }

    seed_prev      = seed_next;
    list_item_prev = list_item;
    list_item      = list_item->next;

  }

  cigar_code_t *cc = cigar_code_new();
  cigar_code_append_new_op(seed_prev->read_end - seed_prev->read_start + 1, 'M', cc);
  seed_prev->info = cc;
  //printf("********************>>||(fill gaps)[%lu-%lu] cigar code: %s\n", seed_prev->read_start, seed_prev->read_end, new_cigar_code_string(cc));
  free(reference);


}

void cal_fill_gaps_2(cal_t *cal,
		     char *query_map,
		     genome_t *genome,
		     metaexons_t *metaexons,
		     avls_list_t *avls_list, 
		     sa_sw_depth_t *sw_depth, 
		     sw_optarg_t *sw_optarg, 
		     sw_multi_output_t *output) {

  int max_size = 2048;
  char *reference = (char *)calloc(max_size, sizeof(char));
  int first, last;
  int min_distance;
  int closed;  
  avl_node_t *node_avl_start;
  avl_node_t *node_avl_end;
  int add_seed, sw_add, sp_found;
  //sw_item_t *sw_item;

  //cal_t *cal = array_list_get(i, meta_alignment->cals_list);
  seed_region_t *seed_prev = NULL, *seed_next;
  linked_list_item_t *list_item, *list_item_prev;
  size_t sp_start, sp_end;
  int distance_aux;
  int distance, sp_type;

  const int FLANK = 8;
  const int MAX_GAPS = 1024;
  const int MIN_GAP = 5;
  //int *gaps_info = (int *)calloc(MAX_GAPS, sizeof(int));
  int total_gaps = 0;
  int set_sw;
  cigar_code_t *cc_sj;

  //Last. Generate cigars and smith-watermans
  list_item = cal->sr_list->first;
  seed_prev = NULL;
  while (list_item != NULL) {
    seed_next = (seed_region_t *)list_item->item;
    int num_match = seed_next->read_end - seed_next->read_start + 1;
    if (seed_prev != NULL) {      
      cc_sj = NULL;
      size_t gap_read = seed_next->read_start - seed_prev->read_end  - 1;
      size_t gap_genome = seed_next->genome_start - seed_prev->genome_end  - 1; 
      int gap = gap_genome - gap_read;

      if (gap >= min_intron) {
	cc_sj = search_splice_junction(sw_optarg,
				       seed_prev, 
				       seed_next,
				       cal->chromosome_id,
				       cal->strand,
				       query_map, 
				       genome, 
				       &sp_start, 
				       &sp_end,
				       &sp_type,
				       &distance);
      }
      
      if (cc_sj) {
	seed_region_t *new_seed = seed_region_new(seed_prev->read_end + 1, seed_next->read_start - 1, 
						  seed_prev->genome_end + 1, seed_next->genome_start - 1, 
						  seed_prev->id + 1, 0, 0);		
	new_seed->info = cc_sj;
	
	linked_list_item_t *new_item = linked_list_item_new(new_seed);	
	list_item_prev->next = new_item;
	new_item->prev = list_item_prev;
	list_item->prev = new_item;
	new_item->next = list_item; 
	cal->sr_list->size++;
	
      } else {
	//Close nt
	if (seed_next->read_start <= seed_prev->read_end || 
	    seed_next->genome_start <= seed_prev->genome_end || 
	    seed_next->read_start - seed_prev->read_end <= MIN_GAP || 
	    seed_next->genome_start - seed_prev->genome_end <= MIN_GAP ) {
	  seed_prev->genome_end   = seed_prev->genome_end - FLANK;
	  seed_next->genome_start = seed_next->genome_start + FLANK;
	  seed_prev->read_end   = seed_prev->read_end - FLANK;
	  seed_next->read_start = seed_next->read_start + FLANK;
	}
	
	if (seed_next->read_start <= seed_prev->read_end) {
	  int dsp = seed_prev->read_end - seed_next->read_start + 1;
	  seed_prev->genome_end   = seed_prev->genome_end - dsp;
	  seed_next->genome_start = seed_next->genome_start + dsp;
	  seed_prev->read_end   = seed_prev->read_end - dsp;
	  seed_next->read_start = seed_next->read_start + dsp;
	}
  
	if (seed_next->genome_start <= seed_prev->genome_end) {
	  int dsp = seed_prev->genome_end - seed_next->genome_start + 1;
	  seed_prev->genome_end   = seed_prev->genome_end - dsp;
	  seed_next->genome_start = seed_next->genome_start + dsp;
	  seed_prev->read_end   = seed_prev->read_end - dsp;
	  seed_next->read_start = seed_next->read_start + dsp;
	}
	
	size_t gap_read = seed_next->read_start - seed_prev->read_end  - 1;
	size_t gap_genome = seed_next->genome_start - seed_prev->genome_end  - 1; 
	
	distance = 0;
	closed = 0;
	add_seed = 0;
	sw_add = 0;
	sp_found = 0;
	
	if (gap_read == gap_genome) {
	  size_t genome_end   = seed_next->genome_start - 1;
	  size_t genome_start = seed_prev->genome_end + 1;
	  int read_end        = seed_next->read_start - 1;
	  int read_start      = seed_prev->read_end + 1;
	  char *query         = &query_map[read_start];
	  
	  first = -1;
	  last = -1;
	  
	  //assert(genome_end - genome_start + 1 < 2048);
	  if (genome_end - genome_start >= max_size) {
	    fprintf(stderr, "Error max reference\n");
	    exit(-1);
	    
	    free(reference);
	    max_size = genome_end - genome_start + 1024;
	    reference = (char *)calloc(max_size, sizeof(char));
	  }
	  
	  genome_read_sequence_by_chr_index(reference, 0, cal->chromosome_id - 1,
					    &genome_start, &genome_end, genome);
	  
	  //printf("[%lu|%i]-GAP-[%i|%lu]\n", genome_start, read_start, read_end, genome_end);
	  for (int k = 0; k < gap_read; k++) {
	    //printf("[q:%c vs r:%c]\n", query[k], reference[k]);
	    if (k >= strlen(query)) {
	      fprintf(stderr, "Error read sequence overflow\n");
	      exit(-1);
	    }
	    if (query[k] != reference[k]) {
	      distance++;
	      if (first == -1) first = k;
	      last = k;
	    }
	  }
	  
	  min_distance = (gap_genome / 3) + 2;
	  //printf("Distance %i <= %i\n", distance, min_distance);
	  if (distance <= min_distance) {
	    closed = 1;
	  }
	}
	
	seed_region_t *new_seed = seed_region_new(seed_prev->read_end + 1, seed_next->read_start - 1, 
						  seed_prev->genome_end + 1, seed_next->genome_start - 1, 
						  seed_prev->id + 1, 0, 0); 	
	
	linked_list_item_t *new_item = linked_list_item_new(new_seed);
	
	list_item_prev->next = new_item;
	new_item->prev = list_item_prev;
	list_item->prev = new_item;
	new_item->next = list_item; 
	cal->sr_list->size++;
	
	if (!closed) {
	  int gap = gap_genome - gap_read; 	
	  //SMITH-WATERMAN
	  
	  char reference[2048];
	  size_t genome_start = seed_prev->genome_end + 1;
	  size_t genome_end   = seed_next->genome_start - 1;
	  int read_start      = seed_prev->read_end + 1;
	  int read_end        = seed_next->read_start - 1;
	  char query[2048];
	  
	  assert(read_start < strlen(query_map));
	  //assert(read_end < strlen(query_map));
	  
	  
	  if (read_start + (read_end - read_start + 1) < strlen(query_map) ||
	      read_end < strlen(query_map) ||
	      genome_end - genome_start + 1 >= 2048 || 
	      read_end - read_start + 1 <= 0) {
	    //fprintf(stderr, "%i - %i + 1 = %i\n", read_end, read_start, read_end - read_start + 1);
	    break;
	  }
	  
	  genome_read_sequence_by_chr_index(reference, 0, cal->chromosome_id - 1,
					    &genome_start, &genome_end, genome);
	  
	  memcpy(query, &query_map[read_start],  read_end - read_start + 1);
	  query[read_end - read_start + 1] = '\0';	  
	  
	  sw_insert_item(query, reference, MIDDLE_SW, new_seed, 
			 sw_optarg, output, sw_depth);	
	} else {
	  cigar_code_t *cc = cigar_code_new();
	  cigar_code_append_new_op(new_seed->read_end - new_seed->read_start + 1, 'M', cc);
	  cc->distance = distance;
	  new_seed->info = cc;	  
	}
      }
    }

    seed_prev      = seed_next;
    list_item_prev = list_item;
    list_item      = list_item->next;
    
  }


  free(reference);

}




void convert_batch_to_str(sa_wf_batch_t *wf_batch) {
  sa_batch_t *sa_batch = wf_batch->mapping_batch;
  sa_rna_input_t *sa_rna = wf_batch->data_input;
  pair_server_input_t *pair_input = sa_rna->pair_input;
  pair_mng_t *pair_mng = pair_input->pair_mng;
  report_optarg_t *report_optarg = pair_input->report_optarg;
  int bam_format = wf_batch->writer_input->bam_format;
  int pair_mode = pair_input->pair_mng->pair_mode;  
  //sa_batch_t *sa_batch = wf_batch->mapping_batch;  
  size_t num_reads, num_mappings;
  //sa_rna_input_t *sa_rna = wf_batch->data;
  //pair_server_input_t *pair_input = sa_rna->pair_input;

  wf_batch->data_output_size = 0;

  num_reads = sa_batch->num_reads;
  if (!num_reads) { return; }
  
  int flag, pnext = 0, tlen = 0;
  char rnext[4] = "*\0";
  char optional_flags[512] = "\0";
  
  fastq_read_t *read;
  array_list_t *read_list = sa_batch->fq_reads;
  alignment_t *alig;
  array_list_t *mapping_list;

  sa_genome3_t *genome = wf_batch->sa_index->genome;
  
  //extern size_t total_reads, unmapped_reads, correct_reads;
  extern size_t total_reads;
  extern size_t reads_no_map;
  extern pthread_mutex_t mutex_sp;
  struct timeval time_free_s, time_free_e;
  extern double time_free_alig, time_free_list, time_free_batch;

  int num_mate_mappings;
  
  read = (fastq_read_t *) array_list_get(0, read_list);

  if (!bam_format) {
    int read_size   = ((strlen(read->id) + 50) + ((read->length + 50)*2) + (1000));
    int buffer_max_size = read_size * (num_reads * MAX_ALIG);
    int buffer_size = 0;

    //Malloc Buffer
    char *buffer = (char *)malloc(sizeof(char)*buffer_max_size);
    buffer[0] = '\0';

    extern basic_statistics_t *basic_st;
    size_t total_reads = num_reads;
    size_t num_mapped_reads = 0;
    size_t total_mappings = 0;
    size_t reads_uniq_mappings = 0;

    if (pair_mode != SINGLE_END_MODE) {
      // PAIR MODE
      int len;
      char *sequence, *quality;
      char *seq;
      alignment_t *alig;
      array_list_t *mate_list;

      for (size_t i = 0; i < num_reads; i++) {
	read = (fastq_read_t *) array_list_get(i, read_list);
	if (i % 2 == 0)  {
	  mate_list = sa_batch->mapping_lists[i+1];
	  num_mate_mappings = array_list_size(mate_list);
	} else {
	  mate_list = mapping_list;
	  num_mate_mappings = num_mappings;
	}
	
	mapping_list = sa_batch->mapping_lists[i];
	num_mappings = array_list_size(mapping_list);


	int mqual = 0;
	if (num_mappings == 1) { 
	  reads_uniq_mappings++; 
	  mqual = 50;
	} else if (num_mappings == 2) {
	  mqual = 3;
	} else if (num_mappings == 3 || 
		   num_mappings == 4) {
	  mqual = 1;
	}

	
	if (num_mappings > 0) {
	  num_mapped_reads++;
	  total_mappings += num_mappings;
	  for (size_t j = 0; j < num_mappings; j++) {
	    alig = (alignment_t *) array_list_get(j, mapping_list);
	    
	    // update alignment
	    alig->secondary_alignment = 0;
	    if (num_mate_mappings != 1) {
	      alig->is_mate_mapped = 0;
	      alig->is_paired_end_mapped = 0;
	      alig->mate_strand = 0;
	    }

	    flag = 0;
	    if (alig->is_paired_end)                              flag += BAM_FPAIRED;
	    if (alig->is_paired_end_mapped)                       flag += BAM_FPROPER_PAIR;
	    if (!alig->is_seq_mapped)                             flag += BAM_FUNMAP;   
	    if ((!alig->is_mate_mapped) && (alig->is_paired_end)) flag += BAM_FMUNMAP;
	    if (alig->mate_strand)                                flag += BAM_FMREVERSE;
	    if (alig->pair_num == 1)	                        flag += BAM_FREAD1;
	    if (alig->pair_num == 2)                              flag += BAM_FREAD2;
	    if (alig->secondary_alignment)                        flag += BAM_FSECONDARY;
	    if (alig->fails_quality_check)                        flag += BAM_FQCFAIL;
	    if (alig->pc_optical_duplicate)                       flag += BAM_FDUP;
	    if (alig->seq_strand) {                               
	      flag += BAM_FREVERSE;
	      //seq = read->revcomp;
	    }

	    char alig_str[read_size];
	    alig_str[0] = '\0';
	    	    
	    sprintf(optional_flags, "AS:%i NM:%i NH:%lu", alig->map_quality, alig->mapq, num_mappings);
	    
	    //num_total_mappings++;
	    sprintf(alig_str, "%s\t%i\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%s\t%s\t%s\n", 
		    read->id,
		    flag,
		    genome->chrom_names[alig->chromosome - 1],
		    alig->position + 1,
		    mqual,
		    alig->cigar,
		    alig->mate_chromosome <= 0 ? "0" : genome->chrom_names[alig->mate_chromosome - 1],
		    alig->mate_position + 1,
		    alig->template_length,
		    alig->sequence,
		    alig->quality,
		    optional_flags
		    );
	    
	    strcpy(&buffer[buffer_size], alig_str);
	    buffer_size += strlen(alig_str);
	    
	    alignment_free(alig);
	    
	  }
	} else {
	  //num_unmapped_reads++;
	  char alig_str[read_size];
	  alig_str[0] = '\0';

	  sprintf(alig_str, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", 
		  read->id,
		  read->sequence,
		  read->quality
		  );

	  assert(buffer_size < buffer_max_size);
	  strcpy(&buffer[buffer_size], alig_str);
	  buffer_size += strlen(alig_str);            

	}
      
	array_list_free(mapping_list, (void *) NULL);
	
      }
    } else {      
      for (size_t i = 0; i < num_reads; i++) {
	read = (fastq_read_t *) array_list_get(i, read_list);
	mapping_list = sa_batch->mapping_lists[i];
	num_mappings = array_list_size(mapping_list);
	
	int mqual = 0;
	if (num_mappings == 1) { 
	  reads_uniq_mappings++; 
	  mqual = 50;
	} else if (num_mappings == 2) {
	  mqual = 3;
	} else if (num_mappings == 3 || 
		   num_mappings == 4) {
	  mqual = 1;
	}

	if (num_mappings > 0) {
	  num_mapped_reads++;
	  total_mappings += num_mappings;
	  for (size_t j = 0; j < num_mappings; j++) {
	    alig = (alignment_t *) array_list_get(j, mapping_list);
	    flag = (alig->seq_strand ? 16 : 0);


	    sprintf(optional_flags, "AS:%i NM:%i NH:%lu", alig->map_quality, alig->mapq, num_mappings);
	
	    char alig_str[read_size];
	    alig_str[0] = '\0';
	
	    sprintf(alig_str, "%s\t%i\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%s\t%s\t%s\n%c", 
		    alig->query_name,
		    flag,
		    genome->chrom_names[alig->chromosome - 1],
		    alig->position + 1,
		    mqual,
		    alig->cigar,
		    rnext,
		    pnext,
		    tlen,
		    alig->seq_strand ? read->revcomp : read->sequence,
		    //alig->sequence,
		    alig->quality,
		    optional_flags,
		    '\0'
		    );
	

	    assert(buffer_size < buffer_max_size);
	    strcpy(&buffer[buffer_size], alig_str);
	    buffer_size += strlen(alig_str);

	    //start_timer(time_free_s);
	    alignment_free(alig);
	    //stop_timer(time_free_s, time_free_e, time_free_alig);

	  }
	} else {
	  char alig_str[read_size];
	  alig_str[0] = '\0';
      
	  sprintf(alig_str, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n%c", 
		  read->id,
		  read->sequence,
		  read->quality,
		  '\0'
		  );
      
	  assert(buffer_size < buffer_max_size);
	  strcpy(&buffer[buffer_size], alig_str);
	  buffer_size += strlen(alig_str);            
	}

	array_list_free(mapping_list, (void *) NULL);
  
      }
    }

    wf_batch->data_output = buffer;
    wf_batch->data_output_size = buffer_size;

    // free memory
    sa_batch_free(sa_batch);  

    basic_statistics_add(total_reads, num_mapped_reads, total_mappings, reads_uniq_mappings, basic_st);

  } else {    
    for (size_t i = 0; i < num_reads; i++) {
      read = (fastq_read_t *) array_list_get(i, read_list);
      mapping_list = sa_batch->mapping_lists[i];
      num_mappings = array_list_size(mapping_list);
      
      alignment_t *alig;
      bam1_t *bam1;
      for (size_t j = 0; j < num_mappings; j++) {
	alig = (alignment_t *) array_list_get(j, mapping_list);	
	if (alig != NULL) {
	  bam1 = convert_to_bam(alig, 33);	
	  alignment_free(alig);	  
	  array_list_replace_at(j, bam1, mapping_list);
	} else {
	  LOG_FATAL_F("alig is NULL, num_items = %lu\n", num_mappings);
	}	
      }
      fastq_read_free(read);
    }

    if (sa_batch->fq_reads) { 
      array_list_free(sa_batch->fq_reads, (void *) NULL);
    }
    
    wf_batch->data_output      = sa_batch->mapping_lists; 
    wf_batch->data_output_size = num_reads;
    free(sa_batch);

  }

  return;

}

int write_to_file(void *data) {
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;
  if (!wf_batch->data_output_size) { return 0; }

  extern pthread_mutex_t mutex_sp;
  struct timeval time_free_s, time_free_e;
  extern double time_write, time_free;

  int bam_format = wf_batch->writer_input->bam_format;

  if (!bam_format) {
    FILE *out_file = (FILE *) wf_batch->writer_input->bam_file;
    fwrite((char *)wf_batch->data_output, sizeof(char), wf_batch->data_output_size, out_file);    
    free(wf_batch->data_output);
  } else {
    bam_file_t *out_file = (bam_file_t *) wf_batch->writer_input->bam_file;
    bam1_t *bam1;
    int num_reads = wf_batch->data_output_size;
    array_list_t **mapping_lists = wf_batch->data_output;
    array_list_t *mapping_list;
    size_t num_mappings;
    for (size_t i = 0; i < num_reads; i++) {
      mapping_list = mapping_lists[i];
      num_mappings = array_list_size(mapping_list);    
      for (size_t j = 0; j < num_mappings; j++) {
	bam1 = array_list_get(j, mapping_list);
	//start_timer(time_free_s);
	bam_fwrite(bam1, out_file);
	//stop_timer(time_free_s, time_free_e, time_write);
	//start_timer(time_free_s);
	bam_destroy1(bam1);
	//stop_timer(time_free_s, time_free_e, time_free);
      }      
      array_list_free(mapping_list, (void *)NULL);      
    }
    free(mapping_lists);
  }

  if (wf_batch) sa_wf_batch_free(wf_batch);

  return 0;

}


extern inline void parse_alignment_data(array_list_t *list, fastq_read_t *fq_read) {
  //start_timer(time_start);
  if (!array_list_size(list)) { return; }
  alignment_data_t *p = array_list_get(0, list);
  array_list_clear(list, (void *)NULL);

  int num_items = p->num_items;

  char cigars_test[num_items][1024];
  size_t actual_read = 0;
  simple_alignment_t *simple_a;
  for (int i = 0; i < num_items; i++) {
    //simple_a = &simple_alignment[i];
    simple_a = &p->simple_alignments_array[i];
    //memcpy(&cigars_test[i], &cigar_buffer[actual_read], simple_a->cigar_len);
    memcpy(&cigars_test[i], &p->cigars_str[actual_read], simple_a->cigar_len);
    cigars_test[i][simple_a->cigar_len] = '\0';
    actual_read += simple_a->cigar_len;
    
    //printf("CIGAR %i: %s\n", i, cigars_test[i]);
    size_t map_len = fq_read->length - simple_a->gap_start - simple_a->gap_end;
    size_t map_genome_len = 0;
    
    cigar_code_t *cc = cigar_code_new_by_string(cigars_test[i]);
    array_list_t *list_aux = array_list_new(5, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    sa_alignment_t *sa_alignment = sa_alignment_new(list_aux);    

    //sa_alignment->c_final = cigar_code_new();

    for (int m = 0; m < array_list_size(cc->ops); m++) {
      cigar_op_t *op = array_list_get(m, cc->ops);
      //cigar_code_append_new_op(op->number, op->name, sa_alignment->c_final);
      if (op->name == 'M' || op->name == 'D' || op->name == 'N') {
	map_genome_len += op->number;
      }
    }

    //printf("-->Simple_a->gap_start = %i\n", simple_a->gap_start);
    //printf("-->Simple_a->gap_start + map_len = %i\n", simple_a->gap_start + map_len);

    if (simple_a->gap_start == 0) {
      sa_alignment->left_close = 1;
    }
    
    if (simple_a->gap_start + map_len == fq_read->length) {
      sa_alignment->right_close = 1;
    }

    //printf("SEED := len_read:%i - gap_read:%i - gap_end:%i = %i, SEED-END = %i\n", fq_read->length, 
    //   simple_a->gap_start, 
    //   simple_a->gap_end, 
    //   map_len, simple_a->gap_start + map_len);

    seed_region_t *s_region = seed_region_new(simple_a->gap_start, 
                                              simple_a->gap_start + map_len - 1,
                                              simple_a->map_start, 
                                              simple_a->map_start + map_genome_len - 1,
                                              0, 0, 0);
    
    //printf("Exit with seed [%i:%i]\n", s_region->read_start, s_region->read_end);    
    linked_list_t *sr_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    linked_list_insert(s_region, sr_list);
    
    cal_t *cal = cal_new(simple_a->map_chromosome, 
                         simple_a->map_strand,
                         simple_a->map_start,
                         simple_a->map_start + map_len - 1,
                         1,
                         sr_list,
                         linked_list_new(COLLECTION_MODE_ASYNCHRONIZED));


    s_region->info = cc;
    cc->distance = simple_a->map_distance;
    
    //cal->info = cc;
    //printf("Cal & Cigar Ok, Insert list\n");    
    //sa_alignment->c_final->distance = simple_a->map_distance;
    
    array_list_insert(cal, sa_alignment->cals_list);
    array_list_insert(sa_alignment, list);
    
    //free(simple_a);
  }

  free(p->simple_alignments_array);
  free(p->cigars_str);
  free(p);

  //stop_timer(time_start, time_end, time_read_proc);

}



int sa_mapped_exact_reads(fastq_read_t *read, 
			  sa_index3_t *sa_index,
			  array_list_t *alignments_list,
			  float match, float mismatch, 
			  float gap_open, float gap_extend) {

  size_t low_p, high_p, suffix_len_p;
  size_t low_n, high_n, suffix_len_n;  
  size_t num_suffixes_p, num_suffixes_n;
  size_t g_start;
  int chrom;

  //******************************************//
  //****           1-Exact Reads          ****//
  //******************************************//
  //****                                  ****//
  //******************************************//

  //Search strand(+)
  num_suffixes_p = search_suffix(&read->sequence[0], sa_index->k_value, 
				 MAX_NUM_SUFFIXES, sa_index, 
				 &low_p, &high_p, &suffix_len_p);
    
  //Search strand(-)
  num_suffixes_n = search_suffix(&read->revcomp[0], sa_index->k_value, 
				 MAX_NUM_SUFFIXES, sa_index, 
				 &low_n, &high_n, &suffix_len_n);
    

  // exact search for exact reads
  if (suffix_len_p == read->length) {
    //Report Exact Maps! (+)
    alignment_t *alignment;
    
    int n_alig = high_p - low_p + 1;
    if (n_alig > MAX_ALIG) {
      n_alig = MAX_ALIG;
    }
    
    size_t suff = low_p;
    for (size_t a = 0; a < n_alig; a++, suff++) {
      chrom = sa_index->CHROM[suff];	
      g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];
      
      //char cigar_str[2048];
      //sprintf(cigar_str, "%i%c", read->length, 'M');     
      alignment = alignment_new();
      alignment->mapq = 0;
      cigar_code_t *cc = cigar_code_new();
      array_list_insert(cigar_op_new(read->length, 'M'), cc->ops);
      int score = get_score(cc, read, match, mismatch, gap_open, gap_extend);

      alignment_init_single_end(strdup(read->id), 
				strdup(read->sequence), 
				strdup(read->quality), 
				0, chrom + 1,
				g_start,
				NULL, 
				1, 
				score, 1, (n_alig > 1),
				0, 0, alignment);
      alignment->alig_data = cc;

      array_list_insert(alignment, alignments_list);
    }
  } 
  
  if (suffix_len_n == read->length) {
    //Report Exact Map! (-)
    alignment_t *alignment;

    int n_alig = high_n - low_n + 1;
    if (n_alig > MAX_ALIG) {
      n_alig = MAX_ALIG;
    }
    
    if (n_alig > 0 && !read->rev_quality) {
      read->rev_quality = str_reverse(read->quality);
    }
    
    size_t suff = low_n;
    for (size_t a = 0; a < n_alig; a++, suff++) {
      chrom = sa_index->CHROM[suff];
      g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];
      
      //char cigar_str[2048];
      //sprintf(cigar_str, "%i%c", read->length, 'M');
      alignment = alignment_new();
      alignment->mapq = 0;
      cigar_code_t *cc = cigar_code_new();
      array_list_insert(cigar_op_new(read->length, 'M'), cc->ops);
      int score = get_score(cc, read, match, mismatch, gap_open, gap_extend);

      alignment_init_single_end(strdup(read->id),
				strdup(read->revcomp), 
				strdup(read->rev_quality),
				1, chrom + 1,
				g_start, 
				NULL,
				1, score, 1, (n_alig > 1),
				0, 0, alignment);
      alignment->alig_data = cc;

      array_list_insert(alignment, alignments_list);
    }
  }
  

  return array_list_size(alignments_list);

}


int sa_generate_cals(fastq_read_t *read, 
		     cal_mng_t *cal_mng_p,
		     cal_mng_t *cal_mng_n,
		     sa_index3_t *sa_index,
		     array_list_t *cals_list) {
  
  //array_list_set_flag(1, sa_alignments_list[r]);      
  //cal_mng_simple_clear(cal_mng_p);//(+)
  //cal_mng_simple_clear(cal_mng_n);//(-)

  array_list_clear(cals_list, (void *)NULL);

  array_list_t *merge_cals;
  int seed_size = sa_index->k_value;
  int seed_inc  = seed_size / 2;
  int read_pos;
  int id_seed = 2;
  cal_mng_t *cal_mng;
  size_t high, low, suffix_len;
  float score;
  int cal_pos;
  seed_region_t *s_prev, *s;
  linked_list_item_t *item;
  size_t num_suffixes;
  size_t g_start;
  int chrom;
  cal_t *cal_prev, *cal_next;
  char *query;
  linked_list_t *cal_list;
  size_t num_cals;

  //===== Seeding Strategy =====//
  int max_read_pos = read->length - seed_size - 1;
  for (int s = 0; s < 2; s++) { //Strand
    read_pos  = 0;
    if (!s) {
      query = read->sequence;
      cal_mng = cal_mng_p;
    } else {
      query = read->revcomp;
      cal_mng = cal_mng_n;
    }
    
    while (read_pos + sa_index->k_value < read->length) {
      num_suffixes = search_suffix(&query[read_pos], sa_index->k_value, 
				   MAX_NUM_SUFFIXES, sa_index, 
				   &low, &high, &suffix_len);
      
	
      if (suffix_len && num_suffixes) {
	//Storage Mappings
	for (size_t suff = low; suff <= high; suff++) {	
	  chrom = sa_index->CHROM[suff];
	  g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom] + 1;
	  
	  //printf("\tSTORE SEED %i:[%lu|%i-%i|%lu]\n", chrom, g_start, read_pos, read_pos + suffix_len - 1, g_start + suffix_len - 1);
	  generate_cals(chrom + 1, s, 
			g_start, g_start + suffix_len - 1, 
			read_pos, read_pos + suffix_len - 1,
			cal_mng->cals_lists[chrom], read->length, id_seed++);
	}
	////////////////////////////////////////
	//|If seed map, extend suffix_len + 1|//
	//|__________________________________|//
	//|---(18)---|--(10)--|              |//
	//|        (28)       |---(18)---|   |//
	//|__________________________________|//
	////////////////////////////////////////
	read_pos += suffix_len + 1;
      } else {
	////////////////////////////////////
	//|  Else, extend seed_size / 2  |//
	//|______________________________|//
	//|---(18)---|                   |//  
	//| (9) |---(18)---|             |//
	//|     | (9) |---(18)---|       |//
	//|______________________________|//
	////////////////////////////////////
	read_pos += seed_inc; 
      }
    } //loop seeds	
    

    if (read->length - 1 != read_pos) {
      read_pos = read->length - sa_index->k_value - 1;
      num_suffixes = search_suffix(&query[read_pos], sa_index->k_value, 
				   MAX_NUM_SUFFIXES, sa_index, 
				   &low, &high, &suffix_len);	      
	
      if (suffix_len && num_suffixes) {
	//Storage Mappings
	for (size_t suff = low; suff <= high; suff++) {	
	  chrom = sa_index->CHROM[suff];
	  g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom] + 1;
	  generate_cals(chrom + 1, s, 
			g_start, g_start + suffix_len - 1, 
			read_pos, read_pos + suffix_len - 1,
			cal_mng->cals_lists[chrom], read->length, id_seed++);
	}
      }
    }
  } //loop strands 
  
    //Merge CALs and select new targets
  for (int st = 0; st < 2; st++) {
    if (!st) { cal_mng = cal_mng_p; }
    else  { cal_mng = cal_mng_n; }
	
    for (unsigned int i = 0; i < cal_mng->num_chroms; i++) {
      cal_list = cal_mng->cals_lists[i];
      num_cals = linked_list_size(cal_list);
      //printf("(%i)%i:Num CAls %i\n", st, i, num_cals);
      if (num_cals) {
	//=====*---------------------*=====//
	merge_cals = array_list_new(100,
				    1.25f,
				    COLLECTION_MODE_ASYNCHRONIZED);
	cal_pos = 1;
	item = cal_list->first;
	cal_prev = item->item;

	array_list_insert(cal_prev, merge_cals); 
	while (cal_pos < num_cals) {

	  item = item->next;
	  cal_next = item->item;
		  
	  s_prev = linked_list_get_last(cal_prev->sr_list);
	  s = linked_list_get_first(cal_next->sr_list);
		  
	  assert(s_prev != NULL);
	  assert(s != NULL);
	      
	  //cal_print(cal_prev);
	  //cal_print(cal_next);
	  //printf("s_prev->id = %i != s->id = %i\n", s_prev->id, s->id);

	  int ok = 0;
	  if (s->read_start > s_prev->read_end) {
	    ok = 1;
	  } else if (s_prev->read_end - s->read_start <= 5) {
	    ok = 1;
	  }

	  if (ok &&
	      cal_prev->chromosome_id == cal_next->chromosome_id && 
	      cal_prev->strand == cal_next->strand && 
	      //abs(s->read_start - s_prev->read_end) <= 5 &&
	      s_prev->id != s->id &&
	      cal_next->start <= (cal_prev->end + max_intron)) {
	    //printf("Merge!! cal_prev->end = %lu, cal_next->start = %lu\n", cal_prev->end, cal_next->start);
	    array_list_insert(cal_next, merge_cals);
	  } else {
	    array_list_insert(merge_cals, cals_list);
	    merge_cals = array_list_new(10,
					1.25f,
					COLLECTION_MODE_ASYNCHRONIZED);
	    array_list_insert(cal_next, merge_cals);
	  }                                                                                  
	  cal_prev = cal_next;
	  cal_pos++;
	}
	array_list_insert(merge_cals, cals_list);	    
	//=====*---------------------*=====//
      }
    }
  } 
  
  return array_list_size(cals_list);

}


int filter_cals_by_score(array_list_t *cals_list, 
			 const int limit_cals, 
			 float *cals_score,
			 fastq_read_t *read,
			 int min_cal_score) {

  int num_targets = array_list_size(cals_list);
  array_list_t *merge_cals;
    
  if (num_targets > limit_cals) {
    //Free other CALs 
    for (int i = num_targets - 1; i >= limit_cals; i--) {
      merge_cals = array_list_remove_at(i, cals_list);
      array_list_free(merge_cals, (void *)NULL);
    }
    num_targets = limit_cals;
  }
  
  //=====*---------------------*=====//
  for (int i = 0; i < num_targets; i++) {
    merge_cals = array_list_get(i, cals_list);
    cals_score[i] = generate_cals_merge_score(merge_cals, read->length);
  }
  //=====*---------------------*=====//    
  
  //=====*---------------------*=====//
  for (int i = 0; i < num_targets - 1; i++) {
    for (int j = i + 1; j < num_targets; j++) {
      if (cals_score[j] > cals_score[i]) {
	array_list_swap(i, j, cals_list);
	float score = cals_score[j];
	cals_score[j] = cals_score[i];
	cals_score[i] = score;
      }
    }
  }
  //=====*---------------------*=====//
  
    
  //Report the first best 10 CALs
  int n_report = num_targets < 5 ? num_targets : 5;
    
  int report_tmp = 1;
  for (int i = 1; i < n_report; i++) {
    if (cals_score[i] < min_cal_score) {
      break;
    }
    report_tmp++;
    //printf("CAL %i: %f\n", i, cals_score[i]);
  }
    
  n_report = report_tmp;
  
  //Free other CALs
  for (int i = num_targets - 1; i >= n_report; i--) {
    merge_cals = array_list_remove_at(i, cals_list);
    array_list_free(merge_cals, (void *)NULL);
  } 
  
  return array_list_size(cals_list);

}


array_list_t *create_list(size_t *valid_items, size_t num_valids, array_list_t *list);

void update_mispaired_pairs(size_t num_items1, size_t num_items2,
			    array_list_t *list1, array_list_t *list2);
void filter_alignments(char report_all,
		       size_t report_n_best, 
		       size_t report_n_hits,
		       int report_best,
		       array_list_t *mapping_list);

void update_mispaired_pair(int pair_num, size_t num_items, array_list_t *list);

void sa_complete_pairs(sa_wf_batch_t *wf_batch) {
  sa_batch_t *batch = wf_batch->mapping_batch;
  sa_rna_input_t *sa_rna = wf_batch->data_input;
  pair_server_input_t *pair_input = sa_rna->pair_input;
  pair_mng_t *pair_mng = pair_input->pair_mng;
  report_optarg_t *report_optarg = pair_input->report_optarg;

  size_t num_items1, num_items2, num_pairs, num_reads = array_list_size(batch->fq_reads);

  int distance;
  //int min_distance = batch->options->pair_min_distance;
  //int max_distance = batch->options->pair_max_distance;
  //int pair_mode = batch->options->pair_mode;

  int min_distance = pair_mng->min_distance;
  int max_distance = pair_mng->max_distance;
  int pair_mode = pair_mng->pair_mode;

  int report_only_paired = report_optarg->only_paired;

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
  int n_best = report_optarg->n_best;
  int n_hits = report_optarg->n_hits;
  int all = report_optarg->all;
  int best = report_optarg->best;

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
	  if (alig1->chromosome <= 0 || alig1->chromosome <= 0) { exit(-1); };

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

cal_mng_t * cal_rna_mng_new(sa_genome3_t *genome) {

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

void cal_rna_mng_free(cal_mng_t *p) {
  if (p) {
    if (p->cals_lists) {
      for (unsigned int i = 0; i < p->num_chroms; i++) {
	if (p->cals_lists[i]) {
	  linked_list_free(p->cals_lists[i], (void *)seed_cal_free);
	}
      }
      free(p->cals_lists);
    }
    free(p);
  }
}

//--------------------------------------------------------------------

int sa_rna_mapper(void *data) {
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;
  sa_batch_t *sa_batch = wf_batch->mapping_batch;
  sa_index3_t *sa_index = (sa_index3_t *) wf_batch->sa_index;
  sa_rna_input_t *sa_rna = wf_batch->data_input;
  genome_t *genome = (genome_t *)sa_rna->genome;
  avls_list_t *avls_list = (avls_list_t *)sa_rna->avls_list;
  metaexons_t *metaexons = (metaexons_t *)sa_rna->metaexons;
  sw_optarg_t *sw_optarg = sa_rna->sw_optarg;
  FILE *file1 = sa_rna->file1;
  FILE *file2 = sa_rna->file2;
  pair_server_input_t *pair_input = sa_rna->pair_input;
  int pair_mode = pair_input->pair_mng->pair_mode;
  report_optarg_t *report_optarg = pair_input->report_optarg;
  
  size_t num_reads = sa_batch->num_reads;
  int max_read_area, min_num_mismatches;
  float max_score;

  MAX_ALIG = sa_rna->max_alig;
  min_score = sa_rna->min_score;
  float gap_open = sw_optarg->gap_open;
  float gap_extend = sw_optarg->gap_extend;
  float match = sw_optarg->subst_matrix['A']['A'];
  float mismatch = sw_optarg->subst_matrix['A']['C'];
  
  // CAL management
  //size_t num_cals;
  //seed_cal_t *cal;
  cal_mng_t *cal_mng;
  //array_list_t *cal_list;
  fastq_read_t *read;  
  //  int saved_pos[2][1024];
  // TODO !!! 20 = min. cal size
  int min_cal_size = sa_rna->cal_optarg->min_cal_size;

  int seed_size = 18;
  //Function old Seeding vars
  int min_cal_score = sa_rna->cal_optarg->min_cal_size;
  int len_seq;
  
  const int MAX_SUFFIXES    = 50;
  const int MIN_SEED_SIZE   = 30;
  const int MIN_COVER       = len_seq * MIN_SEED_SIZE / 100;
  
  double prefix_time, suffix_time;
  int suffix_len0, suffix_len1, num_suffixes, num_suffixes_start, num_suffixes_end;
  size_t low0, high0, low1, high1;
  int read_pos, read_end_pos;
  char reference[2048];
  char *seq;
  char *seq_revcomp;
  char *query;
  int chrom;
  size_t g_start, g_end;
  seed_region_t *seed_region;

  cal_mng_t * cal_mng_p = cal_rna_mng_new(sa_index->genome); //(+)
  cal_mng_t * cal_mng_n = cal_rna_mng_new(sa_index->genome); //(-)

  size_t low_p, high_p, suffix_len_p;
  size_t low_n, high_n, suffix_len_n;  
  size_t low, high, suffix_len;

  size_t num_suffixes_p, num_suffixes_n;
  size_t genome_start, genome_end;
  array_list_t *target_cals = array_list_new(100, 1.25f,
					     COLLECTION_MODE_ASYNCHRONIZED);

  //CALs done! Search CALs with seeds (Double Anchors)
  cal_t *cal;
  linked_list_iterator_t itr;
  linked_list_t *cal_list;
  array_list_t *merge_cals;
  cal_t *cal_prev, *cal_next;
  seed_region_t *seed_first, *seed_last, *seed_first_next;
  size_t num_items, num_cals;
  cal_mng_t *p;
  array_list_t *alignments_list;
  int total_suffix;
  array_list_t *sa_alignments_list[num_reads];
  
  int kk;

  int distance;
  avl_node_t *node_avl_start, *node_avl_end;
  size_t sp_start, sp_end;
  int sp_type;

  sw_multi_output_t *output = sw_multi_output_new(MAX_DEPTH);
  sa_sw_depth_t sw_depth;
  sw_depth.depth = 0;

  sa_alignment_t *sa_alignment_aux;
  int delete_targets[num_reads];
  
  //array_list_t *write_alignments = array_list_new(50, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  seed_region_t *seed_prev, *seed_next, *s_prev;
  int num_sa_alignments;

  const int limit_cals = 200;
  float cals_score[limit_cals];

  char *read_quality_rev = NULL;

  //printf("============================ X =======================\n");

  for (int r = 0; r < num_reads; r++) {
    delete_targets[r] = 0;
    read = array_list_get(r, sa_batch->fq_reads);
    
    //printf("READ %s\n", read->id);

    //Rev-comp
    fastq_read_revcomp(read);
    
    seq = read->sequence;
    seq_revcomp = read->revcomp;

    alignments_list = sa_batch->mapping_lists[r];
    total_suffix = 0;
    len_seq = read->length;

    int tot_mappings = sa_mapped_exact_reads(read, 
					     sa_index,
					     alignments_list,
					     match, mismatch,
					     gap_open, gap_extend);

    if (tot_mappings) {
      continue;
    }

    sa_alignments_list[r] = array_list_new(50, 1.25f,
					   COLLECTION_MODE_ASYNCHRONIZED);

    
    //CALL GENERATE CALS
    int n_report = sa_generate_cals(read, 
				    cal_mng_p,
				    cal_mng_n,
				    sa_index,
				    target_cals);
    
    if (!n_report) { goto cal_mng; }

    n_report = filter_cals_by_score(target_cals,
				    limit_cals, 
				    cals_score,
				    read,
				    min_cal_score);
       
    if (!n_report) { goto cal_mng; }

    int cal_len_prev, cal_len_next, cal_len;
    int fill_type;
    seed_region_t *s_prev_first, *s_prev_last, *s_next_first, *s_next_last;
    cigar_code_t *cc_left, *cc_right;
    int full_left, full_right;
    

    // GENERATE ALIGNMENTS WITH METAEXON
    for (int i = 0; i < n_report; i++) {
      merge_cals = array_list_get(i, target_cals);
      num_cals = array_list_size(merge_cals);

      cal_prev     = array_list_get(0, merge_cals);
      s_prev       = cal_prev->sr_list->first->item; 
      
      cal_next     = array_list_get(num_cals - 1, merge_cals);

      seed_region_t *s_next       = cal_next->sr_list->first->item; 
      sa_alignment_t *sa_alignment = sa_alignment_new(merge_cals);
      array_list_insert(sa_alignment, sa_alignments_list[r]);

      cal_t *selected_cal =NULL;

      int cal_len_prev = cal_prev->end - cal_prev->start + 1;
      int cal_len_next = cal_next->end - cal_next->start + 1;
      selected_cal = cal_len_prev > cal_len_next ? cal_prev : cal_next;

      if (cal_prev->strand) {
	query = seq_revcomp;
      } else {
	query = seq;
      }

      if (s_prev->read_start != 0) {
	//selected_cal = cal_prev;
	sa_alignment->c_left = meta_alignment_fill_extrem_gap(query, 
							      selected_cal,
							      FILL_GAP_LEFT,
							      genome,
							      metaexons, 
							      avls_list);
	//sa_alignment->right_close = 1;  
	//printf("LEFT SEARCH: %s\n", new_cigar_code_string(sa_alignment->c_left));
	if (sa_alignment->c_left) {
	  sa_alignment->left_close = 1;
	}
      }

      if (s_next->read_end != read->length - 1) {
	//selected_cal = cal_next;
	sa_alignment->c_right = meta_alignment_fill_extrem_gap(query, 
							       selected_cal,
							       FILL_GAP_RIGHT,
							       genome,
							       metaexons, 
							       avls_list);	  
	//sa_alignment->left_close = 1;  
	//printf("RIGHT SEARCH: %s\n", new_cigar_code_string(sa_alignment->c_right));
	if (sa_alignment->c_right) {
	  sa_alignment->right_close = 1;
	}
      }

      if (sa_alignment->right_close && sa_alignment->left_close) {
	sa_alignment->complete = 1;
	merge_cals->items[0] = selected_cal;
	merge_cals->size = 1;
      }                
    
      
      if (!sa_alignment->complete) {
	//Close gap righ... Search SP in CALs
	//Write to file	
	//printf("Num cals :::: %i\n", num_cals);
	cal_prev = array_list_get(0, merge_cals);
	for (int c = 1; c < num_cals; c++) {
	  cal_next = array_list_get(c, merge_cals);
	  
	  seed_first = linked_list_get_last(cal_prev->sr_list);
	  seed_last  = linked_list_get_first(cal_next->sr_list);
	  	  
	  int nt = search_simple_splice_junction(seed_first, seed_last,
						 cal_next->chromosome_id,
						 cal_next->strand,
						 query, genome, 
						 &sp_start,
						 &sp_end,
						 &sp_type,
						 &distance);
	  
	  if (!nt) {
	    nt = search_simple_splice_junction_semi_cannonical(seed_first, seed_last,
							       cal_next->chromosome_id, 
							       cal_next->strand, 
							       query, genome, 
							       &sp_start, 
							       &sp_end,
							       &sp_type,
							       &distance);	    
	  }

	  if (nt >= min_intron && 
	      sp_start < sp_end && 
	      seed_first->genome_start < sp_start && 
	      sp_end < seed_last->genome_end) {
	    //printf("distance sp %i, nt = %i\n", distance, nt);
	    sa_alignment->sp_middle_err[sa_alignment->num_sp] = distance;
	    sa_alignment->sp_middle[sa_alignment->num_sp++] = nt;
	    cal_prev->end   = seed_first->genome_end;
	    cal_next->start = seed_last->genome_start;
	  } 
	  
	  //Metaexon actualization...
	  cal_prev = cal_next;
	  
	}
      }
    }

	  
    //Fill gaps CALs
    num_sa_alignments = array_list_size(sa_alignments_list[r]);
    //printf("num_sa_alignments = %i (fill gaps)\n", num_sa_alignments);
    for (int n_sa = 0; n_sa < num_sa_alignments; n_sa++) {
      sa_alignment_t *sa_alignment = array_list_get(n_sa, sa_alignments_list[r]);
      merge_cals = sa_alignment->cals_list;
      
      int num_cals = array_list_size(merge_cals);
      int ok = 1;
      for (int c = 0; c < num_cals; c++) {
	cal = array_list_get(c, merge_cals);
	if (!validate_cal(cal)) { 
	  //Purge CAL. 
	  ok = 0; 
	  break;
	}
      }
      
      if (!ok) { 
	//Purge CAL
	continue;
      }
    
      if (sa_alignment->complete) {
	num_cals = 1;
      } else {
	num_cals = merge_cals->size;
      }
      
      //Close gaps intra CALs
      //----------- generate cigar for each seed ----------//
      int n_sp = 0;
      for (int c = 0; c < num_cals; c++) {
	cal = array_list_get(c, merge_cals);
	if (cal->strand) {
	  query = seq_revcomp;
	} else {
	  query = seq;
	}

	cal_fill_gaps(cal,
		      query,
		      genome,
		      metaexons,
		      avls_list, 
		      &sw_depth, 
		      sw_optarg, 
		      output);
      }
    }
    
    sw_insert_item(NULL, NULL, 0, NULL, 
		   sw_optarg, output, &sw_depth);
    
    //} //End if no alignments found and no metalignments

    //Process metalignments found
    num_sa_alignments = array_list_size(sa_alignments_list[r]);

    //printf("num_sa_alignments = %i (merge cigars)\n", num_sa_alignments);
    for (int n_sa = 0; n_sa < num_sa_alignments; n_sa++) {
      sa_alignment_t *sa_alignment = array_list_get(n_sa, sa_alignments_list[r]);
      merge_cals = sa_alignment->cals_list;
      int num_cals = array_list_size(merge_cals);

      cigar_code_t *cc_new = cigar_code_new();
      cigar_code_t *cc_middle = cigar_code_new(); 
      int n_sp = 0;
      for (int c = 0; c < num_cals; c++) {
	cal = array_list_get(c, merge_cals); 
	linked_list_item_t *item = cal->sr_list->first;
	seed_region_t *sr;
	while (item != NULL) {
	  sr = item->item;
	  cigar_code_t *cc = sr->info;
	  //printf("%i:%lu-%lu\n", cal->chromosome_id, sr->genome_start, sr->genome_end);
	  if (cc) {	      
	    //printf("\t---------------------->[%i-%i] : %s\n", sr->read_start, sr->read_end, new_cigar_code_string(cc));
	    for (int tt = 0; tt < cc->ops->size; tt++) {
	      cigar_op_t *op = array_list_get(tt, cc->ops);
	      cigar_code_append_op(op, cc_middle);
	    }
	    cc_middle->distance += cc->distance;
	    cigar_code_free(cc);
	    sr->info = NULL;
	  } 
	  item = item->next;
	}
	
	//printf("%i < %i\n", n_sp , sa_alignment->num_sp);
	if (n_sp < sa_alignment->num_sp) {
	  cc_middle->distance += sa_alignment->sp_middle_err[n_sp];
	  cigar_code_append_new_op(sa_alignment->sp_middle[n_sp++], 'N', cc_middle);
	}
	
      }

      //if (strcmp("seq.5a", read->id) == 0) {             
      //printf("MIDDLE CIGAR: %s\n", new_cigar_code_string(cc_middle));      
      //}

      //---------------------------------------------//	
      //Left cigar operations
      if (sa_alignment->c_left) {
	//printf("- left : %s\n", new_cigar_code_string(sa_alignment->c_left));
	//Actualization start position
	int dsp = 0;

	cigar_code_t *cc_aux = sa_alignment->c_left;
	for (int tt = 0; tt < cc_aux->ops->size; tt++) {
	  cigar_op_t *op = array_list_get(tt, cc_aux->ops);
	  dsp += op->number;
	  cigar_code_append_op(op, cc_new);
	}
	cc_new->distance += cc_aux->distance;
	cigar_code_free(cc_aux);
	sa_alignment->c_left = NULL;	    
	
	cal = array_list_get(0, merge_cals);	      
	seed_region_t *seed_aux = linked_list_get_first(cal->sr_list);
	seed_aux->genome_start -= dsp;//seed_aux->read_start;
	cal->start -= dsp;//seed_aux->read_start;
	seed_aux->read_start = 0;	
      }
	
      //Insert middle cigar operations
      for (int tt = 0; tt < cc_middle->ops->size; tt++) {
	cigar_op_t *op = array_list_get(tt, cc_middle->ops);
	cigar_code_append_op(op, cc_new);
      }
      cc_new->distance += cc_middle->distance;
      cigar_code_free(cc_middle);
      
      //Right cigar operations
      if (sa_alignment->c_right) {
	//printf("- right : %s\n", new_cigar_code_string(sa_alignment->c_right));
	int dsp = 0;
	for (int tt = 0; tt < sa_alignment->c_right->ops->size; tt++) {
	  cigar_op_t *op = array_list_get(tt, sa_alignment->c_right->ops);
	  //dsp += op->number;
	  //printf("loop %d + %i\n", dsp, op->number);
	  cigar_code_append_op(op, cc_new);
	}
	cc_new->distance += sa_alignment->c_right->distance;
	cigar_code_free(sa_alignment->c_right);
	sa_alignment->c_right = NULL;	  
	
	cal = array_list_get(array_list_size(merge_cals) - 1, merge_cals);	      
	seed_region_t *seed_aux = linked_list_get_last(cal->sr_list);
	//int len_aux = read->length - seed_aux->read_end - 1;
	//printf("increment right in %d = %i\n", dsp, seed_aux->genome_end);
	//seed_aux->genome_end += dsp;
	//printf("NEW! increment right in %d = %i\n", dsp, seed_aux->genome_end);
	//cal->end += dsp;
	seed_aux->read_end = read->length - 1;
      }
      
      sa_alignment->c_final = cc_new;

      //}      
      //printf(" <--------------- CIGAR %s -------------------->\n", new_cigar_code_string(sa_alignment->c_final));
      
      if (cigar_code_validate(read->length, sa_alignment->c_final)) {
	//Cigar generator
	//printf("Data: \n");
	//printf("\trl = %i - num_ops = %i\n", read->length, sa_alignment->c_final->ops->size );
	//char cigar_str[2048] = "\0";
	//for (int k = 0; k < sa_alignment->c_final->ops->size; k++) {
	//cigar_op_t *op = array_list_get(k, sa_alignment->c_final->ops);	      
	//sprintf(cigar_str, "%s%i%c", cigar_str, op->number, op->name);
	//}
	  
	//Report Aligment
	cal = array_list_get(0, merge_cals);	      
	//char *cigar_str = cigar_code_find_and_report_sj(cal->start - 1, sa_alignment->c_final, cal->chromosome_id, 
	//						cal->strand, avls_list, metaexons, genome, read);	

	char *quality_tmp;
	if (cal->strand) {
	  if (!read_quality_rev) {
	    read_quality_rev = str_reverse(read->quality);
	  }
	  quality_tmp = read_quality_rev;
	} else {
	  quality_tmp = read->quality;
	}
	
	//int score = SCORE(read->length, sa_alignment->c_final->distance);	
	int score = get_score (sa_alignment->c_final, read, match, mismatch, gap_open, gap_extend);

	if (score > min_score) {
	  alignment_t *alignment = alignment_new();
	  alignment_init_single_end(strdup(read->id), 
				    strdup(read->sequence), 
				    strdup(quality_tmp),
				    cal->strand, 
				    cal->chromosome_id,
				    cal->start - 1,        
				    NULL, 
				    sa_alignment->c_final->ops->size, 
				    score, 1, 0,
				    0, 0, alignment);

	  alignment->mapq = sa_alignment->c_final->distance;
	  alignment->alig_data = cigar_code_dup(sa_alignment->c_final);

	  //Report SJ	
	  array_list_insert(alignment, alignments_list);
	  sa_alignment->reported = 1;
	}

      }      
    }
    
    //printf("::::::: %i && %i(for)\n", array_list_size(alignments_list), num_sa_alignments);
    if (!array_list_size(alignments_list)) {
      delete_targets[r] = 1;
      num_sa_alignments = array_list_size(sa_alignments_list[r]);
      //array_list_clear(write_alignments, (void *)NULL);
      for (int n_sa = 0; n_sa < num_sa_alignments; n_sa++) {
	sa_alignment_t *sa_alignment = array_list_get(n_sa, sa_alignments_list[r]);
	if (cigar_code_get_num_ops(sa_alignment->c_final) == 0) { 
	  continue;
	}
	sa_alignment_partial_t *sa_partial = sa_alignment_partial_new(sa_alignment,
								      read->length);
	//printf("(%p|%p)%i: SA-PARTIAL.CIGAR: %s\n", sa_partial, sa_partial->cigar, r, sa_partial->cigar);
	array_list_insert(sa_partial, alignments_list);
	//array_list_insert(sa_alignment, write_alignments);      
      }
    }
    
    num_sa_alignments = array_list_size(sa_alignments_list[r]);
    for (int n_sa = 0; n_sa < num_sa_alignments; n_sa++) {
      sa_alignment_t *sa_alignment = array_list_get(n_sa, sa_alignments_list[r]);      
      merge_cals = sa_alignment->cals_list;
      int num_cals = array_list_size(merge_cals);
      cal_t *cal_aux = array_list_get(0, merge_cals);
      assert(cal_aux != NULL);
      //Free sa_alignments
      for (int k = 0; k < num_cals; k++) {
	cal = array_list_get(k, merge_cals);
	linked_list_item_t *item = cal->sr_list->first;
	seed_region_t *sr;	    
	while (item != NULL) {
	  sr = item->item;
	  cigar_code_t *cc = sr->info;
	  if (cc) { 
	    cigar_code_free(cc); 
	  }
	  item = item->next;
	}
      }      
      if (sa_alignment->c_left) {
	cigar_code_free(sa_alignment->c_left);
      }
      if (sa_alignment->c_right) {
	cigar_code_free(sa_alignment->c_right);
      }      
      if (sa_alignment->c_final) {
	array_list_clear(sa_alignment->c_final->ops, (void *)cigar_op_free);
	cigar_code_free(sa_alignment->c_final);	
      }
      array_list_free(merge_cals, (void *)NULL);
      sa_alignment_free(sa_alignment);      
    }

  cal_mng:
    //printf("Clear CALs\n");
    cal_mng_simple_clear(cal_mng_p);//(+)
    cal_mng_simple_clear(cal_mng_n);//(-)

    array_list_free(sa_alignments_list[r], (void *)NULL);
    array_list_clear(target_cals, (void *)NULL);

    free(read_quality_rev);
    read_quality_rev = NULL;
    
  } //End for reads
  

  array_list_t *fq_reads_aux       = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t **mapping_lists_aux = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
  array_list_t *alignments_list_aux;
  int n_reads = 0;
  extern pthread_mutex_t mutex_sp;
  extern size_t total_reads_ph2;

  //printf("NUM READS %i:\n", num_reads);
  if (pair_mode == SINGLE_END_MODE) {
    for (int r = 0; r < num_reads; r++) {
      read = array_list_get(r, sa_batch->fq_reads);
      alignments_list_aux = sa_batch->mapping_lists[r];
      if (delete_targets[r]) {
	//Write to buffer
	pthread_mutex_lock(&mutex_sp);
	sa_file_write_items(SA_PARTIAL_TYPE, read, alignments_list_aux, file1);
	total_reads_ph2++;
	pthread_mutex_unlock(&mutex_sp);
      
	fastq_read_free(read);
	array_list_free(alignments_list_aux, (void *)sa_alignment_partial_free);
      } else {
	array_list_insert(read, fq_reads_aux);
	mapping_lists_aux[n_reads++] = alignments_list_aux;
      }
    }
  } else {
    for (int r = 0; r < num_reads; r++) {
      if (r % 2 == 0) {
	if (delete_targets[r] || delete_targets[r + 1]) {
	  //Write to buffer
	  int type;
	  void *function_callback;
	  
	  if (!delete_targets[r]) {
	    type = SA_ALIGNMENT_TYPE;
	    function_callback = alignment_free;
	  } else {
	    type = SA_PARTIAL_TYPE;
	    function_callback = sa_alignment_partial_free;
	  }
	  read = array_list_get(r, sa_batch->fq_reads);
	  pthread_mutex_lock(&mutex_sp);
	  //printf("P0: %i (%i|%i): \n", type, SA_ALIGNMENT_TYPE, SA_PARTIAL_TYPE);
	  sa_file_write_items(type, read, sa_batch->mapping_lists[r], file1);
	  pthread_mutex_unlock(&mutex_sp);	  
	  fastq_read_free(read);
	  
	  array_list_free(sa_batch->mapping_lists[r], (void *)function_callback);	  	  
	  
	  if (!delete_targets[r + 1]) {
	    type = SA_ALIGNMENT_TYPE;
	    function_callback = alignment_free;
	  } else {
	    type = SA_PARTIAL_TYPE;
	    function_callback = sa_alignment_partial_free;
	  }
	  read = array_list_get(r + 1, sa_batch->fq_reads);
	  pthread_mutex_lock(&mutex_sp);
	  //printf("P1: %i (%i|%i): \n", type, SA_ALIGNMENT_TYPE, SA_PARTIAL_TYPE);
	  sa_file_write_items(type, read, sa_batch->mapping_lists[r + 1], file1);
	  pthread_mutex_unlock(&mutex_sp);
	  fastq_read_free(read);
	
	  pthread_mutex_lock(&mutex_sp);
	  total_reads_ph2 += 2;
	  pthread_mutex_unlock(&mutex_sp);	  
  
	  array_list_free(sa_batch->mapping_lists[r + 1], (void *)function_callback);
	  
	} else {
	  read = array_list_get(r, sa_batch->fq_reads);
	  array_list_insert(read, fq_reads_aux);
	  mapping_lists_aux[n_reads++] = sa_batch->mapping_lists[r];
	  
	  read = array_list_get(r + 1, sa_batch->fq_reads);
	  array_list_insert(read, fq_reads_aux);
	  mapping_lists_aux[n_reads++] = sa_batch->mapping_lists[r + 1];
	}
      }
    }
  }
  
  
  //array_list_free(write_alignments, (void *)NULL);
  array_list_free(sa_batch->fq_reads, (void *)NULL);
  free(sa_batch->mapping_lists);

  sa_batch->fq_reads = fq_reads_aux;
  sa_batch->mapping_lists = mapping_lists_aux;
  sa_batch->num_reads = array_list_size(fq_reads_aux);

  array_list_free(target_cals, (void *)NULL);
  
  cal_mng_simple_free(cal_mng_p);
  cal_mng_simple_free(cal_mng_n);

  sw_multi_output_free(output);

  if (pair_mode != SINGLE_END_MODE) { 
    array_list_t *mapping_list;
    for (int r = 0; r < sa_batch->num_reads; r++) {
      mapping_list = sa_batch->mapping_lists[r];
      size_t n_alig = array_list_size(mapping_list);
      for (int a = 0; a < n_alig; a++) {
	alignment_t *alig = array_list_get(a, mapping_list);
	alig->cigar = new_cigar_code_string(alig->alig_data);
      }
    }

    sa_complete_pairs(wf_batch);

    for (int r = 0; r < sa_batch->num_reads; r++) {
      mapping_list = sa_batch->mapping_lists[r];
      
      //Insert Data in metaexon and free cigars
      size_t n_alig = array_list_size(mapping_list);
      for (int a = 0; a < n_alig; a++) {
	alignment_t *alig = array_list_get(a, mapping_list);
	read = array_list_get(r, sa_batch->fq_reads);
	//printf("%s: \n", read->id);
	if (alig->alig_data == NULL) { printf("\tAlignment data NULL\n"); exit(-1); }

	cigar_code_find_and_report_sj(alig->position, alig->alig_data, alig->chromosome, 
				      alig->seq_strand, avls_list, metaexons, genome, read, 0);
	
	cigar_code_t *c = alig->alig_data;
	array_list_clear(c->ops, (void *)cigar_op_free);
	cigar_code_free(c);
	alig->alig_data = NULL;

      }
    }        

  } else {
    array_list_t *mapping_list;
    for (int r = 0; r < sa_batch->num_reads; r++) {
      mapping_list = sa_batch->mapping_lists[r];
      filter_alignments(report_optarg->all, 
			report_optarg->n_best, 
			report_optarg->n_hits,
			report_optarg->best,
			mapping_list);
      
      //Insert Data in metaexon and free cigars
      size_t n_alig = array_list_size(mapping_list);
      for (int a = 0; a < n_alig; a++) {
	alignment_t *alig = array_list_get(a, mapping_list);
	read = array_list_get(r, sa_batch->fq_reads);
	//printf("%s: \n", read->id);
	if (alig->alig_data == NULL) { printf("\tAlignment data NULL\n"); exit(-1); }

	alig->cigar = cigar_code_find_and_report_sj(alig->position, alig->alig_data, alig->chromosome, 
						    alig->seq_strand, avls_list, metaexons, genome, read, 1);
	
	cigar_code_t *c = alig->alig_data;
	array_list_clear(c->ops, (void *)cigar_op_free);
	cigar_code_free(c);
	alig->alig_data = NULL;

      }
    }        
  } 
  
  
  convert_batch_to_str(wf_batch);
  
  return -1;
  
}


int sa_rna_mapper_last(void *data) {
  array_list_t *sa_list;
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;
  sa_batch_t *sa_batch = wf_batch->mapping_batch;
  size_t num_reads = sa_batch->num_reads;
  sa_index3_t *sa_index = (sa_index3_t *) wf_batch->sa_index;
  sa_rna_input_t *sa_rna = wf_batch->data_input;
  genome_t *genome = (genome_t *)sa_rna->genome;
  avls_list_t *avls_list = (avls_list_t *)sa_rna->avls_list;
  metaexons_t *metaexons = (metaexons_t *)sa_rna->metaexons;
  sw_optarg_t *sw_optarg = sa_rna->sw_optarg;
  FILE *file1 = sa_rna->file1;
  FILE *file2 = sa_rna->file2;
  pair_server_input_t *pair_input = sa_rna->pair_input;
  int pair_mode = pair_input->pair_mng->pair_mode;
  report_optarg_t *report_optarg = pair_input->report_optarg;

  min_score = sa_rna->min_score;
  float gap_open = sw_optarg->gap_open;
  float gap_extend = sw_optarg->gap_extend;
  float match = sw_optarg->subst_matrix['A']['A'];
  float mismatch = sw_optarg->subst_matrix['A']['C'];

  int min_cal_size = sa_rna->cal_optarg->min_cal_size;
  char *query;
  //size_t num_reads = sa_batch->num_reads;
  //array_list_t *sa_list;
  //int min_score = 80;
  array_list_t *alignments_list;
  sw_multi_output_t *output = sw_multi_output_new(MAX_DEPTH);
  sa_sw_depth_t sw_depth;
  extern pthread_mutex_t mutex_sp;
  extern size_t total_sw;

  float cals_score[2048];
  array_list_t *target_cals = array_list_new(100, 1.25f,
					     COLLECTION_MODE_ASYNCHRONIZED);

  sw_depth.depth = 0;

  //metaexons_print(metaexons);
  //printf("||======= SECOND WORKFLOW =========||\n");

  char *read_quality_rev = NULL;
  
  for (int r = 0; r < num_reads; r++) {
    fastq_read_t *read = array_list_get(r, sa_batch->fq_reads);
    sa_list = sa_batch->mapping_lists[r];
    
    if (array_list_get_flag(sa_list) != SA_PARTIAL_TYPE) { continue; }

    fastq_read_revcomp(read);
    parse_alignment_data(sa_list, read);
    
    //if (!array_list_size(sa_list)) { continue; }
    //if (strcmp("seq.5a", read->id) == 0) {       
    //exit(-1);
    //}

    size_t num_items = array_list_size(sa_list);
    
    //extern pthread_mutex_t mutex_sp;
    //pthread_mutex_lock(&mutex_sp);
    //printf("@@@@@@ READ SECOND (%i): %s\n", num_items, read->id);
    //pthread_mutex_unlock(&mutex_sp);
    //printf("seq: %s\n", read->id);

    for (int i = 0; i < num_items; i++) {
      sa_alignment_t *sa_alignment = array_list_get(i, sa_list);
      cal_t *cal = array_list_get(0, sa_alignment->cals_list);
      seed_region_t *region = linked_list_get_first(cal->sr_list);
      
      assert(region);
      
      //if (!strcmp("seq.10a", read->id)) {
      //cal_print(cal);
      //printf("\tSEED: (%i)%i:[%lu|%i - %i|%lu] = %s\n", cal->strand, cal->chromosome_id, region->genome_start, region->read_start, region->read_end, region->genome_end, new_cigar_code_string(sa_alignment->c_final));
      //exit(-1);
      //}
      
      
      if (cal->strand) {
	query = read->revcomp;
      } else {
	query = read->sequence;
      }
      
      cigar_code_t *c_left, *c_right;

      //Close left gap...
      if (region->read_start != 0) {
	sa_alignment->left_close = 1;
	sa_alignment->left_dsp_w2 = 1;
	
	c_left = fill_extrem_gap(query, 
				 cal,
				 FILL_GAP_LEFT,
				 genome,
				 metaexons, 
				 avls_list);

	//if (c_left) {

	//}

	seed_region_t *seed_region_l = seed_region_new(0, region->read_start - 1, 
						       region->genome_start - 1 - region->read_start,
						       region->genome_start - 1, 0, 0, 0);
	linked_list_insert_first(seed_region_l, cal->sr_list);
	seed_region_l->info = c_left;

	if (!c_left && region->read_start <= 30) {
	  pthread_mutex_lock(&mutex_sp);
	  total_sw++;
	  pthread_mutex_unlock(&mutex_sp);

	  char reference_sw[2048];
	  size_t genome_start = region->genome_start - region->read_start;
	  size_t genome_end   = region->genome_start - 1;
	  
	  genome_read_sequence_by_chr_index(reference_sw, 0, cal->chromosome_id - 1,
					    &genome_start, &genome_end, genome);
	  
	  char query_sw[2048];
	  memcpy(query_sw, &query[0], region->read_start);
	  query_sw[region->read_start] = '\0';
	  
	  //printf("(%i)l-Query: %s\n", strlen(query_sw), query_sw);
	  //printf("(%i)l-Refer: %s\n", strlen(reference_sw), reference_sw);
	  
	  sw_insert_item(query_sw, reference_sw, FIRST_SW, seed_region_l,
	  		 sw_optarg, output, &sw_depth);	
	}
	//printf("seed_region_l = %s\n", new_cigar_code_string(seed_region_l->info));
	//if (seed_region_l->info == NULL) {
	//seed_region_l->info = cigar_code_new();
	//cigar_code_append_new_op(100, "S", seed_region_l->info);
      }
      
      //Close right gap...
      if (region->read_end != read->length - 1) {	  
	sa_alignment->right_close = 1;
	sa_alignment->right_dsp_w2 = 1;

	int gap_len = read->length - region->read_end - 1;
	c_right = fill_extrem_gap(query, 
				  cal,
				  FILL_GAP_RIGHT,
				  genome,
				  metaexons, 
				  avls_list);

	//if (c_right) {
	//sa_alignment->right_close = 1;
	//}

	seed_region_t *seed_region_r = seed_region_new(region->read_end + 1, 
						       read->length - 1,
						       //region->read_end + 1 + gap_len, 
						       region->genome_end + 1,
						       region->genome_end + 1 + gap_len, 0, 0, 0);
	seed_region_r->info = c_right;
	linked_list_insert_last(seed_region_r, cal->sr_list);

	if (!c_right && gap_len <= 30) {
	  pthread_mutex_lock(&mutex_sp);
	  total_sw++;
	  pthread_mutex_unlock(&mutex_sp);

	  char reference_sw[2048];
	  size_t genome_start = region->genome_end + 1;
	  size_t genome_end   = genome_start + gap_len;
	  
	  genome_read_sequence_by_chr_index(reference_sw, 0, cal->chromosome_id - 1,
	  				    &genome_start, &genome_end, genome);
	  
	  char query_sw[2048];
	  memcpy(query_sw, &query[region->read_end + 1], gap_len);
	  query_sw[gap_len] = '\0';
	  	  
	  //printf("r-Query: %s\n", query_sw, strlen(query_sw));
	  //printf("r-Refer: %s\n", reference_sw, strlen(reference_sw));
	  
	  sw_insert_item(query_sw, reference_sw, LAST_SW, seed_region_r,
			 sw_optarg, output, &sw_depth);

	}
      }

    }
    //array_list_clear(sa_list, (void *)NULL);
  }

  
  sw_insert_item(NULL, NULL, 0, NULL, 
  		 sw_optarg, output, &sw_depth);
  
  array_list_t *alignments_list_aux;
  for (int r = 0; r < num_reads; r++) {
    fastq_read_t *read = array_list_get(r, sa_batch->fq_reads);
    sa_list = sa_batch->mapping_lists[r];

    if (array_list_get_flag(sa_list) != SA_PARTIAL_TYPE) { continue; }
    
    size_t num_items = array_list_size(sa_list);
    int list_size = num_items == 0 ? 100 : num_items;

    //printf("@@@@@@ READ SECOND (%i): %s\n", num_items, read->id);

    alignments_list_aux = array_list_new(list_size, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    
    //if (num_items) { 
    //printf("Process Read %s (%i items)\n", read->id, num_items);
    //alignments_list = sa_batch->mapping_lists[r];
    
    for (int i = 0; i < num_items; i++) {
      sa_alignment_t *sa_alignment = array_list_get(i, sa_list);
      cal_t *cal = array_list_get(0, sa_alignment->cals_list);
      seed_region_t *region = linked_list_get_first(cal->sr_list);      
      assert(region);
      
      //printf("SEED: %i:[%lu|%i - %i|%lu]\n", cal->chromosome_id, region->genome_start, region->read_start, region->read_end, region->genome_end);
      
      if (sa_alignment->left_dsp_w2 && region->info != NULL) {
	//linked_list_item_t *item = cal->sr_list->first;
	//if (item) {
	//seed_region_t *sr = item->item;
	cigar_code_t *cc_aux = region->info;
	//if (cc_aux) {
	int dsp = 0;	    
	for (int tt = 0; tt < cc_aux->ops->size; tt++) {
	  cigar_op_t *op = array_list_get(tt, cc_aux->ops);
	  if (op->name != 'I') { 
	    dsp += op->number;
	  }
	}
	  
	//cal = array_list_get(0, merge_cals);	      
	//seed_region_t *seed_aux = linked_list_get_first(cal->sr_list);
	region->genome_start -= dsp;//seed_aux->read_start;
	cal->start -= dsp;//seed_aux->read_start;
	region->read_start = 0;	
	//printf("dsp = %i, cal->start = %i \n", dsp, cal->start);
	//}
	//}
      }
      
      cigar_code_t *c_final = cigar_code_new();
      cigar_code_t *c_aux = NULL;
      linked_list_item_t *item = cal->sr_list->first;
      while (item != NULL) {
	seed_region_t *seed_aux = item->item;
	c_aux = seed_aux->info;
	if (c_aux) {
	  for (int c = 0; c < c_aux->ops->size; c++) {
	    cigar_op_t *op = array_list_get(c, c_aux->ops);
	    cigar_code_append_op(op, c_final);
	  }
	  c_final->distance += c_aux->distance;
	  cigar_code_free(c_aux);
	}
	seed_aux->info = NULL;
	item = item->next;
      }

      sa_alignment->c_final = c_final;

      int read_score = cigar_code_score(c_final, read->length);  

      if (cigar_code_validate_(read, sa_alignment->c_final)) {
	//Cigar generator, generate score
	//printf("Data: \n");
	//printf("\trl = %i - num_ops = %i\n", read->length, sa_alignment->c_final->ops->size );
	//char cigar_str[2048] = "\0";
	//for (int k = 0; k < sa_alignment->c_final->ops->size; k++) {
	//cigar_op_t *op = array_list_get(k, sa_alignment->c_final->ops);	      

	//sprintf(cigar_str, "%s%i%c", cigar_str, op->number, op->name);
	//cigar_op_free(op);
	//}
	  
	//printf("Cigar: %s\n", cigar_str);

	//Report Aligment
	
	//if (cal->chromosome_id == 0) {
	//printf("***************** FINAL CIGAR %s(%i) (%i:%lu)*****************\n ", cigar_str, sa_alignment->c_final->distance, cal->chromosome_id, cal->start); 
	//exit(-1);
	//}
	
	cigar_op_t *op = array_list_get(0, sa_alignment->c_final->ops);
	if (op && op->name == 'H') { op->name = 'S'; }

	op = array_list_get(array_list_size(sa_alignment->c_final->ops) - 1, sa_alignment->c_final->ops);
	if (op && op->name == 'H') { op->name = 'S'; }

	//char *cigar_str = cigar_code_find_and_report_sj(cal->start - 1, sa_alignment->c_final, cal->chromosome_id, 
	//						cal->strand, avls_list, metaexons, genome, read);	

	char *quality_tmp;
	if (cal->strand) {
	  if (!read_quality_rev) {
	    read_quality_rev = str_reverse(read->quality);
	  }
	  quality_tmp = read_quality_rev;
	} else {
	  quality_tmp = read->quality;
	}

	//int score = SCORE(read->length, sa_alignment->c_final->distance);	
	int score = get_score (sa_alignment->c_final, read, match, mismatch, gap_open, gap_extend);

	if (score > min_score) {
	  alignment_t *alignment = alignment_new();
	  alignment_init_single_end(strdup(read->id), 
				    strdup(read->sequence), 
				    strdup(quality_tmp),
				    cal->strand, 
				    cal->chromosome_id,
				    cal->start - 1,
				    NULL,
				    sa_alignment->c_final->ops->size, 
				    score,
				    1, 0,
				    0, 0, alignment);
	

	  alignment->mapq = sa_alignment->c_final->distance;
	  alignment->alig_data = cigar_code_dup(sa_alignment->c_final);
	  //Report SJ

	  array_list_insert(alignment, alignments_list_aux);
	}	  
      }

      //-------------- Free sa alignments -------------------------
      array_list_t *merge_cals = sa_alignment->cals_list;
      int num_cals = array_list_size(merge_cals);      
      cal_t *cal_aux = array_list_get(0, merge_cals);
      
      for (int k = 0; k < num_cals; k++) {
	cal = array_list_get(k, merge_cals);
	linked_list_item_t *item = cal->sr_list->first;
	seed_region_t *sr;    
	while (item != NULL) {
	  sr = item->item;
	  cigar_code_t *cc = sr->info;
	  if (cc) {
	    cigar_code_free(cc);	
	  }
	  item = item->next;
	}
      }
      
      array_list_free(merge_cals, (void *)cal_free);
            
      if (sa_alignment->c_final) {
	array_list_clear(sa_alignment->c_final->ops, (void *)cigar_op_free);
	cigar_code_free(sa_alignment->c_final);	
      }
      
      sa_alignment_free(sa_alignment);
      
      //--------------------------------------------------------------------      
      
    }//End items


    size_t num_suffixes_p, num_suffixes_n;
    size_t suffix_len_p, suffix_len_n;
    size_t low_p, high_p, low_n, high_n;
    int chrom;
    size_t g_start;
    char *seq = read->sequence;
    char *seq_revcomp = read->revcomp;
    cigar_code_t *cc_aux;
    
    if (!array_list_size(alignments_list_aux)) {
      //continue;
      //Search strand(+)
      num_suffixes_p = search_suffix(&seq[0], sa_index->k_value, 
				     MAX_NUM_SUFFIXES, sa_index, 
				     &low_p, &high_p, &suffix_len_p);
      
      //Search strand(-)
      num_suffixes_n = search_suffix(&seq_revcomp[0], sa_index->k_value, 
				     MAX_NUM_SUFFIXES, sa_index, 
				     &low_n, &high_n, &suffix_len_n);

      if (suffix_len_p && num_suffixes_p) {
	//Report Exact Maps! (+)
	for (size_t suff = low_p; suff <= high_p; suff++) {
	  chrom = sa_index->CHROM[suff];	
	  g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];
	  
	  cal_t *cal_tmp = cal_simple_new(chrom + 1,
					  0, g_start, g_start + suffix_len_p);
	  
	  seed_region_t *seed_region_aux = seed_region_new(0, suffix_len_p - 1, 
							   g_start, g_start + suffix_len_p, 
							   0, 0, 0);
	  
	  linked_list_insert(seed_region_aux, cal_tmp->sr_list);
	  
	  cc_aux = fill_extrem_gap(seq, 
				   cal_tmp,
				   FILL_GAP_RIGHT,
				   genome,
				   metaexons, 
				   avls_list);
	  
	  if (cc_aux != NULL) {
	    cigar_op_t *c_op_aux = array_list_get(0, cc_aux->ops);
	    if (c_op_aux->name == 'M') { c_op_aux->number += suffix_len_p; }
	    else { array_list_insert(cigar_op_new(suffix_len_p, 'M'), cc_aux->ops); }
	    
	    //printf("=========::======= %s ==========::=======\n", new_cigar_code_string(cc_aux));
	    if (cigar_code_validate_(read, cc_aux)) {	    
	      //char cigar_str[2048] = "\0";
	      //for (int k = 0; k < cc_aux->ops->size; k++) {
	      //cigar_op_t *op = array_list_get(k, cc_aux->ops);	      
	      //if (op->name == 'H') { op->name = 'S'; }
	      //sprintf(cigar_str, "%s%i%c", cigar_str, op->number, op->name);
	      //cigar_op_free(op);
	      //}

	      cigar_op_t *op = array_list_get(0, cc_aux->ops);
	      if (op && op->name == 'H') { op->name = 'S'; }
	      
	      op = array_list_get(array_list_size(cc_aux->ops) - 1, cc_aux->ops);
	      if (op && op->name == 'H') { op->name = 'S'; }

	      //char *cigar_str = cigar_code_find_and_report_sj(cal_tmp->start, cc_aux, cal_tmp->chromosome_id, 
	      //					      cal_tmp->strand, avls_list, metaexons, genome, read);	

	      char *quality_tmp;
	      if (cal_tmp->strand) {
		if (!read_quality_rev) {
		  read_quality_rev = str_reverse(read->quality);
		}
		quality_tmp = read_quality_rev;
	      } else {
		quality_tmp = read->quality;
	      }

	      //int score = SCORE(read->length, cc_aux->distance);	
	      int score = get_score (cc_aux, read, match, mismatch, gap_open, gap_extend);

	      if (score > min_score) {
		alignment_t *alignment = alignment_new();
		alignment_init_single_end(strdup(read->id), 
					  strdup(read->sequence), 
					  strdup(quality_tmp),
					  cal_tmp->strand, 
					  cal_tmp->chromosome_id,
					  cal_tmp->start,
					  //cigar_str,
					  NULL,
					  cc_aux->ops->size, 
					  score,
					  1, 0,
					  0, 0, alignment);
	      
		alignment->mapq = cc_aux->distance;
		alignment->alig_data = cigar_code_dup(cc_aux);

		//Report SJ	      
		array_list_insert(alignment, alignments_list_aux); 
	      }
	    }
	    array_list_clear(cc_aux->ops, (void *)cigar_op_free);
	    cigar_code_free(cc_aux);
	  }
	  cal_simple_free(cal_tmp);
	}
      }
      
      //printf("RESULTS (-):\n");
      if (suffix_len_n && num_suffixes_n) {
	for (size_t suff = low_n; suff <= high_n; suff++) {
	  chrom = sa_index->CHROM[suff];
	  g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];

	  cal_t *cal_tmp = cal_simple_new(chrom + 1,
					  1, g_start, g_start + suffix_len_n);

	  seed_region_t *seed_region_aux = seed_region_new(0, suffix_len_n - 1, 
							   g_start, g_start + suffix_len_n, 
							   0, 0, 0);

	  linked_list_insert(seed_region_aux, cal_tmp->sr_list);

	  cc_aux = fill_extrem_gap(seq_revcomp, 
				   cal_tmp,
				   FILL_GAP_RIGHT,
				   genome,
				   metaexons, 
				   avls_list);

	  if (cc_aux != NULL) {
	    cigar_op_t *c_op_aux = array_list_get(0, cc_aux->ops);
	    if (c_op_aux->name == 'M') { c_op_aux->number += suffix_len_n; }
	    else { array_list_insert(cigar_op_new(suffix_len_n, 'M'), cc_aux->ops); }

	    if (cigar_code_validate_(read, cc_aux)) {	
	      //char cigar_str[2048] = "\0";
	      //for (int k = 0; k < cc_aux->ops->size; k++) {
	      //cigar_op_t *op = array_list_get(k, cc_aux->ops);	      
	      //if (op->name == 'H') { op->name = 'S'; }
	      //sprintf(cigar_str, "%s%i%c", cigar_str, op->number, op->name);
	      //cigar_op_free(op);
	      //}
    
	      //if (cal_tmp->chromosome_id == 0) {
	      //printf("***************** FINAL CIGAR %s (%i)*****************\n ", cigar_str); 
	      //exit(-1);
	      //}

	      cigar_op_t *op = array_list_get(0, cc_aux->ops);
	      if (op && op->name == 'H') { op->name = 'S'; }
	      
	      op = array_list_get(array_list_size(cc_aux->ops) - 1, cc_aux->ops);
	      if (op && op->name == 'H') { op->name = 'S'; }
	      
	      //char *cigar_str = cigar_code_find_and_report_sj(cal_tmp->start, cc_aux, cal_tmp->chromosome_id, 
	      //					      cal_tmp->strand, avls_list, metaexons, genome, read);	

	      char *quality_tmp;
	      if (cal_tmp->strand) {
		if (!read_quality_rev) {
		  read_quality_rev = str_reverse(read->quality);
		}
		quality_tmp = read_quality_rev;
	      } else {
		quality_tmp = read->quality;
	      }

	      //int score = SCORE(read->length, cc_aux->distance);	
	      int score = get_score (cc_aux, read, match, mismatch, gap_open, gap_extend);

	      if (score > min_score) {
		alignment_t *alignment = alignment_new();
		alignment_init_single_end(strdup(read->id), 
					  strdup(read->sequence), 
					  strdup(quality_tmp),
					  cal_tmp->strand, 
					  cal_tmp->chromosome_id,
					  cal_tmp->start,
					  //cigar_str,
					  NULL,
					  cc_aux->ops->size, 
					  score, 
					  1, 0,
					  0, 0, alignment);
	      
		alignment->mapq = cc_aux->distance;
		alignment->alig_data = cigar_code_dup(cc_aux);
		//Report SJ
	      
		array_list_insert(alignment, alignments_list_aux);
	      }
	    }
	    array_list_clear(cc_aux->ops, (void *)cigar_op_free);
	    cigar_code_free(cc_aux);
	  }
	  cal_simple_free(cal_tmp);
	}
      }
    }

    
    if (!array_list_size(alignments_list_aux)) {
      //===== Seeding Strategy =====//
      //extern pthread_mutex_t mutex_sp;
      //pthread_mutex_lock(&mutex_sp);
      //fprintf(stderr, "XXX@@@XXX-3: %s\n", read->id);
      array_list_t *sa_alignments_list = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

      cal_mng_t * cal_mng_p = cal_rna_mng_new(sa_index->genome); //(+)
      cal_mng_t * cal_mng_n = cal_rna_mng_new(sa_index->genome); //(-)
      cal_mng_t * cal_mng, *p;

      linked_list_t *cal_list;
      int num_cals, n_report;

      array_list_t *merge_cals;
      int cal_pos;

      linked_list_item_t *item;
      seed_region_t *s, *s_prev;

      cal_t *cal_prev, *cal_next, *cal;

      int seed_size = 18;
      int seed_inc  = seed_size / 2;
      int read_pos;
      size_t num_suffixes, low, high, suffix_len, len_seq, id_seed;

      int max_read_pos = read->length - seed_size - 1;
      
      for (int s = 0; s < 2; s++) { //Strand
	read_pos  = 0;//seed_size;
	if (!s) {
	  query = seq;
	  cal_mng = cal_mng_p;
	} else {
	  query = seq_revcomp;
	  cal_mng = cal_mng_n;
	}
	
	while (read_pos + sa_index->k_value < read->length) {
	  num_suffixes = search_suffix(&query[read_pos], sa_index->k_value, 
				       MAX_NUM_SUFFIXES, sa_index, 
				       &low, &high, &suffix_len);
	  
	  //printf("(%c)============= Seed (%i + %i), %i ===========\n", s == 0 ? '+' : '-', read_pos, read_pos + sa_index->k_value, suffix_len);
	  
	  if (suffix_len && num_suffixes) {
	    //Storage Mappings
	    for (size_t suff = low; suff <= high; suff++) {	
	      chrom = sa_index->CHROM[suff];
	      g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom] + 1;
	      //printf("\tSTORE SEED %i:[%lu|%i-%i|%lu]\n", chrom, g_start, read_pos, read_pos + suffix_len - 1, g_start + suffix_len - 1);
	      generate_cals(chrom + 1, s, 
			    g_start, g_start + suffix_len - 1, 
			    read_pos, read_pos + suffix_len - 1,
			    cal_mng->cals_lists[chrom], read->length, id_seed++);
	    }
	    ////////////////////////////////////////
	    //|If seed map, extend suffix_len + 1|//
	    //|__________________________________|//
	    //|---(18)---|--(10)--|              |//
	    //|        (28)       |---(18)---|   |//
	    //|__________________________________|//
	    ////////////////////////////////////////
	    read_pos += suffix_len + 1;
	  } else {
	    ////////////////////////////////////
	    //|  Else, extend seed_size / 2  |//
	    //|______________________________|//
	    //|---(18)---|                   |//  
	    //| (9) |---(18)---|             |//
	    //|     | (9) |---(18)---|       |//
	    //|______________________________|//
	    ////////////////////////////////////
	    read_pos += seed_inc; 
	  }
	} //loop seeds	

	
	if (read->length - 1 != read_pos) {
	  read_pos = read->length - sa_index->k_value - 1;
	  num_suffixes = search_suffix(&query[read_pos], sa_index->k_value, 
				       MAX_NUM_SUFFIXES, sa_index, 
				       &low, &high, &suffix_len);	      

	  //printf("L.(%c)============= Seed (%i + %i), %i ===========\n", s == 0 ? '+' : '-', read_pos, read_pos + sa_index->k_value, suffix_len);

	  if (suffix_len && num_suffixes) {
	    //Storage Mappings
	    for (size_t suff = low; suff <= high; suff++) {
	      //printf("\tL.STORE SEED %i:[%lu|%i-%i|%lu]\n", chrom, g_start, read_pos, read_pos + suffix_len - 1, g_start + suffix_len - 1);

	      chrom = sa_index->CHROM[suff];
	      g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom] + 1;
	      generate_cals(chrom + 1, s, 
			    g_start, g_start + suffix_len - 1, 
			    read_pos, read_pos + suffix_len - 1,
			    cal_mng->cals_lists[chrom], read->length, id_seed++);
	    }
	  }
	}
      } //loop strands      

      //printf("==================== CALs result Seeding =====================\n");
      //printf("==(+)==\n");
      //cal_mng_print(cal_mng_p);
      //printf("==(-)==\n");
      //cal_mng_print(cal_mng_n);
      //printf("==============================================================\n");
      
      array_list_clear(target_cals, (void *)NULL);
      //Merge CALs and select new targets
      for (int st = 0; st < 2; st++) {
	if (!st) { p = cal_mng_p; }
	else  { p = cal_mng_n; }
	
	for (unsigned int i = 0; i < p->num_chroms; i++) {
	  cal_list = p->cals_lists[i];
	  num_cals = linked_list_size(cal_list);
	  //printf("(%i)%i:Num CAls %i\n", st, i, num_cals);
	  if (num_cals) {
	    //=====*---------------------*=====//
	    merge_cals = array_list_new(100,
					1.25f,
					COLLECTION_MODE_ASYNCHRONIZED);
	    cal_pos = 1;
	    item = cal_list->first;
	    cal_prev = item->item;

	    array_list_insert(cal_prev, merge_cals); 
	    while (cal_pos < num_cals) {

	      item = item->next;
	      cal_next = item->item;
		  
	      s_prev = linked_list_get_last(cal_prev->sr_list);
	      s = linked_list_get_first(cal_next->sr_list);
		  
	      assert(s_prev != NULL);
	      assert(s != NULL);
	      
	      //cal_print(cal_prev);
	      //cal_print(cal_next);
	      //printf("s_prev->id = %i != s->id = %i\n", s_prev->id, s->id);
	      
	      int ok = 0;
	      if (s->read_start > s_prev->read_end) {
		ok = 1;
	      } else if (s_prev->read_end - s->read_start <= 5) {
		ok = 1;
	      }
	      
	      if (ok &&
		  cal_prev->chromosome_id == cal_next->chromosome_id && 
		  cal_prev->strand == cal_next->strand && 
		  s_prev->id != s->id &&
		  cal_next->start <= (cal_prev->end + max_intron)) {
		//printf("Merge!! cal_prev->end = %lu, cal_next->start = %lu\n", cal_prev->end, cal_next->start);
		array_list_insert(cal_next, merge_cals);
	      } else {
		array_list_insert(merge_cals, target_cals);
		merge_cals = array_list_new(10,
					    1.25f,
					    COLLECTION_MODE_ASYNCHRONIZED);
		array_list_insert(cal_next, merge_cals);
	      }                                                                                  
	      cal_prev = cal_next;
	      cal_pos++;
	    }
	    array_list_insert(merge_cals, target_cals);	    
	    //=====*---------------------*=====//
	  }
	}
      }
      
      
      int num_targets = array_list_size(target_cals);    
      if (!num_targets) { 
	goto free_w2;
      }

      const int limit_cals = 50;

      if (num_targets > limit_cals) {
	//Free other CALs      
	for (int i = limit_cals; i < num_targets; i++) {
	  merge_cals = array_list_get(i, target_cals);
	  array_list_free(merge_cals, (void *)NULL);
	}
	num_targets = limit_cals;
      }

      //=====*---------------------*=====//
      for (int i = 0; i < num_targets; i++) {
	merge_cals = array_list_get(i, target_cals);
	cals_score[i] = generate_cals_merge_score(merge_cals, read->length);
      }
      //=====*---------------------*=====//
    	  
      //=====*---------------------*=====//
      for (int i = 0; i < num_targets - 1; i++) {
	for (int j = i + 1; j < num_targets; j++) {
	  if (cals_score[j] > cals_score[i]) {
	    array_list_swap(i, j, target_cals);
	    float score = cals_score[j];
	    cals_score[j] = cals_score[i];
	    cals_score[i] = score;
	  }
	}
      }
      //=====*---------------------*=====//

      //Report the first best 10 CALs
      n_report = num_targets < 5 ? num_targets : 5;
      
      int report_tmp = 0;
      //int report_tmp = 1;
      for (int i = 0; i < n_report; i++) {
	//printf("CAL %i: %f\n", i, cals_score[i]);
	if (cals_score[i] < min_cal_size) { 
	  break;
	}
	report_tmp++;
      }
    
      n_report = report_tmp;
      
      /*
      printf("Read : %s\n", read->id);
      printf("========================= (%i) =============================\n", n_report);
      for (int i = 0; i < n_report; i++) {
	merge_cals = array_list_get(i, target_cals);
	printf("Merge \%i (CALs %i):\n", i, array_list_size(merge_cals));
	for (int j = 0; j < array_list_size(merge_cals); j++) {
	  cal = array_list_get(j, merge_cals);
	  cal_print(cal);
	}
      }
      printf("======================================================\n");
      */

      //Free other CALs
      for (int i = n_report; i < num_targets; i++) {
	merge_cals = array_list_get(i, target_cals);
	array_list_free(merge_cals, (void *)NULL);
      }     
      
      if (!n_report) { 	goto free_w2; }

      //printf("*************** FINAL CALS ******************\n");
      seed_region_t *seed_first;
      seed_region_t *seed_last;
      linked_list_item_t *list_item, *list_item_prev;
      seed_region_t *seed_prev, *seed_next;
      size_t sp_start, sp_end;
      int distance, sp_type;

      for (int i = 0; i < n_report; i++) {
	merge_cals = array_list_get(i, target_cals);

	int ok = 1;
	for (int c = 0; c < array_list_size(merge_cals); c++) {
	  cal_t *cal = array_list_get(c, merge_cals);
	  if (!validate_cal(cal)) {
	      ok = 0;
	  }
	}

	if (!ok ) { 
	  array_list_free(merge_cals, (void *)NULL);
	  continue;
	}

	sa_alignment_t *sa_alignment = sa_alignment_new(merge_cals);
	array_list_insert(sa_alignment, sa_alignments_list);
	
	//printf("SCORE (%f):\n", cals_score[i]);
	cal_prev = array_list_get(0, merge_cals);
	for (int j = 1; j < array_list_size(merge_cals); j++) {
	  cal_next = array_list_get(j, merge_cals);

	  if (cal_prev->strand) {
	    query = read->revcomp;
	  } else {
	    query = read->sequence;
	  }
	  
	  //Fill gaps prev CAL
	  cal_fill_gaps_2(cal_prev,
			  query,
			  genome,
			  metaexons,
			  avls_list, 
			  &sw_depth, 
			  sw_optarg, 
			  output);
	  
	  //Search big SJ
	  seed_first = linked_list_get_last(cal_prev->sr_list);
	  seed_last  = linked_list_get_first(cal_next->sr_list);	  
	  
	  cigar_code_t *cc_sj = search_splice_junction(sw_optarg,
						       seed_first, 
						       seed_last,
						       cal_next->chromosome_id,
						       cal_next->strand,
						       query, genome, 
						       &sp_start, 
						       &sp_end,
						       &sp_type,
						       &distance);
	  if (cc_sj == NULL) {
	    cc_sj = cigar_code_new();
	    cigar_code_append_new_op(seed_last->read_start - seed_first->read_end - 1, 'M', cc_sj);
	  }

	  sa_alignment->cigar_middle[sa_alignment->num_sp++] = cc_sj;
	  //printf("[%lu-%lu] Cigar %s\n", sp_start, sp_end, new_cigar_code_string(cc_sj));
	  
	  cal_prev = cal_next;

	}
	
	//Fill gaps prev CAL
	cal_fill_gaps_2(cal_prev,
			query,
			genome,
			metaexons,
			avls_list, 
			&sw_depth, 
			sw_optarg, 
			output);

	cal_prev = array_list_get(0, merge_cals);
	seed_region_t *region = linked_list_get_first(cal_prev->sr_list);
	if (region->read_start != 0) { 
	  sa_alignment->left_dsp_w2 = 1;
	  /*
	  cigar_code_t *c_left = fill_extrem_gap(query, 
						 cal_prev,
						 FILL_GAP_LEFT,
						 genome,
						 metaexons, 
						 avls_list);
	  */
	  seed_region_t *seed_region_l = seed_region_new(0, region->read_start - 1, 
							 region->genome_start - 1 - region->read_start,
							 region->genome_start - 1, 0, 0, 0);
	  linked_list_insert_first(seed_region_l, cal_prev->sr_list);
	  //seed_region_l->info = c_left;
	  
	  //c_left = NULL;
	  pthread_mutex_lock(&mutex_sp);
	  total_sw++;
	  pthread_mutex_unlock(&mutex_sp);
	    
	  char reference_sw[2048];
	  size_t genome_start = region->genome_start - region->read_start;
	  size_t genome_end   = region->genome_start - 1;
	    
	  genome_read_sequence_by_chr_index(reference_sw, 0, cal_prev->chromosome_id - 1,
					    &genome_start, &genome_end, genome);
	    
	  char query_sw[2048];
	  memcpy(query_sw, &query[0], region->read_start);
	  query_sw[region->read_start] = '\0';
	    
	  //printf("(%i)l-Query: %s\n", strlen(query_sw), query_sw);
	  //printf("(%i)l-Refer: %s\n", strlen(reference_sw), reference_sw);
	    
	  sw_insert_item(query_sw, reference_sw, FIRST_SW, seed_region_l,
			 sw_optarg, output, &sw_depth);	
	  
	  
	  //printf("seed_region_l = %s\n", new_cigar_code_string(seed_region_l->info));
	  //if (seed_region_l->info == NULL) {
	  //seed_region_l->info = cigar_code_new();
	  //cigar_code_append_new_op(100, "S", seed_region_l->info);	  
	}
      
	//Close right gap...
	cal = array_list_get(array_list_size(merge_cals) - 1, merge_cals);
	region = linked_list_get_last(cal->sr_list); 
	if (region->read_end != read->length - 1) {
	  sa_alignment->right_dsp_w2 = 1;
	  int gap_len = read->length - region->read_end - 1;

	  /*
	  cigar_code_t *c_right = fill_extrem_gap(query,
						  cal,
						  FILL_GAP_RIGHT,
						  genome,
						  metaexons, 
	  					  avls_list);	  
	  */
	  //if (c_right) {
	  //sa_alignment->right_close = 1;
	  //}

	  //c_right = NULL;
	  seed_region_t *seed_region_r = seed_region_new(region->read_end + 1, 
							 read->length - 1,
							 //region->read_end + 1 + gap_len, 
							 region->genome_end + 1,
							 region->genome_end + 1 + gap_len, 0, 0, 0);

	  //seed_region_r->info = c_right;
	  linked_list_insert_last(seed_region_r, cal->sr_list);
	  

	  pthread_mutex_lock(&mutex_sp);
	  total_sw++;
	  pthread_mutex_unlock(&mutex_sp);
	    
	  char reference_sw[2048];
	  size_t genome_start = region->genome_end + 1;
	  size_t genome_end   = genome_start + gap_len;
	    
	  genome_read_sequence_by_chr_index(reference_sw, 0, cal->chromosome_id - 1,
					    &genome_start, &genome_end, genome);
	    
	  char query_sw[2048];
	  memcpy(query_sw, &query[region->read_end + 1], gap_len);
	  query_sw[gap_len] = '\0';
	    
	  //printf("r-Query: %s\n", query_sw, strlen(query_sw));
	  //printf("r-Refer: %s\n", reference_sw, strlen(reference_sw));
	  
	  sw_insert_item(query_sw, reference_sw, LAST_SW, seed_region_r,
			 sw_optarg, output, &sw_depth);
	  
	}
      }
      
      //printf("*************** FINAL CALS ******************\n");
      
      sw_insert_item(NULL, NULL, 0, NULL, 
		     sw_optarg, output, &sw_depth);
      
      
      //Merge cigars      
      for (int i = 0; i < array_list_size(sa_alignments_list); i++) {	
	sa_alignment_t *sa_alignment = array_list_get(i, sa_alignments_list);
	
	sa_alignment->c_final = cigar_code_new();
	merge_cals = sa_alignment->cals_list;
	cal_prev = array_list_get(0, merge_cals);
	cal_next = array_list_get(array_list_size(merge_cals) - 1, merge_cals);

	seed_first = cal_prev->sr_list->first->item;
	seed_last  = cal_next->sr_list->last->item;

	if ((sa_alignment->left_dsp_w2   && seed_first->info == NULL) || 
	    (sa_alignment->right_dsp_w2  && seed_last->info == NULL)) {
	  
	  for (int c = 0; c < array_list_size(merge_cals); c++) {
	    cal_t *cal = array_list_get(c, merge_cals);
	    linked_list_item_t *item = cal->sr_list->first;
	    seed_region_t *sr;
	    while (item != NULL) {
	      sr = item->item;
	      cigar_code_t *cc = sr->info;
	      if (cc) {
		array_list_clear(cc->ops, (void *)cigar_op_free);
		cigar_code_free(cc);
	      }
	      item = item->next;
	    }
	  }

	  array_list_free(merge_cals, (void *)NULL);
	  cigar_code_free(sa_alignment->c_final);

	  if (sa_alignment->num_sp > 0) {
	    for (int x = 0; x < sa_alignment->num_sp; x++) {
	      array_list_clear(((cigar_code_t *)sa_alignment->cigar_middle[x])->ops, (void *)cigar_op_free);
	      cigar_code_free(sa_alignment->cigar_middle[x]);
	    }
	  }
	  
	  sa_alignment_free(sa_alignment);
	  
	  continue; 
	}


	if (sa_alignment->left_dsp_w2) {
	  int dsp = 0;
	  cc_aux = seed_first->info;

	  if (cc_aux == NULL) {
	    //fprintf(stderr, "YES NULL\n");
	  } else {
	    //fprintf(stderr, "NOT NULL\n");
	    for (int tt = 0; tt < cc_aux->ops->size; tt++) {
	      cigar_op_t *op = array_list_get(tt, cc_aux->ops);
	      if (op->name != 'I') { 
		dsp += op->number;
	      }
	    }
	  
	    seed_first->genome_start -= dsp;
	    cal_prev->start -= dsp;
	    seed_first->read_start = 0;		  	  

	  }
	}

	int sp_pos = 0;
			
	for (int c = 0; c < array_list_size(merge_cals); c++) {
	  cal_t *cal = array_list_get(c, merge_cals);
	  linked_list_item_t *item = cal->sr_list->first;
	  seed_region_t *sr;
	  while (item != NULL) {
	    sr = item->item;
	    cigar_code_t *cc = sr->info;
	    //printf("%i:%lu-%lu, %i-%i\n", cal->chromosome_id, sr->genome_start, sr->genome_end, sr->read_start, sr->read_end);
	    if (cc) {	      
	      //printf("\t---------------------->[%i-%i] : %s\n", sr->read_start, sr->read_end, new_cigar_code_string(cc));
	      for (int tt = 0; tt < cc->ops->size; tt++) {
		cigar_op_t *op = array_list_get(tt, cc->ops);
		cigar_code_append_op(op, sa_alignment->c_final);
	      }
	      sa_alignment->c_final->distance += cc->distance;
	      cigar_code_free(cc);
	      //sr->info = NULL;
	    } else {
	      //printf("\t----------------------->NULL, %iM\n", sr->read_end - sr->read_start + 1);
	      cigar_code_append_new_op(sr->read_end - sr->read_start + 1, 'M', sa_alignment->c_final);
	    }
	    item = item->next;
	  }
	 
	  if (c != array_list_size(merge_cals) - 1) {
	    cigar_code_t *cc = sa_alignment->cigar_middle[sp_pos++];
	    if (cc) {
	      //printf("\t----------------------> SJ : %s\n", new_cigar_code_string(cc));
	      for (int tt = 0; tt < cc->ops->size; tt++) {
		cigar_op_t *op = array_list_get(tt, cc->ops);
		cigar_code_append_op(op, sa_alignment->c_final);
	      }
	      sa_alignment->c_final->distance += cc->distance;
	      cigar_code_free(cc);
	    }
	  }
	}
	
	//printf("¡¡¡¡¡¡¡¡¡ALIGNMENT FOUND: %lu: %s (%i)!!!!!!!!!!\n", cal_prev->start, new_cigar_code_string(sa_alignment->c_final), sa_alignment->c_final->distance);
	
	if (cigar_code_validate_(read, sa_alignment->c_final)) {	
	  //char cigar_str[2048] = "\0";
	  //for (int k = 0; k < sa_alignment->c_final->ops->size; k++) {
	  //cigar_op_t *op = array_list_get(k, sa_alignment->c_final->ops);	      
	  //if (op->name == 'H') { op->name = 'S'; }
	  //sprintf(cigar_str, "%s%i%c", cigar_str, op->number, op->name);
	  //cigar_op_free(op);
	  //}

	  cigar_op_t *op = array_list_get(0, sa_alignment->c_final->ops);
	  if (op && op->name == 'H') { op->name = 'S'; }
	  
	  op = array_list_get(array_list_size(sa_alignment->c_final->ops) - 1, sa_alignment->c_final->ops);
	  if (op && op->name == 'H') { op->name = 'S'; }

	  //char *cigar_str = cigar_code_find_and_report_sj(cal_prev->start - 1, sa_alignment->c_final, cal_prev->chromosome_id, 
	  //						  cal_prev->strand, avls_list, metaexons, genome, read);	

	  char *quality_tmp;
	  if (cal_prev->strand) {
	    if (!read_quality_rev) {
	      read_quality_rev = str_reverse(read->quality);
	    }
	    quality_tmp = read_quality_rev;
	  } else {
	    quality_tmp = read->quality;
	  }

	  //int score = SCORE(read->length, sa_alignment->c_final->distance);	
	  int score = get_score (sa_alignment->c_final, read, match, mismatch, gap_open, gap_extend);

	  if (score > min_score) {
	    alignment_t *alignment = alignment_new();
	    alignment_init_single_end(strdup(read->id), 
				      strdup(read->sequence), 
				      strdup(quality_tmp),
				      cal_prev->strand,
				      cal_prev->chromosome_id,
				      cal_prev->start - 1,
				      NULL,
				      sa_alignment->c_final->ops->size, 
				      score,
				      1, 0,
				      0, 0, alignment);
	  
	    alignment->mapq = sa_alignment->c_final->distance;
	    alignment->alig_data = cigar_code_dup(sa_alignment->c_final);

	    //Report SJ	  
	    array_list_insert(alignment, alignments_list_aux);
	  }
	}
	
	//-------------- Free sa alignments -------------------------
	array_list_t *merge_cals = sa_alignment->cals_list;
	//int num_cals = array_list_size(merge_cals);      
	//cal_t *cal_aux = array_list_get(0, merge_cals);	
	
	//if (sa_alignment->c_final) {
	//array_list_clear(sa_alignment->c_final->ops, (void *)cigar_op_free);
	//cigar_code_free(sa_alignment->c_final);	
	//}	

	array_list_free(merge_cals, (void *)NULL);
	
	array_list_clear(sa_alignment->c_final->ops, (void *)cigar_op_free);
	cigar_code_free(sa_alignment->c_final);	
	sa_alignment_free(sa_alignment);
	
      }
      
    free_w2:
      cal_mng_simple_clear(cal_mng_p);
      cal_mng_simple_clear(cal_mng_n);

      cal_mng_simple_free(cal_mng_p);
      cal_mng_simple_free(cal_mng_n);

      array_list_free(sa_alignments_list, (void *)NULL);

    }
    

    array_list_set_flag(array_list_get_flag(sa_batch->mapping_lists[r]), alignments_list_aux);
    array_list_free(sa_batch->mapping_lists[r], (void *)NULL);
    sa_batch->mapping_lists[r] = alignments_list_aux;

    free(read_quality_rev);
    read_quality_rev = NULL;
    
  } //End reads 

  array_list_free(target_cals, (void *)NULL);
  sw_multi_output_free(output);

  if (pair_mode != SINGLE_END_MODE) {
    array_list_t *mapping_list;
    for (int r = 0; r < sa_batch->num_reads; r++) {
      mapping_list = sa_batch->mapping_lists[r];
      if (array_list_get_flag(mapping_list) != SA_PARTIAL_TYPE) { continue; }

      size_t n_alig = array_list_size(mapping_list);
      for (int a = 0; a < n_alig; a++) {
	alignment_t *alig = array_list_get(a, mapping_list);
	alig->cigar = new_cigar_code_string(alig->alig_data);
      }
    }

    for (int r = 0; r < sa_batch->num_reads; r++) {
      mapping_list = sa_batch->mapping_lists[r];      
      size_t n_alig = array_list_size(mapping_list);
      for (int a = 0; a < n_alig; a++) {
	alignment_t *alig = array_list_get(a, mapping_list);
	if (alig->cigar == NULL) {
	  printf("%i - %i, %p\n", array_list_get_flag(mapping_list), SA_PARTIAL_TYPE, alig->alig_data);
	  exit(-1);
	}
      }
    }


    sa_complete_pairs(wf_batch);

    for (int r = 0; r < sa_batch->num_reads; r++) {
      mapping_list = sa_batch->mapping_lists[r];

      //Insert Data in metaexon and free cigars
      size_t n_alig = array_list_size(mapping_list);
      for (int a = 0; a < n_alig; a++) {
	alignment_t *alig = array_list_get(a, mapping_list);
	if (alig->alig_data == NULL) { 
	  alig->alig_data = cigar_code_new_by_string(alig->cigar);
	}

	fastq_read_t *read = array_list_get(r, sa_batch->fq_reads);
	cigar_code_find_and_report_sj(alig->position, alig->alig_data, alig->chromosome, 
				      alig->seq_strand, avls_list, metaexons, genome, read, 0);	

	cigar_code_t *c = alig->alig_data;
	array_list_clear(c->ops, (void *)cigar_op_free);
	cigar_code_free(c);
	alig->alig_data = NULL;
      }
    }

  } else {
    array_list_t *mapping_list;
    for (int r = 0; r < sa_batch->num_reads; r++) {
      mapping_list = sa_batch->mapping_lists[r];
      filter_alignments(report_optarg->all, 
			report_optarg->n_best, 
			report_optarg->n_hits,
			report_optarg->best,
			mapping_list);
      
      //Insert Data in metaexon and free cigars
      size_t n_alig = array_list_size(mapping_list);
      for (int a = 0; a < n_alig; a++) {
	alignment_t *alig = array_list_get(a, mapping_list);
	if (alig->alig_data == NULL) { printf("Alignment NULL\n"); exit(-1); }

	fastq_read_t *read = array_list_get(r, sa_batch->fq_reads);
	alig->cigar = cigar_code_find_and_report_sj(alig->position, alig->alig_data, alig->chromosome, 
						    alig->seq_strand, avls_list, metaexons, genome, read, 1);
	
	cigar_code_t *c = alig->alig_data;
	array_list_clear(c->ops, (void *)cigar_op_free);
	cigar_code_free(c);
	alig->alig_data = NULL;

      }      
    }
  }

  convert_batch_to_str(wf_batch);
  
  return -1;

}


//FOR SEARCH CANNONICAL SPLICE JUNCTION WITHOUT SMITH-WATERMAN ALGORITHM
cigar_code_t* search_splice_junction(sw_optarg_t *sw_optarg,
				     seed_region_t *s_prev, seed_region_t *s_next,
				     int chromosome_id, int strand, 
				     char *sequence, genome_t *genome, 
				     size_t *sp_start, size_t *sp_end,
				     int *sp_type,
				     int *distance) {

  assert(s_prev != NULL);
  assert(s_next != NULL);

  cigar_code_t *cc_final = NULL;
  int seq_len = strlen(sequence);
  int read_start   = s_prev->read_end;
  int read_end     = s_next->read_start;
  int intron_size = 0;
  
  size_t genome_start;
  size_t genome_end;
  
  const int FLANK = 20;
  const int SEQURITY_FLANK = 5;
  const int SEQURITY_FLANK_2 = 10;

  int gap_read = read_end - read_start - 1;

  if (gap_read == 0) {
    gap_read = -1;
  }
  
  //int read_end_aux = s_prev->read_end;
  //int read_start_aux = s_next->read_start;
  //size_t genome_start_aux = s_next->genome_start;
  //size_t genome_end_aux = s_prev->genome_end;

  if (gap_read < 0) {
    gap_read = abs(gap_read) + 5;
    s_prev->read_end     -= gap_read;
    s_next->read_start   += gap_read;

    s_prev->genome_end   -= gap_read;
    s_next->genome_start += gap_read;    
  } else {
    s_prev->read_end     -= SEQURITY_FLANK;
    s_next->read_start   += SEQURITY_FLANK;
    
    s_prev->genome_end   -= SEQURITY_FLANK;
    s_next->genome_start += SEQURITY_FLANK;
  }

  
  read_start = s_prev->read_end;
  read_end = s_next->read_start;
  
  gap_read = read_end - read_start - 1;
  //printf("%i - %i - 1 = %i\n", read_end, read_start, gap_read);
  
  char left_exon[2048];
  char right_exon[2048];
  
  genome_start = s_prev->genome_end + 1;
  genome_end   = s_prev->genome_end + gap_read + FLANK;
  assert(genome_end - genome_start < 2048);
  
  genome_read_sequence_by_chr_index(left_exon, 0, 
				    chromosome_id - 1, 
				    &genome_start, &genome_end, genome);


  genome_start = s_next->genome_start - gap_read - FLANK;
  genome_end   = s_next->genome_start - 1;
  assert(genome_end - genome_start < 2048);

  genome_read_sequence_by_chr_index(right_exon, 0, 
				    chromosome_id - 1, 
				    &genome_start, &genome_end, genome);		  
  
  
  size_t dsp_l, dsp_r, type;
  
  size_t breaks_starts[gap_read];
  int found_starts = 0;
  
  size_t breaks_ends[gap_read];
  int found_ends = 0;
  
  size_t type_starts[gap_read];
  size_t type_ends[gap_read];
  int c_s, c_e;
  int end_search = gap_read + SEQURITY_FLANK;

  //printf("START READ START %i/%lu READ END %i/%lu\n", read_start, s_prev->genome_end, 
  //	 read_end, s_next->genome_start);
  
  //printf("Left  exon: %s\n", left_exon);
  //printf("Right exon: %s\n", right_exon);
  
  //printf("search start!\n");
  // Search step by step (GT)/(AG) 
  for (c_s = 0, c_e = strlen(right_exon) - 1;
       c_s < end_search; c_s++, c_e--) {
    if (left_exon[c_s] == 'G' && left_exon[c_s + 1] == 'T') {
      type_starts[found_starts] = GT_AG_SPLICE;
      breaks_starts[found_starts++] = c_s;
      //printf("S.FOUND GT (%i)\n", c_s);
    } else if (left_exon[c_s] == 'C' && left_exon[c_s + 1] == 'T') {
      type_starts[found_starts] = CT_AC_SPLICE;
      breaks_starts[found_starts++] = c_s;
      //printf("S.FOUND CT (%i)\n", c_s);
    }

    if (right_exon[c_e] == 'G' && right_exon[c_e - 1] == 'A') {
      type_ends[found_ends] = GT_AG_SPLICE;
      breaks_ends[found_ends++] = strlen(right_exon) - c_e - 1;
      //printf("E.FOUND AG (%i)\n", strlen(right_exon) - c_e - 1);
    } else if (right_exon[c_e] == 'C' && right_exon[c_e - 1] == 'A') {
      type_ends[found_ends] = CT_AC_SPLICE;
      breaks_ends[found_ends++] = strlen(right_exon) - c_e - 1;
      //printf("E.FOUND AC (%i)\n", strlen(right_exon) - c_e - 1);
    }
  }

  //Not found any splice junction
  if (found_starts == 0 || found_ends == 0) {
    return NULL;
  } 

  array_list_t *splice_junction = array_list_new(2*(found_starts + found_ends), 
						 1.25f, COLLECTION_MODE_ASYNCHRONIZED);


  //printf("FOUND STARTS %i, FOUND ENDS %i\n", found_starts, found_ends);

  //If found more than one break...
  for (int i = 0; i < found_starts; i++) {
    for (int j = 0; j < found_ends; j++) {
      //printf("%i(%i) + %i(%i)=%i\n", breaks_starts[i], type_starts[i],
      //breaks_ends[j], type_ends[j], breaks_starts[i] + breaks_ends[j]);
      if (type_starts[i] == type_ends[j]) {
	int gap_break = breaks_starts[i] + breaks_ends[j];
	//printf("::(%i, %i), ¿gap_break=%i - gap_read=%i?\n", breaks_starts[i], breaks_ends[j], gap_break, gap_read);
	if (abs(gap_break - gap_read) <= 4) {	  
	  //printf("\tEND:%i-START:%i\n", type_ends[j], type_starts[i]);
	  array_list_insert((void *)breaks_starts[i], splice_junction);
	  array_list_insert((void *)breaks_ends[j], splice_junction);
	  array_list_insert((void *)type_ends[j], splice_junction);
	}
      }
    }
  }


  char q[1024], rl[1024], rr[1024];
  int read_pos = read_start + 1;

  assert(read_pos + gap_read < strlen(sequence));

  memcpy(q, &sequence[read_pos], gap_read);
  q[gap_read] = '\0';

  size_t num_sj  = array_list_size(splice_junction);
  char **query_p = (char **)malloc(sizeof(char *)*num_sj);
  char **ref_p   = (char **)malloc(sizeof(char *)*num_sj);

  int n_sw = 0;

  for (int i = 0; i < num_sj; i += 3) {
    int limit_left  = (size_t)array_list_get(i, splice_junction);
    int limit_right = (size_t)array_list_get(i + 1, splice_junction);      
    
    int gap_r_exon = strlen(right_exon) - limit_right;

    assert(limit_left < 1024);
    memcpy(rl, left_exon, limit_left);
    //printf(">>>>>>> %s (%i)\n", rl, strlen(rl));

    assert(limit_right + limit_left + 1 < 1024);

    memcpy(&rl[limit_left], &right_exon[gap_r_exon], limit_right);
    rl[limit_left + limit_right] = '\0';

    //printf(">>>>>>> %s (%i)\n", rl, strlen(rl));

    //printf("%i - %i (gap_read = %i, gap_genome = %i)\n", limit_left, limit_right, gap_read, limit_left + limit_right);
    //printf("Sequence: %s\n", q);
    //printf("Referenc: %s\n", rl);
    
    query_p[n_sw] = strdup(q);
    ref_p[n_sw++] = strdup(rl);

  }

  sw_multi_output_t *output = sw_multi_output_new(n_sw);
  smith_waterman_mqmr(query_p, ref_p, n_sw,
		      sw_optarg, 1,
		      output);

  float match = sw_optarg->subst_matrix['A']['A'];
  int sw_distance;
  float max_score = 0.0f;
  int sj_select = -1;

  for (int i = 0; i < n_sw; i++) {    
    //printf("==========================SW============================\n");
    //printf("REF: %s\n", output->ref_map_p[i]);
    //printf("SEQ: %s\n", output->query_map_p[i]);
    //printf("======================================================\n");
    
    //printf("/////////// cigar : %s\n", new_cigar_code_string(cigar_code));
    
    float norm_score = NORM_SCORE(output->score_p[i], strlen(query_p[i]), match);
    
    if (norm_score > max_score) {
      max_score = norm_score;
      sj_select = i;
    }
    
    //printf("norm_score = %f\n", norm_score);   
        
  }

  cigar_code_t *cc_sj = NULL;
  if (sj_select >= 0) {
    cc_sj = generate_cigar_code(output->query_map_p[sj_select],
				output->ref_map_p[sj_select],
				strlen(output->ref_map_p[sj_select]),
				output->query_start_p[sj_select], 
				output->ref_start_p[sj_select],
				strlen(query_p[sj_select]), 
				strlen(ref_p[sj_select]),
				&sw_distance, MIDDLE_SW);
  }
  
  
  for (int i = 0; i < n_sw; i++) {    
    free(output->query_map_p[i]);
    free(output->ref_map_p[i]);
    
    free(query_p[i]);
    free(ref_p[i]);

    output->query_map_p[i] = NULL;
    output->ref_map_p[i] = NULL;
  }
  
  free(query_p);
  free(ref_p);

  sw_multi_output_free(output);
  
  if (sj_select < 0) {
    goto exit;
  }

  dsp_l = (size_t)array_list_get((sj_select * 3), splice_junction);
  dsp_r = (size_t)array_list_get((sj_select * 3) + 1, splice_junction);
  type  = (size_t)array_list_get((sj_select * 3) + 2, splice_junction);  
 
  if (type == CT_AC_SPLICE) {
    strand = 1;
  } else {
    strand = 0;
  }  
  
  
  size_t start_splice = s_prev->genome_end   + dsp_l + 1;
  size_t end_splice   = s_next->genome_start - dsp_r - 1;
  
  *sp_start = start_splice;
  *sp_end   = end_splice;
  *sp_type  = type;

  //printf("SP :=> [%i:%lu-%lu]\n", chromosome_id, start_splice, end_splice);
  intron_size = end_splice - start_splice + 1;

  if (intron_size < min_intron) {
    goto exit;
  }

  //Insert SJ in cigar
  cc_final = cigar_code_new();

  int limit_ = dsp_l;
  int limit_tmp = limit_;

  //printf("S.limit_ = %i, cc_sj = %s\n", limit_, new_cigar_code_string(cc_sj));
  for (int c = 0; c < cc_sj->ops->size; c++) {
    cigar_op_t *op = array_list_get(c, cc_sj->ops);
    //printf("op: %i%c, limit_ = %i\n", op->number, op->name, limit_);
    
    if ((limit_ > 0) &&
	(op->name == 'M' || op->name == 'D')) {
      limit_tmp = limit_ - op->number;

      //printf("limit_ = %i, limit_tmp_ = %i\n", limit_, limit_tmp);    
      if (limit_tmp > 0) {
	//printf("1.\n");
	cigar_code_append_new_op(op->number, op->name, cc_final);
      } else if (limit_tmp == 0) {
	//printf("2.\n");
	cigar_code_append_new_op(op->number, op->name, cc_final);
	cigar_code_append_new_op(intron_size, 'N', cc_final);
      } else if (limit_tmp < 0) {
	//printf("3. %i%c, %iN, %i%c\n", abs(limit_), op->name, intron_size, op->number - limit_, op->name);
	cigar_code_append_new_op(limit_, op->name, cc_final);
	cigar_code_append_new_op(intron_size, 'N', cc_final);
	cigar_code_append_new_op(op->number - limit_, op->name, cc_final);
      }
      
      limit_ -= (int)op->number;

    } else {
      cigar_code_append_new_op(op->number, op->name, cc_final);
    }
    
    cigar_op_free(op);

    //printf("FINAL CIGAR: %s\n", new_cigar_code_string(cc_final));

  }

  cigar_code_free(cc_sj);

 exit:
  array_list_free(splice_junction, NULL);

  return cc_final;

}

