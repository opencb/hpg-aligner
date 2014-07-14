#include "sa_rna_mapper.h"

#define MAX_DEPTH 4
#define DEBUG 1

const int MAX_INTRON_SIZE = 500000;
const int MIN_INTRON_SIZE = 40;

void SA_DEBUG(char *msg, ...) {
  if (DEBUG) {
    printf("%s", msg);
  }
}

//Sa Mappings Reader, 2nd round
void *sa_alignments_reader_rna(void *input) {
  sa_wf_input_t *wf_input = (sa_wf_input_t *) input;  
  sa_wf_batch_t *new_wf_batch = NULL;
  sa_wf_batch_t *curr_wf_batch = wf_input->wf_batch;  
  FILE *fd = ((sa_rna_input_t *)curr_wf_batch->data)->file1;

  const int MAX_READS = 100;
  fastq_batch_reader_input_t *fq_reader_input = wf_input->fq_reader_input;
  array_list_t *reads = array_list_new(MAX_READS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  size_t num_reads = 0;
  array_list_t **mapping_lists = (array_list_t **) malloc(MAX_READS * sizeof(array_list_t *));
  size_t num_items;

  //printf("===== Read from %p =====\n", fd);

  while (1) {
    fastq_read_t *fq_read = file_read_fastq_reads(&num_items, fd);

    if (fq_read == NULL) { printf("fq NULL\n"); break; }

    //printf("(num items %i)\nID : %s\nSEQ: %s\nQUA: %s\n", num_items, fq_read->id, fq_read->sequence, fq_read->quality);
    array_list_insert(fq_read, reads);
    
    mapping_lists[num_reads] = array_list_new(50,
					      1.25f, 
					      COLLECTION_MODE_ASYNCHRONIZED);
    
    //printf("Read alignments\n");
    sa_file_read_alignments(num_items, mapping_lists[num_reads],
			    fq_read, fd);

    num_reads++;

    if (num_reads >= MAX_READS) { break; }

  }


  if (num_reads) {
    sa_batch_t *sa_batch = sa_batch_simple_new(reads);
    sa_batch->mapping_lists = mapping_lists;
    sa_batch->num_reads = num_reads;
    new_wf_batch = sa_wf_batch_new(NULL,
				   curr_wf_batch->sa_index,
				   curr_wf_batch->writer_input, 
				   sa_batch,
				   curr_wf_batch->data);
  } else {
    array_list_free(reads, (void *)fastq_read_free);
    free(mapping_lists);
  }

  return new_wf_batch;

}


//Fastq Reader
void *sa_fq_reader_rna(void *input) {
  sa_wf_input_t *wf_input = (sa_wf_input_t *) input;
  
  sa_wf_batch_t *new_wf_batch = NULL;
  sa_wf_batch_t *curr_wf_batch = wf_input->wf_batch;
  
  fastq_batch_reader_input_t *fq_reader_input = wf_input->fq_reader_input;
  array_list_t *reads = array_list_new(fq_reader_input->batch_size, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

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
      fastq_fread_bytes_se(reads, fq_reader_input->batch_size, fq_reader_input->fq_file1);
    } else {
      fastq_fread_bytes_aligner_pe(reads, fq_reader_input->batch_size, 
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
				   curr_wf_batch->data);
  }

  return new_wf_batch;

}

//Fastq Writer
int sa_sam_writer_rna(void *data) {
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;
  
  sa_batch_t *mapping_batch = (sa_batch_t *) wf_batch->mapping_batch;

  int flag, pnext = 0, tlen = 0;
  char rnext[4] = "*\0";
  char optional_flags[512] = "NM:i:3\0";

  fastq_read_t *read;
  array_list_t *read_list = mapping_batch->fq_reads;

  alignment_t *alig;
  array_list_t *mapping_list;
  FILE *out_file = (FILE *) wf_batch->writer_input->bam_file;

  sa_genome3_t *genome = wf_batch->sa_index->genome;

  size_t num_reads, num_mappings;
  num_reads = mapping_batch->num_reads;

  extern size_t total_reads, unmapped_reads, correct_reads;
  extern pthread_mutex_t mutex_sp;

  //pthread_mutex_lock(&mutex_sp);
  //total_reads += num_reads;
  //pthread_mutex_unlock(&mutex_sp);

  for (size_t i = 0; i < num_reads; i++) {
    read = (fastq_read_t *) array_list_get(i, read_list);
    mapping_list = mapping_batch->mapping_lists[i];
    num_mappings = array_list_size(mapping_list);
    //printf("%i.Read %s (num_mappings %i)\n", i, read->id, num_mappings);
    if (num_mappings > 0) {
      //============================= Delete!!! ===============================
      //For debug... Validate Reads
      //--- Extract correct position ---//
      
      int CHROMOSOME;
      size_t START, END;
      
      char *id = read->id;
      int c = 0, len = strlen(id);
      int pos = 0;
      while (c < len && pos < 3) {
	//printf("while pos [%c]\n", id[c]);
	if (id[c++] == '@') { pos++; }
      }
      
      if (pos) {
	//printf("Actual pos %c\n", id[c]);
	char value[128];
	int p = 0;
	while (id[c] != '@') {
	  value[p++] = id[c++];
	}
	c++;
	value[p] = '\0';
	if (strcmp(value, "X") == 0) { CHROMOSOME = 23; }
	else if (strcmp(value, "Y") == 0) { CHROMOSOME = 24; }
	else if (strcmp(value, "MT") == 0) { CHROMOSOME = 25; }
	else { CHROMOSOME = atoi(value); }
      
	//printf("Actual pos %c\n", id[c]);
	p = 0;
	while (id[c] != '@') {
	  //printf("while pos [%c]\n", id[c]);
	  value[p++] = id[c++];
	}
	c++;
	value[p] = '\0';
	START = atol(value);
      
      
	p = 0;
	while (id[c] != '@') {
	  //printf("while pos [%c]\n", id[c]);
	  value[p++] = id[c++];
	}
	c++;
	value[p] = '\0';
	END = atol(value);
      
	//---                          ---//
	int map = 0;
	for (int j = 0; j < array_list_size(mapping_list); j++) {
	  alignment_t *alignment = array_list_get(j, mapping_list);
	  if ((alignment->chromosome - 1) == CHROMOSOME - 1 && 
	      alignment->position >= START && 
	      alignment->position <= END) {
	    map = 1;
	    break;
	  }
	}
      
	//if (map) { 
	  //pthread_mutex_lock(&mutex_sp);
	  //correct_reads++;
	  //pthread_mutex_unlock(&mutex_sp);
	//} else {
	  //pthread_mutex_lock(&mutex_sp);
	  //printf("%s:\n", read->id);
	  //for (int j = 0; j < array_list_size(mapping_list); j++) {
	  //alignment_t *alignment = array_list_get(j, mapping_list);
	  //printf("\t%i:%lu\n", alignment->chromosome, alignment->position);	    
	  //}
	  //pthread_mutex_unlock(&mutex_sp);
	//}
      }

      //============================= Delete!!! ===============================      
      for (size_t j = 0; j < num_mappings; j++) {
	alig = (alignment_t *) array_list_get(j, mapping_list);
	flag = (alig->seq_strand ? 16 : 0);
	//printf("%s\n", alig->query_name);
	//printf("AAAAA %i %s\n",  alig->position, genome->chrom_names[alig->chromosome - 1]);
	//printf("%i %s", alig->position, alig->cigar);
	fprintf(out_file, "%s\t%i\t%s\t%lu\t%i\t%s\t%s\t%lu\t%i\t%s\t%s\t%s\n", 
		alig->query_name,
		flag,
		genome->chrom_names[alig->chromosome - 1],
		alig->position + 1,
		alig->map_quality,
		alig->cigar,
		rnext,
		pnext,
		tlen,
		alig->sequence,
		alig->quality,
		optional_flags
		);

	alignment_free(alig);	 
      }
    } else {
      //printf("###@@@### : %s\n", read->id);
      //pthread_mutex_lock(&mutex_sp);
      //unmapped_reads++;
      //pthread_mutex_unlock(&mutex_sp);

      fprintf(out_file, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", 
	      read->id,
	      read->sequence,
	      read->quality
	      );
    }
    array_list_free(mapping_list, (void *) NULL);
  }

  // free memory
  sa_batch_free(mapping_batch);

  if (wf_batch) sa_wf_batch_free(wf_batch);

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
	  printf(" [%lu-%lu](%i) :\n", cal->start, cal->end, linked_list_size(cal->sr_list));
	  seed_first = linked_list_get_first(cal->sr_list);
	  seed_last = linked_list_get_last(cal->sr_list);

	  list_item = cal->sr_list->first;
	  while (list_item) {
	    seed_first = list_item->item;
	    printf("\t[%lu|%lu-%lu|%lu] ", seed_first->genome_start, seed_first->read_start, seed_first->read_end, 
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
    
      //printf("REF: %s\n", output->ref_map_p[i]);
      //printf("SEQ: %s\n", output->query_map_p[i]);
      //printf("======================================================\n");

      //printf("/////////// cigar : %s\n", new_cigar_code_string(cigar_code));
      ((seed_region_t *)sw_depth->item_ref[i])->info = cigar_code;

      free(sw_depth->q[i]);
      free(sw_depth->r[i]);
    
      free(output->query_map_p[i]);
      free(output->ref_map_p[i]);
      output->query_map_p[i] = NULL;
      output->ref_map_p[i] = NULL;
    } else {
      printf("REF: %s\n", output->ref_map_p[i]);
      printf("SEQ: %s\n", output->query_map_p[i]);      
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
		   int min_intron_size,
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
  const int FLANK = 5;

  //cal_t *cal = array_list_get(i, meta_alignment->cals_list);
  seed_region_t *seed_prev = NULL, *seed_next;
  linked_list_item_t *list_item = cal->sr_list->first, *list_item_prev;
  size_t sp_start, sp_end;
  int distance_aux;

  while (list_item != NULL) {
    seed_next = (seed_region_t *)list_item->item;
    int num_match = seed_next->read_end - seed_next->read_start + 1;
    if (seed_prev != NULL) {
      //Close nt 
      if (seed_next->read_start <= seed_prev->read_end || 
	  seed_next->genome_start <= seed_prev->genome_end) {
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
	  free(reference);
	  max_size = genome_end - genome_start + 1024;
	  reference = (char *)calloc(max_size, sizeof(char));
	}
	genome_read_sequence_by_chr_index(reference, 0, cal->chromosome_id - 1,
					  &genome_start, &genome_end, genome);
	
	//printf("[%lu|%i]-GAP-[%i|%lu]\n", genome_start, read_start, read_end, genome_end);
	for (int k = 0; k < gap_read; k++) {
	  //printf("[q:%c vs r:%c]\n", query[k], reference[k]);
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
	if (gap > min_intron_size) {
	  //printf("Search splice");
	  int sp_type;
	  /*
	  int nt = search_simple_splice_junction_semi_cannonical(seed_prev, seed_next,
								 cal->chromosome_id,
								 cal->strand,
								 query_map, genome, 
								 &sp_start,
								 &sp_end,
								 &sp_type,
								 &distance_aux);
								 */
	  int nt = search_simple_splice_junction(seed_prev, seed_next,
						 cal->chromosome_id, cal->strand, 
						 query_map, genome, 
						 &sp_start, 
						 &sp_end,
						 &sp_type,
						 &distance_aux);
	  
	  if (nt && 
	      seed_prev->genome_start < sp_start && 
	      sp_end < seed_next->genome_end && 
	      sp_start < sp_end) {

	    //printf("FOUND SPLICE CAL Middle fill gaps: (%lu-%lu)%i\n", sp_start, sp_end, nt);
	    allocate_start_node(cal->chromosome_id - 1,
				cal->strand,
				sp_start,
				sp_end,
				sp_start,
				sp_end,
				FROM_READ,
				sp_type,
				NULL, 
				&node_avl_start,
				&node_avl_end, 
				avls_list);		  

	    //assert(seed_prev->genome_start < node_avl_start->position);
	    //assert(node_avl_end->position < seed_next->genome_end);

	    metaexon_insert(cal->strand, cal->chromosome_id - 1,
			    seed_prev->genome_start, sp_start, 40,
			    METAEXON_RIGHT_END, node_avl_start,
			    metaexons);
	    
	    metaexon_insert(cal->strand, cal->chromosome_id - 1,
			    sp_end, seed_next->genome_end, 40,
			    METAEXON_LEFT_END, node_avl_end,
			    metaexons);

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
	//int read_start   = seed_prev->read_end + 1;
	//int read_end     = seed_next->read_start - 1;
	//int genome_start = seed_prev->genome_end + 1;
	//int genome_end   = seed_next->genome_start - 1;
	
	//int diff_read   = abs(read_end - read_start);
	//int diff_genome = abs(genome_end - genome_start);
	//if (diff_read <= 5 ||
	//  diff_genome <= 5) {
	//int diff_min = diff_read < diff_genome ? diff_read : diff_genome;
	//int sub_seed = 5 - diff_min;

	//seed_prev->read_end     -= sub_seed;
	//seed_prev->genome_end   -= sub_seed;
	//seed_next->read_start   += sub_seed;
	//seed_next->genome_start += sub_seed;

	//}
	
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

	  genome_read_sequence_by_chr_index(reference, 0, cal->chromosome_id - 1,
					    &genome_start, &genome_end, genome);
	  char query[2048];

	  assert(read_end - read_start + 1 < 2048);
	  if (read_end - read_start + 1 <= 0) {
	    //fprintf(stderr, "%i - %i + 1 = %i\n", read_end, read_start, read_end - read_start + 1);
	    break;
	  }
	  
	  memcpy(query, &query_map[read_start],  read_end - read_start + 1);
	  query[read_end - read_start + 1] = '\0';
	  
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

    //int num_targets = array_list_size(target_cals);
    //int search_sp = 0;

    /*
    printf("Num targets %i\n", num_targets);
    for (int i = 0; i < num_targets; i++) {
      printf("\tDouble anchors found %i/%i\n", i + 1, num_targets);
      merge_cals = array_list_get(i, target_cals);
      //num_cals = array_list_size(merge_cals);
      
      sa_alignment_t *sa_alignment = sa_alignment_new(merge_cals);
      array_list_insert(sa_alignment, sa_alignments_list[r]);	
      
    }
    */
    
    /*
    //printf("<----------------------- Fill gaps -------------------->\n");
    //Process Results. Are CALs in target_cals??
	int num_targets = array_list_size(target_cals);
	int search_sp = 0;
      
	//printf("Num targets %i\n", num_targets);
	for (int i = 0; i < num_targets; i++) {
	  //printf("Double anchors found %i/%i\n", i + 1, num_targets);
	  merge_cals = array_list_get(i, target_cals);
	  num_cals = array_list_size(merge_cals);
	
	  //Bad CAL??
	  int ok = 1;
	  for (int c = 0; c < num_cals; c++) {
	    cal = array_list_get(c, merge_cals);
	    linked_list_item_t *item_prev = cal->sr_list->first, *item = NULL;
	    seed_region_t *seed_prev, *seed_next;
	    while (item_prev) {
	      item = item_prev->next;
	      if (!item) { break; }
	      else {
		seed_prev = item_prev->item;
		seed_next = item->item;
		if (seed_prev->id == seed_next->id || 
		    seed_next->genome_start < seed_prev->genome_end || 
		    seed_next->read_start < seed_prev->read_end) { ok = 0; break; }
	      }
	      item_prev = item;
	    }
	  }
	
	  if (!ok) { continue; }//goto free; }
	  
	  //printf("==================================\n");
	  //for (int w = 0; w < num_cals; w++) {
	  //cal_t *cal_aux = array_list_get(w, merge_cals);
	  //cal_print(cal_aux);
	  //}
	  //printf("==================================\n");
	  
	  if (num_cals <= 1) {
	    //Fill gaps, one CAL
	    cal = array_list_get(0, merge_cals);
	    int num_seeds = linked_list_size(cal->sr_list);
	    //cal_print(cal);
	    //printf("num_seeds = %i\n", num_seeds);
	    seed_first = linked_list_get_first(cal->sr_list);
	    seed_last = linked_list_get_last(cal->sr_list);
	    //printf("[%i:%lu-%lu]num cals = %i, num_seeds = %i\n", cal->chromosome_id, cal->start, cal->end, num_cals, num_seeds);
	  
	    if (num_seeds < 2) { 
	      //SINGLE ANCHOR - ONE SEED
	      seed_first = linked_list_get_first(cal->sr_list);
	      
	      int seed_size = seed_first->read_end - seed_first->read_start + 1;
	      if (seed_size > MIN_SEED_SIZE) {	    
		if (cal->strand) {
		  query = seq_revcomp;
		} else {
		  query = seq;
		}
		int distance = 0;
		int read_start, read_end;
		int lim;
		
		metaexon_t *metaexon;
		pthread_mutex_lock(&metaexons->mutex[cal->chromosome_id - 1]);
		
		metaexon_search(cal->strand, cal->chromosome_id - 1,
				cal->start, cal->end, &metaexon,
				metaexons);
	      
		//cal_print(cal);
		cigar_code_t *cc;
		if (metaexon != NULL) {
		  if (seed_first->read_start == 0) {
		    int gap_close = len_seq - (seed_first->read_end + 1);
		    cc =search_left_single_anchor(gap_close, 
						  cal,
						  0,
						  metaexon->right_breaks,
						  query,
						  metaexons,
						  genome,
						  avls_list);
		  } else {
		    int gap_close = seed_first->read_start;
		    cc = search_right_single_anchor(gap_close, 
						    cal,
						    0,
						    metaexon->left_breaks,
						    query,
						    metaexons,
						    genome, 
						    avls_list);  
		  }
		
		  if (cc) {
		    char cigar_str[2048];
		    cigar_op_t *op;
		    if (seed_first->read_start == 0) {
		      op = array_list_get(0, cc->ops);
		    } else {
		      op = array_list_get(cc->ops->size - 1, cc->ops);
		    }
		    op->number += (seed_first->read_end - seed_first->read_start + 1);
		  
		    //for (int k = 0; k < cc->ops->size; k++) {
		    //cigar_op_t *op = array_list_get(k, cc->ops);	      
		    // sprintf(cigar_str, "%s%i%c", cigar_str, op->number, op->name);
		    //cigar_op_free(op);
		    //}
		  
		    //cigar_code_free(cc);
		    array_list_t *array_list_tmp = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
		    array_list_insert(cal, array_list_tmp);
		    sa_alignment_aux = sa_alignment_new(array_list_tmp);
		    sa_alignment_aux->c_final = cc;
		    array_list_insert(sa_alignment_aux, sa_alignments_list[r]);		    
		    //printf("Meta-Map CIGAR:%s\n", cigar_str);	  
		    //alignment_t *alignment = alignment_new();
		    //alignment_init_single_end(strdup(read->id), 
		    //			      strdup(query), 
		    //			      strdup(read->quality),
		    //			      cal->strand,
		    //			      cal->chromosome_id,
		    //			      cal->start,
		    //			      strdup(cigar_str), 1, 254, 1, 0,
		    //			      0, 0, alignment);
		    //array_list_insert(alignment, alignments_list);
		    
		  }
		}
		pthread_mutex_unlock(&metaexons->mutex[cal->chromosome_id - 1]);
	      }
	    } else {
	      //DOBLE ANCHOR - TWO SEEDS	    
	      int read_gap = (int)seed_last->read_start - (int)seed_first->read_end - 1;
	      int genome_gap = (int)seed_last->genome_start - (int)seed_first->genome_end - 1;
	      int gaps_dif = abs(genome_gap - read_gap);

	      //printf(" ::::::::::::::::::::::::: read_gap = %i, genome_gap = %i, gaps_dif = %i\n", read_gap, genome_gap, gaps_dif);	  

	      if (read_gap <= 4 || genome_gap <= 4) {
		int dsp = 5;
	      
		int min = dsp*-1;
		if (read_gap < min || genome_gap < min) {
		  dsp = abs(read_gap) < abs(genome_gap) ? abs(read_gap) : abs(genome_gap);
		}
		
		seed_first->read_end    -= dsp;
		seed_first->genome_end  -= dsp;
		seed_last->read_start   += dsp;
		seed_last->genome_start += dsp;
		
		read_gap = (int)seed_last->read_start - (int)seed_first->read_end - 1;
		genome_gap = (int)seed_last->genome_start - (int)seed_first->genome_end - 1;
		
	      }
	      
	      if (cal->strand) {
		query = seq_revcomp + seed_first->read_end + 1;
	      } else {
		query = seq + seed_first->read_end + 1;
	      }

	      //printf("%i - %i\n", read_gap, genome_gap);	  

	      assert(read_gap >= 0);
	      assert(genome_gap >= 0);
	    
	      if (genome_gap == read_gap) {
		//printf("Close");
		//Close gap
		int distance = 0;
		genome_start = seed_first->genome_end + 1;
		genome_end = seed_last->genome_start - 1;
		
		genome_read_sequence_by_chr_index(reference, 0, cal->chromosome_id - 1,
						  &genome_start, &genome_end, genome);    
	
		//printf("genome read reference=> (%i) %i:%lu-%lu\n", cal->strand, cal->chromosome_id, genome_start, genome_end);
		//printf("Ref: %s\n", reference);
		//printf("Que: %s\n", query);
		for (int c = 0; c < read_gap; c++) {
		  if (reference[c] != query[c]) { distance++; }
		}
	      
		//printf("============( %i vs %i )============\n", distance, (read_gap / 5) + 1);

		if (distance <= (read_gap / 5) + 1) {
		  //printf("Read map %i\n", distance);
		  char cigar_str[2048];
		  sprintf(cigar_str, "%luM", len_seq);
		  //printf("chrom %i\n", cal->chromosome_id);
		  alignment_t *alignment = alignment_new();
		  alignment_init_single_end(strdup(read->id), 
					    strdup(read->sequence), 
					    strdup(read->quality),
					    cal->strand, cal->chromosome_id,
					    cal->start,
					    strdup(cigar_str), 1, 254, 1, 0,//??(num_suffixes_n > 1),
					    0, 0, alignment);
		  array_list_insert(alignment, alignments_list);
		
		  metaexon_insert(cal->strand, 
				  cal->chromosome_id - 1,
				  cal->start,
				  cal->start + len_seq,
				  MIN_INTRON_SIZE,
				  METAEXON_RIGHT_END,
				  NULL,
				  metaexons);
		
		
		}
	      } else if (gaps_dif > MIN_INTRON_SIZE) {
		search_sp = 1;
	      } 
	    }
	  } else { //num_cals > 1
	    search_sp = 1;
	  } //num_cals <= 1
	  
	  if (search_sp) {
	    if (num_cals > 1) {
	      cal_prev = array_list_get(0, merge_cals);
	      seed_first = linked_list_get_last(cal_prev->sr_list);
	    
	      cal = array_list_get(array_list_size(merge_cals) - 1, merge_cals);
	      seed_last = linked_list_get_first(cal->sr_list);
	    
	      int read_gap = (int)seed_last->read_start - (int)seed_first->read_end - 1;
	      int genome_gap = (int)seed_last->genome_start - (int)seed_first->genome_end - 1;
	    
	      //cal_print(cal_prev);
	    
	      if (read_gap < 0 || genome_gap < 0) {
		int dsp = 5;
	      
		int min = dsp*-1;
		if (read_gap < min || genome_gap < min) {
		  dsp = abs(read_gap) < abs(genome_gap) ? abs(read_gap) : abs(genome_gap);
		}
	      
		if (dsp > 20) {
		  goto free;
		}
	      
		seed_first->read_end    -= dsp;
		seed_first->genome_end  -= dsp;
		seed_last->read_start   += dsp;
		seed_last->genome_start += dsp;
	      
	      }	  
	    } else {
	      cal_prev = array_list_get(0, merge_cals);
	      seed_first = linked_list_get_first(cal_prev->sr_list);
	    
	      cal = cal_prev;
	      seed_last = linked_list_get_last(cal->sr_list);
	    }
	  
	    if (cal->strand) {
	      query = seq_revcomp;
	    } else {
	      query = seq;
	    }
	  	  
	    //printf("%i:%lu-%lu <-> %lu-%lu\n", cal_prev->chromosome_id, cal->chromosome_id, seed_first->genome_start, 
	    //	   seed_first->genome_end, seed_last->genome_start, seed_last->genome_end);

	    //Map with metaexon
	    int meta_type;
	    int seed_left_len  = seed_first->read_end - seed_first->read_start + 1;
	    int seed_right_len = seed_last->read_end - seed_last->read_start + 1;
	    char cigar_str[2048] = "\0";

	    //cal_print(cal_prev);
	    cigar_code_t *cc = search_double_anchors_cal(query,
							 cal_prev, 
							 cal,
							 metaexons, 
							 genome,
							 read,
							 &meta_type, 
							 avls_list);
	  
	    if (cc != NULL) {
	      cigar_op_t *op = array_list_get(0, cc->ops);
	      op->number += seed_left_len;
	      if (meta_type == META_ALIGNMENT_MIDDLE) {	    
		cigar_op_t *op = array_list_get(cc->ops->size - 1, cc->ops);
		op->number += seed_right_len;
	      }

	      // for (int k = 0; k < cc->ops->size; k++) {
	      //cigar_op_t *op = array_list_get(k, cc->ops);	      
	      //sprintf(cigar_str, "%s%i%c", cigar_str, op->number, op->name);
	      //cigar_op_free(op);
	      //}
	      //cigar_code_free(cc);

	      //printf("Meta-Map CIGAR:%s\n", cigar_str);
	    } else {
	      //cal_print(cal);	    
	      int read_gap = (int)seed_last->read_start - (int)seed_first->read_end - 1;
	      if (read_gap > 18) {
		continue;
	      }
	      
	      //int nt = search_simple_splice_junction_semi_cannonical(seed_first, seed_last,
	      //						     cal->chromosome_id,
	      //						     cal->strand,
	      //						     query, genome, 
	      //						     &sp_start,
	      //						     &sp_end,
	      //						     &sp_type,
	      //						     &distance);
	      
	      int nt = search_simple_splice_junction(seed_first, seed_last,
						     cal->chromosome_id,
						     cal->strand,
						     query, genome, 
						     &sp_start,
						     &sp_end,
						     &sp_type,
						     &distance);	      

	      if (nt < MIN_INTRON_SIZE || sp_start <= seed_first->genome_start || sp_end >= seed_last->genome_end) { 
		continue;
	      }
	    
	      allocate_start_node(cal->chromosome_id - 1,
				  cal->strand,
				  sp_start,
				  sp_end,
				  sp_start,
				  sp_end,
				  FROM_READ,
				  sp_type,
				  NULL, 
				  &node_avl_start,
				  &node_avl_end, 
				  avls_list); 	 
	    
	      assert(seed_first->genome_start < sp_start);
	      assert(sp_end < seed_last->genome_end);

	      metaexon_insert(cal->strand, cal->chromosome_id - 1,
			      seed_first->genome_start, sp_start,
			      MIN_INTRON_SIZE,
			      METAEXON_RIGHT_END, node_avl_start,
			      metaexons);	 
	  
	      metaexon_insert(cal->strand, cal->chromosome_id - 1,
			      sp_end, seed_last->genome_end,
			      MIN_INTRON_SIZE,
			      METAEXON_LEFT_END, node_avl_end,
			      metaexons);

	      seed_left_len  = sp_start - seed_first->genome_start;
	      seed_right_len = seed_last->genome_end - sp_end;
	  
	      //printf("(%lu:%lu) | (start=%lu, end=%lu)------------------------> Search simple splice nt = %i\n", sp_start, sp_end, seed_first->genome_start, seed_last->genome_end, nt);

	      //sprintf(cigar_str, "%iM%iN%iM", seed_left_len,
	      //      nt, seed_right_len);	  
	  
	      //printf("Simple Search CIGAR:%s\n", cigar_str);	  
	      cc = cigar_code_new();
	      cigar_code_append_new_op(seed_left_len, 'M', cc);
	      cigar_code_append_new_op(nt, 'N', cc);
	      cigar_code_append_new_op(seed_right_len, 'M', cc);
	    }

	    array_list_t *array_list_tmp = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
	    if (num_cals > 1) {
	      array_list_insert(cal, array_list_tmp);
	      array_list_insert(cal_prev, array_list_tmp);
	    } else {
	      array_list_insert(cal, array_list_tmp);
	    }

	    sa_alignment_aux = sa_alignment_new(array_list_tmp);
	    sa_alignment_aux->c_final = cc;
	    array_list_insert(sa_alignment_aux, sa_alignments_list[r]);	    
	    //printf("============ Insert to list %i:%lu : %s ==========\n", cal->chromosome_id, cal->start, new_cigar_code_string(cc));
	  }
	  
	free:
	  kk = 1;
	}//for target_cals
	
	for (int m = 0; m < num_targets; m++) {
	  merge_cals = array_list_get(m, target_cals);
	  array_list_free(merge_cals, (void *)NULL);       
	}
	
	array_list_clear(target_cals, (void *)NULL);

	}*/
    //}


int sa_rna_mapper(void *data) {
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;
  sa_batch_t *sa_batch = wf_batch->mapping_batch;
  sa_index3_t *sa_index = (sa_index3_t *) wf_batch->sa_index;
  sa_rna_input_t *sa_rna = wf_batch->data;
  genome_t *genome = (genome_t *)sa_rna->genome;
  avls_list_t *avls_list = (avls_list_t *)sa_rna->avls_list;
  metaexons_t *metaexons = (metaexons_t *)sa_rna->metaexons;
  sw_optarg_t *sw_optarg = sa_rna->sw_optarg;
  FILE *file1 = sa_rna->file1;
  FILE *file2 = sa_rna->file2;

  size_t num_reads = sa_batch->num_reads;
  int max_read_area, min_num_mismatches;
  float max_score;

  // CAL management
  //size_t num_cals;
  //seed_cal_t *cal;
  cal_mng_t *cal_mng;
  //array_list_t *cal_list;
  fastq_read_t *read;  
  //  int saved_pos[2][1024];
  // TODO !!! 20 = min. cal size
  uint min_cal_size = 20;
  int seed_size = 18;
  //Function old Seeding vars
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

  cal_mng_t * cal_mng_p = cal_mng_new(sa_index->genome); //(+)
  cal_mng_t * cal_mng_n = cal_mng_new(sa_index->genome); //(-)

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
  
  array_list_t *write_alignments = array_list_new(50, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  seed_region_t *seed_prev, *seed_next, *s_prev;
  int num_sa_alignments;
  float cals_score[2048];

  //printf("****************** Batch ******************\n");
  for (int r = 0; r < num_reads; r++) {
    delete_targets[r] = 0;
    read = array_list_get(r, sa_batch->fq_reads);
    seq = read->sequence;
    //printf("...== %s ==...\n", read->id);
    seq_revcomp = read->revcomp;
    alignments_list = sa_batch->mapping_lists[r];
    total_suffix = 0;
    len_seq = read->length;
    sa_alignments_list[r] = array_list_new(50, 1.25f,
					   COLLECTION_MODE_ASYNCHRONIZED);
    //printf("READ: %s\n", read->id);
    //******************************************//
    //****           1-Exact Reads          ****//
    //******************************************//
    //****                                  ****//
    //******************************************//
    
    //Search strand(+)
    num_suffixes_p = search_suffix(&seq[0], sa_index->k_value, 
				   MAX_NUM_SUFFIXES, sa_index, 
				   &low_p, &high_p, &suffix_len_p);
    
    //Search strand(-)
    num_suffixes_n = search_suffix(&seq_revcomp[0], sa_index->k_value, 
				   MAX_NUM_SUFFIXES, sa_index, 
				   &low_n, &high_n, &suffix_len_n);
    
    // exact search for exact reads
    if (suffix_len_p == len_seq) {
      //Report Exact Maps! (+)
      alignment_t *alignment;
      for (size_t suff = low_p; suff <= high_p; suff++) {
	chrom = sa_index->CHROM[suff];	
	g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];
	
	char cigar_str[2048];
	sprintf(cigar_str, "%luM", len_seq);     
	alignment = alignment_new();
	alignment_init_single_end(strdup(read->id), 
				  strdup(read->sequence), 
				  strdup(read->quality), 
				  0, chrom + 1,
				  g_start,
				  strdup(cigar_str), 1, 254, 1, (num_suffixes_p > 1),
				  0, 0, alignment);
	array_list_insert(alignment, alignments_list);
	/*
	metaexon_insert(0, 
			chrom,
			g_start,
			g_start + len_seq,
			MIN_INTRON_SIZE,
			METAEXON_RIGHT_END,
			NULL,
			metaexons);      
	*/
      }
    } 

    if (suffix_len_n == len_seq) {
      //Report Exact Map! (-)
      alignment_t *alignment;
      for (size_t suff = low_n; suff <= high_n; suff++) {
	chrom = sa_index->CHROM[suff];
	g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];
	
	char cigar_str[2048];
	sprintf(cigar_str, "%luM", len_seq);

	alignment = alignment_new();
	alignment_init_single_end(strdup(read->id), 
				  strdup(read->revcomp), 
				  strdup(read->quality), 
				  1, chrom + 1,
				  g_start, 
				  strdup(cigar_str), 1, 254, 1, (num_suffixes_n > 1),
				  0, 0, alignment);
	array_list_insert(alignment, alignments_list);
	/*
	metaexon_insert(1, 
			chrom,
			g_start,
			g_start + len_seq,
			MIN_INTRON_SIZE,
			METAEXON_RIGHT_END,
			NULL,
			metaexons);
	*/
      }
    } 

    if (array_list_size(alignments_list)) {
      continue;
    }

    //printf("Search double anchors...\n");
    //No Exact map found! Search with Errors and SJ
    //******************************************//
    //****     1-Two Seeds and close gap    ****//
    //******************************************//    
    //****              Read Xnt            ****//
    //****     |------------------------|   ****//
    //****     |----|              |----|   ****//
    //****      18nt       Gap      18nt    ****//
    //****          |<------------>|        ****//
    //******************************************//
      

    /*
    size_t low_p_2, high_p_2, suffix_len_p_2;
    size_t low_n_2, high_n_2, suffix_len_n_2;
      
    size_t num_suffixes_p_2, num_suffixes_n_2;
      
      
    //Search strand(+)
    num_suffixes_p_2 = search_suffix(&seq[len_seq - sa_index->k_value], sa_index->k_value, 
				     MAX_NUM_SUFFIXES, sa_index, 
				     &low_p_2, &high_p_2, &suffix_len_p_2);
      
    //Search strand(-)
    num_suffixes_n_2 = search_suffix(&seq_revcomp[len_seq - sa_index->k_value], sa_index->k_value, 
				     MAX_NUM_SUFFIXES, sa_index,
				     &low_n_2, &high_n_2, &suffix_len_n_2);
      
    //(+) |--S--|---------->
    if (num_suffixes_p <= MAX_SUFFIXES && suffix_len_p) {
      //printf("1.%i suffix\n", num_suffixes_p);
      total_suffix += num_suffixes_p;
      for (size_t suff = low_p; suff <= high_p; suff++) {	
	chrom = sa_index->CHROM[suff];
	g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom] + 1;
	//printf("(+)1.(%i:%lu)\n", chrom, g_start);
	generate_cals(chrom + 1, 0, g_start, g_start + suffix_len_p - 1, 
		      0, suffix_len_p - 1,
		      cal_mng_p->cals_lists[chrom], len_seq, 0);
      }
    }

    //(+) <-----------|--E--| and extend 
    if (num_suffixes_p_2 <= MAX_SUFFIXES && suffix_len_p_2) {
      //printf("2.%i suffix\n", num_suffixes_p_2);
      total_suffix += num_suffixes_p_2;
      for (size_t suff = low_p_2; suff <= high_p_2; suff++) {	
	chrom = sa_index->CHROM[suff];
	g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom] + 1;
	//printf("(+)2.(%i:%lu)\n", chrom, g_start);
	//Extend seeds
	genome_start = g_start - len_seq - 1;
	genome_end   = g_start - 1;	
	genome_read_sequence_by_chr_index(reference, 0, chrom,
					  &genome_start, &genome_end, genome);
	  
	//printf("Ref: %s\n", reference);
	//printf("Seq: %s\n", seq);
	int s_pos = len_seq - suffix_len_p_2 - 1, nt_extend = 0, r_pos = strlen(reference) - 1;
	while (s_pos >= 0 && seq[s_pos--] == reference[r_pos--]) {
	  //printf("(%i)%c vs (%i)%c\n", s_pos + 1, seq[s_pos + 1], r_pos + 1, reference[r_pos + 1]);
	  nt_extend++;
	}
	//printf("suffix_len_p_2 = %i\n", suffix_len_p_2);
	generate_cals(chrom + 1, 0, g_start - nt_extend, g_start + suffix_len_p_2 - 1, 
		      len_seq - (suffix_len_p_2 + nt_extend), len_seq - 1,
		      cal_mng_p->cals_lists[chrom], len_seq, 100);
      }
    }
      
    //(-) |--S--|---------->
    if (num_suffixes_n <= MAX_SUFFIXES && suffix_len_n) {
      //printf("3.%i suffix\n", num_suffixes_n);
      total_suffix += num_suffixes_n;
      for (size_t suff = low_n; suff <= high_n; suff++) {	
	chrom = sa_index->CHROM[suff];
	g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom] + 1;
	//printf("(-)1.(%i:%lu)\n", chrom, g_start);
	//printf("S----------->\n");
	generate_cals(chrom + 1, 1, g_start, g_start + suffix_len_n - 1, 
		      0, suffix_len_n - 1,
		      cal_mng_n->cals_lists[chrom], len_seq, 0);
      }
    }

    //(-) <-----------|--E--| and rev-extend
    if (num_suffixes_n_2 <= MAX_SUFFIXES && suffix_len_n_2) {
      //printf("4.%i suffix\n", num_suffixes_n_2);
      total_suffix += num_suffixes_n_2;
      for (size_t suff = low_n_2; suff <= high_n_2; suff++) {	
	chrom = sa_index->CHROM[suff];
	g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom] + 1;
	//printf("(-)2.(%i:%lu)\n", chrom, g_start);
	//printf("<-----------E\n");
	  
	//Extend seeds
	genome_start = g_start - len_seq - 1;
	genome_end   = g_start - 1;	
	genome_read_sequence_by_chr_index(reference, 0, chrom,
					  &genome_start, &genome_end, genome);
	  
	int s_pos = len_seq - suffix_len_n_2 - 1, nt_extend = 0, r_pos = strlen(reference) - 1;
	//printf("r-Ref: %s\n", reference);
	//printf("r-Seq: %s\n", seq_revcomp);
	//printf("%i(%c), %c\n", s_pos, seq_revcomp[s_pos], reference[r_pos]);
	  
	while (s_pos >= 0 && seq_revcomp[s_pos--] == reference[r_pos--]) {
	  //printf("(%i)%c vs (%i)%c\n", s_pos + 1, seq_revcomp[s_pos + 1], r_pos + 1, reference[r_pos + 1]);
	  nt_extend++;
	}
	  
	generate_cals(chrom + 1, 1, g_start - nt_extend, g_start + suffix_len_n_2 - 1,
		      len_seq - (suffix_len_n_2 + nt_extend), len_seq - 1,
		      cal_mng_n->cals_lists[chrom], len_seq, 100);
      }
    }
      

    printf("==================== CALs result =====================\n");
    printf("==(+)==\n");
    cal_mng_print(cal_mng_p);
    printf("==(-)==\n");
    cal_mng_print(cal_mng_n);
    printf("======================================================\n");


    array_list_clear(target_cals, (void *)NULL);
    if (total_suffix < 50) {      
      for (int s = 0; s < 2; s++) {
	if (!s) { p = cal_mng_p; }
	else  { p = cal_mng_n; }
	
	for (unsigned int i = 0; i < p->num_chroms; i++) {
	  cal_list = p->cals_lists[i];
	  num_items = linked_list_size(cal_list);
	  if (num_items) {
	    linked_list_iterator_init(cal_list, &itr);	
	    //------------------------------------------------------------------------
	    cal_prev = linked_list_iterator_curr(&itr);
	    while (cal_prev) {
	      seed_first = linked_list_get_first(cal_prev->sr_list);
	      seed_last = linked_list_get_last(cal_prev->sr_list);
	      cal = linked_list_iterator_next(&itr);	  
		
	      if (cal) {
		seed_first_next = linked_list_get_first(cal->sr_list);
		//cal_print(cal_prev);
		//cal_print(cal);
		//printf("%i <= %i && %i != %i\n", cal->start, cal_prev->end, seed_first_next->id, seed_last->id);
		int ok = 0;
		if (seed_first_next->read_start > seed_last->read_end) {
		  ok = 1;
		} else if (seed_last->read_end - seed_first_next->read_start <= 5) {
		  ok = 1;
		}

		if (ok &&
		    cal->start <= (cal_prev->end + MAX_INTRON_SIZE) &&
		    //abs(seed_first_next->read_start - seed_last->read_end) <= 5 && 
		    seed_first_next->id != seed_last->id) {
		  //printf("Yes\n");
		  //Merge CALs! Possible Splice. Insert CAL in targets.
		  merge_cals = array_list_new(5, 1.25f,
					      COLLECTION_MODE_ASYNCHRONIZED);
		    
		  array_list_insert(cal_prev, merge_cals);
		  array_list_insert(cal, merge_cals);		    
		  array_list_insert(merge_cals, target_cals);
		  
		  cal_prev = linked_list_iterator_next(&itr);
		  continue;
		}
	      }

	      if (seed_last->read_end - seed_first->read_start >= MIN_COVER) { 
		//Insert CAL in targets. 
		merge_cals = array_list_new(5, 1.25f,
					    COLLECTION_MODE_ASYNCHRONIZED);
		array_list_insert(cal_prev, merge_cals);
		array_list_insert(merge_cals, target_cals);
	      }	  
	      cal_prev = cal;	  
	    }

	  }
	}
      }      
    } 

    
    //(Filter CALs) INSERT COMMENT REGION HERE
    int num_targets_new = array_list_size(target_cals);
    array_list_t *list_aux = array_list_new(num_targets_new, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    //int delete_cals[num_targets];
    for (int a = 0; a < num_targets_new; a++) {
      merge_cals = array_list_get(a, target_cals);
      
      cal_t *cal_first = array_list_get(0, merge_cals);
      cal_t *cal_last  = array_list_get(array_list_size(merge_cals) - 1, merge_cals);
      
      seed_region_t *seed_first = linked_list_get_first(cal_first->sr_list);
      seed_region_t *seed_last  = linked_list_get_last(cal_last->sr_list);

      //cal_print(cal_first);
      //cal_print(cal_last);
      if (seed_first->read_start == 0 && 
	  seed_last->read_end    == read->length - 1) {
	array_list_insert(merge_cals, list_aux);
	//printf("Insert CAL-> (%i)\n", array_list_size(merge_cals));
      } else {
	//delete_cals[0] = 1;
	array_list_free(merge_cals, (void *)NULL);
      }
    }

    array_list_clear(target_cals, (void *)NULL); 
    //printf("list aux items = %i\n", array_list_size(list_aux));

    for (int a = 0; a < array_list_size(list_aux); a++) {
      merge_cals = array_list_get(a, list_aux);
      array_list_insert(merge_cals, target_cals);
    }
    array_list_free(list_aux, (void *)NULL);
    */

    int n_report;
    if (!array_list_size(alignments_list) &&
	!array_list_size(target_cals)) { 
      //Read no map! Seeding... :/

      //printf("Read Seeding...\n");
      
      for (int a = 0; a < array_list_size(target_cals); a++) {
	merge_cals = array_list_get(a, target_cals);
	array_list_free(merge_cals, (void *)NULL);
      }
      array_list_clear(target_cals, (void *)NULL);
      
      array_list_set_flag(1, sa_alignments_list[r]);
      
      //cal_mng_simple_clear(cal_mng_p);//(+)
      //cal_mng_simple_clear(cal_mng_n);//(-)
      
      int seed_size = sa_index->k_value;
      int seed_inc  = seed_size / 2;
      int read_pos;
      int id_seed   = 2;
      cal_mng_t *cal_mng;
      size_t high, low, suffix_len;
      float score;
      int cal_pos;
      seed_region_t *s_prev, *s;
      linked_list_item_t *item;
	  
      //===== Seeding Strategy =====//
      int max_read_pos = len_seq - seed_size - 1;
      for (int s = 0; s < 2; s++) { //Strand
	read_pos  = seed_size;
	if (!s) {
	  query = seq;
	  cal_mng = cal_mng_p;
	} else {
	  query = seq_revcomp;
	  cal_mng = cal_mng_n;
	}
	    
	while (read_pos <= max_read_pos) {
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
			    cal_mng->cals_lists[chrom], len_seq, id_seed++);
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
	}//loop seeds	
      } //loop strands      
      //===== Seeding Strategy End =====//


      printf("==================== CALs result Seeding =====================\n");
      printf("==(+)==\n");
      cal_mng_print(cal_mng_p);
      printf("==(-)==\n");
      cal_mng_print(cal_mng_n);
      printf("==============================================================\n");


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
		  //abs(s->read_start - s_prev->read_end) <= 5 &&
		  s_prev->id != s->id &&
		  cal_next->start <= (cal_prev->end + MAX_INTRON_SIZE)) {
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
    } //End Seeding, Caling
    
    
    //int num_targets;
    int num_targets = array_list_size(target_cals);	
    const int limit_cals = 200;
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
      cals_score[i] = generate_cals_merge_score(merge_cals, len_seq);
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


    //if (num_targets > 5) { 
    //num_targets = 0;
    //}
    
    //extern pthread_mutex_t mutex_sp;
    //pthread_mutex_lock(&mutex_sp);
    //extern int tot_cals[50];
    //tot_cals[num_targets] += 1;
    //pthread_mutex_unlock(&mutex_sp);    
        
    //printf("Total CALs %i\n", num_targets);
    //Report the first best 10 CALs
    n_report = num_targets < 10 ? num_targets : 10;

    /*
    if (n_report > 5) { 
      extern pthread_mutex_t mutex_sp;
      //pthread_mutex_lock(&mutex_sp);
      //printf("%s\n", read->id);
      for (int z = 0; z < num_targets; z++) {
	merge_cals = array_list_get(z, target_cals);
	cal_t *cal = array_list_get(0, merge_cals);
	//printf("%i. [%i:%lu] score : %f\n", z, cal->chromosome_id, cal->start, cals_score[z]);
      }
      n_report = 0; 
      //exit(-1);
      //pthread_mutex_unlock(&mutex_sp);
    }
    */    
    
    //Free other CALs
    for (int i = n_report; i < num_targets; i++) {
      merge_cals = array_list_get(i, target_cals);
      array_list_free(merge_cals, (void *)NULL);
    }
    
    
    int cal_len_prev, cal_len_next, cal_len;
    int fill_type;
    seed_region_t *s_prev_first, *s_prev_last, *s_next_first, *s_next_last, *s;
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
	  
	  //int nt = search_simple_splice_junction_semi_cannonical(seed_first, seed_last,
	  //							   cal_next->chromosome_id,
	  //							   cal_next->strand,
	  //							   query, genome, 
	  //							   &sp_start,
	  //							   &sp_end,
	  //							   &sp_type,
	  //							   &distance);
	  
	  int nt = search_simple_splice_junction(seed_first, seed_last,
						 cal_next->chromosome_id,
						 cal_next->strand,
						 query, genome, 
						 &sp_start,
						 &sp_end,
						 &sp_type,
						 &distance);
	  
	  //printf("Search %s:\n", read->id);
	  //cal_print(cal_prev);
	  //cal_print(cal_next);

	  //printf("nt = %i\n", nt);

	  if (nt > MIN_INTRON_SIZE && 
	      sp_start < sp_end && 
	      seed_first->genome_start < sp_start && 
	      sp_end < seed_last->genome_end) {

	    sa_alignment->sp_middle[sa_alignment->num_sp++] = nt;
	    cal_prev->end   = seed_first->genome_end;
	    cal_next->start = seed_last->genome_start;
	    //printf("(Start actualization) FOUND SPLICE CAL: (%lu-%lu)%i, %i\n", sp_start, sp_end, nt, sa_alignment->num_sp);

	    allocate_start_node(cal_prev->chromosome_id - 1,
				cal_prev->strand,
				sp_start,
				sp_end,
				sp_start,
				sp_end,
				FROM_READ,
				sp_type,
				NULL, 
				&node_avl_start,
				&node_avl_end, 
				avls_list);		  

	    //printf();

	    metaexon_insert(cal_prev->strand, cal_prev->chromosome_id - 1,
			    seed_first->genome_start, sp_start, 40,
			    METAEXON_RIGHT_END, node_avl_start,
			    metaexons);
	    
	    metaexon_insert(cal_prev->strand, cal_prev->chromosome_id - 1,
			    sp_end, seed_last->genome_end, 40,
			    METAEXON_LEFT_END, node_avl_end,
			    metaexons);
	    
	  } //else {
	    //printf(":( No Found\n");
	  //}
	    //extern pthread_mutex_t mutex_sp;
	    //pthread_mutex_lock(&mutex_sp);
	    //char q_sw[2024];
	    //char r_sw[2024];
	    //info_sp_t* info_sp = sw_reference_splice_junction(cal_prev, cal_next,
	  //						      query, genome,
	  //						      q_sw, r_sw);

	  //sw_insert_item(q_sw, r_sw, SJ_SW, sa_alignment, 
	  //		   sw_optarg, output, &sw_depth);

	    //printf("===========> SW FOR SEARCH INSERT <============\n\n");
	    //pthread_mutex_unlock(&mutex_sp);
	    //exit(-1);
	  //}
	  
	  //Metaexon actualization...
	  cal_prev = cal_next;
	  
	}

      }
    }

    /*
    for (int i = 0; i < n_report; i++) {
      //if (cals_score[i] == 0.0) { continue; }
      
      //sa_alignment_t *sa_alignment = array_list_get(sa_alignments_list[r]);
      merge_cals = array_list_get(i, target_cals);
      num_cals = array_list_size(merge_cals);
      
      if (num_cals > 10) { continue; }
      
      cal_prev     = array_list_get(0, merge_cals);
      s_prev       = cal_prev->sr_list->first->item; 
      //s_prev_last  = cal_prev->sr_list->last->item;
      //cal_len_prev = s_prev_last->read_end - s_prev_first->read_start;
      
      cal_next     = array_list_get(num_cals - 1, merge_cals);
      seed_region_t *s_next       = cal_next->sr_list->first->item; 
      //s_next_last  = cal_next->sr_list->last->item;
      //cal_len_next = s_next_last->read_end - s_next_first->read_start;
      
      //printf("SA ALIGNMENT... %i, %i:%lu\n", i, cal_prev->chromosome_id, cal_prev->start);
      sa_alignment_t *sa_alignment = sa_alignment_new(merge_cals);
      array_list_insert(sa_alignment, sa_alignments_list[r]);

      if (cal_prev->strand) {
	query = seq_revcomp;
      } else {
	query = seq;
      }
      
      //1-Close left and right gaps
      if (s_prev->read_start != 0) {
	sa_alignment->c_left = meta_alignment_fill_extrem_gap(query, 
							      cal_prev,
							      FILL_GAP_LEFT,
							      genome,
							      metaexons, 
							      avls_list);
	//printf("LEFT: (((((((((( %s )))))))))))\n", new_cigar_code_string(sa_alignment->c_left));	
	if (sa_alignment->c_left) {
	  sa_alignment->left_close = 1;
	}
      } else {
	sa_alignment->left_close = 1;	
      }
      
      if (s_next->read_end != len_seq - 1) {
	sa_alignment->c_right = meta_alignment_fill_extrem_gap(query, 
							       cal_next,
							       FILL_GAP_RIGHT,
							       genome,
							       metaexons, 
							       avls_list);
	if (sa_alignment->c_right) {
	  sa_alignment->right_close = 1;
	} 
      } else {
	sa_alignment->right_close = 1;
      }
      
      if (sa_alignment->left_close && sa_alignment->right_close) {
	sa_alignment->complete = 1;
      }
      
    }    
    //=====*---------------------*=====//
    */
	  
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
	  //linked_list_iterator_init(cal->sr_list, &iter);
	  //seed_next = linked_list_iterator_curr(&iter);
	  //while (seed_next != NULL) {
	  //if (seed_prev) {
	  //if ((seed_prev->read_end > seed_next->read_start) || 
	  //(seed_prev->genome_end > seed_next->genome_start)) {
	  //return 0;
	  //}
	  //}
	  //seed_prev = seed_next;
	  //seed_next = linked_list_iterator_next(&iter);
	  //}	  
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

	//printf("fill ..........................\n");
	//cal_print(cal);
	//printf("fill ..........................\n");
	cal_fill_gaps(cal,
		      query,
		      genome,
		      metaexons,
		      avls_list, 
		      MIN_INTRON_SIZE, 
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
      //printf("=========== Alignment %i: ============\n", n_sa);
      
      //printf("-------------------(%i cals)-------------------\n", array_list_size(merge_cals));
      //for (int x = 0; x < array_list_size(merge_cals); x++) {
      //cal = array_list_get(x, merge_cals);
      //cal_print(cal);
      //}
      //printf("-----------------------------------------------\n");
      
      //Generate cigars from sa_alignments seeding
      //if (array_list_get_flag(sa_alignments_list[r]) == 1) {
      //if (num_cals - 1 != sa_alignment->num_sp) {
	
      //}

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
	    cigar_code_free(cc);
	    sr->info = NULL;
	  } 
	  item = item->next;
	}

	//printf("%i < %i\n", n_sp , sa_alignment->num_sp);
	if (n_sp < sa_alignment->num_sp) {
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
	char cigar_str[2048] = "\0";
	for (int k = 0; k < sa_alignment->c_final->ops->size; k++) {
	  cigar_op_t *op = array_list_get(k, sa_alignment->c_final->ops);	      
	  sprintf(cigar_str, "%s%i%c", cigar_str, op->number, op->name);
	  //cigar_op_free(op);
	}
	  
	//printf("XXXXXXXXX-Cigar (%s): %s\n", read->id, cigar_str);
	  
	//Report Aligment
	cal = array_list_get(0, merge_cals);	      
	alignment_t *alignment = alignment_new();
	alignment_init_single_end(strdup(read->id), 
				  strdup(read->sequence), 
				  strdup(read->quality),
				  cal->strand, 
				  cal->chromosome_id,
				  cal->start - 1,				  
				  //strdup("100M"), 1, 254, 1, 0,
				  strdup(cigar_str), 
				  //new_cigar_code_string(sa_alignment->c_final->ops),
				  1, 254, 1, 0,
				  0, 0, alignment);
	  
	array_list_insert(alignment, alignments_list);
	sa_alignment->reported = 1;

      }      
    }



    if (array_list_size(alignments_list)) {
      //Report Other mappings
      num_sa_alignments = array_list_size(sa_alignments_list[r]);      
      for (int n_sa = 0; n_sa < num_sa_alignments; n_sa++) {
	sa_alignment_t *sa_alignment = array_list_get(n_sa, sa_alignments_list[r]);
	//Complete Cigar
	if  (!sa_alignment->reported && 
	     cals_score[n_sa] > 50) {
	  merge_cals = sa_alignment->cals_list;
	  cal = array_list_get(0, merge_cals);

	  seed_prev = linked_list_get_first(cal->sr_list);
	  seed_next = linked_list_get_last(cal->sr_list);
	  
	  //cigar_code_t *c_aux = sa_alignment->c_final;
	  if (seed_prev->read_start != 0) {
	    cigar_op_t *op_new = cigar_op_new(seed_prev->read_start, 'S');
	    array_list_insert_at(0, op_new, sa_alignment->c_final->ops);
	  }
	  
	  if (seed_next->read_end != read->length - 1) {
	    cigar_op_t *op_new = cigar_op_new(read->length - seed_next->read_end - 1, 'S');
	    array_list_insert(op_new, sa_alignment->c_final->ops);
	  }

	  if (cigar_code_validate_(read, sa_alignment->c_final)) {	     
	    char cigar_str[2048] = "\0";
	    for (int k = 0; k < sa_alignment->c_final->ops->size; k++) {
	      cigar_op_t *op = array_list_get(k, sa_alignment->c_final->ops);	      
	      sprintf(cigar_str, "%s%i%c", cigar_str, op->number, op->name);
	    }

	    //printf("XXXX-XXXXX-Cigar (%s): %s\n", read->id, cigar_str);

	    alignment_t *alignment = alignment_new();
	    alignment_init_single_end(strdup(read->id), 
				      strdup(read->sequence), 
				      strdup(read->quality),
				      cal->strand, 
				      cal->chromosome_id,
				      cal->start - 1,
				      strdup(cigar_str),
				      1, 254, 1, 0,
				      0, 0, alignment);
	  
	    array_list_insert(alignment, alignments_list);	  
	  }
	}
      }
    }


    //printf("::::::: %i && %i(for)\n", array_list_size(alignments_list), num_sa_alignments);
    array_list_clear(write_alignments, (void *)NULL);
    for (int n_sa = 0; n_sa < num_sa_alignments; n_sa++) {
      sa_alignment_t *sa_alignment = array_list_get(n_sa, sa_alignments_list[r]);
      if (cigar_code_get_num_ops(sa_alignment->c_final) == 0) { 
	continue;
      }
      array_list_insert(sa_alignment, write_alignments);
    }


    if (!array_list_size(alignments_list) && array_list_size(write_alignments)) {
      //Not Perfect Mappings, sa_alignments will be store in file
      //file_write_sa_alignments(read, sa_alignments_list[r], file1);      
      extern pthread_mutex_t mutex_sp;
      pthread_mutex_lock(&mutex_sp);
      /*
      printf("%s\n", read->id);
      for (int n_sa = 0; n_sa < num_sa_alignments; n_sa++) {
	sa_alignment_t *sa_alignment = array_list_get(n_sa, sa_alignments_list[r]);
	merge_cals = sa_alignment->cals_list;
	int num_cals = array_list_size(merge_cals);
	printf("=========== Alignment %i: ============\n", n_sa);
	
	printf("-------------------(%i cals)-------------------\n", array_list_size(merge_cals));
	for (int x = 0; x < array_list_size(merge_cals); x++) {
	  cal = array_list_get(x, merge_cals);
	  cal_print(cal);
	}
	printf("----------------------(Final Cigar: %s)-------------------------\n", new_cigar_code_string(sa_alignment->c_final));
      }
      */
      //exit(-1);
      //fprintf(stderr, "==========: NOP CIGAR\n");
      //printf("------------> %lu\n", array_list_size(write_alignments));
      sa_file_write_items(read, write_alignments, file1);
      pthread_mutex_unlock(&mutex_sp);
      //Delete reads from batch
      //pthread_mutex_lock(&mutex_sp);
      //printf("XXX-2: .Store new file\n");
      //pthread_mutex_unlock(&mutex_sp);
      
      delete_targets[r] = 1;
    } else if (!array_list_size(alignments_list)) {//array_list_size(sa_alignments_list[r]) && 
      //Write alignments to file
      delete_targets[r] = 1;
      extern pthread_mutex_t mutex_sp;
      pthread_mutex_lock(&mutex_sp);
      sa_file_write_items(read, write_alignments, file1);
      pthread_mutex_unlock(&mutex_sp);

      /*
      for (int n_sa = 0; n_sa < num_sa_alignments; n_sa++) {
	sa_alignment_t *sa_alignment = array_list_get(n_sa, sa_alignments_list[r]);
	cal = array_list_get(0, sa_alignment->cals_list);
	//Complete Cigar
	alignment_t *alignment = alignment_new();
	alignment_init_single_end(strdup(read->id), 
				  strdup(read->sequence), 
				  strdup(read->quality),
				  cal->strand, 
				  cal->chromosome_id,
				  cal->start - 1,
				  strdup("100M"),
				  1, 254, 1, 0,
				  0, 0, alignment);
	
	array_list_insert(alignment, alignments_list);	  
      }
      */
      
      //extern pthread_mutex_t mutex_sp;
      //pthread_mutex_lock(&mutex_sp);
      //printf("%s\n", read->id);
      //pthread_mutex_unlock(&mutex_sp);

      //printf("== %s ==\n", read->id);
      
      //pthread_mutex_lock(&mutex_sp);
      //fprintf(stderr, "XXX@@@XXX-3: %s\n", read->id);
      //pthread_mutex_unlock(&mutex_sp);
      
      //New seeding!...  Bad cigars
      //Search strand(+)
      /*
      num_suffixes_p = search_suffix(&seq[0], sa_index->k_value, 
				     MAX_NUM_SUFFIXES, sa_index, 
				     &low_p, &high_p, &suffix_len_p);
      
      //Search strand(-)
      num_suffixes_n = search_suffix(&seq_revcomp[0], sa_index->k_value, 
				     MAX_NUM_SUFFIXES, sa_index, 
				     &low_n, &high_n, &suffix_len_n);

      printf("%i | %i\n", suffix_len_p, suffix_len_n);

      printf("RESULTS (+):\n");
      if (suffix_len_p && num_suffixes_p) {
	//Report Exact Maps! (+)
	int limit_len = read->length * 0.95;
	printf("%i - %i\n", suffix_len_p, limit_len);
	if (suffix_len_p >= limit_len) {
	  printf("Aliggggggggggggggggggggggggggggggg\n");
	  for (size_t suff = low_p; suff <= high_p; suff++) {
	    chrom = sa_index->CHROM[suff];	
	    g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];
	    
	    char str[2048];
	    sprintf(str, "%i%c", read->length, 'M');
	    alignment_t *alignment = alignment_new();
	    alignment_init_single_end(strdup(read->id), 
				      strdup(read->sequence), 
				      strdup(read->quality),
				      0,
				      chrom + 1,
				      g_start,
				      strdup(str),
				      1, 254, 1, 0,
				      0, 0, alignment);
	  
	    array_list_insert(alignment, alignments_list);	  
	    
	  }
	} else {
	  int len_ref = strlen(&read->sequence[suffix_len_p]);
	  for (size_t suff = low_p; suff <= high_p; suff++) {
	    chrom = sa_index->CHROM[suff];	
	    g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];

	    char reference[2048];
	    size_t genome_start = g_start + suffix_len_p + 1;
	    size_t genome_end = genome_start + len_ref - 1;

	    genome_read_sequence_by_chr_index(reference, 0, chrom,
					      &genome_start, &genome_end, genome);

	    printf("REF: %s\n", reference);
	    printf("SEQ: %s\n", &read->sequence[suffix_len_p]);

	    //int rc = 0;
	    //int num_err = 0;
	    //for (int c = suffix_len_p; c < read->length; c++) {
	    //if (read->sequence[c] != reference[rc++]) {
	    //num_err++;
	    //}
	    //}
	  
	    alig_out_t out;
	    float score = doscadfun(&read->sequence[suffix_len_p], 
				    len_ref, 
				    reference,
				    len_ref,
				    0.1f, //10% Error
				    &out);

	    int max_mismatch = len_ref * 0.1;
	    printf("out.mismatch = %i, max_mismatch = %i\n", out.mismatch, max_mismatch);
	    if (out.st_map_len && 
		out.mismatch < max_mismatch + 1) {
	      int value;
	      int name;
	      int t;
	      cigar_get_op(0, &value, &name, &out.cigar);	    	    
	      char str[2048];
	    
	      if (name == 'M') {
		value += suffix_len_p;
		sprintf(str, "%i%c", value, 'M');
		t = 1;
	      } else {
		sprintf(str, "%i%c", suffix_len_p, 'M');
		t = 0;
	      }

	      for (; t < out.cigar.num_ops; t++) {
		cigar_get_op(t, &value, &name, &out.cigar);
		sprintf(str, "%s%i%c", str, value, name);
	      }
	    
	      printf("str first:%s\n", str);
	    
	      alignment_t *alignment = alignment_new();
	      alignment_init_single_end(strdup(read->id), 
					strdup(read->sequence), 
					strdup(read->quality),
					0,
					chrom + 1,
					g_start,
					strdup(str),
					1, 254, 1, 0,
					0, 0, alignment);
	  
	      array_list_insert(alignment, alignments_list);	  
	    
	    }
	  }
	}
	//printf("(%i) +: %i:%lu, len_cover = %i, err = %i, (cigar=%s, mismatch=%i)\n", suffix_len_p, chrom, g_start, read->length - suffix_len_p + 1, num_err, cigar_to_string(&out.cigar), out.mismatch);
      }	      

      printf("RESULTS (-):\n");
      if (suffix_len_n && num_suffixes_n) {
	int limit_len = read->length * 0.95;
	printf("%i - %i\n", suffix_len_n, limit_len);
	if (suffix_len_n >= limit_len) {
	  printf("Aliggggggggggggggggggggggggggggggg\n");
	  for (size_t suff = low_p; suff <= high_p; suff++) {
	    chrom = sa_index->CHROM[suff];	
	    g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];
	    
	    char str[2048];
	    sprintf(str, "%i%c", read->length, 'M');	    
	    alignment_t *alignment = alignment_new();
	    alignment_init_single_end(strdup(read->id), 
				      strdup(read->sequence), 
				      strdup(read->quality),
				      0,
				      chrom + 1,
				      g_start,
				      strdup(str),
				      1, 254, 1, 0,
				      0, 0, alignment);
	    
	    array_list_insert(alignment, alignments_list);	  

	  }	
	} else {
	  //Report Exact Maps! (+)
	  int len_ref = strlen(&read->sequence[suffix_len_n]);
	  for (size_t suff = low_n; suff <= high_n; suff++) {
	    chrom = sa_index->CHROM[suff];
	    g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];	  

	    char reference[2048];
	    size_t genome_start = g_start + suffix_len_n + 1;
	    size_t genome_end = genome_start + strlen(&read->sequence[suffix_len_n]) - 1;

	    genome_read_sequence_by_chr_index(reference, 0, chrom,
					      &genome_start, &genome_end, genome);

	    printf("REF: %s\n", reference);
	    printf("SEQ: %s\n", &read->revcomp[suffix_len_n]);
	  
	    //int rc = 0;
	    //int num_err = 0;
	    //for (int c = suffix_len_n; c < read->length; c++) {
	    //if (read->revcomp[c] != reference[rc++]) {
	    //num_err++;
	    //}
	    //}

	    alig_out_t out;
	    float score = doscadfun(&read->revcomp[suffix_len_n], 
				    len_ref, 
				    reference,
				    len_ref,
				    0.1f, //10% Error
				    &out);
	
	    int max_mismatch = len_ref * 0.1;
	    printf("out.mismatch = %i, max_mismatch = %i, sccore = %f, st_map_len = %i\n", out.mismatch, max_mismatch, out.score, out.st_map_len);

	    if (out.st_map_len && out.mismatch < max_mismatch + 1) {
	      int value;
	      int name;
	      int t;
	      cigar_get_op(0, &value, &name, &out.cigar);	    	    
	      char str[2048];
	    
	      if (name == 'M') {
		value += suffix_len_n;
		sprintf(str, "%i%c", value, 'M');
		t = 1;
	      } else {
		sprintf(str, "%i%c", suffix_len_n, 'M');
		t = 0;
	      }

	      for (; t < out.cigar.num_ops; t++) {
		cigar_get_op(t, &value, &name, &out.cigar);
		sprintf(str, "%s%i%c", str, value, name);
	      }
	    
	      printf("str first:%s\n", str);
	    
	      alignment_t *alignment = alignment_new();
	      alignment_init_single_end(strdup(read->id), 
					strdup(read->revcomp), 
					strdup(read->quality),
					1,
					chrom + 1,
					g_start,
					strdup(str),
					1, 254, 1, 0,
					0, 0, alignment);
	  
	      array_list_insert(alignment, alignments_list);	  
	    
	    }
	    //printf("(%i) -: %i:%lu, len_cover = %i, err = %i, len_ref=%i, (cigar=%s, mismatch=%i), score=%f\n", suffix_len_n, chrom, g_start, read->length - suffix_len_n + 1, num_err, len_ref, cigar_to_string(&out.cigar), out.mismatch, score);	  
	  }
	}
	}*/
    }

    //else if (!array_list_size(alignments_list) && ) {
      //Nothing found
      //extern pthread_mutex_t mutex_sp;
      //pthread_mutex_lock(&mutex_sp);
      //fprintf(stderr, "xxxxxxxxxx: NOP CIGAR\n");
      //pthread_mutex_unlock(&mutex_sp);
      //...
      //extern pthread_mutex_t mutex_sp;
      //pthread_mutex_lock(&mutex_sp);
      //fprintf(stderr, "@@@@ NOP CIGAR\n");
      //pthread_mutex_unlock(&mutex_sp);
      //} //else {
      //extern pthread_mutex_t mutex_sp;
      //pthread_mutex_lock(&mutex_sp);
      //fprintf(stderr, ":::::::::: XXXXXXXX\n");
      //pthread_mutex_unlock(&mutex_sp);
    //}

    num_sa_alignments = array_list_size(sa_alignments_list[r]);
    for (int n_sa = 0; n_sa < num_sa_alignments; n_sa++) {
      sa_alignment_t *sa_alignment = array_list_get(n_sa, sa_alignments_list[r]);

      merge_cals = sa_alignment->cals_list;
      int num_cals = array_list_size(merge_cals);
      
      cal_t *cal_aux = array_list_get(0, merge_cals);
      assert(cal_aux != NULL);
      //printf("(%p)Free sa_alignments %i/%i (%i:%lu:%lu)\n", cal_aux, n_sa + 1, num_sa_alignments, cal_aux->chromosome_id, cal_aux->start, cal_aux->end);
      
      //Free sa_alignments
      for (int k = 0; k < num_cals; k++) {
	cal = array_list_get(k, merge_cals);
	linked_list_item_t *item = cal->sr_list->first;
	seed_region_t *sr;	    
	while (item != NULL) {
	  sr = item->item;
	  cigar_code_t *cc = sr->info;
	  if (cc) {	      
	    //array_list_clear(cc->ops, (void *)cigar_op_free);
	    cigar_code_free(cc);		
	  }
	  item = item->next;
	}
      }
      
      if (sa_alignment->c_left) {
	//array_list_clear(sa_alignment->c_left->ops, (void *)cigar_op_free);
	cigar_code_free(sa_alignment->c_left);
      }

      if (sa_alignment->c_right) {
	//array_list_clear(sa_alignment->c_right->ops, (void *)cigar_op_free);
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

  } //End for reads
  

  array_list_t *fq_reads_aux       = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t **mapping_lists_aux = (array_list_t **)malloc(sizeof(array_list_t *)*num_reads);
  array_list_t *alignments_list_aux;
  int n_reads = 0;
  for (int r = 0; r < num_reads; r++) {
    read = array_list_get(r, sa_batch->fq_reads);
    alignments_list_aux = sa_batch->mapping_lists[r];
    if (delete_targets[r]) {
      //printf("%i/%i FREE!\n", r, num_reads);
      fastq_read_free(read);
      array_list_free(alignments_list_aux, (void *)NULL);
    } else {
      //printf("%i/%i INSERT! (%i)\n", r, num_reads, array_list_size(alignments_list_aux));
      array_list_insert(read, fq_reads_aux);
      mapping_lists_aux[n_reads++] = alignments_list_aux;
    }
  }

  
  array_list_free(write_alignments, (void *)NULL);
  array_list_free(sa_batch->fq_reads, (void *)NULL);
  free(sa_batch->mapping_lists);

  sa_batch->fq_reads = fq_reads_aux;
  sa_batch->mapping_lists = mapping_lists_aux;
  sa_batch->num_reads = array_list_size(fq_reads_aux);


  array_list_free(target_cals, (void *)NULL);
  
  cal_mng_simple_free(cal_mng_p);
  cal_mng_simple_free(cal_mng_n);

  sw_multi_output_free(output);

  /*
  printf("/////////////////////////////////////////////////////////\n");

  num_reads = array_list_size(sa_batch->fq_reads);
  for (size_t i = 0; i < num_reads; i++) {
    read = (fastq_read_t *) array_list_get(i, sa_batch->fq_reads);
    array_list_t *mapping_list = sa_batch->mapping_lists[i];
    int num_mappings = array_list_size(mapping_list);
    printf("%i.Read %s (num_mappings %i)\n", i, read->id, num_mappings);
  }

  exit(-1);
  */

  return -1;
  
}

int sa_rna_mapper_last(void *data) {
  array_list_t *sa_list;
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;
  sa_batch_t *sa_batch = wf_batch->mapping_batch;
  size_t num_reads = sa_batch->num_reads;
  /*
  for (int r = 0; r < num_reads; r++) {
    fastq_read_t *read = array_list_get(r, sa_batch->fq_reads);
    sa_list = sa_batch->mapping_lists[r];
    array_list_clear(sa_list, (void *)NULL);
  }

  return -1;
  */
  sa_index3_t *sa_index = (sa_index3_t *) wf_batch->sa_index;
  sa_rna_input_t *sa_rna = wf_batch->data;
  genome_t *genome = (genome_t *)sa_rna->genome;
  avls_list_t *avls_list = (avls_list_t *)sa_rna->avls_list;
  metaexons_t *metaexons = (metaexons_t *)sa_rna->metaexons;
  sw_optarg_t *sw_optarg = sa_rna->sw_optarg;
  FILE *file1 = sa_rna->file1;
  FILE *file2 = sa_rna->file2;
  char *query;
  //size_t num_reads = sa_batch->num_reads;
  //array_list_t *sa_list;
  int min_score = 80;
  array_list_t *alignments_list;
  sw_multi_output_t *output = sw_multi_output_new(MAX_DEPTH);
  sa_sw_depth_t sw_depth;
  extern pthread_mutex_t mutex_sp;
  
  sw_depth.depth = 0;

  //metaexons_print(metaexons);

  //printf("||======= SECOND WORKFLOW =========||\n");
  
  for (int r = 0; r < num_reads; r++) {
    fastq_read_t *read = array_list_get(r, sa_batch->fq_reads);
    sa_list = sa_batch->mapping_lists[r];
    //if (!array_list_size(sa_list)) { continue; }

    //if (strcmp("seq.5a", read->id) == 0) {       
    //exit(-1);
    //}

    size_t num_items = array_list_size(sa_list);
    
    //extern pthread_mutex_t mutex_sp;
    //pthread_mutex_lock(&mutex_sp);
    printf("@@@@@@ READ SECOND (%i): %s\n", num_items, read->id);
    //pthread_mutex_unlock(&mutex_sp);
    
    for (int i = 0; i < num_items; i++) {
      sa_alignment_t *sa_alignment = array_list_get(i, sa_list);
      cal_t *cal = array_list_get(0, sa_alignment->cals_list);
      cal_print(cal);
      seed_region_t *region = linked_list_get_first(cal->sr_list);
      
      assert(region);
      region->info = sa_alignment->c_final;
      
      printf("\tSEED: %i:[%lu|%i - %i|%lu] = %s\n", cal->chromosome_id, region->genome_start, region->read_start, region->read_end, region->genome_end, new_cigar_code_string(sa_alignment->c_final));

      if (cal->strand) {
	query = read->revcomp;
      } else {
	query = read->sequence;
      }
      
      cigar_code_t *c_left, *c_right;
      //Close left gap...
      if (region->read_start != 0) {
	sa_alignment->left_close = 1;

	c_left = fill_extrem_gap(query, 
				 cal,
				 FILL_GAP_LEFT,
				 genome,
				 metaexons, 
				 avls_list);

	seed_region_t *seed_region_l = seed_region_new(0, region->read_start - 1, 
						       region->genome_start - 1 - region->read_start,
						       region->genome_start - 1, 0, 0, 0);
	linked_list_insert_first(seed_region_l, cal->sr_list);

	seed_region_l->info = c_left;

	//printf("seed_region_l = %s\n", new_cigar_code_string(seed_region_l->info));

	if (seed_region_l->info == NULL) {
	  //SW...
	  char reference_sw[2048];
	  size_t genome_start = region->genome_start - region->read_start;
	  size_t genome_end   = region->genome_start - 1;

	  genome_read_sequence_by_chr_index(reference_sw, 0, cal->chromosome_id - 1,
					    &genome_start, &genome_end, genome);

	  char query_sw[2048];
	  memcpy(query_sw, &query[0], region->read_start);
	  query_sw[region->read_start] = '\0';
	  
	  //printf("(%i)l-Query: %s\n", query_sw, strlen(query_sw));
	  // printf("(%i)l-Refer: %s\n", reference_sw, strlen(reference_sw));

	  sw_insert_item(query_sw, reference_sw, FIRST_SW, seed_region_l,
			 sw_optarg, output, &sw_depth);
	  
	}
      }
      
      //Close right gap...
      if (region->read_end != read->length - 1) {	  
	int gap_len = read->length - region->read_end - 1;

	c_right = fill_extrem_gap(query, 
				  cal,
				  FILL_GAP_RIGHT,
				  genome,
				  metaexons, 
				  avls_list);

	seed_region_t *seed_region_r = seed_region_new(region->read_end + 1, region->read_end + 1 + gap_len, 
						     region->genome_end + 1,
						     region->genome_end + 1 + gap_len, 0, 0, 0);

	seed_region_r->info = c_right;

	linked_list_insert_last(seed_region_r, cal->sr_list);

	//printf("seed_region_r = %s\n", new_cigar_code_string(seed_region_r->info));

	if (seed_region_r->info == NULL) {
	  //SW...
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
      //printf("\tCAL %i/%i = %i:%lu|%i-%i|%lu, (%s)\n", i + 1, num_items, cal->chromosome_id, region->genome_start, 
      //region->read_start, region->read_end, region->genome_end, new_cigar_code_string(sa_alignment->c_final));      
    }
    //array_list_clear(sa_list, (void *)NULL);
  }

  
  sw_insert_item(NULL, NULL, 0, NULL, 
		 sw_optarg, output, &sw_depth);
  
  array_list_t *alignments_list_aux;
  for (int r = 0; r < num_reads; r++) {
    fastq_read_t *read = array_list_get(r, sa_batch->fq_reads);
    sa_list = sa_batch->mapping_lists[r];
    
    size_t num_items = array_list_size(sa_list);
    int list_size = num_items == 0 ? 100 : num_items;

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
      
      if (sa_alignment->left_close) {
	linked_list_item_t *item = cal->sr_list->first;
	if (item) {
	  seed_region_t *sr = item->item;
	  cigar_code_t *cc_aux = sr->info;
	  if (cc_aux) {
	    int dsp = 0;	    
	    for (int tt = 0; tt < cc_aux->ops->size; tt++) {
	      cigar_op_t *op = array_list_get(tt, cc_aux->ops);
	      dsp += op->number;
	    }
	    
	    //cal = array_list_get(0, merge_cals);	      
	    //seed_region_t *seed_aux = linked_list_get_first(cal->sr_list);
	    sr->genome_start -= dsp;//seed_aux->read_start;
	    cal->start -= dsp;//seed_aux->read_start;
	    sr->read_start = 0;	
	    //printf("dsp = %i, cal->start = %i \n", dsp, cal->start);
	  }
	}
      }
      
      cigar_code_t *c_final = cigar_code_new();
      cigar_code_t *c_aux = NULL;
      linked_list_item_t *item = cal->sr_list->first;
      while (item != NULL) {
	seed_region_t *seed_aux = item->item;
	c_aux = seed_aux->info;
	for (int c = 0; c < c_aux->ops->size; c++) {
	  cigar_op_t *op = array_list_get(c, c_aux->ops);
	  cigar_code_append_op(op, c_final);
	}
	cigar_code_free(c_aux);	
	seed_aux->info = NULL;
	item = item->next;
      }

      sa_alignment->c_final = c_final;

      int read_score = cigar_code_score(c_final, read->length);  

      //if (!cigar_code_validate_(read, sa_alignment->c_final)) {
      //extern pthread_mutex_t mutex_sp;
      //pthread_mutex_lock(&mutex_sp);
      //printf("%s : %s (%i score)\n", read->id, new_cigar_code_string(sa_alignment->c_final), read_score);	
      //pthread_mutex_unlock(&mutex_sp);
      //}

      if (cigar_code_validate_(read, sa_alignment->c_final) && 
	  read_score > min_score) {
	//Cigar generator, generate score
	//printf("Data: \n");
	//printf("\trl = %i - num_ops = %i\n", read->length, sa_alignment->c_final->ops->size );
	char cigar_str[2048] = "\0";
	for (int k = 0; k < sa_alignment->c_final->ops->size; k++) {
	  cigar_op_t *op = array_list_get(k, sa_alignment->c_final->ops);	      
	  if (op->name == 'H') { op->name = 'S'; }
	  sprintf(cigar_str, "%s%i%c", cigar_str, op->number, op->name);
	  //cigar_op_free(op);
	}
	  
	//printf("Cigar: %s\n", cigar_str);

	//Report Aligment
	printf("***************** FINAL CIGAR %s (%i)*****************\n ", cigar_str, read_score); 
	alignment_t *alignment = alignment_new();
	alignment_init_single_end(strdup(read->id), 
				  strdup(read->sequence), 
				  strdup(read->quality),
				  cal->strand, 
				  cal->chromosome_id,
				  cal->start - 1,
				  strdup(cigar_str), 1, 254, 1, 0,
				  0, 0, alignment);
	  
	array_list_insert(alignment, alignments_list_aux);
	  
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

      array_list_clear(merge_cals, (void *)cal_free);

      //if (sa_alignment->c_left) {
      //array_list_clear(sa_alignment->c_left->ops, (void *)cigar_op_free);	
      //cigar_code_free(sa_alignment->c_left);	
      //}      

      //if (sa_alignment->c_right) {
	//array_list_clear(sa_alignment->c_right->ops, (void *)cigar_op_free);
	//cigar_code_free(sa_alignment->c_right);
      //}      

      if (sa_alignment->c_final) {
	array_list_clear(sa_alignment->c_final->ops, (void *)cigar_op_free);
	cigar_code_free(sa_alignment->c_final);	
      }
      
      sa_alignment_free(sa_alignment);

      //--------------------------------------------------------------------      

    }//End items
    //}
      
    size_t num_suffixes_p, num_suffixes_n;
    size_t suffix_len_p, suffix_len_n;
    size_t low_p, high_p, low_n, high_n;
    int chrom;
    size_t g_start;
    char *seq = read->sequence;
    char *seq_revcomp = read->revcomp;
    cigar_code_t *cc_aux;
    
    if (!array_list_size(alignments_list_aux)) {
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
	      alignment_t *alignment = alignment_new();
	      alignment_init_single_end(strdup(read->id), 
					strdup(read->sequence), 
					strdup(read->quality),
					cal_tmp->strand, 
					cal_tmp->chromosome_id,
					cal_tmp->start,
					new_cigar_code_string(cc_aux),//strdup("100M"),
					1, 254, 1, 0,
					0, 0, alignment);
	    
	      array_list_insert(alignment, alignments_list_aux); 
	    } 
	    //}
	    //array_list_clear(cc_aux->ops, (void *)cigar_op_free);
	    //cigar_code_free(cc_aux);	    
	  }
	  //cal_simple_free(cal_tmp);
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
	      alignment_t *alignment = alignment_new();
	      alignment_init_single_end(strdup(read->id), 
					strdup(read->sequence), 
					strdup(read->quality),
					cal_tmp->strand, 
					cal_tmp->chromosome_id,
					cal_tmp->start,
					new_cigar_code_string(cc_aux),
					1, 254, 1, 0,
					0, 0, alignment);
	      
	      array_list_insert(alignment, alignments_list_aux);
	    }
	  }
	}
      }

    }

    /*
    if (!array_list_size(alignments_list_aux)) {
      extern pthread_mutex_t mutex_sp;
      pthread_mutex_lock(&mutex_sp);
      fprintf(stderr, "XXX@@@XXX-3: %s\n", read->id);

      cal_mng_t * cal_mng_p = cal_mng_new(sa_index->genome); //(+)
      cal_mng_t * cal_mng_n = cal_mng_new(sa_index->genome); //(-)
      cal_mng_t * cal_mng;

      int seed_size = 18;
      int seed_inc  = seed_size / 2;
      int read_pos;
      size_t num_suffixes, low, high, suffix_len, len_seq, id_seed;

      //===== Seeding Strategy =====//
      int max_read_pos = read->length - seed_size - 1;
      for (int s = 0; s < 2; s++) { //Strand
	read_pos  = seed_size;
	if (!s) {
	  query = read->sequence;
	  cal_mng = cal_mng_p;
	} else {
	  query = read->revcomp;
	  cal_mng = cal_mng_n;
	}
	    
	while (read_pos <= max_read_pos) {
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
			    cal_mng->cals_lists[chrom], len_seq, id_seed++);
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
	}//loop seeds	
      } //loop strands      
      //===== Seeding Strategy End =====//      

      printf("==================== CALs result Seeding =====================\n");
      printf("==(+)==\n");
      cal_mng_print(cal_mng_p);
      printf("==(-)==\n");
      cal_mng_print(cal_mng_n);
      printf("==============================================================\n");

      cal_mng_t *p;
      array_list_t  *merge_cals;//, *target_cals;
      linked_list_t *cal_list;
      int num_cals;
      int cal_pos;
      linked_list_item_t *item;
      seed_region_t *s_prev, *s;
      cal_t *cal_prev, *cal_next;
      float cals_score[2048];

      array_list_t *target_cals = array_list_new(100, 1.25f,
						 COLLECTION_MODE_ASYNCHRONIZED);
      
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
		  //abs(s->read_start - s_prev->read_end) <= 5 &&
		  s_prev->id != s->id &&
		  cal_next->start <= (cal_prev->end + MAX_INTRON_SIZE)) {
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
      } //End Seeding, Caling      	            
    
    
      //int num_targets;
      int num_targets = array_list_size(target_cals);	
      const int limit_cals = 200;
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
	cals_score[i] = generate_cals_merge_score(merge_cals, len_seq);
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
     

      // GENERATE ALIGNMENTS WITH METAEXON
      for (int i = 0; i < num_targets; i++) {
	merge_cals = array_list_get(i, target_cals);
	num_cals = array_list_size(merge_cals);

	cal_prev     = array_list_get(0, merge_cals);	
	for (int ii = 1; ii < num_cals; ii++) {
	  cal_next     = array_list_get(ii, merge_cals);
	  
	  //s_prev       = cal_prev->sr_list->first->item; 
	  
	  //seed_region_t *s_next       = cal_next->sr_list->first->item; 
	  //sa_alignment_t *sa_alignment = sa_alignment_new(merge_cals);
	  //array_list_insert(sa_alignment, sa_alignments_list[r]);
	  
	  //cal_t *selected_cal =NULL;
	  
	  //int cal_len_prev = cal_prev->end - cal_prev->start + 1;
	  //int cal_len_next = cal_next->end - cal_next->start + 1;
	  //selected_cal = cal_len_prev > cal_len_next ? cal_prev : cal_next;
	
	  if (cal_prev->strand) {
	    query = seq_revcomp;
	  } else {
	    query = seq;
	  }		  
	  
	  printf("========= SEARCH SJ =========\n");
	  cal_print(cal_prev);
	  cal_print(cal_next);

	  char q_sw[2024];
	  char r_sw[2024];
	  info_sp_t* info_sp = sw_reference_splice_junction(cal_prev, cal_next,
							    query, genome,
							    q_sw, r_sw);
	  
	  sw_insert_item(q_sw, r_sw, SJ_SW, NULL, 
			 sw_optarg, output, &sw_depth);
	  
	  printf("========= --------- =========\n");
	  
	  //Metaexon actualization...
	  cal_prev = cal_next;
	  
	}
      }
      
      sw_insert_item(NULL, NULL, 0, NULL, 
		     sw_optarg, output, &sw_depth);
 
      exit(-1);
      pthread_mutex_unlock(&mutex_sp);	      
    }
    */
    /*
    if (!array_list_size(alignments_list_aux)) {
      //extern pthread_mutex_t mutex_sp;
      //pthread_mutex_lock(&mutex_sp);
      //fprintf(stderr, "XXX@@@XXX-3: %s\n", read->id);
      //pthread_mutex_unlock(&mutex_sp);	
      //printf("------------------> %s <----------------\n", read->id);
      //extern pthread_mutex_t mutex_sp;
      //pthread_mutex_lock(&mutex_sp);
      //array_list_t *cals_list_aux = array_list_new(50, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);      
      //New seeding!...  Bad cigars
      //Search strand(+)
      num_suffixes_p = search_suffix(&seq[read->length - 19], sa_index->k_value, 
				     MAX_NUM_SUFFIXES, sa_index, 
				     &low_p, &high_p, &suffix_len_p);
      
      //Search strand(-)
      num_suffixes_n = search_suffix(&seq_revcomp[read->length - 19], sa_index->k_value, 
				     MAX_NUM_SUFFIXES, sa_index, 
				     &low_n, &high_n, &suffix_len_n);
      
      //printf("%i | %i\n", suffix_len_p, suffix_len_n);
      
      //printf("RESULTS (+):\n");
      if (suffix_len_p && num_suffixes_p) {
	//Report Exact Maps! (+)
	for (size_t suff = low_p; suff <= high_p; suff++) {
	  chrom = sa_index->CHROM[suff];	
	  g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];
	  //printf("\t%i:%lu\n", chrom, g_start);
	  
	  cal_t *cal_tmp = cal_simple_new(chrom + 1,
					  0, g_start, g_start + suffix_len_p);
	  
	  //cal_tmp->sr_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);	  
	  seed_region_t *seed_region_aux = seed_region_new(read->length - 19, read->length - 1,
							   g_start, g_start + suffix_len_p, 
							   0, 0, 0);
	  
	  linked_list_insert(seed_region_aux, cal_tmp->sr_list);
	  //array_list_insert(cal_tmp, cals_list_aux);
	  
	  cc_aux = fill_extrem_gap(seq, 
				   cal_tmp,
				   FILL_GAP_LEFT,
				   genome,
				   metaexons, 
				   avls_list);

	  if (cc_aux != NULL) {
	    cigar_op_t *c_op_aux = array_list_get(cc_aux->ops->size - 1, cc_aux->ops);

	    if (c_op_aux->name == 'M') { c_op_aux->number += suffix_len_p; }
	    else { array_list_insert(cigar_op_new(suffix_len_p, 'M'), cc_aux->ops); }

	    alignment_t *alignment = alignment_new();
	    alignment_init_single_end(strdup(read->id), 
				      strdup(read->sequence), 
				      strdup(read->quality),
				      cal_tmp->strand, 
				      cal_tmp->chromosome_id,
				      cal_tmp->start,
				      new_cigar_code_string(cc_aux),//strdup("100M"),
				      1, 254, 1, 0,
				      0, 0, alignment);
	
	    array_list_insert(alignment, alignments_list_aux);	    

	    //extern pthread_mutex_t mutex_sp;
	    //pthread_mutex_lock(&mutex_sp);
	    // fprintf(stderr, "=========::======= %s ==========::=======\n", new_cigar_code_string(cc_aux));
	    //pthread_mutex_unlock(&mutex_sp);
	   
	    //cal_print(cal_tmp);
	    //metaexons_print(metaexons);
	    //exit(-1);
	  }
	}
      }

      //printf("RESULTS (-):\n");
      if (suffix_len_n && num_suffixes_n) {
	for (size_t suff = low_n; suff <= high_n; suff++) {
	  chrom = sa_index->CHROM[suff];
	  g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];
	  //printf("\t%i:%lu\n", chrom, g_start);

	  cal_t *cal_tmp = cal_simple_new(chrom + 1,
					  1, g_start, g_start + suffix_len_n);

	  //cal_tmp->sr_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	  seed_region_t *seed_region_aux = seed_region_new(read->length - 19, read->length - 1,
							   g_start, g_start + suffix_len_n, 
							   0, 0, 0);
	  linked_list_insert(seed_region_aux, cal_tmp->sr_list);
	  //array_list_insert(cal_tmp, cals_list_aux);

	  cc_aux = fill_extrem_gap(seq_revcomp, 
				   cal_tmp,
				   FILL_GAP_LEFT,
				   genome,
				   metaexons, 
				   avls_list);

	  if (cc_aux != NULL) {

	    cigar_op_t *c_op_aux = array_list_get(cc_aux->ops->size - 1, cc_aux->ops);
	    if (c_op_aux->name == 'M') { c_op_aux->number += suffix_len_n; }
	    else { array_list_insert(cigar_op_new(suffix_len_n, 'M'), cc_aux->ops); }

	    //cal_print(cal_tmp);
	    //metaexons_print(metaexons);
	    //exit(-1);
	    alignment_t *alignment = alignment_new();
	    alignment_init_single_end(strdup(read->id), 
				      strdup(read->sequence), 
				      strdup(read->quality),
				      cal_tmp->strand, 
				      cal_tmp->chromosome_id,
				      cal_tmp->start,
				      new_cigar_code_string(cc_aux),//strdup("100M"),strdup("100M"),
				      1, 254, 1, 0,
				      0, 0, alignment);
	    
	    array_list_insert(alignment, alignments_list_aux);
	
	  }
	  
	  //extern pthread_mutex_t mutex_sp;
	  //pthread_mutex_lock(&mutex_sp);
	  //fprintf(stderr, "=========::======= %s ==========::=======\n", new_cigar_code_string(cc_aux));
	  //pthread_mutex_unlock(&mutex_sp);
	  
	}
      }
    }
    */
    
    array_list_free(sa_batch->mapping_lists[r], (void *)NULL);
    sa_batch->mapping_lists[r] = alignments_list_aux;
    
  } //End reads 

  sw_multi_output_free(output);
  
  return -1;

}
