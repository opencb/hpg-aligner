#include "bwt_server.h"

#define SIMPLE_SW 1
#define SP_SW 2

#define SW_NORMAL 0 
#define SW_FINAL 1 

//====================================================================================
// bwt_server_input functions: init
//====================================================================================

void bwt_server_input_init(list_t* read_list_p, unsigned int batch_size, bwt_optarg_t *bwt_optarg_p, 
			   bwt_index_t *bwt_index_p, list_t* write_list_p, unsigned int write_size, 
			   list_t* unmapped_read_list_p, metaexons_t *metaexons, sw_optarg_t *sw_optarg,
			   genome_t *genome, bwt_server_input_t* input_p) {
  input_p->read_list_p = read_list_p;
  input_p->batch_size = batch_size;
  input_p->bwt_optarg_p = bwt_optarg_p;
  input_p->write_list_p = write_list_p;
  input_p->write_size = write_size;
  input_p->bwt_index_p = bwt_index_p;
  input_p->unmapped_read_list_p = unmapped_read_list_p;
  input_p->metaexons = metaexons;
  input_p->sw_optarg = sw_optarg;
  input_p->genome = genome;
}

//====================================================================================
// apply_bwt
//====================================================================================

cal_t *convert_bwt_anchor_to_CAL(bwt_anchor_t *bwt_anchor, size_t read_start, size_t read_end) {
  linked_list_t *linked_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
  seed_region_t *seed_region = seed_region_new(read_start, read_end,
					       bwt_anchor->start, bwt_anchor->end, 0, 0, 0);

  linked_list_insert_first(seed_region, linked_list);

  cal_t *cal = cal_new(bwt_anchor->chromosome + 1, bwt_anchor->strand,
		       bwt_anchor->start, bwt_anchor->end,
		       1, linked_list,
		       linked_list_new(COLLECTION_MODE_ASYNCHRONIZED));

  return cal;

}

size_t bwt_search_pair_anchors(array_list_t *list, unsigned int read_length) {
  bwt_anchor_t *bwt_anchor;
  int max_anchor_length = 0;
  

  bwt_anchor_t *bwt_anchor_back, *bwt_anchor_forw;
  int anchor_length_tmp, anchor_back, anchor_forw;
  int strand = 0, type = 0;
  int found_anchor = 0, found_double_anchor = 0;

  const int MIN_ANCHOR = 25;
  const int MIN_SINGLE_ANCHOR = 40;

  //const int MIN_DOUBLE_ANCHOR = MIN_ANCHOR*2;
  const int MAX_BWT_REGIONS = 50;
  const int MAX_BWT_ANCHOR_DISTANCE = 500000;

  array_list_t *anchor_list_tmp, *forward_anchor_list, *backward_anchor_list;
  cal_t *cal;
  int seed_size, gap_read, gap_genome;

  array_list_t *backward_anchor_list_0 = array_list_new(MAX_BWT_REGIONS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *forward_anchor_list_0 = array_list_new(MAX_BWT_REGIONS, 1.25f , COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *backward_anchor_list_1 = array_list_new(MAX_BWT_REGIONS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *forward_anchor_list_1 = array_list_new(MAX_BWT_REGIONS, 1.25f , COLLECTION_MODE_ASYNCHRONIZED);

  array_list_t *big_anchor_list = array_list_new(MAX_BWT_REGIONS, 1.25f , COLLECTION_MODE_ASYNCHRONIZED);

  //printf("Tot Anchors %i\n", array_list_size(list));
  for (int i = 0; i < array_list_size(list); i++) {
    bwt_anchor = array_list_get(i, list);
    if (bwt_anchor->strand == 1) {
      //printf("(-)bwt anchor %i:%lu-%lu (%i): \n", bwt_anchor->chromosome + 1, bwt_anchor->start, bwt_anchor->end, bwt_anchor->end - bwt_anchor->start + 1);
      if (bwt_anchor->type == FORWARD_ANCHOR) {
	array_list_insert(bwt_anchor, forward_anchor_list_1);
	//printf("FORW\n");
      } else {
	array_list_insert(bwt_anchor, backward_anchor_list_1);
	//printf("BACK\n");
      }
    } else {
      //printf("(+)bwt anchor %i:%lu-%lu (%i): \n", bwt_anchor->chromosome + 1, bwt_anchor->start, bwt_anchor->end, bwt_anchor->end - bwt_anchor->start + 1);
      if (bwt_anchor->type == FORWARD_ANCHOR) {
	array_list_insert(bwt_anchor, forward_anchor_list_0);
	//printf("FORW\n");
      } else {
	array_list_insert(bwt_anchor, backward_anchor_list_0);
	//printf("BACK\n");
      }
    }

    anchor_length_tmp = bwt_anchor->end - bwt_anchor->start + 1;
    if (anchor_length_tmp > MIN_SINGLE_ANCHOR && anchor_length_tmp > max_anchor_length) {
      max_anchor_length = anchor_length_tmp;
      found_anchor = 1;
      strand = bwt_anchor->strand;
      type = bwt_anchor->type;
    }
    
    if (read_length - anchor_length_tmp < 16) {
      array_list_insert(bwt_anchor, big_anchor_list);
    } 
    
  }
  
  array_list_clear(list, NULL);

  if (array_list_size(big_anchor_list) > 0) {
    for (int i = array_list_size(big_anchor_list) - 1; i >= 0; i--) {
      //printf("Insert cal %i\n", i);
      bwt_anchor = array_list_remove_at(i, big_anchor_list);
      size_t seed_size = bwt_anchor->end - bwt_anchor->start;

      if (bwt_anchor->type == FORWARD_ANCHOR) {
	cal = convert_bwt_anchor_to_CAL(bwt_anchor, 0, seed_size);
      } else {
	cal = convert_bwt_anchor_to_CAL(bwt_anchor, read_length - seed_size - 1, read_length - 1);
      }
      
      array_list_insert(cal, list);
    }
    array_list_set_flag(SINGLE_ANCHORS, list);
    
    goto exit;
  }

  for (int type = 1; type >= 0; type--) {
    if (!type) {
      forward_anchor_list = forward_anchor_list_1;
      backward_anchor_list = backward_anchor_list_1;
      //printf("Strand (+): %i-%i\n", array_list_size(forward_anchor_list), array_list_size(backward_anchor_list));
    } else { 
      forward_anchor_list = forward_anchor_list_0;
      backward_anchor_list = backward_anchor_list_0;
      //printf("Strand (-): %i-%i\n", array_list_size(forward_anchor_list), array_list_size(backward_anchor_list));
    }

    int *set_forward  = (int *)calloc(array_list_size(forward_anchor_list),  sizeof(int));
    int *set_backward = (int *)calloc(array_list_size(backward_anchor_list), sizeof(int));

    //Associate Anchors (+)/(-)
    for (int i = 0; i < array_list_size(forward_anchor_list); i++) { 
      if (set_forward[i]) { continue; }
      bwt_anchor_forw = array_list_get(i, forward_anchor_list);
      for (int j = 0; j < array_list_size(backward_anchor_list); j++) { 
	if (set_backward[j]) { continue; }
	bwt_anchor_back = array_list_get(j, backward_anchor_list);
	anchor_forw = (bwt_anchor_forw->end - bwt_anchor_forw->start + 1);
	anchor_back = (bwt_anchor_back->end - bwt_anchor_back->start + 1); 

	anchor_length_tmp = anchor_forw + anchor_back;

	//printf("\tCommpare %i:%lu-%lu with %i:%lu-%lu\n", bwt_anchor_forw->chromosome + 1, 
	//     bwt_anchor_forw->start, bwt_anchor_forw->end, bwt_anchor_back->chromosome + 1, 
	//     bwt_anchor_back->start, bwt_anchor_back->end);
	if (bwt_anchor_forw->chromosome == bwt_anchor_back->chromosome &&
	    abs(bwt_anchor_back->start - bwt_anchor_forw->end) <= MAX_BWT_ANCHOR_DISTANCE && 
	    anchor_forw >= MIN_ANCHOR && anchor_back >= MIN_ANCHOR) {
	  
	  if (bwt_anchor_back->start < bwt_anchor_forw->end) { continue; }
	  
	  gap_read = read_length - (anchor_forw + anchor_back);
	  gap_genome = bwt_anchor_back->start - bwt_anchor_forw->end;

	  //printf("anchor_forw = %i, anchor_back = %i, gap_read = %i, gap_genome = %i\n",
	  //	 anchor_forw, anchor_back, gap_read, gap_genome);
	  	  
	  int apply_flank = 0;
	  if (gap_read < 2 || gap_genome < 2) {
	    int gap;
	    if (gap_read < 0 && gap_genome < 0) {
	      gap = abs(gap_read) > abs(gap_genome) ? abs(gap_read) : abs(gap_genome);
	    } else if (gap_read < 0) {
	      gap = abs(gap_read);
	    } else if (gap_genome < 0) {
	      gap = abs(gap_genome);
	    } else {
	      gap = 2;
	    }
	    
	    int flank  = 5;
	    apply_flank = 1;
	    
	    if (abs(gap) >= flank*2) {
	      //Solve read overlap
	      flank = abs(gap)/2 + flank/2;
	    }
	    //printf("\tgap = %i, flank = %i\n", gap, flank);
	    if (flank >= anchor_forw) {
	      bwt_anchor_forw->end -= anchor_forw/2;	      
	    } else {
	      bwt_anchor_forw->end -= flank;
	    }

	    if (flank >= anchor_back) {
	      bwt_anchor_back->start += anchor_back/2;	    
	    } else {
	      bwt_anchor_back->start += flank;
	    }
	  } 
	  	  
	  cal = convert_bwt_anchor_to_CAL(bwt_anchor_forw, 0, bwt_anchor_forw->end - bwt_anchor_forw->start);
	  //printf("INSERT-1 (%i)[%i:%lu-%lu]\n", cal->strand, cal->chromosome_id, cal->start, cal->end);
	  array_list_insert(cal, list);
	  seed_size = bwt_anchor_back->end - bwt_anchor_back->start + 1;
	  //if (bwt_anchor_forw->end + read_length >= bwt_anchor_back->start) {	    
	  //seed_region_t *seed_region = seed_region_new(read_length - seed_size, read_length - 1,
	  //bwt_anchor_back->start, bwt_anchor_back->end, 1);
	  //cal->end = bwt_anchor_back->end;
	  //linked_list_insert_last(seed_region, cal->sr_list);	
	  //} else {
	  cal = convert_bwt_anchor_to_CAL(bwt_anchor_back, read_length - seed_size, read_length - 1);
	  //printf("INSERT-2 (%i)[%i:%lu-%lu]\n", cal->strand, cal->chromosome_id, cal->start, cal->end);
	  array_list_insert(cal, list);
	  if (array_list_size(list) > 5) { 
	    free(set_backward);
	    free(set_forward);	    
	    goto exit;
	  }

	  array_list_set_flag(DOUBLE_ANCHORS, list);
	  found_double_anchor = 1;
	  set_forward[i]  = 1;
	  set_backward[j] = 1;
	  break;
	}                                                                                                                      
      }         
    }
    free(set_backward);
    free(set_forward);
  }

  if (!found_double_anchor && found_anchor) { 
    //Not Double anchor found but one Yes!!
    if (strand == 1) {
      if (type == FORWARD_ANCHOR) {
	anchor_list_tmp = forward_anchor_list_1;
      } else {
	anchor_list_tmp =  backward_anchor_list_1;
      }
    } else {
      if (type == FORWARD_ANCHOR) {
	anchor_list_tmp =  forward_anchor_list_0;
      } else {
	anchor_list_tmp =  backward_anchor_list_0;
      }
    }

    //printf("LIST SIZE %i\n", array_list_size(anchor_list_tmp));
    for (int i = 0; i < array_list_size(anchor_list_tmp); i++) {
      bwt_anchor = array_list_get(i, anchor_list_tmp);
      size_t seed_size = bwt_anchor->end - bwt_anchor->start;
      //array_list_insert(bwt_anchor_new(bwt_anchor->strand, bwt_anchor->chromosome, 
      //			       bwt_anchor->start, bwt_anchor->end, bwt_anchor->type), anchor_list);
      if (bwt_anchor->type == FORWARD_ANCHOR) {
	//printf("------------------------> start %i\n", 0);
	cal = convert_bwt_anchor_to_CAL(bwt_anchor, 0, seed_size);
      } else {
	//printf("------------------------> start %i\n", read_length - seed_size);
	cal = convert_bwt_anchor_to_CAL(bwt_anchor, read_length - seed_size - 1, read_length - 1);
      }
      array_list_insert(cal, list);
    }
    array_list_set_flag(SINGLE_ANCHORS, list);
  } 

 exit:
  array_list_free(forward_anchor_list_1, (void *)bwt_anchor_free);
  array_list_free(backward_anchor_list_1,  (void *)bwt_anchor_free);
  array_list_free(forward_anchor_list_0,  (void *)bwt_anchor_free);
  array_list_free(backward_anchor_list_0,  (void *)bwt_anchor_free);
  array_list_free(big_anchor_list,  (void *)bwt_anchor_free);

  return array_list_size(list);
  
}

//------------------------------------------------------------------------------------
// DNA
//------------------------------------------------------------------------------------

int apply_bwt(bwt_server_input_t* input, batch_t *batch) {
  //printf("APPLY BWT SERVER...\n");



  mapping_batch_t *mapping_batch = batch->mapping_batch;
  size_t num_reads = array_list_size(mapping_batch->fq_batch);
  size_t num_mappings;
  array_list_t *list;
  size_t *unmapped_indices = mapping_batch->targets;
  size_t num_unmapped = 0;
  size_t num_anchors;

  for (int i = 0; i < num_reads; i++) {
    fastq_read_t *read = array_list_get(i, mapping_batch->fq_batch);
    //printf("BWT: %s\n", read->id);
    list = mapping_batch->mapping_lists[i];    
    array_list_set_flag(0, list);
    num_mappings = bwt_map_inexact_read(read,
					input->bwt_optarg_p,
					input->bwt_index_p,
					list);
    if (array_list_get_flag(list) != 2) { //If flag 2, the read exceded the max number of mappings
      if (array_list_get_flag(list) == 1) {
	if (num_mappings > 0) {
	  num_anchors = bwt_search_pair_anchors(list, read->length);	
	  if (num_anchors == 0) {
	    array_list_set_flag(NOT_ANCHORS, list);
	  } 
	} else {
	  array_list_set_flag(NOT_ANCHORS, list);
	}
	//printf("tot anchors found %i %s\n", num_anchors, read->id);
	unmapped_indices[num_unmapped++] = i;
      } else if (num_mappings <= 0) {
	array_list_set_flag(0, list);
	//printf("Read NO Mapped %i %s\n", num_anchors, read->id);
	unmapped_indices[num_unmapped++] = i;
      }
    } else {
      array_list_set_flag(ALIGNMENTS_EXCEEDED, list);
    }
  }
  // array_list flag: 0 -> Not  BWT Anchors found
  //                  1 -> One  BWT Anchors found
  //                  2 -> Pair BWT Anchors found
  //                  3 -> Alignments found
  //                  4 -> Alignments exceded
  
  mapping_batch->num_targets = num_unmapped;
    
  if (batch->mapping_batch->num_targets > 0) {
    return CAL_STAGE;
  }

  return DNA_POST_PAIR_STAGE;
}

//------------------------------------------------------------------------------------
// RNA
//------------------------------------------------------------------------------------

void fastq_read_revcomp(fastq_read_t *read);

int apply_bwt_rna(bwt_server_input_t* input, batch_t *batch) {

  LOG_DEBUG("========= APPLY BWT RNA =========\n");



  metaexons_t *metaexons = input->metaexons;
  mapping_batch_t *mapping_batch = batch->mapping_batch;
  size_t num_reads = array_list_size(mapping_batch->fq_batch);
  size_t num_mappings;
  array_list_t *list;
  size_t *unmapped_indices = mapping_batch->targets;
  size_t num_unmapped = 0;
  size_t num_anchors;

  extern pthread_mutex_t mutex_sp;
  extern st_bwt_t st_bwt;
  //pthread_mutex_lock(&mutex_sp);
  //extern size_t total_reads;
  //total_reads += num_reads;
  //pthread_mutex_unlock(&mutex_sp);
  
  for (int i = 0; i < num_reads; i++) {
    fastq_read_t *read = array_list_get(i, mapping_batch->fq_batch);

    //Rev-comp
    fastq_read_revcomp(read);

    //printf("BWT: %s\n", read->id);
    list = mapping_batch->mapping_lists[i];    
    
    array_list_set_flag(1, list);

    num_mappings = bwt_map_inexact_read(read,
					input->bwt_optarg_p,
					input->bwt_index_p,
					list);

    if (array_list_get_flag(list) != 2) { //If flag 2, the read exceded the max number of mappings
      if (array_list_get_flag(list) == 1) {
	if (num_mappings > 0) {
	  num_anchors = bwt_search_pair_anchors(list, read->length);	
	  if (num_anchors == 0) {
	    array_list_set_flag(NOT_ANCHORS, list);
	    unmapped_indices[num_unmapped++] = i;
	  }
	} else {
	  array_list_set_flag(NOT_ANCHORS, list);
	  unmapped_indices[num_unmapped++] = i;
	} 
	//printf("tot anchors found %i %s\n", num_anchors, read->id);
      } else if (num_mappings <= 0) {
	array_list_set_flag(0, list);
	//printf("Read NO Mapped %i %s\n", num_anchors, read->id);
	unmapped_indices[num_unmapped++] = i;
      } else {
	//Read Map, Metaexon Actualization
	array_list_set_flag(ALIGNMENTS_FOUND, list);

	pthread_mutex_lock(&mutex_sp);
	st_bwt.map_bwt++;
	pthread_mutex_unlock(&mutex_sp);

	for (int i = 0; i < num_mappings; i++) {
	  alignment_t *alignment = array_list_get(i, list);
	  metaexon_insert(0, alignment->chromosome,
			  alignment->position, alignment->position + read->length, 40,
			  METAEXON_NORMAL, NULL,
			  metaexons);
	  //alignment->alig_data = cigar_code_new_by_string(alignment->cigar);
	}
      }
    } else {
      array_list_set_flag(ALIGNMENTS_EXCEEDED, list);
    }
  
    if (array_list_get_flag(list) == DOUBLE_ANCHORS) {
      //printf("DOUBLE ANCHORS\n");
      for (int j = 0; j < array_list_size(list); j++) {
	//bwt_anchor_t *bwt_anchor_prev = array_list_get(j, list);
	cal_t *cal = array_list_get(j, list);
	metaexon_insert(0/*cal->strand*/, cal->chromosome_id - 1,
			cal->start, cal->end, 40,
			METAEXON_NORMAL, NULL,
			metaexons);
      }    
    } else if (array_list_get_flag(list) == SINGLE_ANCHORS) {
      for (int j = 0; j < array_list_size(list); j++) {
	//bwt_anchor_t *bwt_anchor = array_list_get(j, list);
	cal_t *cal = array_list_get(j, list);
	metaexon_t *metaexon;
	if (metaexon_search(0/*cal->strand*/, cal->chromosome_id - 1,
			    cal->start, cal->end, &metaexon,
			    metaexons)) {
	  metaexon_insert(0/*cal->strand*/, cal->chromosome_id - 1,
			  cal->start, cal->end, 40,
			  METAEXON_NORMAL, NULL,
			  metaexons);
	}
      }
    }
  }
  // array_list flag: 0 -> Not  BWT Anchors found
  //                  1 -> One  BWT Anchors found
  //                  2 -> Pair BWT Anchors found
  //                  3 -> Alignments found
  //                  4 -> Alignments exceded
  
  mapping_batch->num_targets = num_unmapped;

  LOG_DEBUG("========= APPLY BWT RNA END =========\n");
    
  if (batch->mapping_batch->num_targets > 0) {
    return RNA_CAL_STAGE;
  } else {
    return RNA_STAGE;
  }
    
}



//====================================================================================
// apply_bwt
//====================================================================================

int apply_bwt_bs(bwt_server_input_t* input, batch_t *batch) {

}
