#include "region_seeker.h"


void region_seeker_input_init(list_t *unmapped_read_list, cal_optarg_t *cal_optarg, 
			      bwt_optarg_t *bwt_optarg, bwt_index_t *bwt_index, 
			      list_t* region_list, unsigned int region_threads, 
			      unsigned int gpu_enable, int padding_left, int padding_right,
			      genome_t *genome, metaexons_t *metaexons,
			      region_seeker_input_t *input_p) {

  input_p->unmapped_read_list_p = unmapped_read_list;
  input_p->cal_optarg_p = cal_optarg;
  input_p->bwt_optarg_p = bwt_optarg;
  input_p->bwt_index_p = bwt_index;
  input_p->region_list_p = region_list;
  input_p->region_threads = region_threads;
  input_p->gpu_enable = gpu_enable;
  input_p->padding_left = padding_left;
  input_p->padding_right = padding_right;
  input_p->genome = genome;
  input_p->metaexons = metaexons;
}

//Return distance from the read_gap to the reference_gap
int close_double_anchor_gap(size_t start_genome_gap, size_t end_genome_gap,
			    int genome_strand, int genome_chromosome,
			    size_t start_read_gap, size_t end_read_gap, 
			    char *sequence,
			    genome_t *genome) {

  int genome_distance = end_genome_gap - start_genome_gap;
  int read_distance = end_read_gap - start_read_gap;
  char reference[genome_distance + 1];
  char query[read_distance + 1];
  int max_pos = (genome_distance >= read_distance)?genome_distance:read_distance;
  register int i = 0;
  register int distance = 0;

  genome_read_sequence_by_chr_index(reference, 0, 
				    genome_chromosome, &start_genome_gap, &end_genome_gap, genome);

  memcpy(query, sequence + start_read_gap, end_read_gap - start_read_gap);
  query[end_read_gap - start_read_gap] = '\0';

  //printf("CLOSE DOUBLE ANCHOR GAP\n");
  //printf("Reference(%lu-%lu): %s\n", start_genome_gap, end_genome_gap, reference);
  //printf("Query(%lu-%lu): %s\n", start_read_gap, end_read_gap, query);

  for (i = 0; i < max_pos; i++) {
    if (query[i] != reference[i]) { distance++; }
  }
  
  //printf("Distance %i\n", distance);

  return distance;

}

//====================================================================================
// apply_seeding
//====================================================================================

int apply_seeding(region_seeker_input_t* input, batch_t *batch) {
  //printf("APPLY SEEDING...\n");


  //if (time_on) { start_timer(start); }


  mapping_batch_t *mapping_batch = batch->mapping_batch;

  size_t num_mappings;


  int seed_size = input->cal_optarg_p->seed_size;
  size_t min_seed_size = input->cal_optarg_p->min_seed_size;






  size_t num_targets = mapping_batch->num_targets;
  size_t *targets = mapping_batch->targets;
  size_t new_num_targets = 0;
  fastq_read_t *read;

  int min_intron_size = 40;

  int target;
  bwt_anchor_t *bwt_anchor = NULL;
  region_t *region;
  int gap_nt;
  int start_search;
  int end_search;

  // set to zero
  mapping_batch->num_to_do = 0;
  
  //TODO: omp parallel for !!
  /*if (batch->mapping_mode == 1000) {
    for (size_t i = 0; i < num_targets; i++) {
      //printf("Seq (i=%i)(target=%i): %s\n", i, targets[i], read->sequence);
      read = array_list_get(targets[i], mapping_batch->fq_batch);
      num_mappings = bwt_map_exact_seeds_seq(padding_left,
					     padding_right,
					     read->sequence,
					     seed_size,
					     min_seed_size,
					     input->bwt_optarg_p, 
					     input->bwt_index_p, 
					     mapping_batch->mapping_lists[targets[i]],
					     mapping_batch->extra_stage_id[targets[i]]);
      
      //printf("Num mappings %i\n", num_mappings);
      if (num_mappings > 0) {
	array_list_set_flag(2, mapping_batch->mapping_lists[targets[i]]);
	targets[new_num_targets++] = targets[i];
	mapping_batch->num_to_do += num_mappings;
      }
    }
    } else {*/
  
  //size_t new_num_targets = 0;
  //size_t *new_targets = (size_t *)malloc(array_list_size(fq_batch)*sizeof(size_t));
  array_list_t *array_list_aux = array_list_new(256, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  
  //Flag 0: The read has simple anchor or any, and need seeds and normal Cal_Seeker 
  //Flag 1: The read has double anchor and the gap is smaller than MIN_INTRON_SIZE. Cal_Seeker will be make one CAL
  //Flag 2: The read has double anchor but the gap is bigger than MIN_INTRON_SIZE. 
  for (size_t i = 0; i < num_targets; i++) {
    read = array_list_get(targets[i], mapping_batch->fq_batch);    
    //printf("Read Region %s: \n", read->id);
    /* if (array_list_get_flag(mapping_batch->mapping_lists[targets[i]]) == 0 ||
	array_list_get_flag(mapping_batch->mapping_lists[targets[i]]) == 1) {
      array_list_clear(mapping_batch->mapping_lists[targets[i]], bwt_anchor_free);
      continue;
      }
    */
    if (array_list_get_flag(mapping_batch->mapping_lists[targets[i]]) == 0 || 
	array_list_get_flag(mapping_batch->mapping_lists[targets[i]]) == 1) {
      //Flag 0 Case, Not anchors found, Make normal seeds      
      //      printf("***** Normal Case 0. Not anchors found!\n");
      for (int j = array_list_size(mapping_batch->mapping_lists[targets[i]]) - 1; j >= 0; j--) {
	bwt_anchor = array_list_remove_at(j, mapping_batch->mapping_lists[targets[i]]);
	array_list_insert(bwt_anchor, array_list_aux);
      }
      num_mappings = 0;

      num_mappings = bwt_map_exact_seeds_seq(0,
					     0,
					     read->sequence,
					     seed_size,
					     min_seed_size,
					     input->bwt_optarg_p, 
					     input->bwt_index_p, 
					     mapping_batch->mapping_lists[targets[i]],
					     0);
      if (num_mappings > 0) {
	array_list_set_flag(0, mapping_batch->mapping_lists[targets[i]]);
	targets[new_num_targets++] = targets[i];
	//mapping_batch->num_to_do += num_mappings;
      }
    } else if (array_list_get_flag(mapping_batch->mapping_lists[targets[i]]) == 1) {
      //Flag 1 Case, One anchor found, Make displacements seeds                  
      printf("***** Case 1. One anchor found!\n");
      for (int j = array_list_size(mapping_batch->mapping_lists[targets[i]]) - 1; j >= 0; j--) {
	bwt_anchor = array_list_remove_at(j, mapping_batch->mapping_lists[targets[i]]);
	array_list_insert(bwt_anchor, array_list_aux);
      }

      int anchor_nt = bwt_anchor->end - bwt_anchor->start;
      int seed_id = 0;
      int seed_start, seed_end;
      int extra_seed;

      if ((bwt_anchor->type == FORWARD_ANCHOR && bwt_anchor->strand == 0) || 
	  (bwt_anchor->type == BACKWARD_ANCHOR && bwt_anchor->strand == 1 )) {
	start_search = anchor_nt + 1;
	end_search = read->length - 1;
	extra_seed = EXTRA_SEED_END;
      } else {
	start_search = 0;
	end_search = read->length - anchor_nt - 2;
	extra_seed = EXTRA_SEED_START;
      }

      printf("end_start %i - start_search %i = %i >= seed_size %i\n", end_search, start_search, end_search - start_search, seed_size);
      if (end_search - start_search >= seed_size) {
	printf("00 bwt_map_exact_seeds_between_coords --> searching from %i to %i\n", start_search, end_search);
	/*
		num_mappings = bwt_map_exact_seeds_between_coords(start_search,
							  end_search,
							  read->sequence, 
							  seed_size, min_seed_size,
							  input->bwt_optarg_p, 
							  input->bwt_index_p, 
							  mapping_batch->mapping_lists[targets[i]],
							  extra_seed, &seed_id);
	*/
      }

      if (bwt_anchor->type == FORWARD_ANCHOR) {
	seed_id = 0;
	seed_start = 0;
	seed_end = anchor_nt;
      } else {
	seed_id += 1;
	seed_start = read->length - anchor_nt - 1;
	seed_end = read->length - 1;
      }

      for (int j = 0; j < array_list_size(array_list_aux); j++) {
	bwt_anchor_t *bwt_anchor = array_list_get(j, array_list_aux);
	//	printf("\tCreate seed Anchor [%i:%lu|%i-%i|%lu]\n", bwt_anchor->chromosome + 1, bwt_anchor->start, 
	//	       seed_start,seed_end,bwt_anchor->end);
	region = region_bwt_new(bwt_anchor->chromosome + 1,
				bwt_anchor->strand,
				bwt_anchor->start,
				bwt_anchor->end,
				seed_start,
				seed_end,
				read->length,
				seed_id);	  
	array_list_insert(region, mapping_batch->mapping_lists[targets[i]]);
      } 
      array_list_clear(array_list_aux, (void *)bwt_anchor_free); 
      array_list_set_flag(0, mapping_batch->mapping_lists[targets[i]]);
      targets[new_num_targets++] = targets[i];
    } else {
      //Flag 2 Case, Pair of anchors found
      printf("***** Case 2. Double anchor found!\n");
      bwt_anchor_t *bwt_anchor;
      bwt_anchor_t *bwt_anchor_forw, *bwt_anchor_back;

      int read_nt, genome_nt;


      int distance;
      int found = 0;
      region_t *region;
      int seed_id = 0;
      //if (array_list_size(mapping_batch->mapping_lists[targets[i]]) > 2) {
      int *anchors_targets = (int *)calloc(array_list_size(mapping_batch->mapping_lists[targets[i]]), sizeof(int));
      int num = 0;

      //min_intron_size = 0;

      //Search if one anchor is at the same distance from the reference and the read
      for (int b = 0; b < array_list_size(mapping_batch->mapping_lists[targets[i]]); b += 2) {
	bwt_anchor_forw = array_list_get(b, mapping_batch->mapping_lists[targets[i]]);
	bwt_anchor_back = array_list_get(b + 1, mapping_batch->mapping_lists[targets[i]]);
	//printf("FORW=%i:%lu-%lu BACK=%i:%lu-%lu\n", bwt_anchor_forw->chromosome, bwt_anchor_forw->start, bwt_anchor_forw->end,
	//     bwt_anchor_back->chromosome, bwt_anchor_back->start, bwt_anchor_back->end);
	read_nt = read->length - ((bwt_anchor_forw->end - bwt_anchor_forw->start) + (bwt_anchor_back->end - bwt_anchor_back->start));
	genome_nt = bwt_anchor_back->start - bwt_anchor_forw->end;	  
	distance = abs(genome_nt - read_nt);
	//printf("\t%i:Distance %i\n", b, distance);
	if (distance < min_intron_size) {
	  found = 1;
	} else {
	  anchors_targets[num++] = b;
	}
      }

      if (found) {
	//printf("\tFound Exact Case... Delete other anchors\n");
	for (int t = num - 1; t >= 0; t--) {
	  target = anchors_targets[t];
	  //printf("\tDelete %i, %i-->\n", target, target + 1);
	  bwt_anchor = array_list_remove_at(target + 1, mapping_batch->mapping_lists[targets[i]]);
	  bwt_anchor_free(bwt_anchor);
	  
	  bwt_anchor = array_list_remove_at(target, mapping_batch->mapping_lists[targets[i]]);
	  bwt_anchor_free(bwt_anchor);
	}
	array_list_set_flag(1, mapping_batch->mapping_lists[targets[i]]);
      } else {
	//Seeding between anchors
	//printf("\tFound gap between anchors \n");
	array_list_t *anchors_forward = array_list_new(array_list_size(mapping_batch->mapping_lists[targets[i]]),
						       1.25f, COLLECTION_MODE_ASYNCHRONIZED);
	array_list_t *anchors_backward = array_list_new(array_list_size(mapping_batch->mapping_lists[targets[i]]),
							1.25f, COLLECTION_MODE_ASYNCHRONIZED);
	int big_gap = 0;
	int final_anchor_nt = 0;
	int anchor_nt;
	int anchor_type;
	int anchor_strand;
	for (int j = array_list_size(mapping_batch->mapping_lists[targets[i]]) - 1; j >= 0; j -= 2) {
	  bwt_anchor_back = array_list_remove_at(j, mapping_batch->mapping_lists[targets[i]]);
	  array_list_insert(bwt_anchor_back, anchors_backward);

	  bwt_anchor_forw = array_list_remove_at(j - 1, mapping_batch->mapping_lists[targets[i]]);
	  array_list_insert(bwt_anchor_forw, anchors_forward);

	  if (bwt_anchor_forw->strand == 0) {
	    anchor_nt = bwt_anchor_forw->end - bwt_anchor_forw->start;
	    gap_nt = read->length - (anchor_nt + (bwt_anchor_back->end - bwt_anchor_back->start));
	  } else {
	    anchor_nt = bwt_anchor_back->end - bwt_anchor_back->start;
	    gap_nt = read->length - (anchor_nt + (bwt_anchor_forw->end - bwt_anchor_forw->start));
	  }
	  if (gap_nt < 0) { gap_nt = 0; }
	  //printf("Gap nt (%i - %i): %i\n", anchor_nt, bwt_anchor_back->end - bwt_anchor_back->start, gap_nt);
	  if (gap_nt > big_gap) {
	    big_gap = gap_nt;
	    final_anchor_nt = anchor_nt;
	    anchor_type = bwt_anchor_back->type;
	    anchor_strand = bwt_anchor_back->strand;
	  }
	}

	printf("%i, %i\n", big_gap - 2, seed_size);
	if (big_gap - 2 > seed_size) {
	  //if (anchor_type == FORWARD_ANCHOR && anchor_strand == 0 || 
	  //  anchor_type == BACKWARD_ANCHOR && anchor_strand == 1 ) {
	    start_search = final_anchor_nt + 1;
	    end_search = final_anchor_nt + big_gap - 1;
	    //} else {
	    // start_search = final_anchor_nt + big_gap - 1;
	    //end_search = final_anchor_nt + 1;
	    //}
	  
	    //printf("Seeding between anchors... gap=%i\n", big_gap);
	    printf("11 bwt_map_exact_seeds_between_coords --> searching from %i to %i\n", start_search, end_search);
	    /*
	    num_mappings = bwt_map_exact_seeds_between_coords(start_search,
							      end_search,
							      read->sequence, seed_size, min_seed_size,
							      input->bwt_optarg_p, 
							      input->bwt_index_p, 
							      mapping_batch->mapping_lists[targets[i]],
							      EXTRA_SEED_NONE,
							      &seed_id);
	    */
	}

	//printf("Making seeds anchors...\n");
	for (int a = 0; a < array_list_size(anchors_forward); a++) {
	  //Insert the last anchor. (Create new seed)
	  bwt_anchor_forw = array_list_get(a, anchors_forward);
	  bwt_anchor_back = array_list_get(a, anchors_backward);

	  anchor_nt = bwt_anchor_forw->end - bwt_anchor_forw->start;
	  gap_nt = read->length - (anchor_nt + (bwt_anchor_back->end - bwt_anchor_back->start));

	  //printf("\t --> Big Seed: %i, gap_nt: %i, anchor_nt = %i\n", a, gap_nt, anchor_nt);
	  if (gap_nt < 0) {
	    //gap_nt = 0;
	    bwt_anchor_forw->end   += gap_nt;
	    bwt_anchor_back->start -= gap_nt;
	    anchor_nt += gap_nt;
	    gap_nt = 0;
	  } else if (gap_nt == 0) {
	    bwt_anchor_forw->end -= 1;
	    bwt_anchor_back->start += 1;
	    anchor_nt -= 1;
	    gap_nt = 1;	    
	  }


	  region = region_bwt_new(bwt_anchor_forw->chromosome + 1,
				  bwt_anchor_forw->strand,
				  bwt_anchor_forw->start,
				  bwt_anchor_forw->end,
				  0,
				  anchor_nt,
				  read->length,
				  0);
	  //printf("Region: %i-%i\n", region->seq_start, region->seq_end);
	  array_list_insert(region, mapping_batch->mapping_lists[targets[i]]);

	  region = region_bwt_new(bwt_anchor_back->chromosome + 1,
				  bwt_anchor_back->strand,
				  bwt_anchor_back->start,
				  bwt_anchor_back->end,
				  anchor_nt + gap_nt,
				  read->length - 1,
				  read->length,
				  seed_id + 1);
	  //printf("Region: %i-%i\n", region->seq_start, region->seq_end);
	  array_list_insert(region, mapping_batch->mapping_lists[targets[i]]);

	  //printf("\tMaking seeds anchors end, %i seeds\n", array_list_size(mapping_batch->mapping_lists[targets[i]]));

	  bwt_anchor_free(bwt_anchor_back);
	  bwt_anchor_free(bwt_anchor_forw);
	}
	array_list_free(anchors_forward, NULL);
	array_list_free(anchors_backward, NULL);
	//printf("Making seeds anchors end, %i seeds\n", array_list_size(mapping_batch->mapping_lists[targets[i]]));

	array_list_set_flag(2, mapping_batch->mapping_lists[targets[i]]);
      }
      free(anchors_targets);
      targets[new_num_targets++] = targets[i];
    }
  }

  mapping_batch->num_targets = new_num_targets;

  array_list_free(array_list_aux, NULL);

  //if (time_on) { stop_timer(start, end, time); timing_add(time, REGION_SEEKER, timing); }

  //printf("APPLY SEEDING DONE!\n");
  
  return CAL_STAGE;

}

//------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
