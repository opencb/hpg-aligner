#include "preprocess_rna.h"

void preprocess_rna_input_init(size_t max_intron_length, size_t flank_length, 
			       size_t seeds_max_distance, size_t seed_size, 
			       genome_t *genome, preprocess_rna_input_t *preprocess_rna) {

  preprocess_rna->max_intron_length = max_intron_length;
  preprocess_rna->flank_length = flank_length;
  preprocess_rna->seeds_max_distance = seeds_max_distance;
  preprocess_rna->seed_size = seed_size;
  preprocess_rna->genome = genome;
}

preprocess_data_t *preprocess_data_new(size_t data_number, preprocess_data_t *preprocess_data) {
  preprocess_data = (preprocess_data_t *)malloc(sizeof(preprocess_data_t));
  preprocess_data->data_number = data_number;

  preprocess_data->num_cal_targets = (int *)calloc(data_number, sizeof(int));
  preprocess_data->associate_cals = (unsigned char **)malloc(sizeof(unsigned char *)*data_number);
  preprocess_data->cal_targets = (int **)malloc(sizeof(int *)*data_number);
  for (int i = 0; i < data_number; i++) {
    preprocess_data->associate_cals[i] = (unsigned char *)calloc(MAX_RNA_CALS, sizeof(unsigned char));
    preprocess_data->cal_targets[i] = (int *)calloc(MAX_RNA_CALS, sizeof(int));    
  }

  preprocess_data->negative_strand = (unsigned char *)calloc(sizeof(unsigned char), data_number);

  return preprocess_data;
}

void preprocess_data_free(preprocess_data_t *preprocess_data) {
  for (int i = 0; i < preprocess_data->data_number; i++) {
    free(preprocess_data->associate_cals[i]);
    free(preprocess_data->cal_targets[i]);
  }
  free(preprocess_data->num_cal_targets);
  free(preprocess_data->associate_cals);
  free(preprocess_data->cal_targets);
  free(preprocess_data->negative_strand);
  free(preprocess_data);
}


void print_cal(cal_t *cal){
  printf("chr:%i<->strand:(%i) <-> genome:[%lu-%lu] <-> seeds:%i\n", cal->chromosome_id, 
	 cal->strand, cal->start, cal->end, cal->num_seeds);
}

void show_cals_list(array_list_t* cals_list) {
  size_t size_list = array_list_size(cals_list);
  cal_t *cal;
  printf("CALS LIST:\n");
  for (int i = 0; i < size_list; i++) {
    cal = array_list_get(i, cals_list);
    print_cal(cal);
  }
  printf("\n");
}

void show_cals_associate(unsigned char *associate_list, size_t n_cals) {
  size_t i = 0;
  printf("ASOCIATE LIST (%i):\n", n_cals);
  
  for (i = 0; i < n_cals; i++) { 
    printf("%i,", associate_list[i]);
  }
  printf("\n");

  while (i < n_cals) {
    while (associate_list[i] != i && i < n_cals) {
      printf("%i->", associate_list[i]);
      i++;
    }
    printf("%i\n", associate_list[i]);
    i++;
  }  
}

int apply_preprocess_rna(preprocess_rna_input_t* input, batch_t *batch) {
  //printf("APPLY PREPROCESS RNA...\n");
  /*
    struct timeval start, end;
    double time;

    if (time_on) { start_timer(start); }
    
    mapping_batch_t *mapping_batch = batch->mapping_batch;
    int c, p, j;
    size_t num_cals;
    cal_t *cal_prev, *cal_next, *cal_aux;
    array_list_t *cals_list;
    size_t max_intron_size = input->max_intron_length;
    size_t flank_length = input->flank_length;
    size_t seed_size = input->seed_size;
    size_t seeds_max_distance = input->seeds_max_distance;
    size_t num_targets = mapping_batch->num_targets;
    size_t num_reads = array_list_size(mapping_batch->fq_batch);
    size_t start_prev, end_prev, start_next, end_next, start_cal, end_cal;
    size_t len;
    preprocess_data_t *preprocess_data = preprocess_data_new(num_targets, preprocess_data);
    genome_t *genome = input->genome;

    unsigned char **associate_cals = preprocess_data->associate_cals;
    unsigned char *negative_strand = preprocess_data->negative_strand;
    int *cal_targets;
    int *minor_priority_fusion = (int *)malloc(sizeof(int)*MAX_RNA_CALS);
    int *minor_priority_single = (int *)malloc(sizeof(int)*MAX_RNA_CALS);
    int t_f, t_s;
    fastq_read_t *read;
    
    for (int i = 0; i < num_targets; i++) {
      cal_targets = preprocess_data->cal_targets[i];
      cals_list = mapping_batch->mapping_lists[mapping_batch->targets[i]];
      read = array_list_get(mapping_batch->targets[i], mapping_batch->fq_batch);
      len = read->length;
      //printf("%s\n", read->sequence);
      num_cals = array_list_size(cals_list);
      if (!num_cals) { continue; }

      //End Merge Check
      p = 0;
      cal_prev = (cal_t *)array_list_get(0, cals_list);
      negative_strand[i] = negative_strand[i] || cal_prev->strand;
      c = 1;
      //Concatenate CALS
      //printf("Fusion CALs %i...\n", num_cals);
      while (c < num_cals) {
	//printf("\tCAL %i\n", c);
	cal_next = (cal_t *)array_list_get(c, cals_list);
	negative_strand[i] = negative_strand[i] || cal_next->strand;
	if (cal_prev->chromosome_id == cal_next->chromosome_id && 
	    cal_prev->strand == cal_next->strand &&
	    (cal_next->start <= (cal_prev->end + max_intron_size))) {
	  associate_cals[i][p] = c;
	} else {
	  associate_cals[i][p] = c - 1;
	  p = c;
	}
	cal_prev = cal_next;
	c++;
      }
      associate_cals[i][num_cals - 1] = num_cals - 1;
      

      if (len < 90) {
	seed_size = 15;
      } else {
	seed_size = 16;
      }
      
      size_t max_seeds_found = 2 + len/seed_size + ((len - seed_size/2)/seed_size);

      if (len % seed_size > 0) {
	max_seeds_found++;
      }
      
      float seeds_factor = 0.7;
      size_t umbral_seeds = max_seeds_found * seeds_factor;      
      //printf("MAX SEEDS FOUND IN THIS READ %lu(70% => %i)\n", max_seeds_found, umbral_seeds);

      //Search for CALs that are suported by seeds_factor
      int t = 0;
      j = 0;
      t_f = 0;
      t_s = 0;
      while (j < num_cals) {
	cal_aux = (cal_t *)array_list_get(j, cals_list);	
	if (j == associate_cals[i][j]) {
	  if (cal_aux->num_seeds >= umbral_seeds) {
	    cal_targets[t++] = j;
	  } else {
	    minor_priority_single[t_s++] = j;
	  }
	} else {
	  minor_priority_fusion[t_f++] = j;
	  j = associate_cals[i][j];
	}
	j++;
	}

      //if (!t) {
      for (int m = 0; m < t_f; m++) {
	cal_targets[t++] = minor_priority_fusion[m];
      }
      for (int m = 0; m < t_s; m++) {
	cal_targets[t++] = minor_priority_single[m];
	}
	//}
      //printf("%i < %i\n", i, num_targets);
      //preprocess_data->num_cal_targets[i] = t;
      /**************************************************************
    }

    batch->optional_data = (void *)preprocess_data;
    free(minor_priority_single);
    free(minor_priority_fusion);
    //printf("END PREPROCESS RNA\n");
    if (time_on) { stop_timer(start, end, time); timing_add(time, RNA_PREPROCESS, timing); }
    //printf("APPLY PREPROCESS RNA DONE!\n");*/
    return SW_STAGE;
}
