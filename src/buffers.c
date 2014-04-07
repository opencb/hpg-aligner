
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "buffers.h"

//#define MAXLINE 2048

//===================================================================================
batch_t *batch_new(bwt_server_input_t *bwt_input,
                   region_seeker_input_t *region_input,
                   cal_seeker_input_t *cal_input,
                   pair_server_input_t *pair_input,
		   preprocess_rna_input_t *preprocess_rna,
                   sw_server_input_t *sw_input,
                   batch_writer_input_t *writer_input,
		   int mapping_mode,
                   mapping_batch_t *mapping_batch) {

  batch_t *b = (batch_t *) calloc(1, sizeof(batch_t));
  b->bwt_input = bwt_input;
  b->region_input = region_input;
  b->cal_input = cal_input;
  b->pair_input = pair_input;
  b->sw_input = sw_input;
  b->writer_input = writer_input;
  b->mapping_batch = mapping_batch;
  b->mapping_mode = mapping_mode;
  b->preprocess_rna = preprocess_rna;

  return b;
}

void batch_free(batch_t *b) {
  if (b) free(b);
}


//====================================================================================

void region_batch_init(array_list_t **allocate_mapping_p, fastq_batch_t *unmapped_batch_p, region_batch_t *region_batch_p){
  region_batch_p->allocate_mapping_p = allocate_mapping_p;
  region_batch_p->unmapped_batch_p = unmapped_batch_p;
}

void region_batch_free(region_batch_t *region_batch_p){
  for(int i = 0; i < region_batch_p->unmapped_batch_p->num_reads; i++){
    array_list_free(region_batch_p->allocate_mapping_p[i], (void *)region_bwt_free);
  }
  free(region_batch_p->allocate_mapping_p);
  fastq_batch_free(region_batch_p->unmapped_batch_p);
  free(region_batch_p);
  
}

//====================================================================================

sw_batch_t* sw_batch_new(unsigned int num_reads, array_list_t **allocate_cals_p, 
			 fastq_read_t **allocate_reads_p) {
  sw_batch_t* sw_batch_p = (sw_batch_t *)malloc(sizeof(sw_batch_t));

  sw_batch_p->num_reads = num_reads;
  sw_batch_p->allocate_reads_p = allocate_reads_p;
  sw_batch_p->allocate_cals_p = allocate_cals_p;
  
  return sw_batch_p;
}

void sw_batch_free(sw_batch_t *sw_batch_p) {
  
  for(int i = 0; i < sw_batch_p->num_reads; i++){
    array_list_free(sw_batch_p->allocate_cals_p[i], (void *)cal_free);
    fastq_read_free(sw_batch_p->allocate_reads_p[i]);
  }

  free(sw_batch_p->allocate_cals_p);
  free(sw_batch_p->allocate_reads_p);
  free(sw_batch_p);
}
void sw_batch_init(unsigned int num_reads, array_list_t **allocate_cals_p, 
		   fastq_read_t **allocate_reads_p, sw_batch_t *sw_batch_p) {
  sw_batch_p->num_reads = num_reads;
  sw_batch_p->allocate_reads_p = allocate_reads_p;
  sw_batch_p->allocate_cals_p = allocate_cals_p;
}

//====================================================================================
//  write_batch functions
//====================================================================================

write_batch_t* write_batch_new(unsigned int allocate_size, unsigned char flag) {
  write_batch_t* write_batch_p = (write_batch_t*) calloc(1, sizeof(write_batch_t));
  
  write_batch_p->flag = flag;
  write_batch_p->size = 0;
  
  if(flag != MATCH_FLAG){
    write_batch_p->allocated_size = allocate_size;
    write_batch_p->buffer_p = (void *) calloc(allocate_size, sizeof(char));
  }else{
    write_batch_p->allocated_size = allocate_size/sizeof(alignment_t *);
    write_batch_p->buffer_p = (void *) calloc(write_batch_p->allocated_size, sizeof(alignment_t *));
  }
  return write_batch_p;
}

//------------------------------------------------------------------------------------

void write_batch_free(write_batch_t* write_batch_p) {
 if (write_batch_p == NULL) return;
 
 if (write_batch_p->buffer_p != NULL) free(write_batch_p->buffer_p);
 
 free(write_batch_p);
}

//====================================================================================

report_optarg_t *report_optarg_new(int all, int n_best, int n_hits, int only_paired, int best) {
  report_optarg_t *p = (report_optarg_t*) calloc(1, sizeof(report_optarg_t));

  p->all = all;
  p->n_best = n_best;
  p->n_hits = n_hits;
  p->only_paired = only_paired;
  p->best = best;

  return p;
}

//------------------------------------------------------------------------------------

void report_optarg_free(report_optarg_t *p) {
  if (p != NULL)
    free(p);
}

//====================================================================================

pair_mng_t *pair_mng_new(int pair_mode, size_t min_distance, 
			 size_t max_distance, int report_only_paired) {
  pair_mng_t *p = (pair_mng_t*) calloc(1, sizeof(pair_mng_t));

  p->pair_mode = pair_mode;
  p->min_distance = min_distance;
  p->max_distance = max_distance;
  p->report_only_paired = report_only_paired;

  return p;
}

//------------------------------------------------------------------------------------

void pair_mng_free(pair_mng_t *p) {
  if (p != NULL)
    free(p);
}

//====================================================================================

cal_batch_t* cal_batch_new(array_list_t **allocate_mapping, fastq_batch_t *unmapped_batch){
  cal_batch_t* cal_batch = (cal_batch_t *)malloc(sizeof(cal_batch_t));
  
  cal_batch->allocate_mapping = allocate_mapping;
  cal_batch->unmapped_batch = unmapped_batch;
  
  return cal_batch;
}

void cal_batch_free(cal_batch_t *cal_batch){
  for(int i = 0; i < cal_batch->unmapped_batch->num_reads; i++){
    array_list_free(cal_batch->allocate_mapping[i], (void *)region_bwt_free);
  }
  free(cal_batch->allocate_mapping);
  fastq_batch_free(cal_batch->unmapped_batch);
  free(cal_batch);
  
}

//====================================================================================

unsigned int pack_junction(unsigned int chromosome, unsigned int strand, 
			   size_t start, size_t end, size_t junction_id,
			   size_t num_reads, char *type, char* buffer_p){
  int len;
  char str[1024];
  char *chr_p, *p = buffer_p;
  char strand_char[2] = {'+', '-'};

  if (chromosome == 23) { sprintf(str, "%c\0", 'X'); }
  else if (chromosome == 24) { sprintf(str, "%c\0", 'Y'); }
  else if (chromosome == 25) { sprintf(str, "%s\0", "MT"); }
  else { sprintf(str, "%i\0", chromosome); }
 
  len = strlen(str);
  memcpy(p, str, len);
  p += len;
  *p = '\t';
  p++;
  
  sprintf(str, "%lu", start);
  len = strlen(str);
  memcpy(p, str, len); 
  p += len;
  *p = '\t'; 
  p++;
  
  sprintf(str, "%lu", end);
  len = strlen(str);
  memcpy(p, str, len); 
  p += len;
  *p = '\t'; 
  p++;
  
  sprintf(str, "JUNCTION_%lu", junction_id);
  len = strlen(str);
  memcpy(p, str, len); 
  p += len;
  *p = '\t'; 
  p++;
  
  sprintf(str, "%lu", num_reads);
  len = strlen(str);
  memcpy(p, str, len); 
  p += len;
  *p = '\t'; 
  p++;
  
  *p = strand_char[strand]; 
  p++;
  *p = '\t'; 
  p++;
  
  len = strlen(type);
  memcpy(p, type, len);
  p += len;
  *p = '\n';
  p++;

  return (p - buffer_p);
}

//=====================================================================================
//=====================================================================================

mapping_batch_t *mapping_batch_new(array_list_t *fq_batch, pair_mng_t *pair_mng) {

  mapping_batch_t *p = (mapping_batch_t *) calloc(1, sizeof(mapping_batch_t));
  size_t num_reads = array_list_size(fq_batch);

  p->action = BWT_ACTION;
  p->num_targets = 0;
  p->num_extra_targets = 0;
  p->num_allocated_targets = num_reads;
  p->extra_stage_do = 0;

  if (!pair_mng) { 
    p->pair_mng = pair_mng_new(SINGLE_END_MODE, 0, 0, 0); 
  } else {
    p->pair_mng = pair_mng_new(pair_mng->pair_mode, pair_mng->min_distance, 
			       pair_mng->max_distance, pair_mng->report_only_paired); 
  }

  p->num_gaps = 0;
  p->num_sws = 0;
  p->num_ext_sws = 0;

  p->num_to_do = 0;
  p->fq_batch = fq_batch;
  p->targets = (size_t *) calloc(num_reads, sizeof(size_t));
  p->extra_targets = (size_t *) calloc(num_reads, sizeof(size_t));
  p->extra_stage_id = (unsigned char *) calloc(num_reads, sizeof(unsigned char));
  p->mapping_lists = (array_list_t **) calloc(num_reads, sizeof(array_list_t*));

  //for debug. TODO:delete
  p->bwt_mappings = (unsigned char *)calloc(num_reads, sizeof(unsigned char));

  for (size_t i = 0; i < num_reads; i++) {
    p->mapping_lists[i] = array_list_new(500, 
					 1.25f, 
					 COLLECTION_MODE_ASYNCHRONIZED); 
  }
    
  p->histogram_sw = (size_t *)calloc(1024, sizeof(size_t));

  // added by PP for bisulfite
  p->num_targets2 = 0;
  p->num_to_do2 = 0;
  p->targets2 = (size_t *) calloc(num_reads, sizeof(size_t));

  return p;
}

//------------------------------------------------------------------------------------

mapping_batch_t *mapping_batch_new_2(size_t num_reads, 
				     array_list_t *fq_batch, 
				     pair_mng_t *pair_mng) {

  mapping_batch_t *p = (mapping_batch_t *) calloc(1, sizeof(mapping_batch_t));

  p->action = BWT_ACTION;
  p->num_targets = 0;
  p->num_extra_targets = 0;
  p->num_allocated_targets = num_reads;
  p->extra_stage_do = 0;

  if (!pair_mng) { 
    p->pair_mng = pair_mng_new(SINGLE_END_MODE, 0, 0, 0); 
  } else {
    p->pair_mng = pair_mng_new(pair_mng->pair_mode, pair_mng->min_distance, 
			       pair_mng->max_distance, pair_mng->report_only_paired); 
  }

  p->num_gaps = 0;
  p->num_sws = 0;
  p->num_ext_sws = 0;

  p->num_to_do = 0;
  p->fq_batch = fq_batch;
  p->targets = (size_t *) calloc(num_reads, sizeof(size_t));
  p->extra_targets = (size_t *) calloc(num_reads, sizeof(size_t));
  p->extra_stage_id = (unsigned char *) calloc(num_reads, sizeof(unsigned char));
  p->mapping_lists = (array_list_t **) calloc(num_reads, sizeof(array_list_t*));

  /*for (size_t i = 0; i < num_reads; i++) {
    p->mapping_lists[i] = array_list_new(500, 
					 1.25f, 
					 COLLECTION_MODE_ASYNCHRONIZED); 
					 }*/

  //for debug. TODO:delete
  p->bwt_mappings = (unsigned char *)calloc(num_reads, sizeof(unsigned char));

  return p;
}

//------------------------------------------------------------------------------------

/*mapping_batch_t *mapping_batch_new_by_num(size_t num_reads, pair_mng_t *pair_mng) {

  mapping_batch_t *p = (mapping_batch_t *) calloc(1, sizeof(mapping_batch_t));

  p->action = BWT_ACTION;
  p->num_targets = 0;
  p->num_extra_targets = 0;
  p->num_allocated_targets = num_reads;
  p->extra_stage_do = 0;

  if (!pair_mng) { 
    p->pair_mng = pair_mng_new(SINGLE_END_MODE, 0, 0, 0); 
  } else {
    p->pair_mng = pair_mng_new(pair_mng->pair_mode, pair_mng->min_distance, 
			       pair_mng->max_distance, pair_mng->report_only_paired); 
  }

  p->num_to_do = 0;
  p->targets = (size_t *) calloc(num_reads, sizeof(size_t));
  p->extra_targets = (size_t *) calloc(num_reads, sizeof(size_t));
  p->extra_stage_id = (unsigned char *) calloc(num_reads, sizeof(unsigned char));
  p->mapping_lists = (array_list_t **) calloc(num_reads, sizeof(array_list_t*));
  p->old_mapping_lists = (array_list_t **) calloc(num_reads, sizeof(array_list_t*));

  for (size_t i = 0; i < num_reads; i++) {
    p->mapping_lists[i] = array_list_new(10, 
					 1.25f, 
					 COLLECTION_MODE_ASYNCHRONIZED);
  }
    
  return p;
}
*/

//------------------------------------------------------------------------------------

void mapping_batch_free(mapping_batch_t *p) {
  if (p == NULL) return;
  
  if (p->fq_batch) { array_list_free(p->fq_batch, (void *) fastq_read_free); }
  if (p->targets) { free(p->targets); }
  if (p->mapping_lists) { free(p->mapping_lists); }
  if (p->pair_mng) { free(p->pair_mng); }
  if (p->extra_stage_id) { free(p->extra_stage_id); }
  if (p->extra_targets) { free(p->extra_targets); }

  if (p->old_mapping_lists) { free(p->old_mapping_lists); }
  if (p->bwt_mappings) free(p->bwt_mappings);

  // added by PP
  if (p->CT_fq_batch) { array_list_free(p->CT_fq_batch, (void *) fastq_read_free); }
  if (p->CT_rev_fq_batch) { array_list_free(p->CT_rev_fq_batch, (void *) fastq_read_free); }
  if (p->GA_fq_batch) { array_list_free(p->GA_fq_batch, (void *) fastq_read_free); }
  if (p->GA_rev_fq_batch) { array_list_free(p->GA_rev_fq_batch, (void *) fastq_read_free); }
  if (p->mapping_lists2) { free(p->mapping_lists2); }
  if (p->targets2) { free(p->targets2); }
  if (p->bs_status) {free(p->bs_status); }
  
  free(p);
}

//------------------------------------------------------------------------------------
/*
rna_batch_t *rna_batch_new(array_list_t *fq_batch) {
  rna_batch_t *p = (rna_batch_t *) calloc(1, sizeof(rna_batch_t));

  size_t num_reads = array_list_size(fq_batch);

  p->action = BWT_ACTION;
  p->all_targets = 1;
  p->num_targets = 0;
  p->num_allocated_targets = num_reads;
  p->num_mapping_lists = num_reads;

  p->num_done = 0;
  p->num_to_do = 0;

  p->fq_batch = fq_batch;
  p->targets = (size_t *) calloc(num_reads, sizeof(size_t));
  p->mapping_lists = (array_list_t **) calloc(num_reads, sizeof(array_list_t*));
  for (size_t i = 0; i < num_reads; i++) {
    p->mapping_lists[i] = array_list_new(500, 
					 1.25f, 
					 COLLECTION_MODE_ASYNCHRONIZED); 
  }
    
  return p;
}

//------------------------------------------------------------------------------------

void rna_batch_free(rna_batch_t *p) {
  if (p == NULL) return;
  
  if (p->fq_batch != NULL) array_list_free(p->fq_batch, (void *)fastq_read_free);
  if (p->targets != NULL) free(p->targets);
  if (p->mapping_lists != NULL) { free(p->mapping_lists); }
  
  free(p);
}
*/

//------------------------------------------------------------------------------------

buffer_item_t *buffer_item_new() {
  return (buffer_item_t *)malloc(sizeof(buffer_item_t));
}


/*buffer_item_t *buffer_item_complete_new(fastq_read_t *fastq_read, array_list_t *items_list, void *aux_data) {
  buffer_item_t *buffer_item = buffer_item_new();
  buffer_item->read = fastq_read;
  buffer_item->items_list = array_list_new(array_list_size(items_list), 
					   1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  
  for (int i = 0; i < array_list_size(items_list); i++) {
    void *item = array_list_get(i, items_list);
    array_list_insert(item, buffer_item->items_list);
  }

  buffer_item->aux_data = aux_data;

  return buffer_item;
  
}
*/

void alignment_aux_init(alignment_t* alignment, 
			alignment_aux_t *alignment_aux) {
  alignment_aux->seq_strand = alignment->seq_strand;
  alignment_aux->chromosome = alignment->chromosome; 
  alignment_aux->position = alignment->position;
  alignment_aux->num_cigar_operations = alignment->num_cigar_operations;
  alignment_aux->map_quality = alignment->map_quality;
  alignment_aux->optional_fields_length = alignment->optional_fields_length;
  alignment_aux->mapping_len = strlen(alignment->sequence);
  alignment_aux->cigar_len = strlen(alignment->cigar);
}


fastq_read_t *file_read_fastq_reads(size_t *num_items, FILE *fd) {

  size_t sizes_to_read[3], head_len, seq_len;
  int bytes;

  bytes = fread(sizes_to_read, sizeof(size_t), 3, fd);
  if (!bytes) { return NULL; }
  
  head_len   = sizes_to_read[0];
  seq_len    = sizes_to_read[1];
  *num_items = sizes_to_read[2];

  int tot_size = head_len + 2*seq_len;
  char *buffer = (char *)calloc(tot_size + 1, sizeof(char));
  bytes = fread(buffer, sizeof(char), tot_size, fd);

  if (!bytes) {
    free(buffer);
    return NULL; 
  }

  char *id = (char *)calloc(head_len + 1, sizeof(char));
  memcpy(id, buffer, head_len);
  //printf("ID : %s\n", id);

  char *sequence = (char *)calloc(seq_len + 1, sizeof(char));  
  memcpy(sequence, &buffer[head_len], seq_len);
  //printf("SEQ: %s\n", sequence);

  char *quality = (char *)calloc(seq_len + 1, sizeof(char));  
  memcpy(quality, &buffer[head_len + seq_len], seq_len);
  //printf("QUA: %s\n", quality);
  
  fastq_read_t *fq_read = fastq_read_new(id, sequence, quality);

  free(buffer);
  free(id);
  free(sequence);
  free(quality);

  return fq_read;

}

int file_read_cals(size_t num_items, array_list_t *list, 
		   fastq_read_t *fq_read, FILE *fd) {

  if (num_items == 0) { return 0; }

  int bytes;
  bwt_anchor_t bwt_anchors[num_items];
  bytes = fread(bwt_anchors, sizeof(bwt_anchor_t), num_items, fd);
  if (!bytes) { LOG_FATAL("Corrupt file\n"); }
  
  for (int i = 0; i < num_items; i++) {
    //printf("[%i:%lu-%lu]\n", bwt_anchors[i].chromosome, bwt_anchors[i].start, bwt_anchors[i].end);
    size_t seed_size = bwt_anchors[i].end - bwt_anchors[i].start;
    cal_t *cal;
    if (bwt_anchors[i].type == FORWARD_ANCHOR) {
      cal = (cal_t *)convert_bwt_anchor_to_CAL(&bwt_anchors[i], 0, seed_size);
    } else {
      cal = (cal_t *)convert_bwt_anchor_to_CAL(&bwt_anchors[i], fq_read->length - seed_size - 1, fq_read->length - 1);
    }
    //cal_print(cal);
    array_list_insert(cal, list); 
  }  
  
  return 0;

}

int file_read_meta_alignments(size_t num_items, array_list_t *list, 
                              fastq_read_t *fq_read, FILE *fd) {

  if (!num_items) { return 0; }

  simple_alignment_t simple_alignment[num_items];
  simple_alignment_t *simple_a;
  int bytes;

  bytes = fread(simple_alignment, sizeof(simple_alignment_t), num_items, fd);
  if (!bytes) { LOG_FATAL("Corrupt file\n"); }
  
  size_t cigar_tot_len = 0;
  for (int i = 0; i < num_items; i++) {
    simple_a = &simple_alignment[i];
    //printf("ITEM %i: (%i)[%i:%lu] [%i-%i]\n", i, simple_a->map_strand, simple_a->map_chromosome,
    //     simple_a->map_start, simple_a->gap_start, simple_a->gap_end);
    cigar_tot_len += simple_a->cigar_len;
  }
    
  char cigar_buffer[cigar_tot_len];
  bytes = fread(cigar_buffer, sizeof(char), cigar_tot_len, fd);
  if (!bytes) { LOG_FATAL("Corrupt file\n"); }

  char cigars_test[num_items][1024];
  size_t actual_read = 0;
  for (int i = 0; i < num_items; i++) {
    simple_a = &simple_alignment[i];
    memcpy(&cigars_test[i], &cigar_buffer[actual_read], simple_a->cigar_len);
    cigars_test[i][simple_a->cigar_len] = '\0';
    actual_read += simple_a->cigar_len;
    //printf("CIGAR %i: %s\n", i, cigars_test[i]);
    size_t map_len = fq_read->length - simple_a->gap_start - simple_a->gap_end;
    //printf("SEED := len_read:%i - gap_read:%i - gap_end:%i = %i, SEED-END = %i\n", fq_read->length, 
    //     simple_a->gap_start, 
    //     simple_a->gap_end, 
    //     map_len, simple_a->gap_start + map_len);
    seed_region_t *s_region = seed_region_new(simple_a->gap_start, 
                                              simple_a->gap_start + map_len - 1,
                                              simple_a->map_start, 
                                              simple_a->map_start + map_len - 1,
                                              0, 0, 0);
    
    //printf("Exit with seed [%i:%i]\n", s_region->read_start, s_region->read_end);
    
    linked_list_t *sr_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    //s_region->info = cigar_code_new_by_string(cigars_test[i]);
    linked_list_insert(s_region, sr_list);
    
    cal_t *cal = cal_new(simple_a->map_chromosome, 
                         simple_a->map_strand,
                         simple_a->map_start,
                         simple_a->map_start + map_len - 1,
                         1,
                         sr_list,
                         linked_list_new(COLLECTION_MODE_ASYNCHRONIZED));
    cigar_code_t *cc = cigar_code_new_by_string(cigars_test[i]);
    cc->distance = simple_a->map_distance;
    cal->info = cc;

    meta_alignment_t *meta_alignment = meta_alignment_new();    
    for (int m = 0; m < array_list_size(cc->ops); m++) {
      cigar_op_t *op = array_list_get(m, cc->ops);
      cigar_code_append_new_op(op->number, op->name, meta_alignment->cigar_code);
      //array_list_insert(cigar_op, ->ops);
    }
    meta_alignment->cigar_code->distance = simple_a->map_distance;

    array_list_insert(cal, meta_alignment->cals_list);
    array_list_insert(meta_alignment, list);

  }

  return 0;
}

int file_read_alignments(size_t num_items, array_list_t *list, 
			 fastq_read_t *fq_read, FILE *fd) {

  if (!num_items) { return 0; }

  //printf("Read file alignment...\n");

  alignment_aux_t alignments_aux[num_items];
  alignment_aux_t *alignment_a;
  int bytes;

  bytes = fread(alignments_aux, sizeof(alignment_aux_t), num_items, fd);
  if (!bytes) { LOG_FATAL("Corrupt file\n"); }

  size_t cigar_tot_len = 0;
  size_t of_tot_len = 0;
  for (int i = 0; i < num_items; i++) {
    alignment_a = &alignments_aux[i];

    //printf("CIGAR: %i + %i\n", cigar_tot_len,
    //	   alignment_a->cigar_len);
    cigar_tot_len += alignment_a->cigar_len;

    //printf("OF: %i + %i\n", of_tot_len,
    //	   alignment_a->optional_fields_length);
    of_tot_len    +=  alignment_a->optional_fields_length;

  }

  char cigar_buffer[cigar_tot_len];
  bytes = fread(cigar_buffer, sizeof(char), cigar_tot_len, fd);
  if (!bytes) { LOG_FATAL("Corrupt file\n"); }

  uint8_t of_buffer[of_tot_len];
  bytes = fread(of_buffer, sizeof(uint8_t), of_tot_len, fd);
  if (!bytes) { LOG_FATAL("Corrupt file\n"); }

  char cigars_test[num_items][1024];
  size_t pos_cigar = 0, pos_of = 0;

  for (int i = 0; i < num_items; i++) {
    alignment_a = &alignments_aux[i];
    //printf("alignment_a->optional_fields_length = %i\n", alignment_a->optional_fields_length);

    memcpy(&cigars_test[i], &cigar_buffer[pos_cigar], alignment_a->cigar_len);
    cigars_test[i][alignment_a->cigar_len] = '\0';
    pos_cigar += alignment_a->cigar_len;
    
    char op;
    char op_value[1024];
    int c = 0;
    int hc_start = 0, hc_end = 0;

    //printf("CIGAR: %s\n", cigars_test[i]);

    for (int j = 0; j < alignment_a->cigar_len; j++) {
      op = cigars_test[i][j];
      if (op < 58) {
	op_value[c++] = op;
      } else {
	op_value[c] = '\0';
	if (op == 'H') {
	  hc_start = atoi(op_value);
	}
	break;
      }
    }

    c = 0;
    if (cigars_test[i][alignment_a->cigar_len - 1] == 'H') {
      for (int j = alignment_a->cigar_len - 2; j >= 0; j--) {
	op = cigars_test[i][j];
	//printf("Process op= %c\n", op);
	if (op < 58) {
	  op_value[c++] = op;
	} else {
	  op_value[c] = '\0';
	  int len = strlen(op_value);
	  char op_val_aux[len];
	  int pos = len - 1;
	  //printf("(%i) :: %s\n", len, op_value);
	  int t = 0;
	  for (t = 0; t < len; t++) {	    
	    op_val_aux[t] = op_value[pos - t];
	  } 
	  op_val_aux[t] = '\0';
	  //printf("(%i) :: %s\n", t, op_val_aux);
	  hc_end = atoi(op_val_aux);
	  break;
	}
      }
    }

    uint8_t *optional_fields = (uint8_t *)calloc(alignment_a->optional_fields_length, sizeof(uint8_t));

    memcpy(optional_fields, &of_buffer[pos_of], alignment_a->optional_fields_length);
    //optional_fields[alignment_a->optional_fields_length] = '0';
    pos_of += alignment_a->optional_fields_length;

    
    int header_len = strlen(fq_read->id);
    char header_id[header_len + 1];
    get_to_first_blank(fq_read->id, header_len, header_id);
    //char *header_match = (char *)malloc(sizeof(char)*header_len);
    //memcpy(header_match, header_id, header_len);

    int len_read = fq_read->length - (hc_start + hc_end);
    
    //printf("hc_start = %i, hc_end = %i, len_read = %i\n", 
    //	   hc_start, hc_end, len_read);
    char *quality = (char *) calloc (len_read + 1, sizeof(char));
    strncpy(quality, fq_read->quality + hc_start, len_read);
    
    char *query = (char *) calloc (len_read + 1, sizeof(char));
    strncpy(query, fq_read->sequence + hc_start, len_read);

    //Revisar rna_Server get_to_first_blank header copy
    alignment_t *alignment = alignment_new();
    alignment_init_single_end(strdup(header_id),
			      query,
			      quality,
			      alignment_a->seq_strand, 
			      alignment_a->chromosome, 
			      alignment_a->position,
			      strdup(cigars_test[i]),
			      alignment_a->num_cigar_operations,
			      alignment_a->map_quality, 
			      1, 
			      num_items < 1,
			      alignment_a->optional_fields_length,
			      optional_fields, 
			      alignment);

    array_list_insert(alignment, list);
  }  

  return 0;

}

void file_write_alignments(fastq_read_t *fq_read, array_list_t *items, FILE *fd) {
  //size_t head_size = strlen(fq_read->id);
  //size_t seq_size  = fq_read->length;

  size_t num_items = array_list_size(items);
  //printf("Num items %i\n", num_items);
  if (num_items <= 0) { return; }

  int tot_len_cigar = 0, tot_len_of = 0;
  //unsigned char type = ALIGNMENT_TYPE;

  //Write binary file 
  //[type][size head][size seq][num items][HEAD][SEQUENCE][QUALITY][ALIG 0][ALIG n][CIGAR STR][OF STR]
  //fwrite(type, sizeof(unsigned char), 1, fd);
  
  //size_t items_sizes[3] = {head_size, seq_size, num_items};
  //[size head][size seq][num items]
  //fwrite(items_sizes, sizeof(size_t), 3, fd);
  
  //[HEAD][SEQUENCE][QUALITY]
  size_t max_len_cigar = num_items*1024*2;
  char *buffer_cigar = (char *)malloc(sizeof(char)*max_len_cigar);

  size_t max_len_of = num_items*1024*2;
  uint8_t *buffer_of = (uint8_t *)malloc(sizeof(uint8_t)*max_len_of);
 
  //memcpy(buffer, fq_read->id, head_size);
  //memcpy(&buffer[head_size], fq_read->sequence, seq_size);
  //memcpy(&buffer[head_size + seq_size], fq_read->quality, seq_size);
  //fwrite(buffer, sizeof(char), total_size, fd);

  alignment_aux_t alignment_aux[num_items];
  alignment_aux_t *alignment_a;

  memset(alignment_aux, 0, sizeof(alignment_aux_t)*num_items);  
  for (int i = 0; i < num_items; i++) {
    alignment_a = &alignment_aux[i];
    alignment_t *alignment = array_list_get(i, items);
    alignment_aux_init(alignment, 
		       alignment_a);     

    int cigar_len = strlen(alignment->cigar);    
    memcpy(&buffer_cigar[tot_len_cigar], alignment->cigar, cigar_len);
    tot_len_cigar += cigar_len;
    
    if (tot_len_cigar >= max_len_cigar) { 
      max_len_cigar = max_len_cigar * 2;
      buffer_cigar = realloc(buffer_cigar, max_len_cigar); 
    }

    int of_len = alignment->optional_fields_length;
    //printf("ALig of len = %i\n", of_len);
    /*memcpy(&buffer_cigar[tot_len_cigar], alignment->optional_fields, of_len);
    tot_len_cigar += of_len;
    
    if (tot_len_cigar >= max_len_cigar) { 
      max_len_cigar = max_len_cigar * 2;
      buffer_cigar = realloc(buffer_cigar, max_len_cigar); 
    }
    */
    memcpy(&buffer_of[tot_len_of], alignment->optional_fields, of_len);
    tot_len_of += of_len;
    
    if (tot_len_of >= max_len_of) {
      max_len_of = max_len_of * 2;
      buffer_of = realloc(buffer_of, max_len_cigar); 
    }
        
  }

  fwrite(alignment_aux, sizeof(alignment_aux_t), num_items, fd);
  fwrite(buffer_cigar, sizeof(char), tot_len_cigar, fd);  
  fwrite(buffer_of, sizeof(uint8_t), tot_len_of, fd);  

  //free(buffer);
  free(buffer_cigar);  

}

void file_write_meta_alignments(fastq_read_t *fq_read, array_list_t *items, FILE *fd) {
  //size_t head_size = strlen(fq_read->id);
  size_t seq_size  = fq_read->length;
  size_t num_items = array_list_size(items);
  if (!num_items) { return; }

  size_t max_len = num_items * 1024;
  char *cigar_buffer = (char *)calloc(max_len, sizeof(char));
  size_t tot_len = 0;

  simple_alignment_t simple_alignment[num_items];
  simple_alignment_t *simple_a;
  //unsigned char type = MENTA_TYPE;
  //fwrite(type, sizeof(unsigned char), 1, fd);

  memset(simple_alignment, 0, sizeof(simple_alignment_t)*num_items);
  //[type][size head][size seq][num items][HEAD][SEQUENCE][QUALITY][CAL 0][CAL n]
  for (int i = 0; i < num_items; i++) {
    meta_alignment_t *meta_alignment = array_list_get(i, items);
    cal_t *first_cal = array_list_get(0, meta_alignment->cals_list);
    cal_t *last_cal = array_list_get(meta_alignment->cals_list->size - 1, meta_alignment->cals_list);
    seed_region_t *first_seed = linked_list_get_first(first_cal->sr_list);
    seed_region_t *last_seed  = linked_list_get_last(last_cal->sr_list);
    cigar_code_t *cigar_code = meta_alignment->cigar_code;    
    char *cigar_str = new_cigar_code_string(cigar_code);
    int cigar_len = strlen(cigar_str);

    simple_a = &simple_alignment[i];

    if (meta_alignment->cigar_left !=  NULL) {
      //printf("LEFT CIGAR: %s\n", new_cigar_code_string(meta_alignment->cigar_left));
      simple_a->gap_start = 0;
    } else {
      simple_a->gap_start = first_seed->read_start;
    }

    if (meta_alignment->cigar_right !=  NULL) {
      //printf("RIGHT CIGAR: %s\n", new_cigar_code_string(meta_alignment->cigar_right));
      simple_a->gap_end = 0;
    } else {
      simple_a->gap_end = seq_size - last_seed->read_end - 1;
    }

    simple_a->map_strand = first_cal->strand;
    simple_a->map_chromosome = first_cal->chromosome_id;
    simple_a->map_start = first_cal->start;
    simple_a->map_distance = cigar_code->distance;
    simple_a->cigar_len = cigar_len;
    
    //printf(" [%i:%lu] INSERT CIGAR(%i): %s\n", simple_a->map_chromosome, 
    //	   simple_a->map_start, cigar_len, cigar_str);

    memcpy(&cigar_buffer[tot_len], cigar_str, cigar_len);
    tot_len += cigar_len;

    if (tot_len >= max_len) { 
      max_len = max_len * 2;
      cigar_buffer = realloc(cigar_buffer, max_len); 
    }

  }

  fwrite(simple_alignment, sizeof(simple_alignment_t), num_items, fd);
  fwrite(cigar_buffer, sizeof(char), tot_len, fd);  

  free(cigar_buffer);
}

void file_write_cals(fastq_read_t *fq_read, array_list_t *items, FILE *fd) {
  //size_t head_size = strlen(fq_read->id);
  //size_t seq_size  = fq_read->length;
  size_t num_items = array_list_size(items);
  if (!num_items) { return; }
  //Write binary file 
  //[type][size head][size seq][num items][HEAD][SEQUENCE][QUALITY][CAL 0][CAL n]
  //size_t items_sizes[3] = {head_size, seq_size, num_items};

  //[size head][size seq][num items]
  //fwrite(items_sizes, sizeof(size_t), 3, fd);
  
  //[HEAD][SEQUENCE][QUALITY]
  //size_t total_size = head_size + 2*seq_size;
  //char *buffer = (char *)malloc(sizeof(char)*total_size);
  
  //memcpy(buffer, fq_read->id, head_size);
  //memcpy(&buffer[head_size], fq_read->sequence, seq_size);
  //memcpy(&buffer[head_size + seq_size], fq_read->quality, seq_size);
  //fwrite(buffer, sizeof(char), total_size, fd);
  
  bwt_anchor_t bwt_anchor[num_items];
  memset(bwt_anchor, 0, sizeof(bwt_anchor_t)*num_items);

  for (int i = 0; i < num_items; i++) {
    cal_t *cal = array_list_get(i, items);
    //cal_print(cal);
    bwt_anchor[i].strand     = cal->strand;
    bwt_anchor[i].chromosome = cal->chromosome_id - 1;
    bwt_anchor[i].start      = cal->start;
    bwt_anchor[i].end        = cal->end;//cal->start + (cal->end - cal->start);
    seed_region_t *seed = linked_list_get_first(cal->sr_list);
    if (seed->read_start == 0) {
      bwt_anchor[i].type = FORWARD_ANCHOR;
    } else {
      bwt_anchor[i].type = BACKWARD_ANCHOR;
    }
  }

  fwrite(bwt_anchor, sizeof(bwt_anchor_t), num_items, fd);
  
  //free(buffer);  
  
}

void file_write_fastq_read(fastq_read_t *fq_read, size_t num_items, FILE *fd) {
  size_t head_size = strlen(fq_read->id);
  size_t seq_size  = fq_read->length;

  //Write binary file 
  //[type][size head][size seq][num items][HEAD][SEQUENCE][QUALITY][CAL 0][CAL n]
  size_t items_sizes[3] = {head_size, seq_size, num_items};
  //printf("NUM items %i\n", num_items);
  //printf("Insert-id  (%i): %s\n", head_size, fq_read->id);
  //printf("Insert-seq (%i): %s\n", seq_size, fq_read->sequence);
  //printf("Insert-qua (%i): %s\n", seq_size, fq_read->quality);

  //[size head][size seq][num items]
  fwrite(items_sizes, sizeof(size_t), 3, fd);
  //fwrite(&head_size, sizeof(size_t), 1, f_sa);
  //fwrite(&seq_size,  sizeof(size_t), 1, f_sa);
  //fwrite(&num_items, sizeof(size_t), 1, f_sa);
  
  //[HEAD][SEQUENCE][QUALITY]
  size_t total_size = head_size + 2*seq_size;
  char *buffer = (char *)malloc(sizeof(char)*total_size);
  
  memcpy(buffer, fq_read->id, head_size);
  memcpy(&buffer[head_size], fq_read->sequence, seq_size);
  memcpy(&buffer[head_size + seq_size], fq_read->quality, seq_size);

  fwrite(buffer, sizeof(char), total_size, fd);

  free(buffer);

}

void file_write_items(fastq_read_t *fq_read, array_list_t *items, 
		      unsigned char data_type, FILE *fd1, FILE *fd2,
		      int mode) {
  FILE *fd;

  if (mode == 0) {
    fd = fd1;
  } else {
    fd = fd2;
  }

  fwrite(&data_type, sizeof(unsigned char), 1, fd);
  file_write_fastq_read(fq_read, array_list_size(items), fd);

  if (data_type == CAL_TYPE) {
    //printf("======= INSERT CAL ITEMS (%i) ========\n", array_list_size(items));
    file_write_cals(fq_read, items, fd);
  } else if (data_type == META_ALIGNMENT_TYPE) {
    //printf("======= INSERT META ALIGNMENTS ITEMS (%i) ========\n", array_list_size(items));
    file_write_meta_alignments(fq_read, items, fd);    
  } else {
    //printf("======= INSERT ALIGNMENTS ITEMS (%i) ========\n", array_list_size(items));
    file_write_alignments(fq_read, items, fd);    
  }
  
}

//=================================================================
//File SA Functions

void sa_file_write_alignments(fastq_read_t *fq_read, array_list_t *items, FILE *fd) {
  //size_t head_size = strlen(fq_read->id);
  size_t seq_size  = fq_read->length;
  size_t num_items = array_list_size(items);
  if (!num_items) { 
    return;
  } 

  size_t max_len = num_items * 1024;
  char *cigar_buffer = (char *)calloc(max_len, sizeof(char));
  size_t tot_len = 0;

  simple_alignment_t simple_alignment[num_items];
  simple_alignment_t *simple_a;
  //unsigned char type = MENTA_TYPE;
  //fwrite(type, sizeof(unsigned char), 1, fd);

  memset(simple_alignment, 0, sizeof(simple_alignment_t)*num_items);
  //[type][size head][size seq][num items][HEAD][SEQUENCE][QUALITY][CAL 0][CAL n]
  for (int i = 0; i < num_items; i++) {
    sa_alignment_t *sa_alignment = array_list_get(i, items);
    cal_t *first_cal = array_list_get(0, sa_alignment->cals_list);
    cal_t *last_cal = array_list_get(sa_alignment->cals_list->size - 1, sa_alignment->cals_list);
    seed_region_t *first_seed = linked_list_get_first(first_cal->sr_list);
    seed_region_t *last_seed  = linked_list_get_last(last_cal->sr_list);
    cigar_code_t *cigar_code = sa_alignment->c_final;    
    char *cigar_str = new_cigar_code_string(cigar_code);
    if (cigar_str == NULL) { 
      cigar_str = (char *)calloc(1, sizeof(char));
    }
    int cigar_len = strlen(cigar_str);

    simple_a = &simple_alignment[i];

    simple_a->gap_start = first_seed->read_start;
    simple_a->gap_end = seq_size - last_seed->read_end - 1;


    simple_a->map_strand = first_cal->strand;
    simple_a->map_chromosome = first_cal->chromosome_id;
    simple_a->map_start = first_cal->start;
    simple_a->map_distance = cigar_code->distance;
    simple_a->cigar_len = cigar_len;
    
    //printf("==== Write to file ====\n");
    //cal_print(first_cal);
    //cal_print(last_cal);
    // printf("==== Write to file ====\n");
    //printf(" [%i:%lu] INSERT CIGAR(%i): %s\n", simple_a->map_chromosome, 
    //	   simple_a->map_start, cigar_len, cigar_str);

    memcpy(&cigar_buffer[tot_len], cigar_str, cigar_len);
    tot_len += cigar_len;

    if (tot_len >= max_len) { 
      max_len = max_len * 2;
      cigar_buffer = realloc(cigar_buffer, max_len); 
    }

  }

  fwrite(simple_alignment, sizeof(simple_alignment_t), num_items, fd);
  fwrite(cigar_buffer, sizeof(char), tot_len, fd);  

  free(cigar_buffer);

}

//-------------------------------------------------------------------------------

int sa_file_read_alignments(size_t num_items, array_list_t *list, 
			    fastq_read_t *fq_read, FILE *fd) {

  if (!num_items) { return 0; }

  simple_alignment_t simple_alignment[num_items];
  simple_alignment_t *simple_a;
  int bytes;

  bytes = fread(simple_alignment, sizeof(simple_alignment_t), num_items, fd);
  if (!bytes) { LOG_FATAL("Corrupt file\n"); }
  
  size_t cigar_tot_len = 0;
  for (int i = 0; i < num_items; i++) {
    simple_a = &simple_alignment[i];
    //printf("ITEM %i: (%i)[%i:%lu] [%i-%i]\n", i, simple_a->map_strand, simple_a->map_chromosome,
    //	   simple_a->map_start, simple_a->gap_start, simple_a->gap_end);
    cigar_tot_len += simple_a->cigar_len;
  }
    
  char cigar_buffer[cigar_tot_len];
  bytes = fread(cigar_buffer, sizeof(char), cigar_tot_len, fd);
  if (!bytes) { LOG_FATAL("Corrupt file\n"); }

  char cigars_test[num_items][1024];
  size_t actual_read = 0;
  for (int i = 0; i < num_items; i++) {
    simple_a = &simple_alignment[i];
    memcpy(&cigars_test[i], &cigar_buffer[actual_read], simple_a->cigar_len);
    cigars_test[i][simple_a->cigar_len] = '\0';
    actual_read += simple_a->cigar_len;
    //printf("CIGAR %i: %s\n", i, cigars_test[i]);
    size_t map_len = fq_read->length - simple_a->gap_start - simple_a->gap_end;
    size_t map_genome_len = 0;

    cigar_code_t *cc = cigar_code_new_by_string(cigars_test[i]);
    array_list_t *list_aux = array_list_new(5, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    sa_alignment_t *sa_alignment = sa_alignment_new(list_aux);    

    sa_alignment->c_final = cigar_code_new();

    for (int m = 0; m < array_list_size(cc->ops); m++) {
      cigar_op_t *op = array_list_get(m, cc->ops);
      cigar_code_append_new_op(op->number, op->name, sa_alignment->c_final);

      if (op->name == 'M' || op->name == 'D' || op->name == 'N') {
	map_genome_len += op->number;
      }

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


    cc->distance = simple_a->map_distance;
    cal->info = cc;

    //printf("Cal & Cigar Ok, Insert list\n");


    sa_alignment->c_final->distance = simple_a->map_distance;
    
    array_list_insert(cal, sa_alignment->cals_list);
    array_list_insert(sa_alignment, list);

  }

  return 0;

}

//-------------------------------------------------------------------------------

void sa_file_write_items(fastq_read_t *fq_read, array_list_t *items, FILE *fd) {
  file_write_fastq_read(fq_read, array_list_size(items), fd);
  sa_file_write_alignments(fq_read, items, fd);  
}

//=================================================================


/*
void insert_file_item (fastq_read_t *fq_read, array_list_t *items, FILE *f_sa) {
  size_t head_size = strlen(fq_read->id);
  size_t seq_size  = fq_read->length;
  size_t num_items = array_list_size(items);

  //Write binary file 
  //[type][size head][size seq][num items][HEAD][SEQUENCE][QUALITY][CAL 0][CAL n]
  size_t items_sizes[3] = {head_size, seq_size, num_items};
  //printf("NUM items %i\n", num_items);
  //printf("Insert-id  (%i): %s\n", head_size, fq_read->id);
  //printf("Insert-seq (%i): %s\n", seq_size, fq_read->sequence);
  //printf("Insert-qua (%i): %s\n", seq_size, fq_read->quality);

  //[size head][size seq][num items]
  fwrite(items_sizes, sizeof(size_t), 3, f_sa);
  //fwrite(&head_size, sizeof(size_t), 1, f_sa);
  //fwrite(&seq_size,  sizeof(size_t), 1, f_sa);
  //fwrite(&num_items, sizeof(size_t), 1, f_sa);
  
  //[HEAD][SEQUENCE][QUALITY]
  size_t total_size = head_size + 2*seq_size;
  char *buffer = (char *)malloc(sizeof(char)*total_size);
  

  memcpy(buffer, fq_read->id, head_size);
  memcpy(&buffer[head_size], fq_read->sequence, seq_size);
  memcpy(&buffer[head_size + seq_size], fq_read->quality, seq_size);
  fwrite(buffer, sizeof(char), total_size, f_sa);
  
  bwt_anchor_t bwt_anchor[num_items];
  memset(bwt_anchor, 0, sizeof(bwt_anchor_t)*num_items);

  for (int i = 0; i < num_items; i++) {
    cal_t *cal = array_list_get(i, items);
    bwt_anchor[i].strand     = cal->strand;
    bwt_anchor[i].chromosome = cal->chromosome_id - 1;
    bwt_anchor[i].start      = cal->start;
    bwt_anchor[i].end        = cal->start + (cal->end - cal->start + 1);
    seed_region_t *seed = linked_list_get_first(cal->sr_list);
    if (seed->read_start == 0) {
      bwt_anchor[i].type = FORWARD_ANCHOR;
    } else {
      bwt_anchor[i].type = BACKWARD_ANCHOR;
    }
  }

  fwrite(bwt_anchor, sizeof(bwt_anchor_t), num_items, f_sa);
  
  free(buffer);

}

void insert_file_item_2 (fastq_read_t *fq_read, array_list_t *items, FILE *f_hc) {
  size_t head_size = strlen(fq_read->id);
  size_t seq_size  = fq_read->length;
  size_t num_items = array_list_size(items);
  size_t max_len = num_items * 1024;
  char *cigar_buffer = (char *)calloc(max_len, sizeof(char));
  size_t tot_len = 0;

  simple_alignment_t simple_alignment[num_items];
  simple_alignment_t *simple_a;
  memset(simple_alignment, 0, sizeof(simple_alignment_t)*num_items);
  for (int i = 0; i < num_items; i++) {
    meta_alignment_t *meta_alignment = array_list_get(i, items);
    cal_t *first_cal = array_list_get(0, meta_alignment->cals_list);
    cal_t *last_cal = array_list_get(meta_alignment->cals_list->size - 1, meta_alignment->cals_list);
    seed_region_t *first_seed = linked_list_get_first(first_cal->sr_list);
    seed_region_t *last_seed  = linked_list_get_last(last_cal->sr_list);
    cigar_code_t *cigar_code = meta_alignment->cigar_code;    
    char *cigar_str = new_cigar_code_string(cigar_code);
    int cigar_len = strlen(cigar_str);

    simple_a = &simple_alignment[i];

    if (meta_alignment->cigar_left !=  NULL) {
      //printf("LEFT CIGAR: %s\n", new_cigar_code_string(meta_alignment->cigar_left));
      simple_a->gap_start = 0;
    } else {
      simple_a->gap_start = first_seed->read_start;
    }

    if (meta_alignment->cigar_right !=  NULL) {
      //printf("RIGHT CIGAR: %s\n", new_cigar_code_string(meta_alignment->cigar_right));
      simple_a->gap_end = 0;
    } else {
      simple_a->gap_end = seq_size - last_seed->read_end - 1;
    }

    simple_a->map_strand = first_cal->strand;
    simple_a->map_chromosome = first_cal->chromosome_id;
    simple_a->map_start = first_cal->start;
    simple_a->map_distance = cigar_code->distance;
    simple_a->cigar_len = cigar_len;
    
    //printf(" [%i:%lu] INSERT CIGAR(%i): %s\n", simple_a->map_chromosome, 
    //	   simple_a->map_start, cigar_len, cigar_str);

    memcpy(&cigar_buffer[tot_len], cigar_str, cigar_len);
    tot_len += cigar_len;

    if (tot_len >= max_len) { 
      max_len = max_len * 2;
      cigar_buffer = realloc(cigar_buffer, max_len); 
    }

  }

  //Write binary file 
  //[size head][size seq][num items][HEAD][SEQUENCE][QUALITY][CAL 0][CAL n]
  size_t items_sizes[3] = {head_size, seq_size, num_items};

  //[size head][size seq][num items]
  fwrite(items_sizes, sizeof(size_t), 3, f_hc);
  
  //[HEAD][SEQUENCE][QUALITY]
  size_t total_size = head_size + 2*seq_size;
  char *buffer = (char *)malloc(sizeof(char)*total_size);
  
  memcpy(buffer, fq_read->id, head_size);
  memcpy(&buffer[head_size], fq_read->sequence, seq_size);
  memcpy(&buffer[head_size + seq_size], fq_read->quality, seq_size);

  fwrite(buffer, sizeof(char), total_size, f_hc);  
  fwrite(simple_alignment, sizeof(simple_alignment_t), num_items, f_hc);
  fwrite(cigar_buffer, sizeof(char), tot_len, f_hc);  

  free(buffer);
  free(cigar_buffer);
}

/*
void buffer_item_insert_new_item(fastq_read_t *fq_read, 
				 linked_list_t *items_list, 
				 void *data,
				 int type_items,
				 linked_list_t *buffer, 
				 linked_list_t *buffer_hc,
				 int phase) {
  //printf("INSERT TO BUFFER\n");

  buffer_item_t *buffer_item = buffer_item_complete_new(fq_read, items_list, data);
  array_list_set_flag(type_items, buffer_item->items_list);
  
  //Select list to insert
  cal_t *cal_prev;
  cal_t *cal_next;
  int insert_hc = 1;
  int end;

  if (phase == 1) {
    if (type_items == BITEM_SINGLE_ANCHORS) {
      linked_list_insert(buffer_item, buffer_hc);
    }
    return;
  } else {
    linked_list_insert(buffer_item, buffer);
    return;
  }

  if (type_items == BITEM_CALS) {    
    end = 0;
    for (int i = 0; i < array_list_size(items_list); i++) {
      array_list_t *fusion_list = array_list_get(i, items_list);
      for (int j = 0; j < array_list_size(fusion_list); j++) {
	cal_prev = array_list_get(j, fusion_list);
	if (cal_prev == NULL) {
	  //printf("CAL NULL\n");
	  insert_hc = 0;
	  end = 1;
	  break;
	} else {
	  seed_region_t *s_first = linked_list_get_first(cal_prev->sr_list);
	  seed_region_t *s_last = linked_list_get_last(cal_prev->sr_list);
	  
	  assert(s_first);
	  assert(s_last);
	  
	  //printf("s_first->read_start = %i, fq_read->length - s_last->read_end = %i\n",
	  //     s_first->read_start, fq_read->length - s_last->read_end);
	  if (s_first->read_start > 16 ||
	      fq_read->length - s_last->read_end > 16) {
	    //printf("INSERT IN HC = 0\n");
	    insert_hc = 0;
	    end = 1;
	    break;
	  }
	}
      }
      if (end) { break; }
    }
  } if (type_items == BITEM_SINGLE_ANCHORS) {    
    for (int j = 0; j < array_list_size(items_list); j++) {
      cal_prev = array_list_get(j, items_list);
      if (cal_prev == NULL) {
	//printf("CAL NULL\n");
	insert_hc = 0;
	break;
      } else {
	seed_region_t *s_first = linked_list_get_first(cal_prev->sr_list);
	seed_region_t *s_last = linked_list_get_last(cal_prev->sr_list);
	
	assert(s_first);
	assert(s_last);
	
	//printf("s_first->read_start = %i, fq_read->length - s_last->read_end = %i\n",
	//     s_first->read_start, fq_read->length - s_last->read_end);
	if (s_first->read_start > 16 ||
	    fq_read->length - s_last->read_end > 16) {
	  //printf("INSERT IN HC = 0\n");
	  insert_hc = 0;
	  break;
	}
      }
    }  
  } else {
    //printf("Insert buffer meta alignments\n");
    for (int i = 0; i < array_list_size(items_list); i++) {
      meta_alignment_t *meta_alignment = array_list_get(i, items_list);
      cal_prev = array_list_get(0, meta_alignment->cals_list);
      cal_next = array_list_get(array_list_size(meta_alignment->cals_list) - 1, 
				meta_alignment->cals_list);

      if (cal_prev == NULL || cal_next == NULL) {
	insert_hc = 0;
	//printf("NULL\n");
	break;
      }

      seed_region_t *s_first = linked_list_get_first(cal_prev->sr_list);
      seed_region_t *s_last = linked_list_get_last(cal_next->sr_list);

      assert(s_first);
      assert(s_last);

      //cal_print(cal_prev);
      //cal_print(cal_next);

      if (s_first->read_start > 16 ||
	  fq_read->length - s_last->read_end > 16) {
	insert_hc = 0;
	break;
      }
    }
  }

  if (!insert_hc) {
    //printf("INSERT IN HC = 0\n");
    linked_list_insert(buffer_item, buffer);
  } else {
    //printf("INSERT IN HC = 1\n");
    linked_list_insert(buffer_item, buffer_hc);
  }
  
}


void buffer_item_free(buffer_item_t *buffer_item) {
  array_list_free(buffer_item->items_list, NULL);
  free(buffer_item);
  
}
*/

bs_context_t *bs_context_new(size_t num_reads) {
  bs_context_t *p = (bs_context_t*) calloc(1, sizeof(bs_context_t));

  p->context_CpG = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  p->context_CHG = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  p->context_CHH = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  p->context_MUT = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  /*
  p->context_bs_CpG = array_list_bs_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  p->context_bs_CHG = array_list_bs_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  p->context_bs_CHH = array_list_bs_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  p->context_bs_MUT = array_list_bs_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  */
  return p;
}

//------------------------------------------------------------------------------------

void bs_context_free(bs_context_t *p) {
  if (p) {
    // free lists...
    /*
    if (p->context_CpG != NULL) array_list_free(p->context_CpG, NULL);
    if (p->context_CHG != NULL) array_list_free(p->context_CHG, NULL);
    if (p->context_CHH != NULL) array_list_free(p->context_CHH, NULL);
    if (p->context_MUT != NULL) array_list_free(p->context_MUT, NULL);
    */
    /*
    if (p->context_bs_CpG) array_list_bs_free(p->context_bs_CpG, NULL);
    if (p->context_bs_CHG) array_list_bs_free(p->context_bs_CHG, NULL);
    if (p->context_bs_CHH) array_list_bs_free(p->context_bs_CHH, NULL);
    if (p->context_bs_MUT) array_list_bs_free(p->context_bs_MUT, NULL);
    */
    free(p);
  }
}

//------------------------------------------------------------------------------------

void bs_context_init(bs_context_t *bs_context, size_t num_reads) {
  /*
  bs_context->context_CpG = (array_list_t **) calloc(num_reads, sizeof(array_list_t*));
  bs_context->context_CHG = (array_list_t **) calloc(num_reads, sizeof(array_list_t*));
  bs_context->context_CHH = (array_list_t **) calloc(num_reads, sizeof(array_list_t*));
  bs_context->context_MUT = (array_list_t **) calloc(num_reads, sizeof(array_list_t*));

  for (size_t i = 0; i < num_reads; i++) {
    bs_context->context_CpG[i] = array_list_new(500,
						1.25f,
						COLLECTION_MODE_ASYNCHRONIZED);
    bs_context->context_CHG[i] = array_list_new(500,
						1.25f,
						COLLECTION_MODE_ASYNCHRONIZED);
    bs_context->context_CHH[i] = array_list_new(500,
						1.25f,
						COLLECTION_MODE_ASYNCHRONIZED);
    bs_context->context_MUT[i] = array_list_new(500,
						1.25f,
						COLLECTION_MODE_ASYNCHRONIZED);
  }
  */
  
  bs_context->context_CpG = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  bs_context->context_CHG = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  bs_context->context_CHH = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  bs_context->context_MUT = array_list_new(num_reads, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  
  /*
  bs_context->context_CpG = array_list_new(num_reads, 2, COLLECTION_MODE_ASYNCHRONIZED);
  bs_context->context_CHG = array_list_new(num_reads, 2, COLLECTION_MODE_ASYNCHRONIZED);
  bs_context->context_CHH = array_list_new(num_reads, 2, COLLECTION_MODE_ASYNCHRONIZED);
  bs_context->context_MUT = array_list_new(num_reads, 2, COLLECTION_MODE_ASYNCHRONIZED);
  */
  bs_context->CpG_methyl   = 0;
  bs_context->CpG_unmethyl = 0;
  bs_context->CHG_methyl   = 0;
  bs_context->CHG_unmethyl = 0;
  bs_context->CHH_methyl   = 0;
  bs_context->CHH_unmethyl = 0;
  bs_context->MUT_methyl   = 0;
  bs_context->num_bases    = 0;
}

//------------------------------------------------------------------------------------
