#include "workflow_functions.h"

int limit_of_reads = 0;
int n_insert = 0;
int tot_reads2 = 0;

//====================================================================================

wf_input_t *wf_input_new(fastq_batch_reader_input_t *fq_reader_input,
                         batch_t *batch) {

  wf_input_t *wfi = (wf_input_t *) calloc(1, sizeof(wf_input_t));
  wfi->fq_reader_input = fq_reader_input;
  wfi->batch = batch;

  return wfi;
}

wf_input_buffer_t *wf_input_buffer_new(linked_list_t *buffer,
				       batch_t *batch) {
  wf_input_buffer_t *wfi = (wf_input_buffer_t *) calloc(1, sizeof(wf_input_buffer_t));
  wfi->buffer = buffer;
  wfi->batch = batch;

  return wfi;
}

void wf_input_free(wf_input_t *wfi) {
  if (wfi) free(wfi);
}

void wf_input_file_free(wf_input_file_t *wfi) {
  if (wfi) free(wfi);
}

wf_input_file_t *wf_input_file_new(FILE *fd,
				   batch_t *batch) {
  wf_input_file_t *wfi = (wf_input_file_t *) calloc(1, sizeof(wf_input_file_t));
  wfi->file = fd;
  wfi->batch = batch;

  return wfi;
}


//--------------------------------------------------------------------
// workflow producer
//--------------------------------------------------------------------

void *fastq_reader(void *input) {
     struct timeval start, end;
     double time;

     //if (time_on) { start_timer(start); }

     wf_input_t *wf_input = (wf_input_t *) input;
     batch_t *new_batch = NULL;
     batch_t *batch = wf_input->batch;
     fastq_batch_reader_input_t *fq_reader_input = wf_input->fq_reader_input;
     array_list_t *reads = array_list_new(10000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

     if (fq_reader_input->gzip) {
       //Gzip fastq file
       if (fq_reader_input->flags == SINGLE_END_MODE) {
	 fastq_gzread_bytes_se(reads, fq_reader_input->batch_size, fq_reader_input->fq_gzip_file1);
       } else {
	 //printf("Gzip Reader for pair-end not implemented\n");;
	 fastq_gzread_bytes_pe(reads, fq_reader_input->batch_size, fq_reader_input->fq_gzip_file1, fq_reader_input->fq_gzip_file2);
	 //fastq_fread_bytes_aligner_pe(reads, fq_reader_input->batch_size, 
	 //		      fq_reader_input->fq_gzip_file1, fq_reader_input->fq_gzip_file2);
       }
     } else {
       //Fastq file
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
	  mapping_batch_t *mapping_batch = mapping_batch_new(reads, 
							     batch->pair_input->pair_mng);

	  new_batch = batch_new(batch->bwt_input, batch->region_input, batch->cal_input, 
				batch->pair_input, batch->preprocess_rna, batch->sw_input, batch->writer_input, 
				batch->mapping_mode, mapping_batch);
     }

     //if (time_on) { stop_timer(start, end, time); timing_add(time, FASTQ_READER, timing); }

     return new_batch;
}


/*
fastq_read_t *file_fastq_read_new(size_t *num_items, FILE *fd) {

  size_t sizes_to_read[3], head_len, seq_len;

  head_len  = sizes_to_read[0];
  seq_len   = sizes_to_read[1];
  *num_items = sizes_to_read[2];
  
  int bytes = fread(sizes_to_read, sizeof(size_t), 3, fd);
  if (!bytes) { return NULL; }
  
  int tot_size = head_len + 2*seq_len;
  buffer = (char *)calloc(tot_size + 1, sizeof(char));
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

int file_cal_fill(size_t num_items, array_list_t *list, FILE *fd) {

  if (!num_items) { return 0; }
  
  bwt_anchor_t bwt_anchors[num_items];
  bytes = fread(bwt_anchors, sizeof(bwt_anchor_t), num_items, fd);
  if (!bytes) { LOG_FATAL("Corrupt file\n"); }
  
  for (int i = 0; i < num_items; i++) {
    //printf("[%i:%lu-%lu]\n", bwt_anchors[i].chromosome, bwt_anchors[i].start, bwt_anchors[i].end);
    size_t seed_size = bwt_anchors[i].end - bwt_anchors[i].start;
    cal_t *cal;
    if (bwt_anchors[i].type == FORWARD_ANCHOR) {
      cal = convert_bwt_anchor_to_CAL(&bwt_anchors[i], 0, seed_size);
    } else {
      cal = convert_bwt_anchor_to_CAL(&bwt_anchors[i], fq_read->length - seed_size - 1, fq_read->length - 1);
    }
    array_list_insert(cal, list); 
  }  
  
  return 0;

}

int file_meta_alignment_fill(size_t num_items, array_list_t *list, FILE *fd) {

  if (!num_items) { return 0; }

  simple_alignment_t simple_alignment[num_items];
  simple_alignment_t *simple_a;
  
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
					      simple_a->map_start + map_len,
					      0);
    
    //printf("Exit with seed [%i:%i]\n", s_region->read_start, s_region->read_end);
    
    linked_list_t *sr_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    //s_region->info = cigar_code_new_by_string(cigars_test[i]);
    linked_list_insert(s_region, sr_list);
    
    cal_t *cal = cal_new(simple_a->map_chromosome, 
			 simple_a->map_strand,
			 simple_a->map_start,
			 simple_a->map_start + map_len,
			 1,
			 sr_list,
			 linked_list_new(COLLECTION_MODE_ASYNCHRONIZED));
    cal->info = cigar_code_new_by_string(cigars_test[i]);
    
    meta_alignment_t *meta_alignment = meta_alignment_new();
    array_list_insert(cal, meta_alignment->cals_list);
    array_list_insert(meta_alignment, list);
  }

  return 0;
}

int file_alignment_fill(size_t num_items, array_list_t *list, 
			fastq_read_t *fq_read, FILE *fd) {
  if (!num_items) { return 0; }
  
  alignment_aux_t alignments_aux[num_items];
  alignment_aux_t *alignment_a;
  
  bytes = fread(alignments_aux, sizeof(alignment_aux_t), num_items, fd);
  if (!bytes) { LOG_FATAL("Corrupt file\n"); }

  size_t cigar_tot_len = 0;
  for (int i = 0; i < num_items; i++) {
    alignment_a = &simple_alignment[i];
    //printf("ITEM %i: (%i)[%i:%lu] [%i-%i]\n", i, simple_a->map_strand, simple_a->map_chromosome,
    //     simple_a->map_start, simple_a->gap_start, simple_a->gap_end);
    cigar_tot_len += alignmment_a->cigar_len + alignment_a->optional_field_length;
  }

  char cigars_test[num_items][1024];
  char optional_fields[num_items][1024];
  size_t actual_read = 0;
  for (int i = 0; i < num_items; i++) {
    alignment_a = &alignments_aux[i];
    memcpy(&cigars_test[i], &cigar_buffer[actual_read], alignment_a->cigar_len);
    cigars_test[i][alignment_a->cigar_len] = '\0';
    actual_read += simple_a->cigar_len;
    
    char op;
    char op_value[1024];
    int c = 0;
    int hc_start = 0, hc_end;
    for (int j = 0; j < alignment_a->cigar_len; j++) {
      op = cigars_test[j];
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

    if (cigars_test[alignment_a->cigar_len - 1] == 'H') {
      for (int j = alignment_a->cigar_len - 2; j >= 0; j--) {
	op = cigars_test[j];
	if (op < 58) {
	  op_value[c++] = op;
	} else {
	  op_value[c] = '\0';
	  int len = strlen(op_value);
	  char op_val_aux[len];
	  int pos = len - 1;
	  for (int j = 0; j < len; j++) {	    
	    op_val_aux[j] = op_value[pos - j];
	  } 
	  hc_end = atoi(op_val_aux);
	  break;
	}
      }
    }

    memcpy(&optional_fields[i], &cigar_buffer[actual_read], alignment_a->optional_fields_length);
    optional_fields[i][alignment_a->optional_fields_length] = '0';
    actual_read += alignment_a->optional_fields_length;

    int header_len = strlen(fq_read->id);
    char header_id[header_len + 1];
    get_to_first_blank(fq_read->id, header_len, header_id);
    //char *header_match = (char *)malloc(sizeof(char)*header_len);
    //memcpy(header_match, header_id, header_len);

    int len_read = fq_read->length - (hc_start + hc_end);
    char *quality = (char *) calloc (len_read + 1, sizeof(char));
    strncpy(quality, fq_read->quality + hc_start, len_read);
    char *query = (char *) calloc (len_read + 1, sizeof(char));
    strncpy(query, fq_read->query + hc_start, len_read);

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
			      strdup(optional_fields[i]), 
			      alignment);
    
    array_list_insert(alignment, list);
  }  

  return 0;

}
*/
void *file_reader(void *input) {
  wf_input_file_t *wf_input = (wf_input_file_t *) input;
  FILE *fd = wf_input->file;
  batch_t *batch = wf_input->batch;
  int pair_mode = batch->pair_input->pair_mng->pair_mode;

  const int MAX_READS = 100;
  int num_reads = 0;
  batch_t *new_batch = NULL;

  size_t tot_size;
  size_t num_items;
  char *buffer, *id, *sequence, *quality;
  size_t bytes;
  unsigned char type;
  array_list_t *reads = array_list_new(MAX_READS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  mapping_batch_t *mapping_batch = mapping_batch_new_2(MAX_READS, 
						       reads,
						       batch->pair_input->pair_mng);  
  while (1) {
    //[type][size head][size seq][num items]
    bytes = fread(&type, sizeof(unsigned char), 1, fd);
    if (!bytes) { break; }
 
    fastq_read_t *fq_read = file_read_fastq_reads(&num_items, fd);
    if (fq_read == NULL) { break; }
    
    mapping_batch->mapping_lists[num_reads] = array_list_new(50,
							     1.25f, 
							     COLLECTION_MODE_ASYNCHRONIZED);
    //printf("(num items %i)\nID : %s\nSEQ: %s\nQUA: %s\n", num_items, fq_read->id, fq_read->sequence, fq_read->quality);

    array_list_insert(fq_read, reads);
    
    if (type == CAL_TYPE) {
      //printf("\tCal Report\n");
      file_read_cals(num_items, mapping_batch->mapping_lists[num_reads], 
		     fq_read, fd);      
      array_list_set_flag(BITEM_SINGLE_ANCHORS, 
			  mapping_batch->mapping_lists[num_reads]);
    } else if (type == META_ALIGNMENT_TYPE) {
      //printf("\tMeta Alignments Report\n");
      array_list_set_flag(BITEM_META_ALIGNMENTS, 
			  mapping_batch->mapping_lists[num_reads]);
      file_read_meta_alignments(num_items, mapping_batch->mapping_lists[num_reads], 
				fq_read, fd);            
    } else {
      //printf("\tAlignments Report\n");
      file_read_alignments(num_items, mapping_batch->mapping_lists[num_reads], 
			   fq_read, fd);
    }

    /*if (strcmp("@ENST00000496771@ENSG00000000003@processed_transcript@X@99887538@99891686@-1@KNOWN_518_447_1_0_0_0_4:0:0_3:0:0_3/1",
	       fq_read->id) == 0) {
      exit(-1);
      }*/
            
    num_reads++;
    if (num_reads >= MAX_READS) { break; }

  }

  //w2_r += num_reads;
  //printf("W2 Reads: %i\n", w2_r);

  if (num_reads) {
    mapping_batch->num_allocated_targets = num_reads;
    new_batch = batch_new(batch->bwt_input, batch->region_input, batch->cal_input, 
			  batch->pair_input, batch->preprocess_rna, batch->sw_input,
			  batch->writer_input, batch->mapping_mode, mapping_batch); 
  } else {
    //array_list_free(reads, NULL);
    mapping_batch_free(mapping_batch);
  }

  return new_batch;

}

void *file_reader_2(void *input) {
  wf_input_file_t *wf_input = (wf_input_file_t *) input;
  FILE *fd = wf_input->file;
  batch_t *batch = wf_input->batch;
  int pair_mode = batch->pair_input->pair_mng->pair_mode;
  
  const int MAX_READS = 100;
  int num_reads = 0;
  batch_t *new_batch = NULL;

  size_t sizes_to_read[3], head_len, seq_len, num_items;
  size_t tot_size;
  char *buffer, *id, *sequence, *quality;
  size_t bytes;
  unsigned char type;
  array_list_t *reads = array_list_new(MAX_READS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  mapping_batch_t *mapping_batch = mapping_batch_new_2(MAX_READS, 
						       reads,
						       batch->pair_input->pair_mng);
  
  while (1) {
    //[size head][size seq][num items]
    bytes = fread(&type, sizeof(unsigned char), 1, fd);
    if (!bytes) { break; }
 
    //fastq_read_t *fq_read = file_fastq_read_new(&num_items, fd);
    fastq_read_t *fq_read = file_read_fastq_reads(&num_items, fd);
    if (fq_read == NULL) { /*printf("fq NULL\n");*/ break; }
    //printf("(num items %i)\nID : %s\nSEQ: %s\nQUA: %s\n", num_items, fq_read->id, fq_read->sequence, fq_read->quality);

    array_list_insert(fq_read, reads);
    
    mapping_batch->mapping_lists[num_reads] = array_list_new(50,
							     1.25f, 
							     COLLECTION_MODE_ASYNCHRONIZED);
    if (type == CAL_TYPE) {
      //exit(-1);
      //printf("\tCal Report\n");
      file_read_cals(num_items, mapping_batch->mapping_lists[num_reads], 
		     fq_read, fd);      
    } else if (type == META_ALIGNMENT_TYPE) {
      //printf("\tMeta Alignments Report\n");
      file_read_meta_alignments(num_items, mapping_batch->mapping_lists[num_reads], 
				fq_read, fd);            
      array_list_set_flag(BITEM_META_ALIGNMENTS,
			  mapping_batch->mapping_lists[num_reads]);    
    } else {
      //exit(-1);
      //printf("\tAlignments Report\n");
      file_read_alignments(num_items, mapping_batch->mapping_lists[num_reads], 
			   fq_read, fd);
    }

    //printf("W3 file read %i\n", array_list_size(mapping_batch->mapping_lists[num_reads]));
    num_reads++;
    if (num_reads >= MAX_READS) { break; }

  }

  tot_reads2 += num_reads;
  //printf("W3 Reads: %i | %i\n", tot_reads2, num_reads);
  //w3_r += num_reads;
  //printf("W3 Reads: %i\n", w3_r);

  if (num_reads) {
    mapping_batch->num_allocated_targets = num_reads;
    new_batch = batch_new(batch->bwt_input, batch->region_input, batch->cal_input, 
			  batch->pair_input, batch->preprocess_rna, batch->sw_input,
			  batch->writer_input, batch->mapping_mode, mapping_batch); 
  } else {
    mapping_batch_free(mapping_batch);
  }

  return new_batch;

}

//--------------------------------------------------------------------
// workflow consumer
//--------------------------------------------------------------------

int search_hard_clipping(array_list_t *array_list);
void write_mapped_read(array_list_t *array_list, bam_file_t *bam_file);
void write_unmapped_read(fastq_read_t *fq_read, bam_file_t *bam_file);

//--------------------------------------------------------------------

int bam_writer(void *data) {
     struct timeval start, end;
     double time;

     //if (time_on) { start_timer(start); }

     batch_t *batch = (batch_t *) data;
     fastq_read_t *fq_read;
     array_list_t *array_list;
     size_t num_items;

     //bam1_t *bam1;
     //alignment_t *alig;

     mapping_batch_t *mapping_batch = (mapping_batch_t *) batch->mapping_batch;

     batch_writer_input_t *writer_input = batch->writer_input;
     bam_file_t *bam_file = writer_input->bam_file;     
     linked_list_t *linked_list = writer_input->list_p;
     size_t num_reads_b = array_list_size(mapping_batch->fq_batch);
     size_t num_mapped_reads = 0;
     size_t total_mappings = 0;
     unsigned char found_p1 = 0;
     unsigned char found_p2 = 0;
     int i = 0;

     extern size_t bwt_correct;
     extern size_t bwt_error;
     extern pthread_mutex_t bwt_mutex, mutex_sp;

     writer_input->total_batches++;

     extern size_t *histogram_sw;

     extern size_t num_reads_map;
     extern size_t num_reads;
     extern size_t tot_reads;
     
     free(mapping_batch->histogram_sw);
     //
     // DNA/RNA mode
     //
     for (size_t i = 0; i < num_reads_b; i++) {
       num_items = array_list_size(mapping_batch->mapping_lists[i]);
       total_mappings += num_items;
       fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->fq_batch);
       
       // mapped or not mapped ?	 
       if (num_items == 0) {
	 total_mappings++;
	 write_unmapped_read(fq_read, bam_file);
	 if (mapping_batch->mapping_lists[i]) {
	   array_list_free(mapping_batch->mapping_lists[i], NULL);
	 }	 
       } else {
	 num_mapped_reads++;
	 write_mapped_read(mapping_batch->mapping_lists[i], bam_file);
       }
     }
     
     //num_reads     += num_reads_b;
     //num_reads_map += num_mapped_reads;
 
     //fprintf(stderr, "TOTAL READS PROCESS: %lu\n", basic_st->total_reads);
     //if (basic_st->total_reads >= writer_input->limit_print) {
       //LOG_DEBUG_F("TOTAL READS PROCESS: %lu\n", basic_st->total_reads);
       //LOG_DEBUG_F("\tTotal Reads Mapped: %lu(%.2f%)\n", 
       //	   basic_st->num_mapped_reads, 
       //	   (float) (basic_st->num_mapped_reads*100)/(float)(basic_st->total_reads));
       //writer_input->limit_print += 1000000;

       /*
       fprintf(stderr, "TOTAL READS PROCESS: %lu\n", tot_reads);
       printf("\tTotal Reads Mapped: %lu(%.2f%)\n", 
	      basic_st->num_mapped_reads, 
	      (float) (basic_st->num_mapped_reads*100)/(float)(basic_st->total_reads));
       */

       //writer_input->limit_print += 100000;

       //extern size_t TOTAL_SW,
       //TOTAL_READS_PROCESS,
       //TOTAL_READS_SEEDING,
       //TOTAL_READS_SEEDING2,
       //TOTAL_READS_SA;
       

       //fprintf(stderr, "TOTAL READS PROCESS = %lu,\n TOTAL READS SEEDING x1 = %lu,\n TOTAL READS SEEDING x2 = %lu,\n TOTAL SW = %lu,\n TOTAL READS SINGLE ANCHOR FINAL = %lu\n\n",
       //      TOTAL_READS_PROCESS, TOTAL_READS_SEEDING, TOTAL_READS_SEEDING2, TOTAL_SW, TOTAL_READS_SA);
       
     //}
     
     //printf("Batch Write OK!\n");     
     
     if (mapping_batch) {
       mapping_batch_free(mapping_batch);
     }
     
     if (batch) batch_free(batch);
     
     basic_statistics_add(num_reads_b, num_mapped_reads, total_mappings, basic_st);
     
     //if (time_on) { stop_timer(start, end, time); timing_add(time, BAM_WRITER, timing); }
}

//--------------------------------------------------------------------
/*
int search_hard_clipping(array_list_t *array_list) {
  size_t num_items = array_list_size(array_list);
  alignment_t *alig;
  int found = 0;

  for (size_t j = 0; j < num_items; j++) {
    alig = (alignment_t *) array_list_get(j, array_list);
    if (alig->large_hard_clipping) {
      return 0;
    }
  }
  return 0;
}*/

//--------------------------------------------------------------------

void write_mapped_read(array_list_t *array_list, bam_file_t *bam_file) {
  size_t num_items = array_list_size(array_list);
  alignment_t *alig;
  bam1_t *bam1;
  for (size_t j = 0; j < num_items; j++) {
    alig = (alignment_t *) array_list_get(j, array_list);

    //printf("\t******** %i(%i)\n", j, num_items);
    //printf("is null alig->name %i\n", (alig->query_name == NULL));
    //printf("name = %s\n", alig->query_name);
    //printf("read = %s\n", alig->sequence);
    //printf("\t-----> %s\n", alig->cigar);
    LOG_DEBUG("writting bam..\n");
    //alignment_print(alig);

    if (alig != NULL) {
      bam1 = convert_to_bam(alig, 33);
      bam_fwrite(bam1, bam_file);
      bam_destroy1(bam1);	 
      alignment_free(alig);
    } else {
      LOG_FATAL_F("alig is NULL, num_items = %lu\n", num_items)
    }
    //printf("\t**************** %i(%i)\n", j, num_items);
  }
  if (array_list) { array_list_free(array_list, NULL); }
}

//--------------------------------------------------------------------

void write_unmapped_read(fastq_read_t *fq_read, bam_file_t *bam_file) {
  static char aux[1024] = "";
  alignment_t *alig;
  size_t header_len;
  char *id;
  bam1_t *bam1;

  // calculating cigar
  //sprintf(aux, "%luX", fq_read->length);	       
  alig = alignment_new();	       
  //header_len = strlen(fq_read->id);
  //id = (char *) malloc(sizeof(char) * (header_len + 1));
  //get_to_first_blank(fq_read->id, header_len, id);
  //free(fq_read->id);
  
  alignment_init_single_end(strdup(fq_read->id), fq_read->sequence, fq_read->quality,
			    0, -1, -1, /*strdup(aux)*/"", 0, 0, 0, 0, 0, NULL, alig);
  
  bam1 = convert_to_bam(alig, 33);
  bam_fwrite(bam1, bam_file);
  bam_destroy1(bam1);
	       
  alig->sequence = NULL;
  alig->quality = NULL;
  alig->cigar = NULL;
  alignment_free(alig);	       

  //printf("\tWRITE : read %i (%d items): unmapped...done !!\n", i, num_items);
  
}

//--------------------------------------------------------------------
// stage functions
//--------------------------------------------------------------------

int bwt_stage(void *data) {
  batch_t *batch = (batch_t *) data;

  if (batch->mapping_mode == DNA_MODE) {
    return apply_bwt(batch->bwt_input, batch);     
  } else {
    return apply_bwt_rna(batch->bwt_input, batch);     
  }
}

//--------------------------------------------------------------------

int bwt_stage_bs(void *data) {
  batch_t *batch = (batch_t *) data;

  //printf("Init BWT\n");
  return apply_bwt_bs(batch->bwt_input, batch);     
}

//--------------------------------------------------------------------

int seeding_stage(void *data) {
  batch_t *batch = (batch_t *) data;

  return apply_seeding(batch->region_input, batch);
}

//--------------------------------------------------------------------

int seeding_stage_bs(void *data) {
  batch_t *batch = (batch_t *) data;

  //printf("Init seeding\n");
  return apply_seeding_bs(batch->region_input, batch);
}

//--------------------------------------------------------------------

int cal_stage(void *data) {
  batch_t *batch = (batch_t *) data;

  return apply_caling_rna(batch->cal_input, batch);
}

//--------------------------------------------------------------------

int cal_stage_bs(void *data) {
  batch_t *batch = (batch_t *) data;

  //printf("Init CAL\n");
  return apply_caling_bs(batch->cal_input, batch);
}

//--------------------------------------------------------------------

int rna_preprocess_stage(void *data) {
  batch_t *batch = (batch_t *) data;

  return apply_preprocess_rna(batch->preprocess_rna, batch);  
}

//---------------------------------------------------------------------

int pre_pair_stage(void *data) {
  batch_t *batch = (batch_t *) data;

  return apply_pair(batch->pair_input, batch);
}

//---------------------------------------------------------------------

int pre_pair_stage_bs(void *data) {
  batch_t *batch = (batch_t *) data;

  //printf("Init pre_pair\n");
  return apply_pair(batch->pair_input, batch);
}

//--------------------------------------------------------------------

int sw_stage(void *data) {
  batch_t *batch = (batch_t *) data;
     
  if (batch->mapping_mode == RNA_MODE) {
    return apply_sw_rna(batch->sw_input, batch);
  } else {
    return apply_sw(batch->sw_input, batch);
  }
}

//--------------------------------------------------------------------

int sw_stage_bs(void *data) {
  batch_t *batch = (batch_t *) data;

  //printf("Init SW\n");
  return apply_sw_bs(batch->sw_input, batch);
}

//--------------------------------------------------------------------

int rna_last_stage(void *data) {
   batch_t *batch = (batch_t *) data;
   return apply_rna_last(batch->sw_input, batch);
}

//--------------------------------------------------------------------

int rna_last_hc_stage(void *data) {
   batch_t *batch = (batch_t *) data;
   return apply_rna_last_hc(batch->sw_input, batch);
}

//--------------------------------------------------------------------

int post_pair_stage(void *data) {
  batch_t *batch = (batch_t *) data;
  return prepare_alignments(batch->pair_input, batch);
}

//--------------------------------------------------------------------

int post_pair_stage_bs(void *data) {
  batch_t *batch = (batch_t *) data;

  return prepare_alignments_bs(batch->pair_input, batch);
}

//--------------------------------------------------------------------

int bs_status_stage(void *data) {
  batch_t *batch = (batch_t *) data;

  //printf("Init bs_status\n");
  return methylation_status_report(batch->sw_input, batch);
}

//--------------------------------------------------------------------

/*
void *buffer_reader(void *input) {
  wf_input_buffer_t *wf_input = (wf_input_t *) input;

  linked_list_t *buffer = wf_input->buffer;
  batch_t *batch = wf_input->batch;
  buffer_item_t *buffer_item;
  const int MAX_READS = 100;
  int num_reads = 0;
  batch_t *new_batch = NULL;

  if (linked_list_size(buffer) > 0) {
    array_list_t *reads = array_list_new(MAX_READS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    mapping_batch_t *mapping_batch = mapping_batch_new_2(MAX_READS, 
							 reads,
							 batch->pair_input->pair_mng);
    while (num_reads < MAX_READS) {
      buffer_item = linked_list_remove_last(buffer);
      if (buffer_item == NULL) { break; }
      fastq_read_t *read = buffer_item->read;
      array_list_insert(buffer_item->read, reads);
      mapping_batch->mapping_lists[num_reads] = array_list_new(50,
							       1.25f, 
							       COLLECTION_MODE_ASYNCHRONIZED);

      for (int i = 0; i < array_list_size(buffer_item->items_list); i++) {
	void *item = array_list_get(i, buffer_item->items_list);
	array_list_insert(item, mapping_batch->mapping_lists[num_reads]);
      }

      array_list_set_flag(array_list_get_flag(buffer_item->items_list),
			  mapping_batch->mapping_lists[num_reads]);
      num_reads++;
      //printf("TOTAL READS %i\n", num_reads);
      buffer_item_free(buffer_item);
    }
    
    mapping_batch->num_allocated_targets = num_reads;
    new_batch = batch_new(batch->bwt_input, batch->region_input, batch->cal_input, 
			  batch->pair_input, batch->preprocess_rna, batch->sw_input,
			  batch->writer_input, batch->mapping_mode, mapping_batch); 
  }


  return new_batch;
}
*/
