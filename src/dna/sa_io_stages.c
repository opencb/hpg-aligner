#include "sa_io_stages.h"

//====================================================================
// PRODUCER
//====================================================================

//--------------------------------------------------------------------
// sa fq reader
//--------------------------------------------------------------------

size_t fastq_fread_se_ex(array_list_t *reads, size_t num_reads, fastq_file_t *fq_file) {
  size_t count = 0;
  char *p;
  char header1[MAX_READ_ID_LENGTH];
  char sequence[MAX_READ_SEQUENCE_LENGTH];
  char header2[MAX_READ_ID_LENGTH];
  char qualities[MAX_READ_SEQUENCE_LENGTH];
  int header_length, sequence_length, quality_length;
  fastq_read_t *read;
  
  while (count < num_reads && fgets(header1, MAX_READ_ID_LENGTH, fq_file->fd) != NULL) {
    fgets(sequence, MAX_READ_SEQUENCE_LENGTH, fq_file->fd);
    fgets(header2, MAX_READ_ID_LENGTH, fq_file->fd);
    fgets(qualities, MAX_READ_SEQUENCE_LENGTH, fq_file->fd);
    
    header_length = strlen(header1);
    sequence_length = strlen(sequence);
    quality_length = strlen(qualities);
    
    // '\n' char is removed, but '\0' is left
    chomp_at(header1, header_length - 1);
    if ((p = strstr(header1, " ")) != NULL) {
      *p = 0;
    }
    chomp_at(sequence, sequence_length - 1);
    chomp_at(qualities, quality_length - 1);

    read = fastq_read_new(&header1[1], sequence, qualities);
    array_list_insert(read, reads);
    
    count++;
  }
  
  return count;
}

//--------------------------------------------------------------------

void *sa_fq_reader(void *input) {
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
    sa_mapping_batch_t *sa_mapping_batch = sa_mapping_batch_new(reads);
    sa_mapping_batch->bam_format = wf_input->bam_format;

    new_wf_batch = sa_wf_batch_new(curr_wf_batch->options,
				   curr_wf_batch->sa_index,
				   curr_wf_batch->writer_input, 
				   sa_mapping_batch);
  }

  return new_wf_batch;

}

//====================================================================
// CONSUMER
//====================================================================

#ifdef _VERBOSE
int num_dup_reads = 0;
int num_total_dup_reads = 0;
#endif

//--------------------------------------------------------------------
// SAM writer
//--------------------------------------------------------------------

void *write_sam_header(sa_genome3_t *genome, FILE *f) {
  fprintf(f, "@PG\tID:HPG-Aligner\tVN:2.0\n");
  for (int i = 0; i < genome->num_chroms; i++) {
    fprintf(f, "@SQ\tSN:%s\tLN:%lu\n", genome->chrom_names[i], genome->chrom_lengths[i] + 1);
  }
}

//--------------------------------------------------------------------

int sa_sam_writer(void *data) {
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;
  
  sa_mapping_batch_t *mapping_batch = (sa_mapping_batch_t *) wf_batch->mapping_batch;
  if (mapping_batch == NULL) {
    printf("bam_writer1: error, NULL mapping batch\n");
    return 0;
  }

  //  for (int i = 0; i < NUM_COUNTERS; i++) {
  //    counters[i] += mapping_batch->counters[i];
  //  }

  #ifdef _TIMING
  for (int i = 0; i < NUM_TIMING; i++) {
    func_times[i] += mapping_batch->func_times[i];
  }
  #endif

  int flag, pnext = 0, tlen = 0;
  char *rnext = "*", *optional_flags = "NM:i:3";

  fastq_read_t *read;
  array_list_t *read_list = mapping_batch->fq_reads;

  alignment_t *alig;
  array_list_t *mapping_list;
  FILE *out_file = (FILE *) wf_batch->writer_input->bam_file;

  sa_genome3_t *genome = wf_batch->sa_index->genome;

  size_t num_reads, num_mappings;
  num_reads = mapping_batch->num_reads;
  for (size_t i = 0; i < num_reads; i++) {
    read = (fastq_read_t *) array_list_get(i, read_list);
    mapping_list = mapping_batch->mapping_lists[i];
    num_mappings = array_list_size(mapping_list);
    #ifdef _VERBOSE
    if (num_mappings > 1) {
      num_dup_reads++;
      num_total_dup_reads += num_mappings;
    }
    #endif

    if (num_mappings > 0) {
      for (size_t j = 0; j < num_mappings; j++) {
	alig = (alignment_t *) array_list_get(j, mapping_list);
	flag = (alig->seq_strand ? 16 : 0);
	  
	fprintf(out_file, "%s\t%i\t%s\t%lu\t%i\t%s\t%s\t%lu\t%i\t%s\t%s\t%s\n", 
		alig->query_name,
		flag,
		genome->chrom_names[alig->chromosome],
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
      fprintf(out_file, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", 
	      read->id,
	      read->sequence,
	      read->quality
	      );
    }
    array_list_free(mapping_list, (void *) NULL);
  }

  // free memory
  sa_mapping_batch_free(mapping_batch);

  if (wf_batch) sa_wf_batch_free(wf_batch);

  return 0;
}


//--------------------------------------------------------------------
// BAM writer
//--------------------------------------------------------------------

bam_header_t *create_bam_header(sa_genome3_t *genome) {

  bam_header_t *bam_header = (bam_header_t *) calloc(1, sizeof(bam_header_t));

  int num_targets = genome->num_chroms;

  bam_header->n_targets = num_targets;
  bam_header->target_name = (char **) calloc(num_targets, sizeof(char *));
  bam_header->target_len = (uint32_t*) calloc(num_targets, sizeof(uint32_t));
  for (int i = 0; i < num_targets; i++) {
    bam_header->target_name[i] = strdup(genome->chrom_names[i]);
    bam_header->target_len[i] = genome->chrom_lengths[i] + 1;
  }
  bam_header->text = strdup("@PG\tID:HPG-Aligner\tVN:1.0\n");
  bam_header->l_text = strlen(bam_header->text);

  return bam_header;
}

//--------------------------------------------------------------------

int sa_bam_writer(void *data) {
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;
  
  sa_mapping_batch_t *mapping_batch = (sa_mapping_batch_t *) wf_batch->mapping_batch;
  if (mapping_batch == NULL) {
    printf("bam_writer1: error, NULL mapping batch\n");
    return;
  }

  //  for (int i = 0; i < NUM_COUNTERS; i++) {
  //    counters[i] += mapping_batch->counters[i];
  //  }

  #ifdef _TIMING
  for (int i = 0; i < NUM_TIMING; i++) {
    func_times[i] += mapping_batch->func_times[i];
  }
  #endif

  int flag, pnext = 0, tlen = 0;
  char *rnext = "*", *optional_flags = "NM:i:3";

  fastq_read_t *read;
  array_list_t *read_list = mapping_batch->fq_reads;

  bam1_t *bam1;
  alignment_t *alig;
  array_list_t *mapping_list;
  bam_file_t *out_file = wf_batch->writer_input->bam_file;

  sa_genome3_t *genome = wf_batch->sa_index->genome;

  size_t num_reads, num_mappings;
  num_reads = mapping_batch->num_reads;
  for (size_t i = 0; i < num_reads; i++) {
    read = (fastq_read_t *) array_list_get(i, read_list);
    mapping_list = mapping_batch->mapping_lists[i];
    num_mappings = array_list_size(mapping_list);
    #ifdef _VERBOSE
    if (num_mappings > 1) {
      num_dup_reads++;
      num_total_dup_reads += num_mappings;
    }
    #endif

    if (num_mappings > 0) {
      for (size_t j = 0; j < num_mappings; j++) {
	bam1 = (bam1_t *) array_list_get(j, mapping_list);
	bam_fwrite(bam1, out_file);
	bam_destroy1(bam1);
      }
    } else {
      alig = alignment_new();       
      alignment_init_single_end(strdup(read->id), read->sequence, read->quality,
				0, -1, -1, /*strdup(aux)*/"", 0, 0, 0, 0, 0, NULL, alig);
      
      bam1 = convert_to_bam(alig, 33);
      bam_fwrite(bam1, out_file);
        
      // free memory
      bam_destroy1(bam1);
      alig->sequence = NULL;
      alig->quality = NULL;
      alig->cigar = NULL;
      alignment_free(alig);
    }
    array_list_free(mapping_list, (void *) NULL);
  }

  // free memory
  sa_mapping_batch_free(mapping_batch);

  if (wf_batch) sa_wf_batch_free(wf_batch);

  return 0;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
