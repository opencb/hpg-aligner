#include "sa_io_stages.h"

//====================================================================
// PRODUCER
//====================================================================

//--------------------------------------------------------------------
// sa fq reader
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
				   sa_mapping_batch,
				   NULL);
  }

  return new_wf_batch;

}

//----------------------------------------------------------------------
// BAM reader for single-end
//----------------------------------------------------------------------

void *sa_bam_reader_single(void *input) {
  sa_wf_input_t *wf_input = (sa_wf_input_t *) input;
  
  sa_wf_batch_t *new_wf_batch = NULL;
  sa_wf_batch_t *curr_wf_batch = wf_input->wf_batch;
  
  fastq_batch_reader_input_t *fq_reader_input = wf_input->fq_reader_input;
  
  int batch_size = fq_reader_input->batch_size;
  array_list_t *reads = array_list_new(batch_size, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  
  bam_file_t *bam_file = (bam_file_t *) fq_reader_input->fq_file1;
  stats_t *stats = (stats_t *) wf_input->stats;

  bam1_t *bam1;
  
  int size = 0;
  int total_reads = 0;
  
  fastq_read_t *read;
  char *header, *sequence, *quality;

  bam1 = bam_init1();
  while ((bam_read1(bam_file->bam_fd, bam1) > 0) && (size < batch_size) ) {
    // convert bam1_t to fastq_read_t
    total_reads++;
    if (!(bam1->core.flag & BAM_FSECONDARY)) {
      header = strdup(bam1_qname(bam1));
      sequence = calloc(sizeof(char), (int32_t)bam1->core.l_qseq + 1);
      quality = calloc(sizeof(char), (int32_t)bam1->core.l_qseq + 1);
      
      bam1_get_sequence(bam1, sequence);
      if (bam1->core.flag & BAM_FREVERSE){
	revcomp_seq(sequence);
      }
      
      bam1_get_quality(bam1, quality);
      
      size += bam1->core.l_qname + 2 * bam1->core.l_qseq;
      read = fastq_read_new(header, sequence, quality);
      array_list_insert(read, reads);
      free(header);
      free(sequence);
      free(quality);
    } else {
      stats->secondary_reads++;
    }
  } // end of while
  bam_destroy1(bam1);
  
  size_t num_reads = array_list_size(reads);
  
  if (num_reads == 0) {
    array_list_free(reads, (void *)fastq_read_free);
  } else {
    sa_mapping_batch_t *sa_mapping_batch = sa_mapping_batch_new(reads);
    sa_mapping_batch->bam_format = wf_input->bam_format;
    
    new_wf_batch = sa_wf_batch_new(curr_wf_batch->options,
				   curr_wf_batch->sa_index,
				   curr_wf_batch->writer_input, 
				   sa_mapping_batch,
				   NULL);
  }
  
  stats->total_reads+=total_reads;
  
  return new_wf_batch;
}

//----------------------------------------------------------------------
// BAM reader for paired-end)
//----------------------------------------------------------------------

void *sa_bam_reader_pairend(void *input) { 
  sa_wf_input_t *wf_input = (sa_wf_input_t *) input;
  
  sa_wf_batch_t *new_wf_batch = NULL;
  sa_wf_batch_t *curr_wf_batch = wf_input->wf_batch;
  
  fastq_batch_reader_input_t *fq_reader_input = wf_input->fq_reader_input;

  int batch_size = fq_reader_input->batch_size;
  array_list_t *reads = array_list_new(batch_size, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  bam_file_t *bam_file = (bam_file_t *) fq_reader_input->fq_file1;
  bamFile bam_file_aux = bam_open(fq_reader_input->filename1, "r");
  bam_file_t *fnomapped = (bam_file_t *)wf_input->data;
  stats_t *stats = (stats_t *)wf_input->stats;

  bam1_t *bam1; // mate 1
  bam1 = bam_init1();
  
  bam1_t *bam2; // mate 2
  bam2 = bam_init1();

  bam_iter_t iter;
  
  int ret, found, size = 0, total_reads = 0;
  char *header, *sequence, *quality;
  fastq_read_t *read;
    
  while ((bam_read1(bam_file->bam_fd, bam1) > 0) && (size < batch_size)) {

    if ((bam1->core.flag & BAM_FREAD1) 
	&& !(bam1->core.flag & BAM_FMUNMAP) 
	&& !(bam1->core.flag & BAM_FSECONDARY)) {

      // bam_fetch modified function
      iter = bam_iter_query(wf_input->idx, (int) bam1->core.mtid,
			    (int) bam1->core.mpos, (int) bam1->core.mpos + 1);
      if (iter == NULL){
	printf("Wrong fail in the index or in bam file, pair not found");	
      }

      found = 0;
      bam_seek(bam_file_aux, 0, SEEK_SET);
      while ((ret = bam_iter_read(bam_file_aux, iter, bam2)) >= 0) {
	if ((!strcmp(bam1_qname(bam1), bam1_qname(bam2))) && (bam1->core.mtid == bam2->core.tid)) {
	  found = 1;
	  break;
	}
      }
      bam_iter_destroy(iter);
      
      if (!found) {
	printf("\n Corrupt Bam, missing a mate, there is a single read \n");
	printf("\n The single read is:  %s \n", bam1_qname(bam1));
	printf("\n ***Check your Bam***\n\n");
	stats->total_reads++;
	stats->alone_reads++;
      } else {
	// convert bam1_t to fastq_read_t
	header = strdup(bam1_qname(bam1));
	sequence = calloc(sizeof(char), (int32_t) bam1->core.l_qseq + 1);
	quality = calloc(sizeof(char), (int32_t) bam1->core.l_qseq + 1);

	bam1_get_sequence(bam1, sequence);

	if (bam1->core.flag & BAM_FREVERSE){
	  revcomp_seq(sequence);
	}

	bam1_get_quality(bam1, quality);

	size += bam1->core.l_qname + 2 * bam1->core.l_qseq;
	read = fastq_read_new(header, sequence, quality);
	array_list_insert(read, reads);

	free(header);
	free(sequence);
	free(quality);

	header = strdup(bam1_qname(bam2));
	sequence = calloc(sizeof(char), (int32_t) bam2->core.l_qseq + 1);
	quality = calloc(sizeof(char), (int32_t) bam2->core.l_qseq + 1);

	bam1_get_sequence(bam2, sequence);
	if (bam2->core.flag  & BAM_FREVERSE) {
	  revcomp_seq(sequence);
	}

	bam1_get_quality(bam2, quality);

	size += bam2->core.l_qname + 2 * bam2->core.l_qseq;
	read = fastq_read_new(header, sequence, quality);
	array_list_insert(read, reads);

	free(header);
	free(sequence);
	free(quality);

	total_reads += 2;
      }
    } else if (((bam1->core.flag & BAM_FUNMAP) || 
		(bam1->core.flag & BAM_FMUNMAP)) && 
	       !(bam1->core.flag & BAM_FSECONDARY)) { 
      // Unmapped or mate, save it into a tmp file
      total_reads++;
      if (bam_fwrite(bam1, fnomapped) < 0) {
	bam_destroy1(bam1);
	bam_destroy1(bam2);
	printf("Fail in file fnomap");
	exit(-1);
      }
    } else if (bam1->core.flag & BAM_FSECONDARY) {
      stats->secondary_reads++;
      total_reads++;
    }
  } // end of while
  bam_destroy1(bam1);
  bam_destroy1(bam2);
  
  bam_close(bam_file_aux);
  
  size_t num_reads = array_list_size(reads);
  
  if (num_reads == 0) {
    array_list_free(reads, (void *) fastq_read_free);
  } else {
    sa_mapping_batch_t *sa_mapping_batch = sa_mapping_batch_new(reads);
    sa_mapping_batch->bam_format = wf_input->bam_format;
    
    new_wf_batch = sa_wf_batch_new(curr_wf_batch->options,
				   curr_wf_batch->sa_index,
				   curr_wf_batch->writer_input, 
				   sa_mapping_batch,
				   NULL);
  }
  
  stats->total_reads+=total_reads;
  
  return new_wf_batch;
}

//----------------------------------------------------------------------
// BAM reader for tmp unmapped reads
//----------------------------------------------------------------------

void *sa_bam_reader_unmapped(void *input) {

  sa_wf_input_t *wf_input = (sa_wf_input_t *) input;

  sa_wf_batch_t *new_wf_batch = NULL;
  sa_wf_batch_t *curr_wf_batch = wf_input->wf_batch;

  fastq_batch_reader_input_t *fq_reader_input = wf_input->fq_reader_input;

  int batch_size = fq_reader_input->batch_size;
  array_list_t *reads = array_list_new(batch_size, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  bam_file_t *fnomap = (bam_file_t *)wf_input->data;
  stats_t *stats = (stats_t *)wf_input->stats;

  bam1_t *bam1; // mate 1
  bam1 = bam_init1();

  bam1_t *bam2; // mate 2
  bam2 = bam_init1();

  int size = 0, found = 0;
  char *header, *sequence, *quality;
  fastq_read_t *read;

  char *id2 = "Hola:)";

  while ((bam_read1(fnomap->bam_fd, bam1) > 0) && (size < batch_size)) {

    char *id1 = bam1_qname(bam1);
    if (!strcmp(id1, id2)) {
      // if is same were sending
    } else if (bam1->core.flag & BAM_FUNMAP) {
      found = 0;
      while (bam_read1(fnomap->bam_fd, bam2) > 0) {
	id2 = bam1_qname(bam2);
	if (!strcmp(id1, id2)) { 
	  // need are same
	  found = 1;
	  break;
	} else {
	  printf("\n Corrupt Bam, missing a mate, there is a single read \n");
	  printf("\n The single read is:  %s \n", bam1_qname(bam1));
	  printf("\n ***Check your Bam***\n\n");
	  bam1 = bam2;
	  stats->alone_reads++;
	}
      }
    } else {
      found = 0;
      while (bam_read1(fnomap->bam_fd, bam2) > 0) {
	id2 = bam1_qname(bam2);
	if (!strcmp(id1, id2) && (bam2->core.flag & BAM_FUNMAP)) {
	  found = 1;
	  break;
	}
      }
    }
    if (!found) {
      printf("\n Corrupt Bam, missing a mate, there is a single read \n");
      printf("\n The single read is:  %s \n", bam1_qname(bam1));
      printf("\n ***Check your Bam***\n\n");
      stats->alone_reads++;
    } else {
      // ckeck the mates
      if ((bam1->core.flag & BAM_FREAD2) || (bam2->core.flag & BAM_FREAD1)) {

	bam1_t *bamaux; // mate 1
	bamaux = bam1;
	bam1 = bam2;
	bam2 = bamaux;
      }

      // convert bam1_t to fastq_read_t
      header = strdup(bam1_qname(bam1));
      sequence = calloc(sizeof(char), (int32_t) bam1->core.l_qseq + 1);
      quality = calloc(sizeof(char), (int32_t) bam1->core.l_qseq + 1);

      bam1_get_sequence(bam1, sequence);

      //check if read is reverse
      if (bam1->core.flag & BAM_FREVERSE){
	revcomp_seq(sequence);
      }

      bam1_get_quality(bam1, quality);

      size += bam1->core.l_qname + 2 * bam1->core.l_qseq;
      read = fastq_read_new(header, sequence, quality);
      array_list_insert(read, reads);

      free(header);
      free(sequence);
      free(quality);

      header = strdup(bam1_qname(bam2));
      sequence = calloc(sizeof(char), (int32_t) bam2->core.l_qseq + 1);
      quality = calloc(sizeof(char), (int32_t) bam2->core.l_qseq + 1);

      bam1_get_sequence(bam2, sequence);

      //check if read is reverse
      if (bam2->core.flag & BAM_FREVERSE){
	revcomp_seq(sequence);
      }

      bam1_get_quality(bam2, quality);

      size += bam2->core.l_qname + 2 * bam2->core.l_qseq;
      read = fastq_read_new(header, sequence, quality);
      array_list_insert(read, reads);

      free(header);
      free(sequence);
      free(quality);
    }

  } // end of while
  bam_destroy1(bam1);
  bam_destroy1(bam2);

  size_t num_reads = array_list_size(reads);

  if (num_reads == 0) {
    array_list_free(reads, (void *)fastq_read_free);
  } else {
    sa_mapping_batch_t *sa_mapping_batch = sa_mapping_batch_new(reads);
    sa_mapping_batch->bam_format = wf_input->bam_format;
    
    new_wf_batch = sa_wf_batch_new(curr_wf_batch->options,
				   curr_wf_batch->sa_index,
				   curr_wf_batch->writer_input,
				   sa_mapping_batch,
				   NULL);
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

size_t num_mapped_reads = 0;
size_t num_unmapped_reads = 0;
size_t num_total_mappings = 0;
size_t num_multihit_reads = 0;
size_t num_unmapped_reads_by_invalid_cal = 0;
size_t num_unmapped_reads_by_cigar_length = 0;

//--------------------------------------------------------------------
// SAM writer
//--------------------------------------------------------------------

void write_sam_header(options_t *options, sa_genome3_t *genome, FILE *f) {
  fprintf(f, "@HD\tVN:1.4\tSO:unsorted\n");
  fprintf(f, "@PG\tID:HPG-Aligner\tVN:%s\tCL:%s\n", HPG_ALIGNER_VERSION, options->cmdline);
  for (unsigned short int i = 0; i < genome->num_chroms; i++) {
    fprintf(f, "@SQ\tSN:%s\tLN:%lu\n", genome->chrom_names[i], genome->chrom_lengths[i]);
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

  #ifdef _TIMING
  for (int i = 0; i < NUM_TIMING; i++) {
    func_times[i] += mapping_batch->func_times[i];
  }
  #endif

  int num_mismatches, num_cigar_ops;
  size_t flag, pnext = 0, tlen = 0;
  char *cigar_string, *cigar_M_string, *rnext = "*";

  fastq_read_t *read;
  array_list_t *read_list = mapping_batch->fq_reads;

  array_list_t *mapping_list;
  FILE *out_file = (FILE *) wf_batch->writer_input->bam_file;

  sa_genome3_t *genome = wf_batch->sa_index->genome;

  size_t num_reads, num_mappings;
  num_reads = mapping_batch->num_reads;

  if (mapping_batch->options->pair_mode != SINGLE_END_MODE) {
    // PAIR MODE
    int len;
    char *sequence, *quality;

    char *opt_fields;
    alignment_t *alig;
  

    for (size_t i = 0; i < num_reads; i++) {
      read = (fastq_read_t *) array_list_get(i, read_list);

      mapping_list = mapping_batch->mapping_lists[i];
      num_mappings = array_list_size(mapping_list);
      num_total_mappings += num_mappings;

      #ifdef _VERBOSE
      if (num_mappings > 1) {
	num_dup_reads++;
	num_total_dup_reads += num_mappings;
      }
      #endif
      
      if (num_mappings > 0) {
	num_mapped_reads++;
	if (num_mappings > 1) {
	  num_multihit_reads++;
	}
	for (size_t j = 0; j < num_mappings; j++) {
	  alig = (alignment_t *) array_list_get(j, mapping_list);

	  // decoy management
	  if (genome->chrom_flags[alig->chromosome] == DECOY_FLAG) {
	    if (num_mappings == 1) {
	      fprintf(out_file, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tXD:Z:%s\n", 
		      read->id, alig->sequence, alig->quality, genome->chrom_names[alig->chromosome]);
	    }
	    // free alignment and continue
	    alignment_free(alig);	 
	    continue;
	  }

	  if (alig->optional_fields) {
	    opt_fields = (char *)alig->optional_fields;
	  } else {
	    opt_fields = NULL;
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
	  if (alig->seq_strand)                                 flag += BAM_FREVERSE;

	  fprintf(out_file, "%s\t%lu\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%s\t%s\t%s\n", 
		  read->id,
		  flag,
		  genome->chrom_names[alig->chromosome],
		  alig->position + 1,
		  (num_mappings > 1 ? 0 : alig->mapq),
		  alig->cigar,
		  (alig->chromosome == alig->mate_chromosome ? "=" : genome->chrom_names[alig->mate_chromosome]),
		  alig->mate_position + 1,
		  alig->template_length,
		  alig->sequence,
		  alig->quality,
		  (opt_fields == NULL ? "" : opt_fields)
		  );

	  // free memory
	  alignment_free(alig);	 
	} // end for num_mappings
      } else {
	num_unmapped_reads++;

	if (read->adapter) {
	  len = read->length + abs(read->adapter_length);
	  sequence = (char *) malloc(len + 1);
	  quality = (char *) malloc(len + 1);

	  if (read->adapter_length < 0) {
	    strcpy(quality, read->adapter_quality);
	    strcat(quality, read->quality);
	  } else {
	    strcpy(quality, read->quality);
	    strcat(quality, read->adapter_quality);
	  }
	  
	  if ((read->adapter_strand == 0 && read->adapter_length < 0) || 
	      (read->adapter_strand == 1 && read->adapter_length > 0)) {
	    strcpy(sequence, read->adapter);
	    strcat(sequence, read->sequence);
	  } else {
	    strcpy(sequence, read->sequence);
	    strcat(sequence, read->adapter);
	  }

	  sequence[len] = 0; 
	  quality[len] = 0; 
	} else {
	  sequence = read->sequence;
	  quality = read->quality;
	}

	fprintf(out_file, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", 
		read->id,
		sequence,
		quality
		);

	if (read->adapter) {
	  free(sequence);
	  free(quality);
	}
      }
      array_list_free(mapping_list, (void *) NULL);
    }
  } else {
    // SINGLE MODE
    int len;
    char *seq;
    seed_cal_t *cal;

    cigar_t *cigar;
    char *sequence, *revcomp, *quality;

    for (size_t i = 0; i < num_reads; i++) {
      read = (fastq_read_t *) array_list_get(i, read_list);
      mapping_list = mapping_batch->mapping_lists[i];
      num_mappings = array_list_size(mapping_list);
      num_total_mappings += num_mappings;

      #ifdef _VERBOSE
      if (num_mappings > 1) {
	num_dup_reads++;
	num_total_dup_reads += num_mappings;
      }
      #endif
      
      if (num_mappings > 0) {
	num_mapped_reads++;
	if (num_mappings > 1) {
	  num_multihit_reads++;
	}

	for (size_t j = 0; j < num_mappings; j++) {
	  cal = (seed_cal_t *) array_list_get(j, mapping_list);

	  if (read->adapter) {
	    // sequences and cigar
	    len = read->length + abs(read->adapter_length);
	    sequence = (char *) malloc(len + 1);
	    revcomp = (char *) malloc(len + 1);
	    quality = (char *) malloc(len + 1);
	    cigar = cigar_new_empty();

	    if (read->adapter_length < 0) {
	      strcpy(quality, read->adapter_quality);
	      strcat(quality, read->quality);
	    } else {
	      strcpy(quality, read->quality);
	      strcat(quality, read->adapter_quality);
	    }
	    
	    if ( (cal->strand == 1 && 
		  ((read->adapter_strand == 0 && read->adapter_length > 0) || 
		   (read->adapter_strand == 1 && read->adapter_length < 0)))
		 ||
		 (cal->strand == 0 && 
		  ((read->adapter_strand == 0 && read->adapter_length < 0) ||
		   (read->adapter_strand == 1 && read->adapter_length > 0))) ) {
	      strcpy(sequence, read->adapter);
	      strcat(sequence, read->sequence);
	      strcpy(revcomp, read->adapter_revcomp);
	      strcat(revcomp, read->revcomp);
	      
	      cigar_append_op(abs(read->adapter_length), 'S', cigar);
	      cigar_concat(&cal->cigar, cigar);
	    } else {
	      strcpy(sequence, read->sequence);
	      strcat(sequence, read->adapter);
	      strcpy(revcomp, read->revcomp);
	      strcat(revcomp, read->adapter_revcomp);
	      
	      cigar_concat(&cal->cigar, cigar);
	      cigar_append_op(read->adapter_length, 'S', cigar);
	    }
	    sequence[len] = 0; 
	    revcomp[len] = 0; 
	    quality[len] = 0; 
	  } else {
	    // sequences and cigar
	    sequence = read->sequence;
	    revcomp = read->revcomp;
	    quality = read->quality;
	    cigar = &cal->cigar;
	  }

	  // decoy management
	  if (genome->chrom_flags[cal->chromosome_id] == DECOY_FLAG) {
	    if (num_mappings == 1) {
	      fprintf(out_file, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\tXD:Z:%s\n", 
		      read->id, sequence, quality, genome->chrom_names[cal->chromosome_id]);
	    }
	    // go to free memory
	    goto free_memory1;
	  }

	  if (cal->strand) {
	    flag = 16;
	    seq = revcomp;
	  } else {
	    flag = 0;
	    seq = sequence;
	  }

	  cigar_string = cigar_to_string(cigar);
	  cigar_M_string = cigar_to_M_string(&num_mismatches, &num_cigar_ops, cigar);
	  if (num_mappings > 1) {
	    cal->mapq = 0;
	  }
	  fprintf(out_file, "%s\t%lu\t%s\t%lu\t%i\t%s\t%s\t%lu\t%lu\t%s\t%s\tAS:i:%i\tNM:i:%i\n", 
		  read->id,
		  flag,
		  genome->chrom_names[cal->chromosome_id],
		  cal->start + 1,
		  (num_mappings == 1 ? cal->mapq : 0),
		  cigar_M_string,
		  rnext,
		  pnext,
		  tlen,
		  seq,
		  quality,
		  (int) cal->score,
		  num_mismatches
		  );

	  // free memory
	  free(cigar_M_string);
	  free(cigar_string);
	free_memory1:
	  seed_cal_free(cal);	 
	  if (read->adapter) {
	    free(sequence);
	    free(revcomp);
	    free(quality);
	    cigar_free(cigar);
	  }
	}
      } else {
	num_unmapped_reads++;

	if (read->adapter) {
	  // sequences and cigar
	  len = read->length + abs(read->adapter_length);
	  sequence = (char *) malloc(len + 1);
	  quality = (char *) malloc(len + 1);

	  if (read->adapter_length < 0) {
	    strcpy(quality, read->adapter_quality);
	    strcat(quality, read->quality);
	  } else {
	    strcpy(quality, read->quality);
	    strcat(quality, read->adapter_quality);
	  }
	  
	  if ((read->adapter_strand == 0 && read->adapter_length < 0) || 
	      (read->adapter_strand == 1 && read->adapter_length > 0)) {
	    strcpy(sequence, read->adapter);
	    strcat(sequence, read->sequence);
	  } else {
	    strcpy(sequence, read->sequence);
	    strcat(sequence, read->adapter);
	  }

	  sequence[len] = 0; 
	  quality[len] = 0; 
	} else {
	  // sequences
	  sequence = read->sequence;
	  quality = read->quality;
	}
	
	fprintf(out_file, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", 
		read->id,
		sequence,
		quality
		);

	if (read->adapter) {
	  free(sequence);
	  free(quality);
	}
      }
      
      array_list_free(mapping_list, (void *) NULL);
    } // end for num_reads
  }

  // free memory
  sa_mapping_batch_free(mapping_batch);

  if (wf_batch) sa_wf_batch_free(wf_batch);

  return 0;
}


//--------------------------------------------------------------------
// BAM writer
//--------------------------------------------------------------------

bam_header_t *create_bam_header(options_t *options, sa_genome3_t *genome) {

  bam_header_t *bam_header = (bam_header_t *) calloc(1, sizeof(bam_header_t));

  int num_targets = genome->num_chroms;

  bam_header->n_targets = num_targets;
  bam_header->target_name = (char **) calloc(num_targets, sizeof(char *));
  bam_header->target_len = (uint32_t*) calloc(num_targets, sizeof(uint32_t));
  for (int i = 0; i < num_targets; i++) {
    bam_header->target_name[i] = strdup(genome->chrom_names[i]);
    bam_header->target_len[i] = genome->chrom_lengths[i];
  }

  char pg[1024];
  sprintf(pg, "@HD\tVN:1.4\tSO:unsorted\n");
  sprintf(pg, "@PG\tID:HPG-Aligner\tVN:%s\tCL:%s\n", HPG_ALIGNER_VERSION, options->cmdline);
  bam_header->text = strdup(pg);
  bam_header->l_text = strlen(bam_header->text);

  return bam_header;
}

//--------------------------------------------------------------------

int sa_bam_writer(void *data) {
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;
  
  sa_mapping_batch_t *mapping_batch = (sa_mapping_batch_t *) wf_batch->mapping_batch;
  if (mapping_batch == NULL) {
    printf("bam_writer1: error, NULL mapping batch\n");
    return 0;
  }

  #ifdef _TIMING
  for (int i = 0; i < NUM_TIMING; i++) {
    func_times[i] += mapping_batch->func_times[i];
  }
  #endif

  int len;
  char *sequence, *quality;

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
    num_total_mappings += num_mappings;

    #ifdef _VERBOSE
    if (num_mappings > 1) {
      num_dup_reads++;
      num_total_dup_reads += num_mappings;
    }
    #endif

    if (num_mappings > 0) {
      num_mapped_reads++;
      if (num_mappings > 1) {
	num_multihit_reads++;
      }
      for (size_t j = 0; j < num_mappings; j++) {
	alig = (alignment_t *) array_list_get(j, mapping_list);
	
	// decoy management
	if (genome->chrom_flags[alig->chromosome] == DECOY_FLAG) {
	  if (num_mappings == 1) {
	    alignment_t *aux_alig = alignment_new();       
	    alignment_init_single_end(strdup(read->id), alig->sequence, alig->quality,
				      0, -1, -1, "", 0, 0, 0, 0, 0, NULL, aux_alig);
	    bam1_t *aux_bam1 = convert_to_bam(aux_alig, 33);
	    bam_fwrite(aux_bam1, out_file);
	    // free memory
	    bam_destroy1(aux_bam1);
	    alignment_free(aux_alig);
	    alig->sequence = NULL;
	    alig->quality = NULL;
	  }
	  // free alignment and continue
	  alignment_free(alig);	 
	  continue;
	}

	// update alignment
	if (num_mappings > 1) {
	  alig->map_quality = 0;
	} else {
	  alig->map_quality = alig->mapq;
	}

	bam1 = convert_to_bam(alig, 33);
	bam_fwrite(bam1, out_file);
	bam_destroy1(bam1);
	alignment_free(alig);
      }
    } else {
      num_unmapped_reads++;

      if (read->adapter) {
	// sequences and cigar
	len = read->length + abs(read->adapter_length);
	sequence = (char *) malloc(len + 1);
	quality = (char *) malloc(len + 1);

	if (read->adapter_length < 0) {
	  strcpy(quality, read->adapter_quality);
	  strcat(quality, read->quality);
	} else {
	  strcpy(quality, read->quality);
	  strcat(quality, read->adapter_quality);
	}
	
	if ((read->adapter_strand == 0 && read->adapter_length < 0) || 
	    (read->adapter_strand == 1 && read->adapter_length > 0)) {
	  strcpy(sequence, read->adapter);
	  strcat(sequence, read->sequence);
	} else {
	  strcpy(sequence, read->sequence);
	  strcat(sequence, read->adapter);
	}
	sequence[len] = 0; 
	quality[len] = 0; 
      } else {
	// sequences
	sequence = read->sequence;
	quality = read->quality;
      }
      
      alig = alignment_new();       
      alignment_init_single_end(strdup(read->id), sequence, quality,
				0, -1, -1, "", 0, 0, 0, 0, 0, NULL, alig);
      
      bam1 = convert_to_bam(alig, 33);
      bam_fwrite(bam1, out_file);
        
      // free memory
      bam_destroy1(bam1);
      alig->sequence = NULL;
      alig->quality = NULL;
      alig->cigar = NULL;
      alignment_free(alig);
      if (read->adapter) {
	free(sequence);
	free(quality);
      }
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
