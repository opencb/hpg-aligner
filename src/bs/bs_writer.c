#include "bs_writer.h"
#include "methylation.h"

int bs_writer(void *data) {

  struct timeval start, end;
  double time;
  
  if (time_on) { start_timer(start); }
  
  batch_t *batch = (batch_t *) data;
  fastq_read_t *fq_read;
  array_list_t *array_list;
  size_t num_items;
  
  mapping_batch_t *mapping_batch = (mapping_batch_t *) batch->mapping_batch;
  bs_context_t *bs_context = (mapping_batch->bs_context);
  /*
  {
    LOG_DEBUG("before reverting 2...\n");
    array_list_t **mapping_lists = mapping_batch->mapping_lists2;
    size_t num_reads = array_list_size(mapping_batch->fq_batch);
    for (size_t i = 0; i < num_reads; i++) {
      num_items = array_list_size(mapping_lists[i]);
      // mapped or not mapped ?	 
      if (num_items > 0) {
	for (size_t j = 0; j < num_items; j++) {
	  alignment_t *alig = (alignment_t *) array_list_get(j, mapping_lists[i]);
	  alignment_print(alig);
	}
      }
    }
  }
  */
  
  // set the sequences of the mapping to the original
  revert_mappings_seqs(mapping_batch->mapping_lists, mapping_batch->mapping_lists2, mapping_batch->fq_batch);
  
  batch_writer_input_t *writer_input = batch->writer_input;
  bam_file_t *bam_file = writer_input->bam_file;     
  linked_list_t *linked_list = writer_input->list_p;
  size_t num_reads = array_list_size(mapping_batch->fq_batch);
  size_t num_mapped_reads = 0;
  size_t total_mappings = 0;
  unsigned char found_p1 = 0;
  unsigned char found_p2 = 0;
  int i = 0;
  
  extern size_t bwt_correct;
  extern size_t bwt_error;
  extern pthread_mutex_t bwt_mutex;
  
  writer_input->total_batches++;
  
  array_list_t **mapping_lists;
  int *found = (int *) calloc(num_reads, sizeof(int));
  metil_file_t *metil_file = writer_input->metil_file;

  array_list_t *bs_stat = mapping_batch->bs_status;
  char *bs_seq;

  // process mapping_lists and mapping_lists2
  for (int k = 0; k < 2; k++) {
    //printf("\t-----> mapping_lists %i\n", k);
    
    mapping_lists = (k == 0) ? mapping_batch->mapping_lists : mapping_batch->mapping_lists2;
    
    for (size_t i = 0; i < num_reads; i++) {
      num_items = array_list_size(mapping_lists[i]);
      total_mappings += num_items;
      fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->fq_batch);

      // mapped or not mapped ?	 
      if (num_items == 0) {
	if (mapping_lists[i]) {
	  array_list_free(mapping_lists[i], NULL);
	}
      } else {
	found[i] = 1;
	LOG_DEBUG_F("printing alignment (bs_id = %i)...\n", k);
	write_mapped_read(mapping_lists[i], bam_file);
      }
    }
  }
  
  for (size_t i = 0; i < num_reads; i++) {
    if (found[i]) {
      num_mapped_reads++; //
    } else {
      total_mappings++;
      fq_read = (fastq_read_t *) array_list_get(i, mapping_batch->fq_batch);
      write_unmapped_read(fq_read, bam_file);
    }
  }

  // to use with the postprocess
  //write_bs_context(metil_file, bs_context);

  
  /*
  //write_context_bs(metil_file, bs_context);
  printf("4 FINAL       \t%3lu\tmapp               \t%3lu\tno map (discard) \t%3lu\n",
	 num_reads, num_mapped_reads, num_reads - num_mapped_reads);
  */
  
  if (basic_st->total_reads >= writer_input->limit_print) {
    LOG_DEBUG_F("TOTAL READS PROCESS: %lu\n", basic_st->total_reads);
    LOG_DEBUG_F("\tTotal Reads Mapped: %lu(%.2f%)\n", 
		basic_st->num_mapped_reads, 
		(float) (basic_st->num_mapped_reads*100)/(float)(basic_st->total_reads));
    writer_input->limit_print += 1000000;
  }
  
  //printf("Batch Write OK!\n");     
  
  if (mapping_batch) mapping_batch_free(mapping_batch);
  
  if (batch) batch_free(batch);
  
  if (found) free(found);

  basic_statistics_add(num_reads, num_mapped_reads, total_mappings, basic_st);
  
  if (time_on) { stop_timer(start, end, time); timing_add(time, BAM_WRITER, timing); }
}
