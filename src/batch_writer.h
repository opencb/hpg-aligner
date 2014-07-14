#ifndef BATCH_WRITER_H
#define BATCH_WRITER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "buffers.h"
#include "timing.h"

#include "commons/commons.h"
#include "commons/system_utils.h"

#include "containers/list.h"
#include "containers/cprops/hashtable.h"

#include "bioformats/fastq/fastq_file.h"
#include "bioformats/fastq/fastq_batch.h"
#include "bioformats/bam/bam_file.h"

#include "buffers.h"
#include "timing.h"

#include "bs/methylation.h"

//====================================================================================
//  structures and prototypes
//====================================================================================

struct batch_writer_input {
  char* match_filename;
  char* mismatch_filename;

  char* splice_exact_filename;
  char* splice_extend_filename;
  
  //  char* header_filename;
  genome_t* genome;

  linked_list_t* list_p;

  // internal
  bam_file_t *bam_file;
  size_t total_batches;
  size_t total_reads;
  size_t total_mappings;
  size_t num_mapped_reads;
  size_t limit_print;

  // for methylation only
  metil_file_t *metil_file;
};

//------------------------------------------------------------------------------------

void batch_writer_input_init(char* match_filename, char* splice_exact_filename, 
			     char* splice_extend_filename, linked_list_t* list_p, 
			     genome_t* genome, batch_writer_input_t* input);

//====================================================================================

bam_header_t *create_bam_header_by_genome(genome_t *genome);

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

#endif // BATCH_WRITER_H
