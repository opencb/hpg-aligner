/*
 * rna_server.c
 *
 *  Created on: Jan 20, 2012
 *  Last modified: Dec 20, 2013
 *  Author: Hector Martinez
 */

#ifndef RNA_ALIGNER_MPI_H
#define RNA_ALIGNER_MPI_H

#include <mpi.h>

#include "batch_writer.h"
#include "error.h"
#include "timing.h"
#include "buffers.h"
#include "bwt_server.h"
#include "region_seeker.h"
#include "cal_seeker.h"
#include "sw_server.h"
#include "pair_server.h"

#include "options.h"
#include "statistics.h"
#include "workflow_functions.h"
#include "rna/rna_server.h"

//For SA Mapping
#include "sa/sa_index3.h"
#include "sa/sa_search.h"
#include "rna/sa_rna_mapper.h"

#include "dna/sa_dna_commons.h"
#include "dna/doscadfun.h"
#include "dna/sa_io_stages.h"
#include "dna/sa_mapper_stage.h"

#include "containers/linked_list.h"
#include "containers/array_list.h"

#include "bioformats/fastq/fastq_file.h"
#include "bioformats/fastq/fastq_batch_reader.h"
#include "aligners/sw/smith_waterman.h"

void rna_aligner_mpi(options_t *options, int argc, char *argv[]);
//void rna_aligner_mpi_work_stealing(options_t *options, int argc, char *argv[]);


size_t MPI_fastq_fread_bytes_se(array_list_t *reads, size_t bytes, fastq_file_t *fq_file);
void* SA_MPI_producer(void *input);

void MPI_consumer(void *batch);
void MPI_output(void *data);


#endif
