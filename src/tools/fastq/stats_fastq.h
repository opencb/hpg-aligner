#ifndef STATS_FASTQ_H
#define STATS_FASTQ_H

/*
 * stats_fastq.h
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

#include <stdlib.h>
#include <string.h>

//#include "argtable2.h"
//#include "libconfig.h"
//#include "commons/log.h"
//#include "commons/system_utils.h"
//#include "commons/file_utils.h"

#include "commons/workflow_scheduler.h"
#include "containers/array_list.h"
#include "containers/khash.h"

#include "bioformats/fastq/fastq_file.h"
#include "bioformats/fastq/fastq_stats.h"

#include "stats_options.h"

//--------------------------------------------------------------------
// stats counters
//--------------------------------------------------------------------

KHASH_MAP_INIT_INT(32, int)

typedef struct stats_counters {
  size_t num_reads;
  int phred;

  int filter_on;
  size_t num_passed;
  size_t num_failed;

  int min_length;
  int max_length;
  size_t acc_length;
  float mean_length;

  float acc_quality;
  float mean_quality;

  size_t num_As;
  size_t num_Cs;
  size_t num_Ts;
  size_t num_Gs;
  size_t num_Ns;
  
  int kmers_on;
  kmer_t kmers[NUM_KMERS];

  khash_t(32) *kh_length_histogram;
  khash_t(32) *kh_quality_histogram;
  khash_t(32) *kh_gc_histogram;

  khash_t(32) *kh_count_quality_per_nt;
  khash_t(32) *kh_acc_quality_per_nt;
  
  khash_t(32) *kh_num_As_per_nt;
  khash_t(32) *kh_num_Cs_per_nt;
  khash_t(32) *kh_num_Ts_per_nt;
  khash_t(32) *kh_num_Gs_per_nt;
  khash_t(32) *kh_num_Ns_per_nt;

} stats_counters_t;

//------------------------------------------------------------------------

void stats_fastq(stats_options_t *opts);

//------------------------------------------------------------------------

#endif // end of STATS_FASTQ_H

//------------------------------------------------------------------------
//------------------------------------------------------------------------
