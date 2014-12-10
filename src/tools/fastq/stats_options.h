#ifndef STATS_OPTIONS_H
#define STATS_OPTIONS_H

/*
 * stats_options.h
 *
 *  Created on: May 22, 2013
 *      Author: jtarraga
 */

#include <stdlib.h>
#include <string.h>

#include "argtable/argtable2.h"
#include "config/libconfig.h"
#include "commons/log.h"
#include "commons/system_utils.h"
#include "commons/file_utils.h"

#include "bioformats/fastq/fastq_stats.h"

#include "commons_fastq.h"

//------------------------------------------------------------------------

#define NUM_STATS_OPTIONS    15

//------------------------------------------------------------------------

typedef struct stats_options { 
  int kmers_on;
  int filter_on;
  int log_level;
  int verbose;
  int help;
  int num_threads;
  int batch_size;
  int quality_encoding_value;
  
  int max_read_length;
  int min_read_length;
  int max_N;
  int max_read_quality;
  int min_read_quality;
  int left_length;
  int max_left_quality;
  int min_left_quality;
  int right_length;
  int max_right_quality;
  int min_right_quality;
  int max_out_of_quality;

  char *quality_encoding_name;
  char *read_length_range;
  char *read_quality_range;
  char *left_quality_range;
  char *right_quality_range;
  
  char* in_filename;
  char* out_dirname;

  char *exec_name;
  char *command_name;
} stats_options_t;

//------------------------------------------------------------------------

stats_options_t *stats_options_new(char *exec_name, char *command_nane);

stats_options_t *stats_options_parse(char *exec_name, char *command_nane,
			 int argc, char **argv);

void stats_options_free(stats_options_t *opts);

void stats_options_validate(stats_options_t *opts);

void stats_options_display(stats_options_t *opts);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif
