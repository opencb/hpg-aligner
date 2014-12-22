#ifndef FILTER_OPTIONS_H
#define FILTER_OPTIONS_H

/*
 * filter_options.h
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

#define NUM_FILTER_OPTIONS    13

//------------------------------------------------------------------------

typedef struct filter_options { 
  int log_level;
  int verbose;
  int help;
  int num_threads;
  int batch_size;
  
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

  char *read_length_range;
  char *read_quality_range;
  char *left_quality_range;
  char *right_quality_range;
  
  char* in_filename;
  char* out_dirname;

  char *exec_name;
  char *command_name;
} filter_options_t;

//------------------------------------------------------------------------

filter_options_t *filter_options_new(char *exec_name, char *command_nane);

filter_options_t *filter_options_parse(char *exec_name, char *command_nane,
			 int argc, char **argv);

void filter_options_free(filter_options_t *opts);

void filter_options_validate(filter_options_t *opts);

void filter_options_display(filter_options_t *opts);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif
