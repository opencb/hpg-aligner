#ifndef EDIT_OPTIONS_H
#define EDIT_OPTIONS_H

/*
 * edit_options.h
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

#define NUM_EDIT_OPTIONS    12

//------------------------------------------------------------------------

typedef struct edit_options { 
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
} edit_options_t;

//------------------------------------------------------------------------

edit_options_t *edit_options_new(char *exec_name, char *command_nane);

edit_options_t *edit_options_parse(char *exec_name, char *command_nane,
			 int argc, char **argv);

void edit_options_free(edit_options_t *opts);

void edit_options_validate(edit_options_t *opts);

void edit_options_display(edit_options_t *opts);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif
