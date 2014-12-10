#ifndef OPTIONS_H
#define OPTIONS_H

/*
 * options.h
 *
 *  Created on: May 22, 2013
 *      Author: jtarraga
 */

#include <stdlib.h>
#include <string.h>

#include "commons/argtable/argtable2.h"
#include "commons/config/libconfig.h"
#include "commons/log.h"
#include "commons/system_utils.h"
#include "commons/file_utils.h"

//------------------------------------------------------------------------

#define QUALITY_PHRED33_VALUE  33
#define QUALITY_PHRED33_NAME  "phred33"

#define QUALITY_PHRED64_VALUE  64
#define QUALITY_PHRED64_NAME  "phred64"


//------------------------------------------------------------------------

#define NO_VALUE       -1
#define NUM_OPTIONS    17

//------------------------------------------------------------------------

typedef struct options { 
  int kmers_on;
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
} options_t;

//------------------------------------------------------------------------

options_t *options_new(char *exec_name, char *command_nane);

options_t *options_parse(char *exec_name, char *command_nane,
			 int argc, char **argv);

void options_free(options_t *opts);

void options_validate(options_t *opts);

void options_display(options_t *opts);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif
