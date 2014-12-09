#ifndef SORT_OPTIONS_H
#define SORT_OPTIONS_H

/*
 * sort_options.h
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

#include <stdlib.h>
#include <string.h>

#include "argtable/argtable2.h"
#include "config/libconfig.h"
#include "commons/log.h"
#include "commons/system_utils.h"
#include "commons/file_utils.h"

//============================ DEFAULT VALUES ============================

//------------------------------------------------------------------------

#define NUM_SORT_OPTIONS	5

//------------------------------------------------------------------------

typedef struct sort_options { 
  int help;

  size_t max_memory;
  char *criteria;

  char* in_filename;
  char* out_dirname;

  char *exec_name;
  char *command_name;
} sort_options_t;

//------------------------------------------------------------------------

sort_options_t *sort_options_new(char *exec_name, char *command_nane);

sort_options_t *sort_options_parse(char *exec_name, char *command_nane,
				     int argc, char **argv);

void sort_options_free(sort_options_t *opts);

void sort_options_validate(sort_options_t *opts);

void sort_options_display(sort_options_t *opts);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif
