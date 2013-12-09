#ifndef STATS_OPTIONS_H
#define STATS_OPTIONS_H

/*
 * stats_options.h
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

#include <stdlib.h>
#include <string.h>

#include "commons/argtable/argtable2.h"
#include "commons/config/libconfig.h"
#include "commons/log.h"
#include "commons/system_utils.h"
#include "commons/file_utils.h"

#include "bioformats/db/db_utils.h"

#include "commons_bam.h"

//============================ DEFAULT VALUES ============================

//------------------------------------------------------------------------

#define NUM_STATS_OPTIONS	10

//------------------------------------------------------------------------

typedef struct stats_options { 
  int log_level;
  int verbose;
  int help;
  int num_threads;
  int batch_size;
  int db_on;

  khash_t(stats_chunks) *hash;
  sqlite3 *db;
  region_table_t *region_table;

  char* in_filename;
  char* out_dirname;
  char* out_dbname;
  char* gff_region_filename;
  char* region_list;

  char *exec_name;
  char *command_name;
} stats_options_t;

//------------------------------------------------------------------------

stats_options_t *stats_options_new(char *exec_name, char *command_nane);
void stats_options_free(stats_options_t *opts);


stats_options_t *stats_options_parse(char *exec_name, char *command_nane,
				     int argc, char **argv);

void stats_options_validate(stats_options_t *opts);

void stats_options_display(stats_options_t *opts);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

#endif
