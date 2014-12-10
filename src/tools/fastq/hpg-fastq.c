/*
 * hpg-fastq.c
 *
 *  Created on: Mar 8, 2013
 *      Author: jtarraga
 */

#include "stats_options.h"
#include "filter_options.h"
#include "edit_options.h"

#include "stats_fastq.h"
#include "filter_fastq.h"
#include "edit_fastq.h"

//------------------------------------------------------------------------

void usage(char *exec_name) {
    printf("Program: %s (High-performance tools for handling FastQ files)\n", exec_name);
    printf("\n");
    printf("Usage: %s <command> [options]\n", exec_name);
    printf("\n");
    printf("Command: stats\t\tstatistics summary\n");
    printf("         filter\t\tfilter a FastQ file by using advanced criteria\n");
    printf("         edit\t\tedit a FastQ file according the specified options\n");
    printf("\n");    
    printf("For more information about a certain command, type %s <command> --help\n", exec_name);
    exit(-1);
}

//------------------------------------------------------------------------
//                    M A I N     F U N C T I O N
//------------------------------------------------------------------------

int main (int argc, char *argv[]) {
  // init logs, then after parsing the command-line
  // logs will be re-set according to the command-line
  log_level = LOG_FATAL_LEVEL;
  log_verbose = 1;
  log_file = NULL;

  char *exec_name = argv[0];  
  if (argc == 1 || !strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) {
    usage(exec_name);
  }

  char *command_name = argv[1];  

  argc--;
  argv++;

  if (strcmp(command_name, "stats") == 0) {

    //--------------------------------------------------------------------
    //                  S T A T S     C O M M A N D
    //--------------------------------------------------------------------
    
    // parse, validate and display stats options
    stats_options_t *opts = stats_options_parse(exec_name, command_name, 
						 argc, argv);
    stats_options_validate(opts);
    stats_options_display(opts);
    
    // now, we can set logs according to the command-line
    init_log_custom(opts->log_level, 1, "hpg-fastq.log", "w");

    // run stats
    stats_fastq(opts);

    // free memory
    stats_options_free(opts);

  } else if (strcmp(command_name, "filter" ) == 0) {

    //--------------------------------------------------------------------
    //                  F I L T E R     C O M M A N D
    //--------------------------------------------------------------------

    // parse, validate and display filter options
    filter_options_t *opts = filter_options_parse(exec_name, command_name, 
						  argc, argv);
    filter_options_validate(opts);
    filter_options_display(opts);
 
    // now, we can set logs according to the command-line
    init_log_custom(opts->log_level, 1, "hpg-fastq.log", "w");

    // run filter
    filter_fastq(opts);

    // free memory
    filter_options_free(opts);

  } else if (strcmp(command_name, "edit" ) == 0) {

    //--------------------------------------------------------------------
    //                  E D I T     C O M M A N D
    //--------------------------------------------------------------------

    // parse, validate and display edit options
    edit_options_t *opts = edit_options_parse(exec_name, command_name, 
					      argc, argv);
    edit_options_validate(opts);
    edit_options_display(opts);
    
    // now, we can set logs according to the command-line
    init_log_custom(opts->log_level, 1, "hpg-fastq.log", "w");

    // run filter
    edit_fastq(opts);

    // free memory
    edit_options_free(opts);

  } else {

    //--------------------------------------------------------------------
    //               U N K N O W N     C O M M A N D
    //--------------------------------------------------------------------

    usage(exec_name);
  }


  LOG_INFO("Done !\n");
}

//------------------------------------------------------------------------
//------------------------------------------------------------------------
