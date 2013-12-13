/*
 * main.c
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

//#include "containers/khash.h"

//#include "bioformats/features/region/region_table.h"
//#include "bioformats/features/region/region_table_utils.h"
//#include "bioformats/bam/bam_db.h"
//#include "bioformats/bam/bam_stats.h"
//#include "bioformats/bam/bam_stats_report.h"
//#include "bioformats/bam/bam_filter.h"

//#include "bioformats/db/db_utils.h"

#include "stats_bam.h"
#include "filter_bam.h"
#include "index_options.h"
#include "sort_options.h"

//------------------------------------------------------------------------

int recalibrate_bam(int argc, char **argv);

//------------------------------------------------------------------------

void usage(char *exec_name) {
    printf("Program: %s (High-performance tools for handling BAM files)\n", exec_name);
    printf("Version: 1.0.0\n");
    printf("\n");
    printf("Usage: %s <command> [options]\n", exec_name);
    printf("\n");
    printf("Command: stats\t\tstatistics summary\n");
    printf("         filter\t\tfilter a BAM file by using advanced criteria\n");
    printf("         recalibrate\tbase quality recalibrate from a BAM file\n");
    printf("         realign\tlocal realign from BAM file\n");
    printf("         index\t\tindex a BAM file (using the samtools 0.1.18)\n");
    printf("         sort\t\tsort a BAM file (using the samtools 0.1.18)\n");

    //    printf("         compare\tcompare two BAM files\n");
    //    printf("         realignment\trealign locally a BAM file\n");
    //    printf("         recalibrate\trecalibrate a BAM file\n");
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
    init_log_custom(opts->log_level, 1, "hpg-bam.log", "w");

    // run stats
    stats_bam(opts);

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

    // run filter
    filter_bam(opts);

    // free memory
    filter_options_free(opts);

    /*    
    // create the region table
    region_table_t *region_table = build_region_table(opts->region_list,
						      opts->gff_region_filename);
    
    // set parameters
    int by_num_errors, min_num_errors, max_num_errors;
    if (opts->min_num_errors != -1 ||
	opts->max_num_errors != -1) {
      by_num_errors = 1;
      min_num_errors = opts->min_num_errors; 
      max_num_errors = opts->max_num_errors; 
    } else {
      by_num_errors = 0;
    }
    int by_quality, min_quality, max_quality;
    if (opts->min_quality != -1 ||
	opts->max_quality != -1) {
      by_quality = 1;
      min_quality = opts->min_quality; 
      max_quality = opts->max_quality; 
    } else {
      by_quality = 0;
    }
    int by_length, min_length, max_length;
    if (opts->min_length != -1 ||
	opts->max_length != -1) {
      by_length = 1;
      min_length = opts->min_length; 
      max_length = opts->max_length; 
    } else {
      by_length = 0;
    }

    bam_filter_input_t* input = new_bam_filter_input(opts->in_filename, opts->out_dirname,
						     opts->mapped, opts->unmapped, 
						     opts->proper_pairs, opts->unique, 
						     by_num_errors, min_num_errors, max_num_errors,
						     by_quality, min_quality, max_quality,
						     by_length, min_length, max_length,	
						     region_table,
						     opts->num_threads, opts->batch_size);
    
    // run filter
    filter_bam(input);

    // free memory
    free_bam_filter_input(input);
    if (region_table) free_table(region_table);
    free_filter_options(opts);
    */

  } else if (strcmp(command_name, "recalibrate" ) == 0) {

    //--------------------------------------------------------------------
    //           R E C A L I B R A T E     C O M M A N D
    //--------------------------------------------------------------------

    printf("***********************\n");
    // recalibrate
    recalibrate_bam(argc, argv);

  } else if (strcmp(command_name, "realign" ) == 0) {

    //--------------------------------------------------------------------
    //           R E A L I G N M E N T     C O M M A N D
    //--------------------------------------------------------------------

    printf("***********************\n");
    // local realignment
    alig_bam(argc, argv);

  } else if (strcmp(command_name, "index" ) == 0) {

    //--------------------------------------------------------------------
    //                  I N D E X     C O M M A N D
    //--------------------------------------------------------------------

    // parse, validate and display index options
    index_options_t *opts = index_options_parse(exec_name, command_name, 
					  argc, argv);
    index_options_validate(opts);
    index_options_display(opts);

    // set parameters
    char idx_filename[strlen(opts->out_dirname) + strlen(opts->in_filename) + 10];
    char *p = strrchr(opts->in_filename, '/');
    if (p) {
      sprintf(idx_filename, "%s/%s.bai", opts->out_dirname, p + 1);
    } else {
      sprintf(idx_filename, "%s/%s.bai", opts->out_dirname, opts->in_filename);
    }

    // run index
    bamFile bf = bam_open(opts->in_filename, "r");
    bam_index_t *idx = bam_index_core(bf);
    bam_close(bf);
    
    FILE *idxf = fopen(idx_filename, "wb");
    if (idxf == NULL) {
      LOG_FATAL_F("Could not open file %s", idx_filename);
    }
    bam_index_save(idx, idxf);
    bam_index_destroy(idx);
    fclose(idxf);

    printf("Index created in %s\n", idx_filename);

    // free memory
    index_options_free(opts);

  } else if (strcmp(command_name, "sort" ) == 0) {

    //--------------------------------------------------------------------
    //                  S O R T     C O M M A N D
    //--------------------------------------------------------------------

    // parse, validate and display sort options
    sort_options_t *opts = sort_options_parse(exec_name, command_name, 
					      argc, argv);
    sort_options_validate(opts);
    sort_options_display(opts);

    // set parameters
    int is_by_qname = 0, is_stdout = 0;
    char *sorted_filename[strlen(opts->in_filename)];
    char *path[strlen(opts->out_dirname) + strlen(opts->in_filename) + 100];

    strcpy(sorted_filename, opts->in_filename);
    char *ext = strstr(sorted_filename, ".bam");
    if (ext) {
      *ext = 0;
    }
    sprintf(path, "%s/%s.sorted", opts->out_dirname, sorted_filename);

    if (strcmp("name", opts->criteria) == 0) {
      is_by_qname = 1;
    }

    // run sort
    bam_sort_core_ext(is_by_qname, opts->in_filename, path, opts->max_memory, is_stdout);

    printf("Sorted BAM file in %s.bam\n", path);

    // free memory
    sort_options_free(opts);
  } else {

    //--------------------------------------------------------------------
    //                  U N K N O W N     C O M M A N D
    //--------------------------------------------------------------------

    usage(exec_name);
  }
  printf("Done !\n");
}

//------------------------------------------------------------------------

//region_table_t *build_region_table(char *bam_filename, 
//				   char *by_string, char *by_gff_file);

//extern int bam_index_build2(const char *fn, const char *_fnidx);

//------------------------------------------------------------------------


//------------------------------------------------------------------------
//------------------------------------------------------------------------
