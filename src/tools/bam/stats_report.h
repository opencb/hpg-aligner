#ifndef STATS_REPORT_H
#define STATS_REPORT_H

/*
 * stats_report.h
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

//#include "commons/workflow_scheduler.h"
//#include "containers/array_list.h"
//#include "bioformats/bam-sam/bam_file.h"
//#include "bioformats/features/region/region_table.h"

#include "stats_bam.h"

//------------------------------------------------------------------------

#define MAX_SERIES_GRAPH	10
#define MAX_LINE_GRAPH  	1024

//------------------------------------------------------------------------

#define NUM_ERRORS_STATS 20
#define QUALITY_STATS   256

//------------------------------------------------------------------------

typedef struct stats_counters {

  // from filter
  size_t num_passed;
  size_t num_failed;

  // global statistics
  size_t ref_length;
  int num_sequences;
  int single_end;

  size_t num_reads;
  size_t num_unique_alignments;
  size_t num_mapped_reads;
  size_t num_unmapped_reads;
  size_t num_mapped_reads_1;
  size_t num_unmapped_reads_1;
  size_t num_mapped_reads_2;
  size_t num_unmapped_reads_2;

  size_t min_alignment_length;
  size_t max_alignment_length;

  // stats per strand
  size_t num_unique_alignments_strand[2];
  size_t num_mapped_reads_strand[2];

  // errors stats
  size_t num_indels;
  size_t indels_acc;
  size_t num_errors[NUM_ERRORS_STATS + 1];

  // nucleotide content
  size_t num_nucleotides;
  size_t num_As;
  size_t num_Cs;
  size_t num_Ts;
  size_t num_Gs;  
  size_t num_Ns;  
  size_t GC_content[100];

  // insert
  size_t min_insert_size;
  size_t max_insert_size;
  size_t insert_size_acc;

  // quality
  size_t min_quality;
  size_t max_quality;
  size_t quality_acc;
  size_t quality[QUALITY_STATS];

  // coverage (depth): global and per chromosome
  size_t unmapped_nts;
  double depth;
  char** sequence_labels;
  //int **sequence_depths_per_nt;
  uint16_t **sequence_depths_per_nt;
  size_t* sequence_lengths;
  double* depth_per_sequence;

  //  khash_t(32) *gc_hash;
} stats_counters_t;

//------------------------------------------------------------------------

stats_counters_t *stats_counters_new();
void stats_counters_free(stats_counters_t *p);

//------------------------------------------------------------------------

/**
* @brief Structure for report graphs parameters
* 
* Structure containing parameters for graphics representation of report graphs
*/
typedef struct report_graph {
    int x_autoscale;					/**< Autoscale flag for X axis. */
    int x_start;					/**< X axis start coordinate. */
    int x_end;						/**< X axis end coordinate. */
    int y_autoscale;					/**< Autoscale flag for Y axis. */
    int y_start;					/**< Y axis start coordinate. */
    int y_end;						/**< Y axis end coordinate. */
    int lmargin;					/**< Left margin. */
    int rmargin;					/**< Right margin. */	
    int tmargin;					/**< Top margin. */
    int bmargin;					/**< Bottom margin. */
    int num_y_columns;					/**< Number of columns in the Y axis. */
    int x_column;					/**< Number of column in data file with X axis data. */
    int y_columns[MAX_SERIES_GRAPH];			/**< Numbers of columns in data file that are represented in Y axis. */
    char y_titles[MAX_SERIES_GRAPH][MAX_LINE_GRAPH];			/**< Titles of series represented in Y axis. */
    char title[MAX_LINE_GRAPH];					/**< Title of the graphic. */
    char xlabel[MAX_LINE_GRAPH];					/**< Label of X axis. */	
    char ylabel[MAX_LINE_GRAPH];					/**< Label of Y axis. */
    char type[MAX_LINE_GRAPH]; 					/**< Type of graph. */
} report_graph_t;


//------------------------------------------------------------------------

void report_stats(char *bam_filename, char *outdir, 
		  void *db, stats_counters_t *stats);

//------------------------------------------------------------------------

#endif // STATS_REPORT_H

//------------------------------------------------------------------------
//------------------------------------------------------------------------
