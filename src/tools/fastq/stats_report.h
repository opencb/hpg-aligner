#ifndef STATS_REPORT_H
#define STATS_REPORT_H

/*
 * stats_report.h
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

#include "stats_fastq.h"
#include "stats_options.h"

//--------------------------------------------------------------------
// stats report
//--------------------------------------------------------------------

void stats_report(stats_counters_t *counters, stats_options_t *opts);
//void stats_report(char *in_filename, stats_counters_t *counters, char *out_dir);

//--------------------------------------------------------------------
// report graph
//--------------------------------------------------------------------

#define MAX_SERIES_GRAPH	10
#define MAX_LINE_GRAPH  	1024

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
    char title[MAX_LINE_GRAPH];					/**< Title of the graphic. */
    char xlabel[MAX_LINE_GRAPH];					/**< Label of X axis. */	
    char ylabel[MAX_LINE_GRAPH];					/**< Label of Y axis. */
    char type[MAX_LINE_GRAPH]; 					/**< Type of graph. */
    int y_columns[MAX_SERIES_GRAPH];			/**< Numbers of columns in data file that are represented in Y axis. */
    char y_titles[MAX_SERIES_GRAPH][MAX_LINE_GRAPH];			/**< Titles of series represented in Y axis. */
} report_graph_t;

//------------------------------------------------------------------------

#endif // STATS_REPORT_H

//------------------------------------------------------------------------
//------------------------------------------------------------------------
