/*
 * stats_report.c
 *
 *  Created on: Feb 19, 2013
 *      Author: jtarraga
 */

#include "bioformats/db/db_utils.h"

#include "stats_report.h"

//--------------------------------------------------------------------
// stats counters
//--------------------------------------------------------------------

stats_counters_t *stats_counters_new() {
  stats_counters_t *p = (stats_counters_t *) calloc(1, sizeof(stats_counters_t));

  // filter
  p->num_passed = 0;
  p->num_failed = 0;

  // global statistics
  p->ref_length = 0;
  p->num_sequences = 0;
  
  p->single_end = 1;
  
  p->num_reads = 0;
  p->num_unique_alignments = 0;
  p->num_mapped_reads = 0;
  p->num_unmapped_reads = 0;

  p->min_alignment_length = 999999999;
  p->max_alignment_length = 0;

  // stats per strand
  p->num_unique_alignments_strand[0] = 0;
  p->num_unique_alignments_strand[1] = 0;
  p->num_mapped_reads_strand[0] = 0;
  p->num_mapped_reads_strand[1] = 0;

  // errors stats
  p->num_indels = 0;
  p->indels_acc = 0;
  for (int i = 0; i <= NUM_ERRORS_STATS ; i++) {
    p->num_errors[i] = 0;
  }

  // nucleotide content
  p->num_As = 0;
  p->num_Cs = 0;
  p->num_Gs = 0;
  p->num_Ts = 0;
  p->num_Ns = 0;
  for (int i = 0; i < 100 ; i++) {
    p->GC_content[i] = 0;
  }

  // insert
  p->min_insert_size = 99999999;
  p->max_insert_size = 0;
  p->insert_size_acc = 0;

  // quality
  p->min_quality = 999999;
  p->max_quality = 0;
  p->quality_acc = 0;
  for (int i = 0; i < QUALITY_STATS ; i++) {
    p->quality[i] = 0;
  }

  // coverage (depth): global and per chromosome
  p->unmapped_nts = 0;
  p->sequence_labels = NULL;
  p->sequence_depths_per_nt = NULL;
  p->sequence_lengths = NULL;
  p->depth_per_sequence = NULL;
  p->depth = 0.0f;

  return p;
}

void stats_counters_free(stats_counters_t *p) {
  if (p) {
    
    if (p->sequence_labels) {
      for(int i = 0; i < p->num_sequences; i++) {
	if (p->sequence_labels[i]) free(p->sequence_labels[i]);
      }
      free(p->sequence_labels);
    } 
    
    if (p->sequence_depths_per_nt) {
      for(int i = 0; i < p->num_sequences; i++) {
	if (p->sequence_depths_per_nt[i]) free(p->sequence_depths_per_nt[i]);
      }
      free(p->sequence_depths_per_nt);
    }
    
    if (p->sequence_lengths) {
      free(p->sequence_lengths);
    }

    if (p->depth_per_sequence) {
      free(p->depth_per_sequence);
    }

    free(p);
  }
}

//--------------------------------------------------------------------

void init_report_graph(report_graph_t* graph) {
    graph->x_autoscale = 1;
    graph->x_start = 1;
    graph->x_end = 100;
    graph->y_autoscale = 1;
    graph->y_start = 0;
    graph->y_end = 100;
    graph->lmargin = 10;
    graph->rmargin = 4;
    graph->tmargin = 3;
    graph->bmargin = 4;
    *graph->title = 0;
    *graph->xlabel = 0;
    *graph->ylabel = 0;
    *graph->type = 0;
    graph->x_column = 0;
    graph->num_y_columns = 1;
    graph->y_columns[0] = 1;
}

//--------------------------------------------------------------------

void generate_gnuplot_image(report_graph_t *graph, char *data_filename, char *prefix) {
    // lines specifying input data and output graph are declared and filled
    char line[2048];

    char gnuplot_filename[strlen(prefix) + 100];
    sprintf(gnuplot_filename, "%s.gnuplot", prefix);

    // open the file for writing gnuplot lines
    FILE* f = fopen(gnuplot_filename, "w");
    
    if (f == NULL) {
      LOG_FATAL("Impossible save file for BAM report");
    }

    fprintf(f, "set output '%s.png'\n", prefix);
    fprintf(f, "set terminal png nocrop enhanced font arial 10 size 640,360\n");
    fprintf(f, "set ylabel '%s'\n", graph->ylabel);
    fprintf(f, "set xlabel '%s'\n", graph->xlabel);
    fprintf(f, "set ytics border in scale 1,0.5 mirror norotate  offset character 0, 0, 0\n");
    fprintf(f, "set title '%s'\n", graph->title);

    if (graph->x_autoscale == 1) {
      fprintf(f, "%s", "set autoscale x\n");
    } else {
      fprintf(f, "set xrange [ %i : %i ] noreverse nowriteback\n", graph->x_start, graph->x_end);
    }

    if (graph->y_autoscale == 1) {
      fprintf(f, "%s", "set autoscale y\n");
    } else {
      fprintf(f, "set yrange [ %i : %i ] noreverse nowriteback\n", graph->x_start, graph->x_end);
    }

    fprintf(f, "set lmargin '%i'\n", graph->lmargin);
    fprintf(f, "set rmargin '%i'\n", graph->rmargin);
    fprintf(f, "set tmargin '%i'\n", graph->tmargin);
    fprintf(f, "set bmargin '%i'\n", graph->bmargin);

    sprintf(line, "%s", "plot ");

    for (int i = 0; i < graph->num_y_columns; i++) {
      sprintf(line, "%s%s '%s' using %i:%i title '%s' with %s", 
	      line, (i == 0 ? "" : ", "), data_filename, 
	      graph->x_column, graph->y_columns[i], 
	      graph->y_titles[i], graph->type);
    }
    fprintf(f, "%s\n", line);

    fclose(f);    
    
    // build the command line by calling gnuplot followed by is instruction file
    // execute command line: gnuplot filename.gnuplot
    sprintf(line, "gnuplot %s;", gnuplot_filename);
    int res = system(line);
}

//--------------------------------------------------------------------

void report_summary(char *prefix, stats_counters_t *output) {

  FILE *f; 
  char path[strlen(prefix) + 100];
  sprintf(path, "%s.summary.txt", prefix);

  if ( (f = fopen(path, "w")) == NULL) {
    LOG_FATAL_F("Impossible save stats summary (%s)for the BAM report", path);
  }

  fprintf(f, "\n");
  fprintf(f, "===============================================================\n");
  fprintf(f, "=            S T A T I S T I C S     S U M M A R Y            =\n");
  fprintf(f, "===============================================================\n");

  fprintf(f, "\nFilter:\n");
  if (output->num_passed != 0 || output->num_passed != 0) {
    fprintf(f, "\tEnabled\n");
    fprintf(f, "\tNumber of alignments in file  : %lu\n", output->num_passed + output->num_failed);
    fprintf(f, "\tNumber of processed alignments: %lu (%0.2f %%)\n", output->num_reads, 
	    100.0f * output->num_reads / (output->num_passed + output->num_failed));
  } else {
    fprintf(f, "\tDisabled\n");
    fprintf(f, "\tNumber of alignments in file  : %lu\n", output->num_reads);
  }

  fprintf(f, "\nCoverage:\n");
  fprintf(f, "\tReference length       : %lu\n", output->ref_length);
  size_t mapped_nts = output->ref_length - output->unmapped_nts;;
  fprintf(f, "\tMapped reference length: %lu (%0.2f %%)\n", 
	 mapped_nts, 100.0 * mapped_nts / output->ref_length);
  fprintf(f, "\n");
  fprintf(f, "\tAverage coverage (whole reference length)      : %0.2f\n", output->depth);
  fprintf(f, "\tAverage coverage (only mapped reference length): %0.2f\n", 
	 1.0f * output->depth * output->ref_length/ mapped_nts);

  fprintf(f, "\nAlignment (and read) counters:\n");
  fprintf(f, "\tNumber of alignments       : %lu\n", output->num_reads);
  fprintf(f, "\tNumber of unique alignments: %lu (%0.2f %%)\n", 
	 output->num_unique_alignments, 100.0 * output->num_unique_alignments / output->num_reads);

  fprintf(f, "\n\tMinimum alignment length: %lu\n", output->min_alignment_length);
  fprintf(f, "\tMaximum alignment length: %lu\n", output->max_alignment_length);
  fprintf(f, "\tAverage alignment length: %0.2f\n", 1.0f * output->num_nucleotides / output->num_reads);

  if (output->single_end) {
    fprintf(f, "\n\tMode: single\n");
  } else {
    fprintf(f, "\n\tMode: pair\n");
  }
  fprintf(f, "\tNumber of mapped reads     : %lu (%0.2f %%)\n", 
	 output->num_mapped_reads, 100.0 * output->num_mapped_reads / output->num_reads);
  fprintf(f, "\tNumber of unmapped reads   : %lu (%0.2f %%)\n",
	 output->num_unmapped_reads, 100.0 * output->num_unmapped_reads / output->num_reads);
  if (!output->single_end) {
    fprintf(f, "\tNumber of mapped reads in pair #1  : %lu (%0.2f %%)\n", 
	   output->num_mapped_reads_1, 100.0 * output->num_mapped_reads_1 / output->num_reads);
    fprintf(f, "\tNumber of unmapped reads in pair #1: %lu (%0.2f %%)\n", 
	   output->num_unmapped_reads_1, 100.0 * output->num_unmapped_reads_1 / output->num_reads);
    fprintf(f, "\tNumber of mapped reads in pair #2  : %lu (%0.2f %%)\n", 
	   output->num_mapped_reads_2, 100.0 * output->num_mapped_reads_2 / output->num_reads);
    fprintf(f, "\tNumber of unmapped reads in pair #2: %lu (%0.2f %%)\n", 
	   output->num_unmapped_reads_2, 100.0 * output->num_unmapped_reads_2 / output->num_reads);

    fprintf(f, "\n\tInsert size:\n");
    fprintf(f, "\t\tMinimum insert size: %lu\n", output->min_insert_size);
    fprintf(f, "\t\tMaximum insert size: %lu\n", output->max_insert_size);
    fprintf(f, "\t\tAverage insert size: %0.2f\n", 1.0f * output->insert_size_acc / output->num_reads);
  }

  fprintf(f, "\nNumber of errors (mismatches, indels):\n");
  for (int i = 0; i < NUM_ERRORS_STATS; i++) {
  fprintf(f, "\tNumber of reads with %i errors: %lu (%0.2f %%)\n", 
	  i, output->num_errors[i], 100.0 * output->num_errors[i] / output->num_mapped_reads);
  }
  fprintf(f, "\tNumber of reads with more than %i errors: %lu (%0.2f %%)\n", 
	  NUM_ERRORS_STATS, output->num_errors[NUM_ERRORS_STATS], 100.0 * output->num_errors[NUM_ERRORS_STATS] / output->num_mapped_reads);
  /*
  fprintf(f, "\tNumber of reads with 0 errors: %lu (%0.2f %)\n", 
	 output->num_0_errors, 100.0 * output->num_0_errors / output->num_mapped_reads);
  fprintf(f, "\tNumber of reads with 1 ~ 10 errors: %lu (%0.2f %)\n", 
	 output->num_1_10_errors, 100.0 * output->num_1_10_errors / output->num_mapped_reads);
  fprintf(f, "\tNumber of reads with more than 10 errors: %lu (%0.2f %)\n", 
	 output->num_11_errors, 100.0 * output->num_11_errors / output->num_mapped_reads);
  */
  fprintf(f, "\n\tNumber of indels: %lu\n", 
	 output->num_indels);
  fprintf(f, "\tAverage indel size: %0.2f\n", 
	 1.0f * output->indels_acc / output->num_indels);
  fprintf(f, "\tAverage indels per read: %0.2f\n", 
	 1.0f * output->num_indels / output->num_mapped_reads);

  fprintf(f, "\nPer strand:\n");
  fprintf(f, "Strand\tNum. unique alignments\tNum. mapped reads\n");
  fprintf(f, "---------------------------------------------------------------------------------------\n");
  for (int i = 0; i < 2; i++) {
    fprintf(f, "%c\t%lu (%0.2f %%)\t%lu (%0.2f %%)\n",
	   (i == 0 ? '+' : '-'), 
	   output->num_unique_alignments_strand[i], 100.0 * output->num_unique_alignments_strand[i] / output->num_mapped_reads,
	   output->num_mapped_reads_strand[i], 100.0 * output->num_mapped_reads_strand[i] / output->num_mapped_reads);
  }

  fprintf(f, "\nAlignment quality:\n");
  fprintf(f, "\tMinimum alignment quality: %lu\n", output->min_quality);
  fprintf(f, "\tMaximum alignment quality: %lu\n", output->max_quality);
  fprintf(f, "\tAverage alignment quality: %0.2f\n", 1.0f * output->quality_acc / output->num_mapped_reads);

  fprintf(f, "\nNucleotide content:\n");
  fprintf(f, "\tNumber of A's: %lu (%0.2f %%)\n", 
	 output->num_As, 100.0 * output->num_As / output->num_nucleotides);
  fprintf(f, "\tNumber of C's: %lu (%0.2f %%)\n", 
	 output->num_Cs, 100.0 * output->num_Cs / output->num_nucleotides);
  fprintf(f, "\tNumber of T's: %lu (%0.2f %%)\n", 
	 output->num_Ts, 100.0 * output->num_Ts / output->num_nucleotides); 
  fprintf(f, "\tNumber of G's: %lu (%0.2f %%)\n", 
	 output->num_Gs, 100.0 * output->num_Gs / output->num_nucleotides);
  fprintf(f, "\tNumber of N's: %lu (%0.2f %%)\n", 
	 output->num_Ns, 100.0 * output->num_Ns / output->num_nucleotides);
  fprintf(f, "\tGC percentage: %0.2f %%\n", 
	 100.0 * (output->num_Cs + output->num_Gs) / output->num_nucleotides);

  fprintf(f, "\nPer sequences (chromosomes) statistics:\n");
  fprintf(f, "Number of sequences: %i\n", output->num_sequences);
  fprintf(f, "Name\tLength\tCoverage\n");
  fprintf(f, "--------------------------------------\n");
  for (int i = 0; i < output->num_sequences; i++) {
    fprintf(f, "%s\t%lu\t%0.2f\n", output->sequence_labels[i], output->sequence_lengths[i], output->depth_per_sequence[i]);
  }

  if (f != stdout) fclose(f);
}

//--------------------------------------------------------------------

void report_summary_sqlite3(sqlite3 *db, stats_counters_t *output) {

  sqlite3_stmt* stmt;
  char aux[256], aux1[256], aux2[256];

  prepare_statement_global_stats(db, &stmt);
  
  char* errorMessage;
  sqlite3_exec(db, "BEGIN TRANSACTION", NULL, NULL, &errorMessage);

  sprintf(aux, "%lu", output->ref_length);
  insert_statement_global_stats("REF_LENGTH", "Reference length", aux, stmt, db);

  // ref. length
  size_t mapped_nts = output->ref_length - output->unmapped_nts;;
  sprintf(aux, "%lu", mapped_nts);
  insert_statement_global_stats("MAPPED_REF_LENGTH", "Mapped reference length", aux, stmt, db);
  sprintf(aux, "%0.2f", 100.0 * mapped_nts / output->ref_length);
  insert_statement_global_stats("MAPPED_REF_LENGTH_PERC", "Mapped reference length (%)", aux, stmt, db);

  // coverage
  sprintf(aux, "%0.2f", output->depth);
  insert_statement_global_stats("AVG_COVERAGE_WHOLE_REF", "Average coverage (whole reference length)", aux, stmt, db);

  sprintf(aux, "%0.2f", 1.0f * output->depth * output->ref_length/ mapped_nts);
  insert_statement_global_stats("AVG_COVERAGE_MAPPED_REF", "Average coverage (only mapped reference length)", aux, stmt, db);

  // num. alignments, unique
  sprintf(aux, "%lu", output->num_reads);
  insert_statement_global_stats("NUM_ALIGNMENTS", "Number of alignments", aux, stmt, db);

  sprintf(aux, "%lu", output->num_unique_alignments);
  insert_statement_global_stats("NUM_UNIQUE_ALIGNMENTS", "Number of unique alignments", aux, stmt, db);

  sprintf(aux, "%02.f", 100.0 * output->num_unique_alignments / output->num_reads);
  insert_statement_global_stats("NUM_UNIQUE_ALIGNMENTS_PERC", "Number of unique alignments (%)", aux, stmt, db);

  // alignment length
  sprintf(aux, "%lu", output->min_alignment_length);
  insert_statement_global_stats("MIN_ALIGNMENT_LENGTH", "Minimum alignment length", aux, stmt, db);

  sprintf(aux, "%lu", output->max_alignment_length);
  insert_statement_global_stats("MAX_ALIGNMENT_LENGTH", "Maximum alignment length", aux, stmt, db);

  sprintf(aux, "%0.2f", 1.0f * output->num_nucleotides / output->num_reads);
  insert_statement_global_stats("AVG_ALIGNMENT_LENGTH", "Average alignment length", aux, stmt, db);

  // mode
  insert_statement_global_stats("MODE", "Mode", (output->single_end ? "single" : "pair"), stmt, db);

			   /*
  fprintf(f, "\tNumber of mapped reads     : %lu (%0.2f %)\n", 
	 output->num_mapped_reads, 100.0 * output->num_mapped_reads / output->num_reads);
  fprintf(f, "\tNumber of unmapped reads   : %lu (%0.2f %)\n",
	 output->num_unmapped_reads, 100.0 * output->num_unmapped_reads / output->num_reads);
  if (!output->single_end) {
    fprintf(f, "\tNumber of mapped reads in pair #1  : %lu (%0.2f %)\n", 
	   output->num_mapped_reads_1, 100.0 * output->num_mapped_reads_1 / output->num_reads);
    fprintf(f, "\tNumber of unmapped reads in pair #1: %lu (%0.2f %)\n", 
	   output->num_unmapped_reads_1, 100.0 * output->num_unmapped_reads_1 / output->num_reads);
    fprintf(f, "\tNumber of mapped reads in pair #2  : %lu (%0.2f %)\n", 
	   output->num_mapped_reads_2, 100.0 * output->num_mapped_reads_2 / output->num_reads);
    fprintf(f, "\tNumber of unmapped reads in pair #2: %lu (%0.2f %)\n", 
	   output->num_unmapped_reads_2, 100.0 * output->num_unmapped_reads_2 / output->num_reads);
    fprintf(f, "\n");
			   */

  // insert
  sprintf(aux, "%lu", output->min_insert_size);
  insert_statement_global_stats("MIN_INSERT_LENGTH", "Minimum insert length", aux, stmt, db);

  sprintf(aux, "%lu", output->max_insert_size);
  insert_statement_global_stats("MAX_INSERT_LENGTH", "Maximum insert length", aux, stmt, db);

  sprintf(aux, "%0.2f", 1.0f * output->insert_size_acc / output->num_reads);
  insert_statement_global_stats("AVG_INSERT_LENGTH", "Average insert length", aux, stmt, db);

			  
  // errors
  for (int i = 0; i < NUM_ERRORS_STATS; i++) {
    sprintf(aux, "%lu", output->num_errors[i]);
    sprintf(aux1, "NUM_READS_WITH_%i_ERRORS", i);
    sprintf(aux2, "Number of reads with %i errors", i);
    insert_statement_global_stats(aux1, aux2, aux, stmt, db);

    sprintf(aux, "%0.2f", 100.0 * output->num_errors[i] / output->num_mapped_reads);
    sprintf(aux1, "NUM_READS_WITH_%i_ERRORS_PERC", i);
    sprintf(aux2, "Number of reads with %i errors (%%)", i);
    insert_statement_global_stats(aux1, aux2, aux, stmt, db);
  }

  sprintf(aux, "%lu", output->num_errors[NUM_ERRORS_STATS]);
  sprintf(aux1, "NUM_READS_WITH_MORE_%i_ERRORS", NUM_ERRORS_STATS);
  sprintf(aux2, "Number of reads with more than %i errors", NUM_ERRORS_STATS);
  insert_statement_global_stats(aux1, aux2, aux, stmt, db);

  sprintf(aux, "%0.2f", 100.0 * output->num_errors[NUM_ERRORS_STATS] / output->num_mapped_reads);
  sprintf(aux1, "NUM_READS_WITH_MORE_%i_ERRORS_PERC", NUM_ERRORS_STATS);
  sprintf(aux2, "Number of reads with %i errors (%%)", NUM_ERRORS_STATS);
  insert_statement_global_stats(aux1, aux2, aux, stmt, db);

  // indels
  sprintf(aux, "%lu", output->num_indels);
  insert_statement_global_stats("NUM_INDELS", "Number of indels", aux, stmt, db);
  
  sprintf(aux, "%0.2f", 1.0f * output->indels_acc / output->num_indels);
  insert_statement_global_stats("AVG_INDELS_LENGTH", "Average indels length", aux, stmt, db);

  sprintf(aux, "%0.2f", 1.0f * output->num_indels / output->num_mapped_reads);
  insert_statement_global_stats("AVG_INDELS_PER_READ", "Average indels per read", aux, stmt, db);
  
  // per forward strand
  sprintf(aux, "%lu", output->num_unique_alignments_strand[0]);
  insert_statement_global_stats("NUM_UNIQUE_ALIGNMENTS_FORWARD_STRAND", "Number of unique alignments in forward strand", aux, stmt, db);

  sprintf(aux, "%0.2f", 100.0 * output->num_unique_alignments_strand[0] / output->num_mapped_reads);
  insert_statement_global_stats("NUM_UNIQUE_ALIGNMENTS_FORWARD_STRAND_PERC", "Number of unique alignments in forward strand (%)", aux, stmt, db);

  sprintf(aux, "%lu", output->num_mapped_reads_strand[0]);
  insert_statement_global_stats("NUM_ALIGNMENTS_FORWARD_STRAND", "Number of alignments in forward strand", aux, stmt, db);

  sprintf(aux, "%0.2f", 100.0 * output->num_mapped_reads_strand[0] / output->num_mapped_reads);
  insert_statement_global_stats("NUM_ALIGNMENTS_FORWARD_STRAND_PERC", "Number of alignments in forward strand (%)", aux, stmt, db);

  // per reverse strand
  sprintf(aux, "%lu", output->num_unique_alignments_strand[1]);
  insert_statement_global_stats("NUM_UNIQUE_ALIGNMENTS_REVERSE_STRAND", "Number of unique alignments in reverse strand", aux, stmt, db);

  sprintf(aux, "%0.2f", 100.0 * output->num_unique_alignments_strand[1] / output->num_mapped_reads);
  insert_statement_global_stats("NUM_UNIQUE_ALIGNMENTS_REVERSE_STRAND_PERC", "Number of unique alignments in reverse strand (%)", aux, stmt, db);

  sprintf(aux, "%lu", output->num_mapped_reads_strand[1]);
  insert_statement_global_stats("NUM_ALIGNMENTS_REVERSE_STRAND", "Number of alignments in reverse strand", aux, stmt, db);

  sprintf(aux, "%0.2f", 100.0 * output->num_mapped_reads_strand[1] / output->num_mapped_reads);
  insert_statement_global_stats("NUM_ALIGNMENTS_REVERSE_STRAND_PERC", "Number of alignments in reverse strand (%)", aux, stmt, db);

  // alignment quality
  sprintf(aux, "%lu", output->min_quality);
  insert_statement_global_stats("MIN_ALIGNMENT_QUALITY", "Minimum alignment quality", aux, stmt, db);

  sprintf(aux, "%lu", output->max_quality);
  insert_statement_global_stats("MAX_ALIGNMENT_QUALITY", "Maximum alignment quality", aux, stmt, db);

  sprintf(aux, "%0.2f", 1.0f * output->quality_acc / output->num_mapped_reads);
  insert_statement_global_stats("AVG_ALIGNMENT_QUALITY", "Average alignment quality", aux, stmt, db);

  // nucleotide content
  sprintf(aux, "%lu", output->num_As);
  insert_statement_global_stats("NUM_A", "Number of A", aux, stmt, db);

  sprintf(aux, "%0.2f", 100.0 * output->num_As / output->num_nucleotides);
  insert_statement_global_stats("NUM_A_PERC", "Number of A (%)", aux, stmt, db);

  sprintf(aux, "%lu", output->num_Cs);
  insert_statement_global_stats("NUM_C", "Number of C", aux, stmt, db);

  sprintf(aux, "%0.2f", 100.0 * output->num_Cs / output->num_nucleotides);
  insert_statement_global_stats("NUM_C_PERC", "Number of C (%)", aux, stmt, db);

  sprintf(aux, "%lu", output->num_Ts);
  insert_statement_global_stats("NUM_T", "Number of T", aux, stmt, db);

  sprintf(aux, "%0.2f", 100.0 * output->num_Ts / output->num_nucleotides);
  insert_statement_global_stats("NUM_T_PERC", "Number of T (%)", aux, stmt, db);

  sprintf(aux, "%lu", output->num_Gs);
  insert_statement_global_stats("NUM_G", "Number of G", aux, stmt, db);

  sprintf(aux, "%0.2f", 100.0 * output->num_Gs / output->num_nucleotides);
  insert_statement_global_stats("NUM_G_PERC", "Number of G (%)", aux, stmt, db);

  sprintf(aux, "%lu", output->num_Ns);
  insert_statement_global_stats("NUM_N", "Number of N", aux, stmt, db);

  sprintf(aux, "%0.2f", 100.0 * output->num_Ns / output->num_nucleotides);
  insert_statement_global_stats("NUM_N_PERC", "Number of N (%)", aux, stmt, db);

  sprintf(aux, "%0.2f", 100.0 * (output->num_Cs + output->num_Gs) / output->num_nucleotides);
  insert_statement_global_stats("NUM_GC_PERC", "GC (%)", aux, stmt, db);

  // per sequences (chromosomes) statistics
  sprintf(aux, "%i", output->num_sequences);
  insert_statement_global_stats("NUM_CHR", "Number of chromosomes", aux, stmt, db);
  
  int len = output->num_sequences;
  for (int i = 0; i < output->num_sequences; i++) {
    len += strlen(output->sequence_labels[i]);
  }
  char chr_list[len];
  for (int i = 0; i < output->num_sequences; i++) {
    sprintf(chr_list, "%s%s%s", chr_list, (i > 0 ? "," : ""), output->sequence_labels[i]);
  }
  chr_list[len] = '\0';
  insert_statement_global_stats("CHR_LIST", "List of chromosomes", chr_list, stmt, db);
  
  for (int i = 0; i < output->num_sequences; i++) {
    sprintf(aux, "%lu", output->sequence_lengths[i]);
    sprintf(aux1, "CHR_%s_LENGTH", output->sequence_labels[i]);
    sprintf(aux2, "Chromosome %s length", output->sequence_labels[i]);
    insert_statement_global_stats(aux1, aux2, aux, stmt, db);

    sprintf(aux, "%0.2f", output->depth_per_sequence[i]);
    sprintf(aux1, "CHR_%s_COVERAGE", output->sequence_labels[i]);
    sprintf(aux2, "Chromosome %s coverage", output->sequence_labels[i]);
    insert_statement_global_stats(aux1, aux2, aux, stmt, db);
  }

  sqlite3_exec(db, "COMMIT TRANSACTION", NULL, NULL, &errorMessage);

  finalize_statement_global_stats(stmt);
}

//--------------------------------------------------------------------

void report_coverage(char *prefix, stats_counters_t *output) {

  FILE *f; 
  int name_length = strlen(prefix) + 100;
  char data_filename[name_length];
  char img_prefix[name_length];

  int num_cols = 100;

  // coverage histogram
  int hist[num_cols];
  memset(hist, 0, num_cols * sizeof(int));

  size_t index, len;
  for (int i = 0; i < output->num_sequences; i++) {
    len = output->sequence_lengths[i];
    for (int j = 0; j < len; j++) {
      index = output->sequence_depths_per_nt[i][j];
      if (index >= num_cols) {
	index = num_cols - 1;
      }
      hist[index]++;
    }
  }

  // coverage histogram data
  sprintf(data_filename, "%s.coverage.histogram.data", prefix);
  if ( (f = fopen(data_filename, "w")) == NULL) {
    LOG_FATAL("Impossible save stats summary for the BAM report");
  }

  for (int i = 0; i < num_cols; i++) {
    fprintf(f, "%i\t%i\n", i, hist[i]);
  }
  fclose(f);

  // coverage histogram image
  sprintf(img_prefix, "%s.coverage.histogram", prefix);
  report_graph_t graph;
  init_report_graph(&graph);
  
  strcpy(graph.title, "Coverage Histogram");
  strcpy(graph.xlabel, "Coverage");
  strcpy(graph.ylabel, "Number of genomic locations");
  strcpy(graph.type, "boxes");
  graph.x_autoscale = 0;
  graph.x_start = 0;
  graph.x_end = num_cols;
  graph.y_autoscale = 1;
  graph.x_column = 1;
  graph.num_y_columns = 1;
  graph.y_columns[0] = 2;
  
  generate_gnuplot_image(&graph, data_filename, img_prefix);

  // coverage fraction
  double acc_hist[num_cols];
  memset(acc_hist, 0, num_cols * sizeof(double));

  size_t acc = 0;
  for (int i = num_cols - 1; i > 0; i--) {
    acc += hist[i];
    acc_hist[i] = 100.0f * acc / output->ref_length;
  }

  // coverage fraction data
  sprintf(data_filename, "%s.coverage.fraction.data", prefix);
  if ( (f = fopen(data_filename, "w")) == NULL) {
    LOG_FATAL("Impossible save stats summary for the BAM report");
  }

  for (int i = 0; i < num_cols; i++) {
    fprintf(f, "%i\t%0.2f\n", i, acc_hist[i]);
  }
  fclose(f);

  // coverage histogram image
  sprintf(img_prefix, "%s.coverage.fraction", prefix);
  init_report_graph(&graph);
  
  strcpy(graph.title, "Genome Fraction Coverage");
  strcpy(graph.xlabel, "Coverage");
  strcpy(graph.ylabel, "Fraction of refence (%)");
  strcpy(graph.type, "boxes");
  graph.x_autoscale = 0;
  graph.x_start = 0;
  graph.x_end = num_cols;
  graph.y_autoscale = 1;
  graph.x_column = 1;
  graph.num_y_columns = 1;
  graph.y_columns[0] = 2;
  
  generate_gnuplot_image(&graph, data_filename, img_prefix);
}

//--------------------------------------------------------------------

void report_num_errors(char *prefix, stats_counters_t *output) {

  FILE *f; 
  int name_length = strlen(prefix) + 100;
  char data_filename[name_length];
  char img_prefix[name_length];

  int num_cols = NUM_ERRORS_STATS + 1;

  // num. errors histogram
  int hist[num_cols];
  memset(hist, 0, num_cols * sizeof(int));

  // num. errors histogram data
  sprintf(data_filename, "%s.num.errors.histogram.data", prefix);
  if ( (f = fopen(data_filename, "w")) == NULL) {
    LOG_FATAL("Impossible save stats summary for the BAM report");
  }

  for (int i = 0; i < num_cols; i++) {
    fprintf(f, "%i\t%lu\n", i, output->num_errors[i]);
  }
  fclose(f);

  // num. errors histogram image
  sprintf(img_prefix, "%s.num.errors.histogram", prefix);
  report_graph_t graph;
  init_report_graph(&graph);
  
  strcpy(graph.title, "Num. Errors Histogram");
  strcpy(graph.xlabel, "Num. errors");
  strcpy(graph.ylabel, "Number of reads");
  strcpy(graph.type, "boxes");
  graph.x_autoscale = 0;
  graph.x_start = 0;
  graph.x_end = num_cols;
  graph.y_autoscale = 1;
  graph.x_column = 1;
  graph.num_y_columns = 1;
  graph.y_columns[0] = 2;
  
  generate_gnuplot_image(&graph, data_filename, img_prefix);
}

//--------------------------------------------------------------------

void report_quality(char *prefix, stats_counters_t *output) {

  FILE *f; 
  int name_length = strlen(prefix) + 100;
  char data_filename[name_length];
  char img_prefix[name_length];

  int num_cols = QUALITY_STATS;

  // quality histogram
  int hist[num_cols];
  memset(hist, 0, num_cols * sizeof(int));

  // quality histogram data
  sprintf(data_filename, "%s.quality.histogram.data", prefix);
  if ( (f = fopen(data_filename, "w")) == NULL) {
    LOG_FATAL("Impossible save stats summary for the BAM report");
  }

  for (int i = 0; i < num_cols; i++) {
    fprintf(f, "%i\t%lu\n", i, output->quality[i]);
  }
  fclose(f);

  // quality histogram image
  sprintf(img_prefix, "%s.quality.histogram", prefix);
  report_graph_t graph;
  init_report_graph(&graph);
  
  strcpy(graph.title, "Quality Histogram");
  strcpy(graph.xlabel, "Quality");
  strcpy(graph.ylabel, "Number of reads");
  strcpy(graph.type, "boxes");
  graph.x_autoscale = 0;
  graph.x_start = 0;
  graph.x_end = num_cols;
  graph.y_autoscale = 1;
  graph.x_column = 1;
  graph.num_y_columns = 1;
  graph.y_columns[0] = 2;
  
  generate_gnuplot_image(&graph, data_filename, img_prefix);
}

//--------------------------------------------------------------------

void report_gc_content(char *prefix, stats_counters_t *output) {

  FILE *f; 
  int name_length = strlen(prefix) + 100;
  char data_filename[name_length];
  char img_prefix[name_length];

  int num_cols = 100;

  // GC content histogram
  int hist[num_cols];
  memset(hist, 0, num_cols * sizeof(int));

  // GC content histogram data
  sprintf(data_filename, "%s.gc.content.histogram.data", prefix);
  if ( (f = fopen(data_filename, "w")) == NULL) {
    LOG_FATAL("Impossible save stats summary for the BAM report");
  }

  for (int i = 0; i < num_cols; i++) {
    fprintf(f, "%i\t%lu\n", i + 1, output->GC_content[i]);
  }
  fclose(f);

  // GC content histogram image
  sprintf(img_prefix, "%s.gc.content.histogram", prefix);
  report_graph_t graph;
  init_report_graph(&graph);
  
  strcpy(graph.title, "GC Content Distribution");
  strcpy(graph.xlabel, "GC content (%)");
  strcpy(graph.ylabel, "Number of reads");
  strcpy(graph.type, "lines");
  graph.x_autoscale = 0;
  graph.x_start = 0;
  graph.x_end = num_cols;
  graph.y_autoscale = 1;
  graph.x_column = 1;
  graph.num_y_columns = 1;
  graph.y_columns[0] = 2;
  
  generate_gnuplot_image(&graph, data_filename, img_prefix);
}

//--------------------------------------------------------------------

void report_stats(char *bam_filename, char *outdir, 
		  void *db, stats_counters_t *stats) {

  int len = strlen(bam_filename) + strlen(outdir) + 100;

  char filename[len], prefix[len], *p;
  sprintf(filename, "%s", ((p = strrchr(bam_filename, '/')) ? (p+1) : bam_filename));
  sprintf(prefix, "%s/%s", outdir, filename);

  report_summary(prefix, stats);

  if (db) {
    report_summary_sqlite3((sqlite3 *)db, stats);
  }
  
  report_coverage(prefix, stats);
  report_num_errors(prefix, stats);
  report_quality(prefix, stats);
  report_gc_content(prefix, stats);

  //  report_insert(prefix, stats);
  //  report_nucleotides(prefix, stats);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
/*


*/
