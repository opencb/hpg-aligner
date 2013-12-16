#ifndef STATISTICS_H
#define STATISTICS_H

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "timing.h"

typedef struct statistics{
  int num_sections;
  int num_subsections;
  pthread_mutex_t mutex;
  unsigned int *num_values;
  char **section_sublabels_p;
  char **section_labels_p;
  size_t **values_p;
}statistics_t;

typedef struct basic_statistics {
  size_t total_reads;
  size_t num_mapped_reads;
  size_t total_mappings;
  size_t total_sp;
  size_t uniq_sp;
  pthread_mutex_t mutex;
} basic_statistics_t;

typedef struct cal_st {
  size_t max_cals;
  size_t no_cals;
  pthread_mutex_t mutex;
} cal_st_t;

extern cal_st_t cal_st;
extern char statistics_on;
extern statistics_t *statistics_p;
extern basic_statistics_t *basic_st;

statistics_t* statistics_new(char** section_labels_p, char** section_sublabels_p, unsigned int *num_values, int num_sections, int num_subsection);

void statistics_free(statistics_t* statistics_p);

void statistics_set(unsigned int section, unsigned int subsection, size_t value, statistics_t* statistics_p);

void statistics_add(unsigned int section, unsigned int subsection, size_t value, statistics_t* statistics_p);

void statistics_display(statistics_t* statistics_p);

void basic_statistics_display(basic_statistics_t *statistics, int rna_mode, float alig_time, float load_time);

void basic_statistics_init(size_t total_reads, size_t num_mapped_reads, size_t total_mappings, basic_statistics_t *statistics);

void basic_statistics_sp_init(size_t total_sp, size_t uniq_sp, basic_statistics_t *statistics);

void timing_and_statistics_display(statistics_t* statistics_p, 
				   timing_t* timing_p);

void basic_statistics_add(size_t total_reads, size_t num_mapped_reads, size_t total_mappings, basic_statistics_t *basic);

basic_statistics_t *basic_statistics_new();

#endif // end of if TIMING
