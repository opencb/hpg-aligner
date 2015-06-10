#ifndef TIMING_H
#define TIMING_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h>


typedef struct timing {
  int num_sections;
  char** section_labels;
  double* section_times;
  pthread_mutex_t time_mutex;
} timing_t;


timing_t* timing_new(char** section_labels, int num_sections);
void timing_free(timing_t* timing_p);
void timing_display(timing_t* timing_p);
void timing_add(double time, size_t section_id, timing_t *t_p);

extern char time_on;
extern timing_t *timing;

// for debugging
extern double bwt_time[100], seeding_time[100], cal_time[100], sw_time[100];
extern int thr_bwt_items[100], thr_seeding_items[100], thr_cal_items[100], thr_sw_items[100];
extern int thr_batches[100];

#endif // end of if TIMING






