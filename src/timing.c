#include "timing.h"

//-----------------------------------------------------

double bwt_time[100], seeding_time[100], cal_time[100], sw_time[100];
int thr_bwt_items[100], thr_seeding_items[100], thr_cal_items[100], thr_sw_items[100];
int thr_batches[100];

//-----------------------------------------------------

timing_t* timing_new(char** section_labels, int num_sections) {

  // for debugging
  timing_t* t_p = (timing_t*) calloc(1, sizeof(timing_t));

  t_p->num_sections = num_sections;
  t_p->section_labels = (char**) calloc(num_sections, sizeof(char*));
  t_p->section_times = (double*) calloc(num_sections, sizeof(double));

  for(int i = 0 ; i < num_sections ; i++) {
    t_p->section_labels[i] = strdup(section_labels[i]);
  }
  pthread_mutex_init(&t_p->time_mutex, NULL);
  return t_p;
} 

//-----------------------------------------------------

void timing_free(timing_t* t_p) {
  if (t_p == NULL) return;

  for(int i = 0 ; i < t_p->num_sections ; i++) {
    if (t_p->section_labels[i] != NULL) { free(t_p->section_labels[i]); } 
  }

  if (t_p->section_labels != NULL) { free(t_p->section_labels); }
  if (t_p->section_times != NULL) { free(t_p->section_times); }
  
  free(t_p);
} 

//-----------------------------------------------------

void timing_add(double time, size_t section_id, timing_t *t_p) {
  pthread_mutex_lock(&t_p->time_mutex);
  t_p->section_times[section_id] += time;
  pthread_mutex_unlock(&t_p->time_mutex);
}

//-----------------------------------------------------

void timing_display(timing_t* t_p) {
  if (t_p == NULL) return;
  printf("\n");
  printf("===========================================================\n");
  printf("=                    T i m i n g                          =\n");
  printf("===========================================================\n");
  for(int i=0 ; i < t_p->num_sections ; i++) {
    if(i == t_p->num_sections - 1) {
      printf("-----------------------------------------------------------\n");
    }
    printf("%s\t: %4.04f sec\n", t_p->section_labels[i], t_p->section_times[i] / 1000000);
  }
  printf("===========================================================\n");
  printf("===========================================================\n");
  printf("\n");
}

//-----------------------------------------------------

