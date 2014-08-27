#ifndef WORKFLOW_SCHEDULER_SA_H
#define WORKFLOW_SCHEDULER_SA_H

#include <sys/syscall.h>
#include <sched.h>
#include <stdio.h>
#include <stdlib.h>
//#include <omp.h>
#include <unistd.h>
#include <pthread.h>


#include "commons/commons.h"
#include "containers/array_list.h"
//#include "extrae_user_events.h"

//----------------------------------------------------------------------------------------

#define WORKFLOW_STATUS_RUNNING  1
#define WORKFLOW_STATUS_FINISHED 2

//----------------------------------------------------------------------------------------
// work_item
//----------------------------------------------------------------------------------------

typedef struct work_item_SA {
  int stage_id;
  void *data;

  void *context;
} work_item_SA_t;

work_item_SA_t *work_item_SA_new(int stage_id, void *data);
void work_item_SA_free(work_item_SA_t *item);

typedef int (*workflow_stage_function_SA_t) (void *data);

typedef void* (*workflow_producer_function_SA_t) (void *data);
typedef int (*workflow_consumer_function_SA_t) (void *data);

//----------------------------------------------------------------------------------------
// workflow
//----------------------------------------------------------------------------------------

typedef struct workflow_SA workflow_SA_t;

struct workflow_SA {
  int num_threads;
  int max_num_work_items;
  int num_stages;
  int completed_producer;
  int num_processing_items;
  int num_pending_items;     
  int running_producer;
  int running_consumer;
  
  int complete_extra_stage;
  
  pthread_cond_t producer_cond;
  pthread_cond_t consumer_cond;
  pthread_cond_t workers_cond;
  
  pthread_mutex_t producer_mutex;
  pthread_mutex_t consumer_mutex;     
  
  pthread_mutex_t main_mutex;
  
  double workflow_time;
  double producer_time;
  double consumer_time;
  double *stage_times;

  pthread_mutex_t *stage_times_mutex;     

  array_list_t **pending_items;
  array_list_t *completed_items;
  
  workflow_stage_function_SA_t *stage_functions;
  char** stage_labels;
  
  workflow_producer_function_SA_t *producer_function;
  char* producer_label;
  
  
  workflow_consumer_function_SA_t *consumer_function;
  char* consumer_label;
  //int (*status_function)(workflow_t *);
};

workflow_SA_t *workflow_SA_new();
void workflow_SA_free(workflow_SA_t *wf);

void workflow_set_stages_SA(int num_stages, workflow_stage_function_SA_t *functions, 
			    char **labels, workflow_SA_t *wf);
void workflow_set_producer_SA(workflow_producer_function_SA_t *function, 
			     char *label, workflow_SA_t *wf);
void workflow_set_consumer_SA(workflow_consumer_function_SA_t *function, 
			      char *label, workflow_SA_t *wf);

int workflow_get_num_items_SA(workflow_SA_t *wf);
int workflow_get_num_items_at_SA(int stage_id, workflow_SA_t *wf);
int workflow_get_num_completed_items_SA(workflow_SA_t *wf);
int workflow_is_producer_finished_SA(workflow_SA_t *wf);

void workflow_insert_item_SA(void *data, workflow_SA_t *wf);
void workflow_insert_item_at_SA(int stage_id, void *data, workflow_SA_t *wf);
void workflow_insert_stage_item_at_SA(void *data, int new_stage, workflow_SA_t *wf);
void *workflow_remove_item_SA(workflow_SA_t *wf);
void *workflow_remove_item_at_SA(int stage_id, workflow_SA_t *wf);

int workflow_get_status_SA(workflow_SA_t *wf);
int workflow_get_simple_status_SA(workflow_SA_t *wf);
void workflow_producer_finished_SA(workflow_SA_t *wf);
void workflow_change_status_function_SA(workflow_SA_t *wf);

void workflow_run_SA(void *input, workflow_SA_t *wf);
void workflow_run_with_SA(int num_threads, void *input, workflow_SA_t *wf);

void workflow_display_timing_SA(workflow_SA_t *wf);

//----------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------

#endif
