#include "statistics.h"

//---------------------------------------------------------------------------------------------------------------------

statistics_t* statistics_new(char** section_labels_p, char** section_sublabels_p, unsigned int *num_values, int num_sections, int num_subsections) {
  statistics_t* s_p = (statistics_t*)calloc(1, sizeof(statistics_t));
  
  pthread_mutex_init(&(s_p->mutex), NULL);
  
  s_p->num_sections = num_sections;
  s_p->num_subsections = num_subsections;
  s_p->section_labels_p = (char**) calloc(num_sections, sizeof(char*));
  s_p->section_sublabels_p = (char**) calloc(num_subsections, sizeof(char*));
  s_p->values_p   = (size_t**) calloc(num_sections, sizeof(size_t*));
  s_p->num_values = (unsigned int *)calloc(num_sections, sizeof(unsigned int));
  int i;

  for (i = 0; i < num_sections; i++) {
    s_p->section_labels_p[i] = strdup(section_labels_p[i]);
    s_p->values_p[i] = (size_t*)calloc(num_values[i], sizeof(size_t));
    s_p->num_values[i] = num_values[i];
  }

  for ( i = 0; i < num_subsections; i++) {
     s_p->section_sublabels_p[i] = strdup(section_sublabels_p[i]);
  }
  return s_p;
} 

//----------------------------------------------------------------------------------------------------------------------

void statistics_free(statistics_t* statistics_p) {
    for (int i = 0; i < statistics_p->num_sections; i++) {
      if (statistics_p->section_labels_p[i] != NULL) { free(statistics_p->section_labels_p[i]); } 
    }
    for (int i = 0; i < statistics_p->num_subsections; i++) {
      if (statistics_p->section_sublabels_p[i] != NULL) { free(statistics_p->section_sublabels_p[i]); } 
    }
}

//----------------------------------------------------------------------------------------------------------------------

void statistics_set(unsigned int section, unsigned int subsection, size_t value, statistics_t* statistics_p) {
  statistics_p->values_p[section][subsection] = value;
}

//----------------------------------------------------------------------------------------------------------------------

void statistics_add(unsigned int section, unsigned int subsection, size_t value, statistics_t* statistics_p) {
  pthread_mutex_lock(&statistics_p->mutex);
  statistics_p->values_p[section][subsection] += value;
  pthread_mutex_unlock(&statistics_p->mutex);
}

//---------------------------------------------------------------------------------------------------------------------

void statistics_display(statistics_t* statistics_p) {
  if (statistics_p == NULL) return;

  printf("\n");
  printf("===========================================================\n");
  printf("=                  S t a t i s t i c s                    =\n");
  printf("===========================================================\n");
  int i, j, s, start = 0, end = 0;
  for (i = 0; i < statistics_p->num_sections; i++) {
      printf("-----------------------------------------------------------\n");
      printf("%s\n", statistics_p->section_labels_p[i]);
      printf("-----------------------------------------------------------\n");
      end += statistics_p->num_values[i];
      s = 0;
      for (j = start; j < end; j++) {
	printf("\t%s:%-7lu\n", statistics_p->section_sublabels_p[j], statistics_p->values_p[i][s]);
	s++;
      }
      //printf("\n");
      start = j;
  }
  
  printf("===========================================================\n");
  printf("===========================================================\n");
  printf("\n");
}

//-------------------------------------------------------------------------------------------------------------------

void timing_and_statistics_display(statistics_t* statistics_p, timing_t *t_p) {

  /*  if (t_p == NULL) return;

  printf("\n");
  printf("=============================================================================\n");
  printf("=                  O P E R A T I O N S   N U M B E R                        =\n");
  printf("=============================================================================\n");
  int i, j;

  for (i = 3 ; i < t_p->num_sections - 3; i++) {
    t_p->section_times_p[i] = 0;
    for (j = 0; j < t_p->num_threads_p[i]; j++) {
      if (t_p->section_times_p[i] < t_p->thread_times_p[i][j]) {
	t_p->section_times_p[i] = t_p->thread_times_p[i][j];
      }
    }

    if ( i ==  (t_p->num_sections - 4)) {
      printf("%s\t: %lu/%4.04f = %4.04f(op/sec)\n", t_p->section_labels_p[i], statistics_p->values_p[i - 2][4], t_p->section_times_p[i]/1000000, 
	     (statistics_p->values_p[i - 2][4] / (t_p->section_times_p[i]/1000000)) );
    } else {
      printf("%s\t: %lu/%4.04f = %4.04f(op/sec)\n", t_p->section_labels_p[i], statistics_p->values_p[i- 2][1], t_p->section_times_p[i]/1000000, 
	     (statistics_p->values_p[i][1] / (t_p->section_times_p[i]/1000000)) );
    }
  }
  printf("=============================================================================\n");
  printf("=============================================================================\n");
  printf("\n");
  */
}

//---------------------------------------------------------------------------------------

void basic_statistics_display(basic_statistics_t *statistics, int rna_mode, float alig_time, float load_time) {
  size_t total_reads = statistics->total_reads;
  size_t num_mapped_reads = statistics->num_mapped_reads;
  size_t total_mappings = statistics->total_mappings;

  size_t total_sp = statistics->total_sp;
  size_t uniq_sp = statistics->uniq_sp;

  /*  printf("-------------------------------------------------\n");
  printf("-                MAPPING STATISTICS             -\n");
  printf("-------------------------------------------------\n");
  printf("Num. total reads   : %lu\n", total_reads);
  printf("    Mapped reads   : %lu\t(%0.2f %)\n", num_mapped_reads, num_mapped_reads * 100.0 / total_reads);
  printf("    Unmapped reads : %lu\t(%0.2f %)\n", total_reads - num_mapped_reads, (total_reads - num_mapped_reads) * 100.0 / total_reads);
  printf("Num. total mappings: %lu\n", total_mappings);
  if (rna_mode) {
    printf("Total Splice Junctions: %lu\n", total_sp);
    printf("    Differents Splice Junctions: %lu\n", uniq_sp);
  }
  printf("-------------------------------------------------\n");*/
  printf("+--------------------------------------------------------------------------------------+\n");
  printf("|                                     GLOBAL STATISTICS                                |\n");
  printf("+--------------------------------------------------------------------------------------+\n");
  printf("| Loading Time (s)  : %-65.2f", load_time);
  printf("|\n");
  printf("| Alignment Time (s): %-65.2f", alig_time);
  printf("|\n");
  printf("| Total Time (s)    : %-65.2f", load_time + alig_time);
  printf("|\n");
  printf("========================================================================================\n");
  printf("| Total Reads Processed: %-62llu", total_reads);
  printf("|\n");
  printf("+-------------------------------------------+------------------------------------------+\n");
  printf("| Reads Mapped: %-18llu  %6.2f", num_mapped_reads, num_mapped_reads * 100.0 / total_reads);
  printf("% | ");
  printf(" Reads Unmapped: %-14llu  %6.2f", total_reads - num_mapped_reads, (total_reads - num_mapped_reads) * 100.0 / total_reads);
  printf("%  |\n");
  if (rna_mode) {
    /*  printf("========================================================================================\n");
    printf("| Total Splice Junctions: %-61llu", total_sp);
    printf("|\n");
    printf("+--------------------------------------------------------------------------------------+\n");
    printf("| Differents Splice Junctions: %-56llu", uniq_sp);
    printf("|\n");
    printf("+-------------------------------------------+------------------------------------------+\n");
    printf("| Cannonical: %-20llu  %6.2f", uniq_sp, 100.0);
    printf("%  | Semi-Cannonical: %-14llu  %6.2f", 0, 0);
    printf("%  |\n");*/
    printf("+-------------------------------------------+------------------------------------------+\n");
  } else {
    printf("+-------------------------------------------+------------------------------------------+\n");
  }

}

//---------------------------------------------------------------------------------------

void basic_statistics_init(size_t total_reads, size_t num_mapped_reads, size_t total_mappings, basic_statistics_t *statistics){
  statistics->total_reads =  total_reads;
  statistics->num_mapped_reads = num_mapped_reads;
  statistics->total_mappings = total_mappings;
}

//---------------------------------------------------------------------------------------

void basic_statistics_sp_init(size_t total_sp, size_t uniq_sp, basic_statistics_t *statistics){
  statistics->total_sp = total_sp;
  statistics->uniq_sp = uniq_sp;
}

//------------------------------------------------------------------------------------------


basic_statistics_t *basic_statistics_new() {
  basic_statistics_t *basic = (basic_statistics_t *)malloc(sizeof(basic_statistics_t));
  pthread_mutex_init(&(basic->mutex), NULL);
  basic->total_reads = 0;
  basic->num_mapped_reads = 0;
  basic->total_mappings = 0;
  basic->total_sp = 0;
  basic->uniq_sp = 0;
  return basic;
}


//-------------------------------------------------------------------------------------------

void basic_statistics_add(size_t total_reads, size_t num_mapped_reads, size_t total_mappings, basic_statistics_t *basic) {
  pthread_mutex_lock(&basic->mutex);
  basic->total_reads += total_reads;
  basic->num_mapped_reads += num_mapped_reads;
  basic->total_mappings += total_mappings;
  pthread_mutex_unlock(&basic->mutex);
}
