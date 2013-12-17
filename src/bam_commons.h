
#ifndef BAM_COMMONS_H
#define BAM_COMMONS_H

#include <pthread.h>
#include <sys/time.h>

#define OPTIMAL_THREADS_FOR_COMPUTE_CAPABILITY_20 512

#define PHRED33  33
#define PHRED64  64

#define LIST_OF_CHROMOSOMES {"1", "2", "3", "4", "5", "6", "7, "8, "9, "10", "11", "12", "13", "14", "15", "16", "17, "18, "19, "20", "21", "22+X", "22+Y", "22+MT"}
#define DEFAULT_NUM_OF_CHROMOSOMES 	25
#define MAX_NUM_OF_CHROMOSOMES 		100
#define ALL_CHROMOSOMES  		9999

#define MAX_ALIGNMENT_LENGTH 1000

/* **************************************************************
 *    		Global structures and variables   		*
 * *************************************************************/

extern int number_of_batchs;
extern int bam_reader_alive;

extern double sort_time;
extern struct timeval t1_sort, t2_sort;

extern double convert_time;
extern struct timeval t1_convert, t2_convert;

extern double filter_time;
extern struct timeval t1_filter, t2_filter;

extern double reader_server_time;
extern struct timeval t1_reader_server, t2_reader_server;

extern double qc_calc_server_time;
extern struct timeval t1_qc_calc_server, t2_qc_calc_server;

extern double cpus_server_time;
extern struct timeval t1_cpus_server, t2_cpus_server;

extern double results_server_time;
extern struct timeval t1_results_server, t2_results_server;

extern double gpus_standby_time, cpus_standby_time, results_standby_time;
extern struct timeval t1_active_reader, t1_active_gpus, t1_active_cpus, t1_active_results;

extern int num_alignments;
extern pthread_mutex_t read_count_lock;

extern unsigned int nts_with_coverage;
extern unsigned long mean_coverage;

extern int num_of_chromosomes;

#endif /* BAM_COMMONS_H */
