#include "dna/dna_aligner.h"
#include "rna/rna_aligner.h"

#include "build-index/index_builder.h"


//--------------------------------------------------------------------
// constants
//--------------------------------------------------------------------

#define OPTIONS 		       29
#define MIN_ARGC  			5
#define NUM_SECTIONS_STATISTICS 	5
#define NUM_SECTIONS_STATISTICS_SB     21

//--------------------------------------------------------------------
// global variables for timing and capacity meausures
//--------------------------------------------------------------------

//====================================================================
int min_intron, max_intron;
//====================================================================

basic_statistics_t *basic_st;

pthread_mutex_t mutex_sp;

FILE *fd_log;
size_t junction_id;

size_t total_reads = 0;
size_t reads_no_map = 0;

size_t total_sw = 0;

unsigned char mute;

double time_write = 0;
double time_free = 0;
double time_free_batch = 0;

double time_timer0 = 0;
double time_timer1 = 0;
double time_timer2 = 0;
double time_timer3 = 0;

double time_read_fq   = 0;
double time_read_fq_process   = 0;
double time_read_alig = 0;
double time_read_proc = 0;

char convert_ASCII[128];

size_t search_calls = 0;
size_t insert_calls = 0;
double time_search = 0.0;
double time_insert = 0.0;
pthread_mutex_t mutex_calls;

size_t fd_read_bytes = 0;
size_t fd_total_bytes = 0;

size_t total_reads_ph2 = 0;
size_t reads_ph2 = 0;

int redirect_stdout = 0;
int gziped_fileds = 0;

st_bwt_t st_bwt;
int w1_end;
int w2_end;
int w3_end;

size_t total_reads_w2, total_reads_w3;
size_t reads_w2, reads_w3;

//--------------------------------------------------------------------
//--------------------------------------------------------------------
