#include "dna/dna_aligner.h"
#include "rna/rna_aligner.h"
#include "bs/bs_aligner.h"
#include "build-index/index_builder.h"


double emboss_matrix_t = 0.0f, emboss_tracking_t = 0.0f;
double sse_matrix_t = 0.0f, sse_tracking_t = 0.0f;
double sse1_matrix_t = 0.0f, sse1_tracking_t = 0.0f;
double avx_matrix_t = 0.0f, avx_tracking_t = 0.0f;
double avx1_matrix_t = 0.0f, avx1_tracking_t = 0.0f;

//--------------------------------------------------------------------
// constants
//--------------------------------------------------------------------

#define OPTIONS 			30
#define MIN_ARGC  			5
#define NUM_SECTIONS_STATISTICS 	5
#define NUM_SECTIONS_STATISTICS_SB	21

//--------------------------------------------------------------------
// global variables for timing and capacity meausures
//--------------------------------------------------------------------

basic_statistics_t *basic_st;

pthread_mutex_t mutex_sp;

FILE *fd_log;
size_t junction_id;

//--------------------------------------------------------------------
// main parameters support
//--------------------------------------------------------------------
int main(int argc, char* argv[]) {
  pthread_mutex_init(&mutex_sp, NULL);
  
  //memset(tot_cals, 0, sizeof(int)*50);
  //const char HEADER_FILE[1024] = "Human_NCBI37.hbam\0";

  basic_st = basic_statistics_new();

  // init logs, after parsing the command-line
  // logs will be re-set according to the command-line
  log_level = LOG_FATAL_LEVEL;
  log_verbose = 1;
  log_file = NULL;

  if (argc <= 1) {
    LOG_FATAL("Missing command.\nValid commands are:\n\tdna: to map DNA sequences\n\trna: to map RNA sequences\n\tbuild-sa-index: to create the genome SA index (suffix array).\n\tbuild-bwt-index: to create the genome BWT index.\nUse -h or --help to display hpg-aligner options.\n");
  }

  if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
    usage_cli(DNA_MODE);
  }

  char *command = argv[1];  

  // We need to consume command: {dna | rna | bs | build-index}
  argc -= 1;
  argv += 1;
  
  if(strcmp(command, "dna") != 0 && 
     strcmp(command, "rna") != 0 &&
     strcmp(command, "bs" ) != 0 && 
     strcmp(command, "build-sa-index") != 0 &&
     strcmp(command, "build-bwt-index") != 0) {
    LOG_FATAL("Command unknown.\nValid commands are:\n\tdna: to map DNA sequences\n\trna: to map RNA sequences\n\tbs: to map BS sequences\n\tbuild-sa-index: to create the genome sa index.\n\tbuild-bwt-index: to create the genome bwt index.\nUse -h or --help to display hpg-aligner options.\n");

  }

  if (!strcmp(command, "build-bwt-index") || 
      !strcmp(command, "build-sa-index")) {
      run_index_builder(argc, argv, command);
  }

  // parsing options
  options_t *options = parse_options(argc, argv);

  // now, we can set logs according to the command-line
  init_log_custom(options->log_level, 1, "hpg-aligner.log", "w");
  LOG_DEBUG_F("Command Mode: %s\n", command);

  
  if(strcmp(command, "dna") == 0) { 
    // DNA command
    validate_options(options);
    dna_aligner(options);
  } else if (strcmp(command, "rna") == 0)  { 
    // RNA command
    validate_options(options);
    rna_aligner(options);
  }

  //} else if (strcmp(command, "bs") == 0) {
  // BS commnad
  //printf("Run BS mode...");
  //run_bs_aligner(genome1, genome2, genome, bwt_index1, bwt_index2,
  //                   bwt_optarg, cal_optarg, pair_mng, report_optarg, options);
  //}   

  /*
  if (options->fast_mode) {  
    size_t mapped_reads = total_reads - unmapped_reads;
    printf("TOTAL READS    : %lu\n", total_reads);
    printf("    UNMAPPED READS : %lu (%f%%)\n", unmapped_reads, ((float)unmapped_reads * 100)/((float)total_reads));
    printf("    MAPPED READS   : %lu (%f%%)\n", mapped_reads, ((float)mapped_reads * 100)/((float)total_reads));
    printf("-------------------------------\n");
    printf("    CORRECT READS    : %lu (%f%%)\n", correct_reads,  ((float)correct_reads * 100)/((float)total_reads));
    printf("    INCORRECT READS  : %lu (%f%%)\n", mapped_reads - correct_reads,  ((float)(mapped_reads - correct_reads) * 100)/((float)total_reads));
    

    //printf("Seeds with one Error %i (%0.2f)\n", seeds_1err, seeds_1err*100/total_reads);
  }
  */

  options_free(options);
  
  return 0;

}


//--------------------------------------------------------------------
