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
//#define REQUIRED 1
//#define NO_REQUIRED 0

//--------------------------------------------------------------------
// global variables for log functions
//--------------------------------------------------------------------

//int log_level = DEFAULT_LOG_LEVEL;
//bool verbose = true;

//--------------------------------------------------------------------
// global variables for timing and capacity meausures
//--------------------------------------------------------------------

//char time_on = 0;
//char statistics_on = 0;
//timing_t* timing = NULL;
basic_statistics_t *basic_st;
//cal_st_t cal_st;
//double kl_time;

//pthread_cond_t cond_sp;
pthread_mutex_t mutex_sp;

//size_t bwt_correct = 0;
//size_t bwt_error = 0;
//size_t seeding_reads = 0;
//pthread_mutex_t bwt_mutex;


//struct timeval time_start_alig, time_end_alig;
//double time_alig;

//size_t tot_reads_in = 0;
//size_t tot_reads_out = 0;

FILE *fd_log;
size_t junction_id;

// timing
//double main_time;
//size_t TOTAL_SW, 
//TOTAL_READS_PROCESS,
//TOTAL_READS_SEEDING,
//TOTAL_READS_SEEDING2,
//TOTAL_READS_SA;

//--------------------------------------------------------------------
// main parameters support
//--------------------------------------------------------------------
int main(int argc, char* argv[]) {

  pthread_mutex_init(&mutex_sp, NULL);
  
  const char HEADER_FILE[1024] = "Human_NCBI37.hbam\0";
  basic_st = basic_statistics_new();

  // init logs, after parsing the command-line
  // logs will be re-set according to the command-line
  log_level = LOG_FATAL_LEVEL;
  log_verbose = 1;
  log_file = NULL;

  if (argc <= 1) {
    LOG_FATAL("Missing command.\nValid commands are:\n\tdna: to map DNA sequences\n\trna: to map RNA sequences\n\tbs: to map BS sequences\n\tbuild-index: to create the genome index.\nUse -h or --help to display hpg-aligner options.");
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
    LOG_FATAL("Command unknown.\nValid commands are:\n\tdna: to map DNA sequences\n\trna: to map RNA sequences\n\tbs: to map BS sequences\n\tbuild-index: to create the genome index.\nUse -h or --help to display hpg-aligner options.");    
  }

  // parsing options
  options_t *options = parse_options(argc, argv);

  // now, we can set logs according to the command-line
  init_log_custom(options->log_level, 1, "hpg-aligner.log", "w");

  LOG_DEBUG_F("Command Mode: %s\n", command);

  if(strcmp(command, "dna") == 0) { 
    // DNA command
    validate_options(options, DNA_MODE);
    dna_aligner(options);
  } else if (strcmp(command, "rna") == 0)  { 
    // RNA command
    validate_options(options, RNA_MODE);
    rna_aligner(options);
  } else if (strcmp(command, "bs") == 0) {
    // BS commnad
    printf("Run BS mode...");
    //run_bs_aligner(genome1, genome2, genome, bwt_index1, bwt_index2,
    //		   bwt_optarg, cal_optarg, pair_mng, report_optarg, options);
  } else if (!strcmp(command, "build-bwt-index")) {
    // BWT-INDEX command
    if (options->bs_index == 0) {
      run_index_builder(options->genome_filename, options->bwt_dirname, 
			options->index_ratio, false, "ACGT");
      LOG_DEBUG("Done !!\n");
      exit(0);
    } else { // bisulphite index generation
      //printf("\nBisulphite index generation\n");
      //run_index_builder(options->genome_filename, options->bwt_dirname, 
      //		options->index_ratio, false, "ACGT");
      //exit(0);
    }
  } else if (!strcmp(command, "build-sa-index")) {
    // SA-INDEX command
    printf("SA Index...");
  }

  return 0;
}

//--------------------------------------------------------------------
