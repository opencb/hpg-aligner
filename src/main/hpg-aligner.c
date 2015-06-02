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

extern int redirect_stdout;
extern int gziped_fileds;

extern char convert_ASCII[128];

extern st_bwt_t st_bwt;

//--------------------------------------------------------------------

void display_main_help() {
  printf("%s <command>\n", HPG_ALIGNER_BIN);
  printf("\n");
  printf("Commands:\n");
  printf("\tdna: to map DNA sequences\n");
  printf("\trna: to map RNA sequences\n");
  printf("\tbuild-sa-index: to create the genome SA index (suffix array).\n");
  printf("\tbuild-bwt-index: to create the genome BWT index (only available for RNA mapping).\n");
  printf("\n");
  printf("Use -h or --help to display hpg-aligner options.\n");
  printf("Use -v or --version to display hpg-aligner version.\n");
}

//--------------------------------------------------------------------
// main parameters support
//--------------------------------------------------------------------

int main(int argc, char* argv[]) {
  redirect_stdout = 0;
  gziped_fileds = 0;

  if (!isatty(fileno(stdout))) {
    redirect_stdout = 1;
  }

  st_bwt.multi_alig = 0;
  st_bwt.single_alig = 0;
  st_bwt.total_reads = 0;
  st_bwt.map_bwt = 0;

  st_bwt.map_w1 = 0;
  st_bwt.map_w2 = 0;
  st_bwt.map_w3 = 0;

  st_bwt.tot_sj = 0;
  st_bwt.dif_sj = 0;
  st_bwt.cannonical_sj = 0;
  st_bwt.semi_cannonical_sj = 0;

  pthread_mutex_init(&mutex_sp, NULL);

  basic_st = basic_statistics_new();

  // init logs, after parsing the command-line
  // logs will be re-set according to the command-line
  log_level = LOG_FATAL_LEVEL;
  log_verbose = 1;
  log_file = NULL;

  if (argc <= 1) {
    printf("Error: missing command.\n");
    display_main_help();
    exit(-1);
  }

  if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
    display_main_help();
    exit(0);
  }

  if (strcmp(argv[1], "-v") == 0 || strcmp(argv[1], "--version") == 0) {
    display_version();
    exit(0);
  }

  char *command = argv[1];  

  argc -= 1;
  argv += 1;
  
  if(strcmp(command, "dna") != 0 && 
     strcmp(command, "rna") != 0 &&
     strcmp(command, "build-sa-index") != 0 &&
     strcmp(command, "build-bwt-index") != 0) {
    printf("Error: unknown command (%s)\n", command);
    display_main_help();
    exit(-1);
  }


  if (!strcmp(command, "build-bwt-index") || 
      !strcmp(command, "build-sa-index")) {
    run_index_builder(argc, argv, command);
    exit(0);
  }

  // parsing options
  options_t *options = parse_options(argc, argv);

  if (options->version) {
    display_version();
    options_free(options);
    exit(0);
  }

  options->cmdline = create_cmdline(argc, argv);

  if (options->adapter) {
    options->adapter_revcomp = strdup(options->adapter);
    seq_reverse_complementary(options->adapter_revcomp, strlen(options->adapter_revcomp));
  }

  // now, we can set logs according to the command-line
  init_log_custom(options->log_level, 1, "hpg-aligner.log", "w");
  LOG_DEBUG_F("Command Mode: %s\n", command);

  //convert ASCII fill 
  convert_ASCII['a'] = 'T';
  convert_ASCII['A'] = 'T';

  convert_ASCII['c'] = 'G';
  convert_ASCII['C'] = 'G';

  convert_ASCII['g'] = 'C';
  convert_ASCII['G'] = 'C';

  convert_ASCII['t'] = 'a';
  convert_ASCII['T'] = 'A';

  convert_ASCII['n'] = 'N';
  convert_ASCII['N'] = 'N';

  
  if(strcmp(command, "dna") == 0) { 
    // DNA command
    validate_options(options);
    dna_aligner(options);
  } else if (strcmp(command, "rna") == 0)  { 
    // RNA command
    validate_options(options);
    rna_aligner(options);
  }
  options_free(options);

  return 0;

}
