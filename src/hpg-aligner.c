#include <stdio.h>
#include <stdlib.h>
#include <options.h>
#include <multi/multi_aligner.h>

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

#define RNA_MODE                         1
#define DNA_MODE                         0

//--------------------------------------------------------------------
// main parameters support
//--------------------------------------------------------------------
int main(int argc, char* argv[]) {


  if (argc <= 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
    usage_cli(DNA_MODE);
  }

  // parsing options
  options_t *options = parse_options(argc, argv);

  validate_options(options);
  
  hpg_multialigner_main(options, argc, argv);
  
  options_free(options);
  
  return 0;

}
