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
int min_intron, max_intron;

int main(int argc, char* argv[]) {
    
  hpg_multialigner_main(argc, argv);  
  
  return 0;
  
}
