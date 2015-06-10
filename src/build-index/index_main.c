#include <stdio.h>
#include "index_builder.h"


int main(int argc, char* argv[]) {

  if (argc <= 1) {
    fprintf(stderr, "Missing command.\nValid commands are:\n\tbuild-simple-index: to create only compress genome file\n\tbuild-sa-index: to create the genome SA index (suffix array).\n\tbuild-bwt-index: to create the genome BWT index.\nUse -h or --help to display hpg-aligner options.\n");
    exit(-1);
  }

  
  char *command = argv[1];  
  
  // We need to consume command: {dna | rna | bs | build-index}
  argc -= 1;
  argv += 1;
  
  if(strcmp(command, "build-sa-index") != 0 &&
     strcmp(command, "build-simple-index") != 0 &&
     strcmp(command, "build-bwt-index") != 0) {
    fprintf(stderr, "Missing command.\nValid commands are:\n\tbuild-simple-index: to create only compress genome file\n\tbuild-sa-index: to create the genome SA index (suffix array).\n\tbuild-bwt-index: to create the genome BWT index.\nUse -h or --help to display hpg-aligner options.\n");
    exit(-1);
  }

  run_index_builder(argc, argv, command);
  
}
