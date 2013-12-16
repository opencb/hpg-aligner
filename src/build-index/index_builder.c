#include "index_builder.h"


//------------------------------------------------------------------------------------

void check_index_builder(char *genome_filename, char *bwt_dirname, int bwt_ratio) {

     if (!exists(genome_filename)) {
	  printf("Reference genome does not exist.\n\n");
	  help_index_builder();
     }

     if (!exists(bwt_dirname)) {
	  printf("BWT index directory does not exist.\n\n");
	  help_index_builder();
     }
	  
     if (bwt_ratio <= 0) {
	  printf("Invalid BWT index ratio. It must be greater than 0.\n\n");
	  help_index_builder();
     }
}

//------------------------------------------------------------------------------------

void help_index_builder() {
     printf("./hpg-aligner build-index -[i|--bwt-index=<directory>] [-g|--ref-genome=<file>] [--index-ratio=<int>] [-bs|--bisulphite-index]\n");
     printf("-i, --bwt-index=<directory>\tBWT directory name\n");
     printf("-g, --ref-genome=<file>\t\tReference genome\n");
     printf("-r, --index-ratio=<int>\t\tBWT index ratio of compression\n");
     printf("-bs, --bisulphite-index\t\tIndicates the generation of index for bisulphite case\n");
     exit(0);
}

//------------------------------------------------------------------------------------

void run_index_builder(char *genome_filename, char *bwt_dirname, 
		       int bwt_ratio, bool duplicate_strand, char *bases) {

     check_index_builder(genome_filename, bwt_dirname, bwt_ratio);

     char binary_filename[strlen(bwt_dirname) + 128];
     sprintf(binary_filename, "%s/dna_compression.bin", bwt_dirname);

     LOG_DEBUG("Compressing reference genome...\n");
     generate_codes(binary_filename, genome_filename);
     LOG_DEBUG("...done !\n");

     LOG_DEBUG("Building BWT index...\n");
     bwt_generate_index_files(genome_filename, bwt_dirname, bwt_ratio, duplicate_strand, bases);
     LOG_DEBUG("...done !\n");
}

//------------------------------------------------------------------------------------

