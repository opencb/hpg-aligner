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
     printf("./hpg-aligner {build-bwt-index | build-sa-index} -[i|--bwt-index=<directory>] [-g|--ref-genome=<file>] [--index-ratio=<int>] [-bs|--bisulphite-index]\n");
     printf("-i, --index=<directory>\tBWT directory name\n");
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

void run_index_builder_bs(char *genome_filename, char *bwt_dirname, 
			  int bwt_ratio, bool duplicate_strand, char *bases) {
      /** **************************************************************************	*
       * 										*
       * Generates the genome transform from the input and builds the index		*
       * 										*
       * The genome transformed are stored in the directory give by the user,		*
       * and the index are stored in subfolders				       		*
       * 										*
       * ***************************************************************************	*/
  /*
      run_index_builder(options->genome_filename, options->bwt_dirname, options->index_ratio, false, "ACGT");

      // generate binary code for original genome
      char binary_filename[strlen(options->bwt_dirname) + 128];
      sprintf(binary_filename, "%s/dna_compression.bin", options->bwt_dirname);
      generate_codes(binary_filename, options->genome_filename);

      char bs_dir1[256];
      sprintf(bs_dir1, "%s/AGT_index", options->bwt_dirname);
      //if (is_directory(bs_dir1) == 0) {
      create_directory(bs_dir1);
      //}

      LOG_DEBUG("Generation of AGT index\n");
      char genome_1[256];
      sprintf(genome_1, "%s/AGT_genome.fa", options->bwt_dirname);
      char gen1[256];
      sprintf(gen1, "sed 's/C/T/g' %s > %s",options->genome_filename, genome_1);
      system(gen1);

      run_index_builder(genome_1, bs_dir1, options->index_ratio, false, "AGT");
      LOG_DEBUG("AGT index Done !!\n");

      LOG_DEBUG("Generation of ACT index\n");
      char bs_dir2[256];
      sprintf(bs_dir2, "%s/ACT_index", options->bwt_dirname);
      //if (is_directory(bs_dir2) == 0) {
      create_directory(bs_dir2);
      //}

      char genome_2[256];
      sprintf(genome_2, "%s/ACT_genome.fa", options->bwt_dirname);
      char gen2[256];
      sprintf(gen2, "sed 's/G/A/g' %s > %s",options->genome_filename, genome_2);
      system(gen2);

      run_index_builder(genome_2, bs_dir2, options->index_ratio, false, "ACT");
      LOG_DEBUG("ACT index Done !!\n");
  */
}
