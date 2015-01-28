#include "index_builder.h"


//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------

index_options_t *index_options_new() {
  index_options_t *options = (index_options_t *)malloc(sizeof(index_options_t));

  options->mode = 0;

  options->version = 0;
  options->help = 0;
  options->index_ratio = 0;

  options->ref_genome = NULL;
  options->alt_filename = NULL;
  options->decoy_genome = NULL;
  options->index_filename = NULL;

  return options;
}

//------------------------------------------------------------------------------------

void index_options_free(index_options_t *options) {
  if (options) {
    if (options->ref_genome) { free(options->ref_genome); }
    if (options->alt_filename) { free(options->alt_filename); }
    if (options->decoy_genome) { free(options->decoy_genome); }
    if (options->index_filename) { free(options->index_filename); }
    free(options);
  }
}

//------------------------------------------------------------------------------------

void** argtable_index_options_new(int mode) {
  int num_options = NUM_INDEX_OPTIONS;
  if (mode == BWT_INDEX) { 
    num_options += NUM_INDEX_BWT_OPTIONS; 
  }
    
  // NUM_OPTIONS +1 to allocate end structure
  void **argtable = (void**)malloc((num_options + 1) * sizeof(void*));	

  int count = 0;
  argtable[count++] = arg_file1("i", "index", NULL, "Index directory name");
  argtable[count++] = arg_file1("g", "ref-genome", NULL, "Reference genome (FASTA format)");

  if (mode == BWT_INDEX) {
    argtable[count++] = arg_int0("r", "index-ratio", NULL, "BWT index compression ratio. Default: 8");
  } else {
    argtable[count++] = arg_file0("a", "alternative-map", NULL, "Alternative mapping filename. This two-columns file contains the alternative sequence names with their corresponding chromosome names (only for SA index)");
    argtable[count++] = arg_file0("d", "decoy-genome", NULL, "Decoy genome in FASTA format (only for SA index)");
  }

  argtable[count++] = arg_lit0("v", "version", "Display version");
  argtable[count++] = arg_lit0("h", "help", "Help option");

  argtable[num_options] = arg_end(count);
  
  return argtable;
}

void argtable_index_options_free(void **argtable, int num_options) {
  if(argtable != NULL) {
    arg_freetable(argtable, num_options);        // struct end must also be freed
    free(argtable);
  }
}

//------------------------------------------------------------------------------------

index_options_t *read_CLI_index_options(void **argtable, index_options_t *options, int mode) {        

  int count = -1;
  if (((struct arg_file*)argtable[++count])->count) { options->index_filename = strdup(*(((struct arg_file*)argtable[count])->filename)); }
  if (((struct arg_file*)argtable[++count])->count) { options->ref_genome = strdup(*(((struct arg_file*)argtable[count])->filename)); }
  if (mode == BWT_INDEX) {
    if (((struct arg_int*)argtable[++count])->count) { 
      options->index_ratio = *(((struct arg_int*)argtable[count])->ival); 
    } else {
      options->index_ratio = BWT_RATIO_DEFAULT;
    }
  } else {
    if (((struct arg_file*)argtable[++count])->count) { options->alt_filename = strdup(*(((struct arg_file*)argtable[count])->filename)); }
    if (((struct arg_file*)argtable[++count])->count) { options->decoy_genome = strdup(*(((struct arg_file*)argtable[count])->filename)); }
  }

  if (((struct arg_int*)argtable[++count])->count) { options->version = ((struct arg_int*)argtable[count])->count; }
  if (((struct arg_int*)argtable[++count])->count) { options->help = ((struct arg_int*)argtable[count])->count; }

  return options;
}


//------------------------------------------------------------------------------------

void usage_index(void **argtable, int mode) {
  printf("\nUsage:\n\t%s %s <options>\n", HPG_ALIGNER_BIN, 
	 (mode == BWT_INDEX ? "build-bwt-index" : "build-sa-index"));

  //arg_print_syntaxv(stdout, argtable, "\n");
  printf("\nOptions:\n");
  arg_print_glossary(stdout, argtable, "\t%-50s\t%s\n");

  exit(0);
}

//------------------------------------------------------------------------------------

index_options_t *parse_index_options(int argc, char **argv) {
  int mode = SA_INDEX, num_options = NUM_INDEX_OPTIONS;

  if (strcmp(argv[0], "build-bwt-index") == 0) {
    mode = BWT_INDEX;
    num_options += NUM_INDEX_BWT_OPTIONS;
  } else if (strcmp(argv[0], "build-sa-index") == 0) {
    mode = SA_INDEX;
  } 

  void **argtable = argtable_index_options_new(mode);
  index_options_t *options = index_options_new();
  if (argc < 2) {
    usage_index(argtable, mode);
    exit(-1);
  } else {
    int num_errors = arg_parse(argc, argv, argtable);
    if (num_errors > 0) {
      fprintf(stdout, "\nErrors:\n");
      // struct end is always allocated in the last position               
      arg_print_errors(stdout, argtable[num_options], "\t");
      usage_index(argtable, mode);
      exit(-1);
    } else {
      options = read_CLI_index_options(argtable, options, mode);
      if(options->help) {
        usage_index(argtable, mode);
        argtable_index_options_free(argtable, num_options);
        index_options_free(options);
        exit(0);
      }
      if (options->version) {
	display_version();
        argtable_index_options_free(argtable, num_options);
        index_options_free(options);
        exit(0);
      }
    }
  }
  
  argtable_index_options_free(argtable, num_options);

  options->mode = mode;
  
  return options;
}

//------------------------------------------------------------------------------------

void validate_index_options(index_options_t *options, int mode) {
  if (!exists(options->ref_genome)) {
    fprintf(stdout, "\nError: Your reference genome (%s) does not exist.\n", 
	    options->ref_genome);
    exit(-1);
  }
  
  if (!exists(options->index_filename)) {
    fprintf(stdout, "\nError: Your index directory (%s) does not exist.\n", 
	    options->index_filename);
    exit(-1);
  }

  if (options->alt_filename && !exists(options->alt_filename)) {
    fprintf(stdout, "\nError: Your alternative mapping file (%s) does not exist.\n", 
	    options->alt_filename);
    exit(-1);
  }

  if (options->decoy_genome && !exists(options->decoy_genome)) {
    fprintf(stdout, "\nError: Your decoy genome (%s) does not exist.\n", 
	    options->decoy_genome);
    exit(-1);
  }
  
  if (mode == BWT_INDEX && options->index_ratio <= 0) {
    fprintf(stdout, "\nError: Your compression ratio (%i) is invalid. It must be greater than 0.\n", 
	    options->index_ratio);
    exit(-1);
  }
}

//------------------------------------------------------------------------------------

void run_index_builder(int argc, char **argv, char *mode_str) {
  int mode = BWT_INDEX;

  if (!strcmp(mode_str, "build-bwt-index")) {
    mode = BWT_INDEX;
  } else if (!strcmp(mode_str, "build-sa-index")) {
    mode = SA_INDEX;
  }

  index_options_t *options = parse_index_options(argc, argv);

  argtable_index_options_new(mode);
  validate_index_options(options, mode);
  
  options->cmdline = create_cmdline(argc, argv);

  index_options_display(options);

  if (mode == SA_INDEX) {
    const uint prefix_value = 18;
    char binary_filename[strlen(options->index_filename) + 128];
    char *final_genome;
    sprintf(binary_filename, "%s/dna_compression.bin", options->index_filename);
    printf("Generating SA Index...\n");
    if (options->decoy_genome) {
      final_genome = calloc(strlen(options->index_filename) + 128, sizeof(char));
      sprintf(final_genome, "%s/tmp.concat.genomes.fa", options->index_filename);
      merge_genomes(options->ref_genome, options->decoy_genome, final_genome);
    } else {
      final_genome = strdup(options->ref_genome);
    }
    if (options->alt_filename) {
      sa_index3_build_k18_alt(final_genome, options->alt_filename, 
			      prefix_value, options->index_filename);
    } else {
      sa_index3_build_k18(final_genome, prefix_value, options->index_filename);
    }
    if (options->decoy_genome) {
      sa_index3_set_decoy(options->decoy_genome, options->index_filename);
      remove(final_genome);
    }
    free(final_genome);
    generate_codes(binary_filename, options->ref_genome);
    printf("SA Index generated!\n");

    LOG_DEBUG("Compressing reference genome...\n");
    generate_codes(binary_filename, options->ref_genome);
    LOG_DEBUG("...done !\n");

  } else {
    char binary_filename[strlen(options->index_filename) + 128];
    sprintf(binary_filename, "%s/dna_compression.bin", options->index_filename);
    
    LOG_DEBUG("Compressing reference genome...\n");
    generate_codes(binary_filename, options->ref_genome);
    LOG_DEBUG("...done !\n");
    
    LOG_DEBUG("Building BWT index...\n");
    bwt_generate_index_files(options->ref_genome, options->index_filename, options->index_ratio, false, "ACGT");
    LOG_DEBUG("...done !\n");    
  }
}

//--------------------------------------------------------------------

void index_options_display(index_options_t *options) {
     printf("\n");
     printf("+===============================================================+\n");
     printf("|             HPG ALIGNER PARAMETERS CONFIGURATION              |\n");
     printf("+===============================================================+\n");
     printf("HPG Aligner version: %s\n", HPG_ALIGNER_VERSION);
     printf("Command line: %s\n", options->cmdline);
     printf("\n");
     printf("General parameters\n");
     printf("\tReference genome: %s\n", options->ref_genome);
     if (options->mode == SA_INDEX) {
       printf("\tAlternative sequence names: %s\n", (options->alt_filename ? options->alt_filename : "None"));
       printf("\tDecoy genome: %s\n", (options->decoy_genome ? options->decoy_genome : "None"));
     }
     printf("\t%s index directory name: %s\n", (options->mode == SA_INDEX ? "SA" : "BWT"),
	    options->index_filename);
     if (options->mode == BWT_INDEX) {
       printf("\tCompression ration: %i\n", options->index_ratio);
     }
     printf("+===============================================================+\n");     
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
