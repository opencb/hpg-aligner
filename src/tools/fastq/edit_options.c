#include "edit_options.h"

//------------------------------------------------------------------------

void usage_edit_options(edit_options_t *opts);

void **new_argtable_edit_options();
edit_options_t *read_cli_edit_options(void **argtable, edit_options_t *opts);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

edit_options_t *edit_options_new(char *exec_name, char *command_name) {
  edit_options_t *opts = (edit_options_t*) calloc (1, sizeof(edit_options_t));
  
  opts->kmers_on = 0;
  opts->filter_on = 0;
  opts->log_level = 0;
  opts->verbose = 0;
  opts->help = 0;
  opts->num_threads = 2;
  opts->batch_size = 10000;

  opts->quality_encoding_value = 0;

  // filter or trim edit_options
  opts->max_read_length = -1;
  opts->min_read_length = -1;
  opts->max_N = -1;
  opts->max_read_quality = -1;
  opts->min_read_quality = -1;
  opts->left_length = -1;
  opts->max_left_quality = -1;
  opts->min_left_quality = -1;
  opts->right_length = -1;
  opts->max_right_quality = -1;
  opts->min_right_quality = -1;
  opts->max_out_of_quality = -1;
 
  opts->quality_encoding_name = NULL;
  opts->read_length_range = NULL;
  opts->read_quality_range = NULL;
  opts->left_quality_range = NULL;
  opts->right_quality_range = NULL;

  // filename edit_options
  opts->in_filename = NULL;
  opts->out_dirname = NULL;

  opts->exec_name = strdup(exec_name);
  opts->command_name = strdup(command_name);

  return opts;
}

//------------------------------------------------------------------------

edit_options_t *edit_options_parse(char *exec_name, char *command_name,
			 int argc, char **argv) {
  void **argtable = new_argtable_edit_options();
  
  edit_options_t *opts = edit_options_new(exec_name, command_name);
  if (argc < 2) {
    usage_argtable(exec_name, command_name, argtable);
  } else {  
    int num_errors = arg_parse(argc, argv, argtable);
    
    // show help
    if (((struct arg_int*) argtable[0])->count) {
      usage_argtable(exec_name, command_name, argtable);
    }
        
    if (num_errors > 0) {
      arg_print_errors(stdout, argtable[NUM_EDIT_OPTIONS], exec_name);
      usage_argtable(exec_name, command_name, argtable);
    } else {
      opts = read_cli_edit_options(argtable, opts);
      if (opts->help) {
	usage_argtable(exec_name, command_name, argtable);
      }
    }
  }
  
  free_argtable(NUM_EDIT_OPTIONS + 1, argtable);

  return opts;
}

//------------------------------------------------------------------------

void edit_options_free(edit_options_t *opts) {
  if (opts == NULL) { return; }
  
  if (opts->quality_encoding_name) { free(opts->quality_encoding_name); }
  if (opts->read_length_range) { free(opts->read_length_range); }
  if (opts->read_quality_range) { free(opts->read_quality_range); }
  if (opts->left_quality_range) { free(opts->left_quality_range); }
  if (opts->right_quality_range) { free(opts->right_quality_range); }

  if (opts->in_filename) { free(opts->in_filename); }
  if (opts->out_dirname) { free(opts->out_dirname); }
  
  if (opts->exec_name) { free(opts->exec_name); }
  if (opts->command_name) { free(opts->command_name); }
  
  free(opts);
}

//------------------------------------------------------------------------

void edit_options_validate(edit_options_t *opts) {
  // input filename
  if (! exists(opts->in_filename)) {
    printf("\nError: Input file name not found !\n");
    usage_edit_options(opts);
  }

  // output dirname
  if (! exists(opts->out_dirname)) {
    opts->out_dirname = strdup(".");
  }

  // quality encoding
  if (opts->quality_encoding_name) {
    if (strcmp(opts->quality_encoding_name, QUALITY_PHRED33_NAME) == 0) {
      opts->quality_encoding_value = QUALITY_PHRED33_VALUE;
    } else if (strcmp(opts->quality_encoding_name, QUALITY_PHRED64_NAME) == 0) {
      opts->quality_encoding_value = QUALITY_PHRED64_VALUE;
    } else {
      printf("\nError: Invalid quality encoding value (%s). Valid values: %s, %s\n",
	     opts->quality_encoding_name, QUALITY_PHRED33_NAME, QUALITY_PHRED64_NAME);
      usage_edit_options(opts);
    }
  } else {
    opts->quality_encoding_name = strdup(QUALITY_PHRED33_NAME);
    opts->quality_encoding_value = QUALITY_PHRED33_VALUE;
  }

  // read length range
  if (!parse_range(&opts->min_read_length, &opts->max_read_length, 
		   opts->read_length_range, "read length range")) {
    usage_edit_options(opts);
  }

  // read quality range
  if (!parse_range(&opts->min_read_quality, &opts->max_read_quality, 
		   opts->read_quality_range, "read quality range")) {
    usage_edit_options(opts);
  }

  // left quality range
  if (!parse_range(&opts->min_left_quality, &opts->max_left_quality, 
		   opts->left_quality_range, "left quality range")) {
    usage_edit_options(opts);
  }

  // right quality range
  if (!parse_range(&opts->min_right_quality, &opts->max_right_quality, 
		   opts->right_quality_range, "right quality range")) {
    usage_edit_options(opts);
  }
}

//------------------------------------------------------------------------

void edit_options_display(edit_options_t *opts) {
  printf("PARAMETERS CONFIGURATION\n");
  printf("=================================================\n");
  printf("Command name : %s\n", opts->command_name);
  printf("\n");
  printf("Main edit_options\n");
  printf("\tFastQ input filename : %s\n", opts->in_filename);
  printf("\tOutput dirname       : %s\n", opts->out_dirname);

  int edit_count = 0;
  printf("\nEdit options\n");
  if (opts->left_length != NO_VALUE && opts->left_quality_range) {
    edit_count++;
    printf("\tTrim left length         : %i nucleotides\n", opts->left_length);
    printf("\tTrim left quality range  : %s\n", opts->left_quality_range);
  }
  if (opts->right_length != NO_VALUE && opts->right_quality_range) {
    edit_count++;
    printf("\tTrim right length        : %i nucleotides\n", opts->right_length);
    printf("\tTrim right quality range : %s\n", opts->right_quality_range);
  }
  if (edit_count == 0) {
    printf("\tNone.\n\n");
  }

  int filter_count = 0;
  printf("\nFilter options\n");
  if (opts->read_length_range) {
    filter_count++;
    printf("\tRead length range   : %s\n", opts->read_length_range);
  }
  if (opts->read_quality_range) {
    filter_count++;
    printf("\tRead quality range  : %s\n", opts->read_quality_range);
  }

  if (opts->max_N != NO_VALUE) {
    filter_count++;
    printf("\tMax. number of Ns   : %i\n", opts->max_N);
  }
  if (opts->max_out_of_quality != NO_VALUE && opts->read_quality_range) {
    filter_count++;
    printf("\tMax. out of quality : %i nucletotides\n", opts->max_out_of_quality);
  }
  
  if (filter_count == 0) {
    printf("\tNone.\n\n");
    opts->filter_on = 0;
  } else {
    opts->filter_on = 1;
  }
  /*  
  printf("Report options\n");
  printf("\tLog level: %d\n",  (int) opts->log_level);
  printf("\tVerbose  : %d\n",  (int) opts->verbose);
  printf("\n");
  */
  printf("\nArchitecture options\n");
  printf("\tNum. threads: %d\n",  (int) opts->num_threads);
  printf("\tBatch size  : %d alignments\n",  (int) opts->batch_size);
  printf("=================================================\n");

  if (edit_count == 0) {
    printf("\n\nNothing to edit, no edit options specified !\n\n");
    exit(-1);
  }
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

edit_options_t *read_cli_edit_options(void **argtable, edit_options_t *opts) {	
  if (((struct arg_int*)argtable[0])->count) { opts->help = ((struct arg_int*)argtable[0])->count; }
  if (((struct arg_file*)argtable[1])->count) { opts->in_filename = strdup(*(((struct arg_file*)argtable[1])->filename)); }
  if (((struct arg_file*)argtable[2])->count) { opts->out_dirname = strdup(*(((struct arg_file*)argtable[2])->filename)); }
  //  if (((struct arg_int*)argtable[3])->count) { opts->log_level = *(((struct arg_int*)argtable[3])->ival); }
  //  if (((struct arg_int*)argtable[4])->count) { opts->verbose = *(((struct arg_int*)argtable[4])->ival); }
  if (((struct arg_int*)argtable[3])->count) { opts->num_threads = *(((struct arg_int*)argtable[3])->ival); }
  if (((struct arg_int*)argtable[4])->count) { opts->batch_size = *(((struct arg_int*)argtable[4])->ival); }
  if (((struct arg_str*)argtable[5])->count) { opts->read_length_range = strdup(*(((struct arg_str*)argtable[5])->sval)); }
  if (((struct arg_str*)argtable[6])->count) { opts->read_quality_range = strdup(*(((struct arg_str*)argtable[6])->sval)); }
  if (((struct arg_int*)argtable[7])->count) { opts->left_length = *(((struct arg_int*)argtable[7])->ival); }
  if (((struct arg_str*)argtable[8])->count) { opts->left_quality_range = strdup(*(((struct arg_str*)argtable[8])->sval)); }
  if (((struct arg_int*)argtable[9])->count) { opts->right_length = *(((struct arg_int*)argtable[9])->ival); }
  if (((struct arg_str*)argtable[10])->count) { opts->right_quality_range = strdup(*(((struct arg_str*)argtable[10])->sval)); }
  if (((struct arg_int*)argtable[11])->count) { opts->max_N = *(((struct arg_int*)argtable[11])->ival); }
  if (((struct arg_int*)argtable[12])->count) { opts->max_out_of_quality = *(((struct arg_int*)argtable[12])->ival); }

  return opts;
}

//--------------------------------------------------------------------

void usage_edit_options(edit_options_t *opts) {
  void **argtable = new_argtable_edit_options();
  usage_argtable(opts->exec_name, opts->command_name, argtable);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

void** new_argtable_edit_options() {
  void **argtable = (void**)malloc((NUM_EDIT_OPTIONS + 1) * sizeof(void*));
  
  // NOTICE that order cannot be changed as is accessed by index in other functions
  argtable[0] = arg_lit0("h", "help", "Help option");
  argtable[1] = arg_file0("f", "fastq-file", NULL, "Input file name (FastQ format)");
  argtable[2] = arg_file0("o", "outdir", NULL, "Output directory name");
  //  argtable[3] = arg_int0(NULL, "log-level", NULL, "Log debug level");
  //  argtable[4] = arg_int0("v", "verbose", NULL, "Verbose");
  argtable[3] = arg_int0(NULL, "num-threads", NULL, "Number of threads");
  argtable[4] = arg_int0(NULL, "batch-size", NULL, "Batch size (in number of alignments)");
  argtable[5] = arg_str0(NULL, "read-length-range", NULL, "Read length range, eg. 80,110");
  argtable[6] = arg_str0(NULL, "read-quality-range", NULL, "Read quality range, eg. 20,40");
  argtable[7] = arg_int0(NULL, "left-length", NULL, "Number of leftmost nucleotides to take into account to trim");
  argtable[8] = arg_str0(NULL, "left-quality-range", NULL, "Quality range for the leftmost nucleotides, eg. 15,45");
  argtable[9] = arg_int0(NULL, "right-length", NULL, "Number of rightmost nucleotides to take into account to trim");
  argtable[10] = arg_str0(NULL, "right-quality-range", NULL, "Quality range for the rightmost nucleotides, eg. 10,60");
  argtable[11] = arg_int0(NULL, "max-N", NULL, "Maximum number of Ns in the sequences");
  argtable[12] = arg_int0(NULL, "max-out-of-quality", NULL, "Maximum number of nucleotides out of the read quality range");
  
  argtable[NUM_EDIT_OPTIONS] = arg_end(20);
  
  return argtable;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
