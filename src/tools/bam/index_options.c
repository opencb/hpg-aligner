#include "index_options.h"

//------------------------------------------------------------------------

void usage_index_options(index_options_t *opts);

void **new_argtable_index_options();
index_options_t *read_cli_index_options(void **argtable, index_options_t *opts);

extern void free_argtable(int num_options, void **argtable);
extern void usage_argtable(char *exec_name, char *command_name, void **argtable);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

index_options_t *index_options_new(char *exec_name, char *command_name) {
  index_options_t *opts = (index_options_t*) calloc (1, sizeof(index_options_t));
  
  opts->help = 0;
  opts->in_filename = NULL;

  opts->exec_name = strdup(exec_name);
  opts->command_name = strdup(command_name);

  return opts;
}

//------------------------------------------------------------------------

index_options_t *index_options_parse(char *exec_name, char *command_name,
				     int argc, char **argv) {
  void **argtable = new_argtable_index_options();
  
  index_options_t *opts = index_options_new(exec_name, command_name);
  if (argc < 2) {
    usage_argtable(exec_name, command_name, argtable);
  } else {  
    int num_errors = arg_parse(argc, argv, argtable);
    
    // show help
    if (((struct arg_int*) argtable[0])->count) {
      usage_argtable(exec_name, command_name, argtable);
    }
        
    if (num_errors > 0) {
      arg_print_errors(stdout, argtable[NUM_INDEX_OPTIONS], exec_name);
      usage_argtable(exec_name, command_name, argtable);
    } else {
      opts = read_cli_index_options(argtable, opts);
      if (opts->help) {
	usage_argtable(exec_name, command_name, argtable);
      }
    }
  }

  free_argtable(NUM_INDEX_OPTIONS + 1, argtable);

  return opts;
}

//------------------------------------------------------------------------

void index_options_free(index_options_t *opts) {
  if (opts == NULL) { return; }
  
  if (opts->in_filename) { free(opts->in_filename); }
  
  if (opts->exec_name) { free(opts->exec_name); }
  if (opts->command_name) { free(opts->command_name); }
  
  free(opts);
}

//------------------------------------------------------------------------

void index_options_validate(index_options_t *opts) {
  if (! exists(opts->in_filename)) {
    printf("\nError: Input file name not found !\n\n");
    usage_index_options(opts);
  }

  if (! exists(opts->out_dirname)) {
    opts->out_dirname = strdup(".");
  }
}

//------------------------------------------------------------------------

void index_options_display(index_options_t *opts) {
  printf("PARAMETERS CONFIGURATION\n");
  printf("=================================================\n");
  printf("Main options\n");
  printf("\tBAM input filename  : %s\n", opts->in_filename);
  printf("\tOutput dirname      : %s\n", opts->out_dirname);
  printf("\n");
  /*
  printf("Report options\n");
  printf("\tLog level: %d\n",  (int) opts->log_level);
  printf("\tVerbose  : %d\n",  (int) opts->verbose);
  printf("\n");

  printf("Architecture options\n");
  printf("\tNum. threads: %d\n",  (int) opts->num_threads);
  printf("\tBatch size  : %d alignments\n",  (int) opts->batch_size);
  */
  printf("=================================================\n");
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

index_options_t *read_cli_index_options(void **argtable, index_options_t *opts) {	
  if (((struct arg_int*)argtable[0])->count) { opts->help = ((struct arg_int*)argtable[0])->count; }
  if (((struct arg_file*)argtable[1])->count) { opts->in_filename = strdup(*(((struct arg_file*)argtable[1])->filename)); }
  if (((struct arg_file*)argtable[2])->count) { opts->out_dirname = strdup(*(((struct arg_file*)argtable[2])->filename)); }
  
  return opts;
}

//--------------------------------------------------------------------

void usage_index_options(index_options_t *opts) {
  void **argtable = new_argtable_index_options();
  usage_argtable(opts->exec_name, opts->command_name, argtable);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

void** new_argtable_index_options() {
  void **argtable = (void**)malloc((NUM_INDEX_OPTIONS + 1) * sizeof(void*));
  
  // NOTICE that order cannot be changed as is accessed by index in other functions
  argtable[0] = arg_lit0("h", "help", "Help option");
  argtable[1] = arg_file0("b", "bam-file", NULL, "Input file name (BAM format)");
  argtable[2] = arg_file0("o", "outdir", NULL, "Output file name (BAM format)");
  
  argtable[NUM_INDEX_OPTIONS] = arg_end(20);
  
  return argtable;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------



/*
int read_config_file(const char *filename, options_t *opts) {
	if (filename == NULL || opts == NULL) {
		return -1;
	}

	config_t *config = (config_t*) calloc (1, sizeof(config_t));
	int ret_code = config_read_file(config, filename);
	if (ret_code == CONFIG_FALSE) {
		LOG_ERROR_F("Configuration file error: %s\n", config_error_text(config));
		return -1;
	}

	const char *tmp_string;
	long tmp_int;

	config_destroy(config);
	free(config);

	return ret_code;
}
*/
