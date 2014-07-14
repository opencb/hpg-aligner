#include "sort_options.h"

//------------------------------------------------------------------------

void usage_sort_options(sort_options_t *opts);

void **new_argtable_sort_options();
sort_options_t *read_cli_sort_options(void **argtable, sort_options_t *opts);

extern void free_argtable(int num_options, void **argtable);
extern void usage_argtable(char *exec_name, char *command_name, void **argtable);

//------------------------------------------------------------------------
//------------------------------------------------------------------------

sort_options_t *sort_options_new(char *exec_name, char *command_name) {
  sort_options_t *opts = (sort_options_t*) calloc (1, sizeof(sort_options_t));
  
  opts->help = 0;

  opts->max_memory = 500000000;
  opts->criteria = NULL; //strdup("coord");

  opts->in_filename = NULL;
  opts->out_dirname = NULL;

  opts->exec_name = strdup(exec_name);
  opts->command_name = strdup(command_name);

  return opts;
}

//------------------------------------------------------------------------

sort_options_t *sort_options_parse(char *exec_name, char *command_name,
				   int argc, char **argv) {
  void **argtable = new_argtable_sort_options();
  
  sort_options_t *opts = sort_options_new(exec_name, command_name);
  if (argc < 2) {
    usage_argtable(exec_name, command_name, argtable);
  } else {  
    int num_errors = arg_parse(argc, argv, argtable);
    
    // show help
    if (((struct arg_int*) argtable[0])->count) {
      usage_argtable(exec_name, command_name, argtable);
    }
        
    if (num_errors > 0) {
      arg_print_errors(stdout, argtable[NUM_SORT_OPTIONS], exec_name);
      usage_argtable(exec_name, command_name, argtable);
    } else {
      opts = read_cli_sort_options(argtable, opts);
      if (opts->help) {
	usage_argtable(exec_name, command_name, argtable);
      }
    }
  }

  free_argtable(NUM_SORT_OPTIONS + 1, argtable);

  return opts;
}

//------------------------------------------------------------------------

void sort_options_free(sort_options_t *opts) {
  if (opts == NULL) { return; }
  
  if (opts->criteria) { free(opts->criteria); }

  if (opts->in_filename) { free(opts->in_filename); }
  if (opts->out_dirname) { free(opts->out_dirname); }
  
  if (opts->exec_name) { free(opts->exec_name); }
  if (opts->command_name) { free(opts->command_name); }
  
  free(opts);
}

//------------------------------------------------------------------------

void sort_options_validate(sort_options_t *opts) {
  if (! exists(opts->in_filename)) {
    printf("\nError: Input file name not found !\n\n");
    usage_sort_options(opts);
  }

  if (! exists(opts->out_dirname)) {
    opts->out_dirname = strdup(".");
  }

  if (opts->max_memory < 0) {
    printf("\nError: Invalid maximum memory (%i), it must be greater than 0 !\n\n",
	   opts->max_memory);
    usage_sort_options(opts);
  }

  if (!opts->criteria) {
    opts->criteria = strdup("coord");
  }

  if (strcmp("name", opts->criteria)  != 0 &&
      strcmp("coord", opts->criteria) != 0    ) {
    printf("\nError: Invalid criteria by sorting (%s), valid values are 'coord' to sort by chromosomal coordinates (default value), and 'name' to sort by read names\n", opts->criteria);
    usage_sort_options(opts);
  }
}

//------------------------------------------------------------------------

void sort_options_display(sort_options_t *opts) {
  printf("PARAMETERS CONFIGURATION\n");
  printf("=================================================\n");
  printf("Main options\n");
  printf("\tBAM input filename  : %s\n", opts->in_filename);
  printf("\tOutput dirname      : %s\n", opts->out_dirname);
  printf("\n");

  printf("Sort criteria\n");
  if (strcmp("name", opts->criteria) == 0) {
    printf("\tby read names\n");
  } else {
    printf("\tby chromosomal coordinates\n");
  }
  printf("\n");

  printf("Architecture options\n");
  printf("\tMax. memory: %lu\n", opts->max_memory);
  printf("=================================================\n");
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

sort_options_t *read_cli_sort_options(void **argtable, sort_options_t *opts) {	
  if (((struct arg_int*)argtable[0])->count) { opts->help = ((struct arg_int*)argtable[0])->count; }
  if (((struct arg_file*)argtable[1])->count) { opts->in_filename = strdup(*(((struct arg_file*)argtable[1])->filename)); }
  if (((struct arg_file*)argtable[2])->count) { opts->out_dirname = strdup(*(((struct arg_file*)argtable[2])->filename)); }
  if (((struct arg_int*)argtable[3])->count) { opts->max_memory = *(((struct arg_int*)argtable[3])->ival); }
  if (((struct arg_str*)argtable[4])->count) { opts->criteria = strdup(*(((struct arg_str*)argtable[4])->sval)); }
  
  return opts;
}

//--------------------------------------------------------------------

void usage_sort_options(sort_options_t *opts) {
  void **argtable = new_argtable_sort_options();
  usage_argtable(opts->exec_name, opts->command_name, argtable);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

void** new_argtable_sort_options() {
  void **argtable = (void**)malloc((NUM_SORT_OPTIONS + 1) * sizeof(void*));
  
  // NOTICE that order cannot be changed as is accessed by index in other functions
  argtable[0] = arg_lit0("h", "help", "Help option");
  argtable[1] = arg_file0("b", "bam-file", NULL, "Input file name (BAM format)");
  argtable[2] = arg_file0("o", "outdir", NULL, "Output file name (BAM format)");
  argtable[3] = arg_int0(NULL, "max-memory", NULL, "Approximately the maximum required memory [500000000]");
  argtable[4] = arg_str0("c", "criteria", NULL, "Sorting criteria: 'coord' to sort by chromosomal coordinates, and 'name' by read names [coord]");
  
  argtable[NUM_SORT_OPTIONS] = arg_end(20);
  
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
