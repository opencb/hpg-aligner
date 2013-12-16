#include "options.h"

const char DEFAULT_OUTPUT_NAME[30] = "hpg-aligner_output";
//const char SPLICE_EXACT_FILENAME[30]   = "exact_junctions.bed";
//const char SPLICE_EXTEND_FILENAME[30]  = "extend_junctions.bed";
//const char INDEX_NAME[30]  = "index";

options_t *options_new(void) {
  options_t *options = (options_t*) calloc (1, sizeof(options_t));
  size_t num_cores = 0;
  
  //======================= COMMON OPTIONS ====================
  options->in_filename = NULL;
  options->in_filename2 = NULL;
  options->report_all =  0;
  options->log_level = LOG_INFO_LEVEL;
  options->output_name = strdup(DEFAULT_OUTPUT_NAME);
  options->num_gpu_threads = DEFAULT_GPU_THREADS;
  //GET Number System Cores
  //----------------------------------------------
  if (num_cores = get_optimal_cpu_num_threads()) {
    options->num_cpu_threads = num_cores;
  }else {
    options->num_cpu_threads = DEFAULT_CPU_THREADS;
  }
  //----------------------------------------------
  options->max_intron_length = DEFAULT_MAX_INTRON_LENGTH;
  options->num_seeds = DEFAULT_NUM_SEEDS;
  options->min_num_seeds_in_cal = DEFAULT_MIN_NUM_SEEDS_IN_CAL;
  options->cal_seeker_errors = DEFAULT_CAL_SEEKER_ERRORS;
  options->write_size = DEFAULT_WRITE_BATCH_SIZE;
  options->min_seed_padding_left = DEFAULT_MIN_SEED_PADDING_LEFT;
  options->min_seed_padding_right = DEFAULT_MIN_SEED_PADDING_RIGHT;

  options->min_score = DEFAULT_MIN_SCORE;

  options->match = DEFAULT_SW_MATCH;
  options->mismatch = DEFAULT_SW_MISMATCH;
  options->gap_open = DEFAULT_SW_GAP_OPEN;
  options->gap_extend = DEFAULT_SW_GAP_EXTEND;
  options->min_intron_length = DEFAULT_MIN_INTRON_LENGTH;
  options->pair_mode = DEFAULT_PAIR_MODE;
  options->pair_min_distance = DEFAULT_PAIR_MIN_DISTANCE;
  options->pair_max_distance = DEFAULT_PAIR_MAX_DISTANCE;
  options->timming = 0;
  options->statistics = 0;
  options->report_n_best = 0;
  options->report_n_hits = 0;
  options->report_best = 0;
  options->report_only_paired = 0;
  options->workflow_enable = 1;
  options->bwt_set = 0;
  options->reg_set = 0;
  options->cal_set = 0;
  options->sw_set = 0;
  options->filter_read_mappings = DEFAULT_FILTER_READ_MAPPINGS;
  options->filter_seed_mappings = DEFAULT_FILTER_SEED_MAPPINGS;
  //=========================================================

  options->min_cal_size = 0; 
  options->seeds_max_distance = 0;
  options->batch_size = 0;
  options->min_seed_size = 0;
  options->seed_size = 0;
  options->flank_length = 0;

  //new variables for bisulphite case in index generation
  options->bs_index = 0;
  return options;
}

void validate_options(options_t *options, char *mode) {
  int value_dir = exists(options->output_name);
  int DEFAULT_READ_BATCH_SIZE;
  int DEFAULT_SEED_SIZE;
  int DEFAULT_FLANK_LENGTH;
  int DEFAULT_MIN_SEED_SIZE;
  int DEFAULT_MIN_CAL_SIZE;
  int DEFAULT_SEEDS_MAX_DISTANCE;

  if (strcmp("dna", mode) == 0) {
    strcpy(options->mode, "DNA");
    DEFAULT_READ_BATCH_SIZE = 20000;
    DEFAULT_SEED_SIZE	= 20;
    DEFAULT_FLANK_LENGTH = 5;
    DEFAULT_MIN_SEED_SIZE = 16;
    DEFAULT_MIN_CAL_SIZE = 30;
    DEFAULT_SEEDS_MAX_DISTANCE = 100;
  }else if (strcmp("bs", mode) == 0) {
    strcpy(options->mode, "BISULFITE");
    DEFAULT_READ_BATCH_SIZE = 20000;
    DEFAULT_SEED_SIZE	= 20;
    DEFAULT_FLANK_LENGTH = 5;
    DEFAULT_MIN_SEED_SIZE = 16;
    DEFAULT_MIN_CAL_SIZE = 30;
    DEFAULT_SEEDS_MAX_DISTANCE = 100;
  }else if (strcmp("rna", mode) == 0) {
    strcpy(options->mode, "RNA");
    DEFAULT_READ_BATCH_SIZE = 200000;
    DEFAULT_SEED_SIZE = 16;
    DEFAULT_FLANK_LENGTH = 30;
    DEFAULT_MIN_SEED_SIZE = 16;
    DEFAULT_MIN_CAL_SIZE = 20;
    DEFAULT_SEEDS_MAX_DISTANCE = 60;
    options->pair_max_distance = DEFAULT_PAIR_MAX_DISTANCE + options->max_intron_length;
  }

  if (strcmp("dna", mode) == 0 || strcmp("rna", mode) == 0) {
    if (!value_dir) {
      create_directory(options->output_name);
    }
    
    if (!options->in_filename) {
      printf("Not filename input found. Please, insert it with option '-f FILENAME'.\n");
      usage_cli();
    }

    
    if (!options->bwt_dirname) {
      printf("Not BWT index input found. Please, insert it with option '-i DIRNAME'.\n");
      usage_cli();
    }
  }

  if (!options->min_cal_size) {
    options->min_cal_size = DEFAULT_MIN_CAL_SIZE;
  }
 
  if (!options->seeds_max_distance) {
    options->seeds_max_distance = DEFAULT_SEEDS_MAX_DISTANCE;
  }
  
  if (!options->batch_size) {
    options->batch_size = DEFAULT_READ_BATCH_SIZE;
  }
   
  if (options->min_seed_size <= 10) {
    options->min_seed_size = DEFAULT_MIN_SEED_SIZE;
  }
  
  if (options->seed_size <= 10) {
    options->seed_size = DEFAULT_SEED_SIZE;
  }

  if (!options->flank_length) {
    options->flank_length = DEFAULT_FLANK_LENGTH;
  }


  if (options->report_best) {
    options->report_all = 0;
    options->report_n_hits = 0;
    options->report_n_best = 0;
  }else if (options->report_n_best) {
    options->report_all = 0;
    options->report_n_hits = 0;   
    options->report_best = 0; 
  } else if (options->report_n_hits) {
    options->report_all = 0;
    options->report_n_best = 0;
    options->report_best = 0;
  } else if (options->report_all) {
    options->report_n_best = 0;
    options->report_n_hits = 0;
    options->report_best = 0;
  } else {
    options->report_best = 1;
    options->report_n_best = 0;
    options->report_n_hits = 0;    
    options->report_all = 0;
  }
}


void options_free(options_t *options) {
     if(options == NULL) { return; }

     if (options->in_filename  != NULL)	{ free(options->in_filename); }
     if (options->in_filename2  != NULL) { free(options->in_filename2); }
     if (options->bwt_dirname  != NULL)	{ free(options->bwt_dirname); }     
     if (options->genome_filename  != NULL) { free(options->genome_filename); }
     if (options->output_name  != NULL)	{ free(options->output_name); }
     if (options->prefix_name != NULL) { free(options->prefix_name); }
     if (options->transcriptome_filename != NULL) { free(options->transcriptome_filename); }

     free(options);
}


void options_display(options_t *options) {
     char* in_filename = strdup(options->in_filename);
     char* in_filename2 = NULL;
     if (options->in_filename2 != NULL) {
	  in_filename2 = strdup(options->in_filename2);
     }
     char* bwt_dirname =  strdup(options->bwt_dirname);
     char* genome_filename =  NULL;
     if (options->genome_filename != NULL) {
	  genome_filename =  strdup(options->genome_filename);
     }
     unsigned int  report_all = (unsigned int)options->report_all;
     unsigned int  report_n_best = (unsigned int)options->report_n_best;
     unsigned int  report_n_hits = (unsigned int)options->report_n_hits;
     unsigned int  report_only_paired = (unsigned int)options->report_only_paired;
     unsigned int  report_best = (unsigned int)options->report_best;
          
     char* output_name =  strdup(options->output_name);
     unsigned int num_gpu_threads =  (unsigned int)options->num_gpu_threads;
     unsigned int num_cpu_threads =  (unsigned int)options->num_cpu_threads;
     unsigned int cal_seeker_errors =  (unsigned int)options->cal_seeker_errors; 
     unsigned int min_cal_size =  (unsigned int)options->min_cal_size; 
     unsigned int seeds_max_distance =  (unsigned int)options->seeds_max_distance; 
     unsigned int batch_size =  (unsigned int)options->batch_size; 
     unsigned int write_size =  (unsigned int)options->write_size;  
     unsigned int min_seed_size =  (unsigned int)options->min_seed_size;
     unsigned int seed_size =  (unsigned int)options->seed_size;
     unsigned int num_seeds =  (unsigned int)options->num_seeds;
     int min_num_seeds_in_cal =  (int)options->min_num_seeds_in_cal;
     unsigned int max_intron_length =  (unsigned int)options->max_intron_length;
     unsigned int flank_length =  (unsigned int)options->flank_length;
     unsigned int pair_mode =  (unsigned int)options->pair_mode;
     unsigned int pair_min_distance =  (unsigned int)options->pair_min_distance;
     unsigned int pair_max_distance =  (unsigned int)options->pair_max_distance;
     unsigned int min_intron_length =  (unsigned int)options->min_intron_length;
     //unsigned int gpu_process = (unsigned int)options->gpu_process;

     int min_score    =  (int)options->min_score;
     float match      =  (float)options->match;
     float mismatch   =  (float)options->mismatch;
     float gap_open   =  (float)options->gap_open;
     float gap_extend =  (float)options->gap_extend;

     printf("\n");
     printf("+--------------------------------------------------------------------------------------+\n");
     printf("|                               PARAMETERS CONFIGURATION                               |\n");
     printf("+--------------------------------------------------------------------------------------+\n");
     //     printf("Num gpu threads %d\n", num_gpu_threads);
     //     printf("GPU Process: %s\n",  gpu_process == 0 ? "Disable":"Enable");
     printf("General parameters\n");
     printf("\tMode: %s\n", options->mode);
     if (in_filename2) {
       printf("\tInput FastQ filename, pair #1: %s\n", in_filename);
       printf("\tInput FastQ filename, pair #2: %s\n", in_filename2);
     } else {
       printf("\tInput FastQ filename: %s\n", in_filename);
     }
     printf("\tBWT index directory name: %s\n", bwt_dirname);
     printf("\tOutput directory name: %s\n", output_name);
     printf("\n");
     printf("Architecture parameters\n");
     printf("\tNumber of cpu threads %d\n",  num_cpu_threads);
     //printf("CAL seeker errors: %d\n",  cal_seeker_errors);
     printf("\tBatch size: %d bytes\n",  batch_size);
     //     printf("\tWrite size: %d bytes\n",  write_size);
     printf("\n");
     printf("Report parameters\n");
     printf("\tReport all hits: %s\n",  report_all == 0 ? "Disable":"Enable");
     printf("\tReport n best hits: %d\n",  report_n_best);
     printf("\tReport n hits: %d\n",  report_n_hits);
     printf("\tReport best hits: %s\n",  report_best == 0 ? "Disable":"Enable");
     printf("\tReport unpaired reads: %s\n",  report_only_paired == 0 ? "Enable":"Disable");
     printf("\n");
     printf("Seeding and CAL parameters\n");
     printf("\tNumber of seeds: %d\n",  num_seeds);
     if (seed_size) {
       printf("\tSeed size: %d\n",  seed_size);
       printf("\tMin seed size: %d\n",  min_seed_size);
     }
     else {
       printf("\tSeeds optimus autoconf\n");
     }
     printf("\tMin CAL size: %d\n",  min_cal_size);
     if (min_num_seeds_in_cal < 0) {
       printf("\tMin. number of seeds in a CAL: the maximun\n");
     } else {
       printf("\tMin. number of seeds in a CAL: %d\n",  min_num_seeds_in_cal);
     }
     printf("\tSeeds max distance: %d\n",  seeds_max_distance);
     printf("\tFlank length: %d\n", flank_length);
     printf("\n");
     printf("Mapping filters\n");
     printf("\tFor reads: %d mappings maximum, otherwise discarded\n", options->filter_read_mappings);
     printf("\tFor seeds: %d mappings maximum, otherwise discarded\n", options->filter_seed_mappings);
     printf("\n");
     printf("Pair-mode parameters\n");
     printf("\tPair mode: %d\n", pair_mode);
     printf("\tMin. distance: %d\n", pair_min_distance);
     printf("\tMax. distance: %d\n", pair_max_distance);
     printf("\n");
     printf("Smith-Waterman parameters\n");
     printf("\tMatch      : %0.4f\n", match);
     printf("\tMismatch   : %0.4f\n", mismatch);
     printf("\tGap open   : %0.4f\n", gap_open);
     printf("\tGap extend : %0.4f\n", gap_extend);
     printf("\n");

     if (strcmp(options->mode, "RNA") == 0) {
       printf("RNA parameters\n");
       printf("\tMax intron length: %d\n", max_intron_length);
       printf("\tMin intron length: %d\n", min_intron_length);
       printf("\tMin score        : %d\n", min_score);
     }
     printf("+--------------------------------------------------------------------------------------+\n");
     
     free(in_filename);
     if (in_filename2 != NULL) free(in_filename2);
     free(bwt_dirname);
     free(genome_filename);
     free(output_name);
}

//--------------------------------------------------------------------

void** argtable_options_new(void) {
     void **argtable = (void**)malloc((NUM_OPTIONS + 1) * sizeof(void*));	// NUM_OPTIONS +1 to allocate end structure
     // NOTICE that order cannot be changed as is accessed by index in other functions
     argtable[0] = arg_file0("f", "fq,fastq", NULL, "Reads file input");
     argtable[1] = arg_file0("i", "bwt-index", NULL, "BWT directory name");
     argtable[2] = arg_int0("l", "log-level", NULL, "Log debug level");
     argtable[3] = arg_lit0(NULL, "report-all", "Report all alignments");
     argtable[4] = arg_file0("o", "outdir", NULL, "Output directory");
     argtable[5] = arg_int0(NULL, "gpu-threads", NULL, "Number of GPU Threads");
     argtable[6] = arg_int0(NULL, "cpu-threads", NULL, "Number of CPU Threads");
     argtable[7] = arg_int0("r", "index-ratio", NULL, "BWT index compression ratio");
     argtable[8] = arg_int0(NULL, "cal-seeker-errors", NULL, "Number of errors in CAL Seeker");
     argtable[9] = arg_int0(NULL, "min-cal-size", NULL, "Minimum CAL size");
     argtable[10] = arg_int0(NULL, "max-distance-seeds", NULL, "Maximum distance between seeds");
     argtable[11] = arg_int0(NULL, "read-batch-size", NULL, "Batch Size");
     argtable[12] = arg_int0(NULL, "write-batch-size", NULL, "Write Size");
     argtable[13] = arg_int0(NULL, "num-cal-seekers", NULL, "Number of CAL Seekers");
     argtable[14] = arg_int0(NULL, "extra-seed-left-padding", NULL, "Nucleotides padding of left min seed");
     argtable[15] = arg_int0(NULL, "extra-seed-right-padding", NULL, "Nucleotides padding of right min seed");
     argtable[16] = arg_file0(NULL, "transcriptome-file", NULL, "Transcriptome file to help search splice junctions");
     argtable[17] = arg_int0(NULL, "seed-size", NULL, "Number of nucleotides in a seed");
     argtable[18] = arg_int0(NULL, "min-seed-size", NULL, "Minimum number of nucleotides in a seed");
     argtable[19] = arg_int0(NULL, "cal-flank-size", NULL, "Flank length for CALs");
     argtable[20] = arg_dbl0(NULL, "sw-match", NULL, "Match value for Smith-Waterman algorithm");
     argtable[21] = arg_dbl0(NULL, "sw-mismatch", NULL, "Mismatch value for Smith-Waterman algorithm");
     argtable[22] = arg_dbl0(NULL, "sw-gap-open", NULL, "Gap open penalty for Smith-Waterman algorithm");
     argtable[23] = arg_dbl0(NULL, "sw-gap-extend", NULL, "Gap extend penalty for Smith-Waterman algorithm");

     argtable[24] = arg_int0(NULL, "min-score", NULL, "Minimum score for valid mappings");

     argtable[25] = arg_int0(NULL, "max-intron-size", NULL, "Maximum intron size");
     argtable[26] = arg_int0(NULL, "min-intron-size", NULL, "Minimum intron size");
     argtable[27] = arg_lit0("t", "time", "Timming mode active");
     argtable[28] = arg_lit0("s", "stats", "Statistics mode active");
     argtable[29] = arg_lit0("h", "help", "Help option");
     argtable[30] = arg_str0(NULL, "prefix", NULL, "File prefix name");
     argtable[31] = arg_file0("g", "ref-genome", NULL, "Reference genome");
     argtable[32] = arg_file0("j", "fq2,fastq2", NULL, "Reads file input #2 (for paired mode)");
     argtable[33] = arg_int0(NULL, "paired-mode", NULL, "Pair mode: 0 = single-end, 1 = paired-end, 2 = mate-pair [Default 0]");
     argtable[34] = arg_int0(NULL, "paired-min-distance", NULL, "Minimum distance between pairs");
     argtable[35] = arg_int0(NULL, "paired-max-distance", NULL, "Maximum distance between pairs");
     argtable[36] = arg_int0(NULL, "report-n-best", NULL, "Report the <n> best alignments");
     argtable[37] = arg_int0(NULL, "report-n-hits", NULL, "Report <n> hits");
     argtable[38] = arg_int0(NULL, "num-seeds", NULL, "Number of seeds per read");
     argtable[39] = arg_int0(NULL, "min-num-seeds", NULL, "Minimum number of seeds to create a CAL (if -1, the maxixum will be taken)");
     argtable[40] = arg_lit0(NULL, "workflow-disable", "Disable Workflow Pipeline");
     argtable[41] = arg_lit0(NULL, "report-only-paired", "Report only the paired reads");
     argtable[42] = arg_int0(NULL, "filter-read-mappings", NULL, "Reads that map in more than <n> locations are discarded");
     argtable[43] = arg_int0(NULL, "filter-seed-mappings", NULL, "Seeds that map in more than <n> locations are discarded");
     argtable[44] = arg_lit0(NULL, "report-best", "Report all alignments with best score");
     argtable[45] = arg_lit0(NULL, "bs-index", "Indicate the use of bisulphite generation of the index");

     argtable[NUM_OPTIONS] = arg_end(20);
     
     return argtable;
}


void argtable_options_free(void **argtable) {
     if(argtable != NULL) {
	  arg_freetable(argtable, NUM_OPTIONS + 1);	// struct end must also be freed
	  free(argtable);
     }
}


int read_config_file(const char *filename, options_t *options) {
	if (filename == NULL || options == NULL) {
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

	/*if(config_lookup_string(config, "app.outdir", &tmp_string)) { options->output_directory = strdup(tmp_string); }
	if(config_lookup_int(config, "app.cpu-num-threads", &tmp_int)) { options->cpu_num_threads = (int)tmp_int; }
	*/

	config_destroy(config);
	free(config);
//	free(tmp_string);

	return ret_code;
}


/**
 * @brief Initializes an options_t structure from argtable parsed CLI with default values. Notice that options are order dependent.
 * @return A new options_t structure initialized with default values.
 *
 * Initializes the only default options from options_t.
 */
options_t *read_CLI_options(void **argtable, options_t *options) {	
  if (((struct arg_file*)argtable[0])->count) { options->in_filename = strdup(*(((struct arg_file*)argtable[0])->filename)); }
  if (((struct arg_file*)argtable[1])->count) { options->bwt_dirname = strdup(*(((struct arg_file*)argtable[1])->filename)); }
  if (((struct arg_file*)argtable[2])->count) { options->log_level = *(((struct arg_int*)argtable[2])->ival); }
  if (((struct arg_file*)argtable[3])->count) { options->report_all = (((struct arg_int *)argtable[3])->count); }
  if (((struct arg_file*)argtable[4])->count) { free(options->output_name); options->output_name = strdup(*(((struct arg_file*)argtable[4])->filename)); }  
  if (((struct arg_int*)argtable[5])->count) { options->num_gpu_threads = *(((struct arg_int*)argtable[5])->ival); }
  if (((struct arg_int*)argtable[6])->count) { options->num_cpu_threads = *(((struct arg_int*)argtable[6])->ival); }
  if (((struct arg_int*)argtable[7])->count) { options->index_ratio = *(((struct arg_int*)argtable[7])->ival); }
  if (((struct arg_int*)argtable[8])->count) { options->cal_seeker_errors = *(((struct arg_int*)argtable[8])->ival); }
  if (((struct arg_int*)argtable[9])->count) { options->min_cal_size = *(((struct arg_int*)argtable[9])->ival); }
  if (((struct arg_int*)argtable[10])->count) { options->seeds_max_distance = *(((struct arg_int*)argtable[10])->ival); }
  if (((struct arg_int*)argtable[11])->count) { options->batch_size = *(((struct arg_int*)argtable[11])->ival); }
  if (((struct arg_int*)argtable[12])->count) { options->write_size = *(((struct arg_int*)argtable[12])->ival); }
  if (((struct arg_int*)argtable[13])->count) { options->cal_set = 1; options->num_cal_seekers = *(((struct arg_int*)argtable[13])->ival); }
  if (((struct arg_int*)argtable[14])->count) { options->min_seed_padding_left = *(((struct arg_int*)argtable[14])->ival); }
  if (((struct arg_int*)argtable[15])->count) { options->min_seed_padding_right = *(((struct arg_int*)argtable[15])->ival); }
  if (((struct arg_file*)argtable[16])->count) { options->transcriptome_filename = strdup(*(((struct arg_file*)argtable[16])->filename)); }
  if (((struct arg_int*)argtable[17])->count) { options->seed_size = *(((struct arg_int*)argtable[17])->ival); }
  if (((struct arg_int*)argtable[18])->count) { options->min_seed_size = *(((struct arg_int*)argtable[18])->ival); }
  if (((struct arg_int*)argtable[19])->count) { options->flank_length = *((struct arg_int*)argtable[19])->ival; }
  if (((struct arg_dbl*)argtable[20])->count) { options->match = *((struct arg_dbl*)argtable[20])->dval; }
  if (((struct arg_dbl*)argtable[21])->count) { options->mismatch = *(((struct arg_dbl*)argtable[21])->dval); }
  if (((struct arg_dbl*)argtable[22])->count) { options->gap_open = *(((struct arg_dbl*)argtable[22])->dval); }
  if (((struct arg_dbl*)argtable[23])->count) { options->gap_extend = *(((struct arg_dbl*)argtable[23])->dval); }

  if (((struct arg_int*)argtable[24])->count) { options->min_score = *(((struct arg_int*)argtable[24])->ival); }

  if (((struct arg_int*)argtable[25])->count) { options->max_intron_length = *(((struct arg_int*)argtable[25])->ival); }
  if (((struct arg_int*)argtable[26])->count) { options->min_intron_length = *(((struct arg_int*)argtable[26])->ival); }
  if (((struct arg_int*)argtable[27])->count) { options->timming = ((struct arg_int*)argtable[27])->count; }
  if (((struct arg_int*)argtable[28])->count) { options->statistics = ((struct arg_int*)argtable[28])->count; }
  if (((struct arg_int*)argtable[29])->count) { options->help = ((struct arg_int*)argtable[29])->count; }
  if (((struct arg_str*)argtable[30])->count) { options->prefix_name = strdup(*(((struct arg_str*)argtable[30])->sval)); }
  if (((struct arg_file*)argtable[31])->count) { options->genome_filename = strdup(*(((struct arg_file*)argtable[31])->filename)); }
  if (((struct arg_file*)argtable[32])->count) { options->in_filename2 = strdup(*(((struct arg_file*)argtable[32])->filename)); }
  if (((struct arg_int*)argtable[33])->count) { options->pair_mode = *(((struct arg_int*)argtable[33])->ival); }
  if (((struct arg_int*)argtable[34])->count) { options->pair_min_distance = *(((struct arg_int*)argtable[34])->ival); }
  if (((struct arg_int*)argtable[35])->count) { options->pair_max_distance = *(((struct arg_int*)argtable[35])->ival); }
  if (((struct arg_int*)argtable[36])->count) { options->report_n_best = *(((struct arg_int*)argtable[36])->ival); }
  if (((struct arg_int*)argtable[37])->count) { options->report_n_hits = *(((struct arg_int*)argtable[37])->ival); }
  if (((struct arg_int*)argtable[38])->count) { options->num_seeds = *(((struct arg_int*)argtable[38])->ival); }
  if (((struct arg_int*)argtable[39])->count) { options->min_num_seeds_in_cal = *(((struct arg_int*)argtable[39])->ival); }
  if (((struct arg_int*)argtable[40])->count) { options->workflow_enable = (((struct arg_int *)argtable[40])->count) == 1 ? 0 : 1; }
  if (((struct arg_int*)argtable[41])->count) { options->report_only_paired = (((struct arg_int*)argtable[41])->count); }
  if (((struct arg_int*)argtable[42])->count) { options->filter_read_mappings = *(((struct arg_int*)argtable[42])->ival); }
  if (((struct arg_int*)argtable[43])->count) { options->filter_seed_mappings = *(((struct arg_int*)argtable[43])->ival); }
  if (((struct arg_int*)argtable[44])->count) { options->report_best = (((struct arg_int*)argtable[44])->count); }


  // new value

  if (((struct arg_int*)argtable[45])->count) { options->bs_index = (((struct arg_int*)argtable[45])->count); }

  return options;
}


options_t *parse_options(int argc, char **argv) {
  void **argtable = argtable_options_new();
  //	struct arg_end *end = arg_end(10);
  //	void **argtable = argtable_options_get(argtable_options, end);
  
  options_t *options = options_new();
  if (argc < 2) {
    usage(argtable);
    exit(-1);
  }else {

    int num_errors = arg_parse(argc, argv, argtable);

    if (((struct arg_int*)argtable[29])->count) {
      usage(argtable);
	argtable_options_free(argtable);
	options_free(options);
	exit(0);
    }
        
    if (num_errors > 0) {
      arg_print_errors(stdout, argtable[NUM_OPTIONS], "hpg-aligner");	// struct end is always allocated in the last position
      usage(argtable);
      exit(-1);
    }else {
      options = read_CLI_options(argtable, options);
      if(options->help) {
	usage(argtable);
	argtable_options_free(argtable);
	options_free(options);
	exit(0);
      }
      // Check if 'help' option has been provided.
    }
    
  }
  //	exit:
  argtable_options_free(argtable);
  //	free(end);
  //	free(argtable_options);

  // in previous versions, min. flank length was 20
  //  if (options->flank_length < 5) {
  //    options->flank_length = 5;
  //  }

  return options;
}

void usage(void **argtable) {
  printf("Usage:\n./hpg-aligner {dna | rna | bs | build-index}");
  arg_print_syntaxv(stdout, argtable, "\n");
  arg_print_glossary(stdout, argtable, "%-50s\t%s\n");
}

void usage_cli() {
  void **argtable = argtable_options_new();
  usage(argtable);
  exit(0);
}
