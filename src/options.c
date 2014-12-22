#include "options.h"

const char DEFAULT_OUTPUT_NAME[29] = "hpg_aligner_output";

//const char SPLICE_EXACT_FILENAME[30]   = "exact_junctions.bed";
//const char SPLICE_EXTEND_FILENAME[30]  = "extend_junctions.bed";
//const char INDEX_NAME[30]  = "index";

//--------------------------------------------------------------------

options_t *options_new(void) {
  options_t *options = (options_t*) calloc (1, sizeof(options_t));
  
  //======================= COMMON OPTIONS ====================
  options->version = 0;
  options->in_filename = NULL;
  options->in_filename2 = NULL;
  options->report_all =  0;
  options->log_level = LOG_ERROR_LEVEL;
  options->output_name = strdup(DEFAULT_OUTPUT_NAME);
  options->num_gpu_threads = DEFAULT_GPU_THREADS;

  options->realignment = 0;
  options->recalibration = 0;
  options->adapter = NULL;
  options->input_format = FASTQ_FORMAT;
  options->output_format = SAM_FORMAT;

  //GET Number System Cores
  //----------------------------------------------
  options->num_cpu_threads = get_optimal_cpu_num_threads();
  /*
  if (num_cores == get_optimal_cpu_num_threads()) {
    options->num_cpu_threads = num_cores;
  } else {
    options->num_cpu_threads = DEFAULT_CPU_THREADS;
  }
  */
  //----------------------------------------------
  options->max_intron_length = DEFAULT_MAX_INTRON_LENGTH;
  options->num_seeds = 0;
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
  options->fast_mode = 1;

  options->adapter_length = 0;

  options->adapter_length = 0;
  options->set_bam_format = 0;
  options->set_cal = 0;

  //new variables for bisulphite case in index generation
  options->bs_index = 0;

  options->cmdline = NULL;

  return options;
}

//--------------------------------------------------------------------

void validate_options(options_t *options) {
  int value_dir = exists(options->output_name);
  int DEFAULT_READ_BATCH_SIZE = 200000;
  int DEFAULT_SEED_SIZE = 18;
  int DEFAULT_FLANK_LENGTH = 5;
  int DEFAULT_NUM_SEEDS = 20;
  int DEFAULT_MIN_SEED_SIZE = 16;
  int DEFAULT_MIN_CAL_SIZE = 20;
  int DEFAULT_SEEDS_MAX_DISTANCE = 100;

  int mode = options->mode;

  if (mode == DNA_MODE) {
    strcpy(options->str_mode, "DNA");
    DEFAULT_READ_BATCH_SIZE = DEFAULT_DNA_READ_BATCH_SIZE;
    DEFAULT_NUM_SEEDS	= DEFAULT_DNA_NUM_SEEDS;
    DEFAULT_SEED_SIZE	= 18;
    DEFAULT_FLANK_LENGTH = 5;
    DEFAULT_MIN_SEED_SIZE = 16;
    DEFAULT_MIN_CAL_SIZE = DEFAULT_DNA_MIN_CAL_SIZE;
    DEFAULT_SEEDS_MAX_DISTANCE = 100;
  } else if (mode == RNA_MODE) {
    strcpy(options->str_mode, "RNA");
    DEFAULT_READ_BATCH_SIZE = DEFAULT_RNA_READ_BATCH_SIZE;;
    DEFAULT_NUM_SEEDS	= DEFAULT_RNA_NUM_SEEDS;
    DEFAULT_SEED_SIZE = DEFAULT_RNA_SEED_SIZE;
    DEFAULT_FLANK_LENGTH = 30;
    DEFAULT_MIN_SEED_SIZE = 16;
    DEFAULT_MIN_CAL_SIZE = DEFAULT_RNA_MIN_CAL_SIZE;
    DEFAULT_SEEDS_MAX_DISTANCE = 60;
    options->pair_max_distance = DEFAULT_PAIR_MAX_DISTANCE + options->max_intron_length;
  }

  //printf("Value %i mode => %i == (%i|%i)\n", value_dir, mode, DNA_MODE, RNA_MODE);
  if (mode == DNA_MODE || mode == RNA_MODE) {
    if (!value_dir) {
      //printf("Create directory %s\n", options->output_name);
      create_directory(options->output_name);
    }
    
    if (!options->in_filename) {
      printf("Filename input is missing. Please, insert it with option '-f FILENAME'.\n");
      usage_cli(mode);
    }

    
    if (!options->bwt_dirname) {
      printf("Index directory is missing. Please, insert it with option '-i DIRNAME'.\n");
      usage_cli(mode);
    }
  }

  if (!options->num_seeds) {
    options->num_seeds = DEFAULT_NUM_SEEDS;
  }    

  if (!options->min_cal_size || options->min_cal_size <= 20) {
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

  //if (!options->fast_mode) {
  //options->bam_format = 1;
  //}
}

//--------------------------------------------------------------------

void options_free(options_t *options) {
     if(!options) { return; }

     if (options->in_filename)	{ free(options->in_filename); }
     if (options->in_filename2) { free(options->in_filename2); }
     if (options->bwt_dirname)	{ free(options->bwt_dirname); }     
     if (options->genome_filename) { free(options->genome_filename); }
     if (options->output_name)	{ free(options->output_name); }
     if (options->prefix_name) { free(options->prefix_name); }
     if (options->adapter) { free(options->adapter); }

     if (options->mode == RNA_MODE) {
       if (options->transcriptome_filename) { free(options->transcriptome_filename); }
     }

     if (options->cmdline)  { free(options->cmdline); }

     free(options);
}

//--------------------------------------------------------------------

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

     unsigned int num_cpu_threads =  (unsigned int)options->num_cpu_threads;

     unsigned int min_cal_size =  (unsigned int)options->min_cal_size; 

     unsigned int batch_size =  (unsigned int)options->batch_size; 


     unsigned int seed_size =  (unsigned int)options->seed_size;
     unsigned int num_seeds =  (unsigned int)options->num_seeds;

     unsigned int max_intron_length =  (unsigned int)options->max_intron_length;

     unsigned int pair_mode =  (unsigned int)options->pair_mode;
     unsigned int pair_min_distance =  (unsigned int)options->pair_min_distance;
     unsigned int pair_max_distance =  (unsigned int)options->pair_max_distance;
     unsigned int min_intron_length =  (unsigned int)options->min_intron_length;
     unsigned int fast_mode =   (unsigned int)options->fast_mode;
     char* adapter =  NULL;
     if (options->adapter) {
       adapter = strdup(options->adapter);
     }

     int input_format = options->input_format;

     //unsigned int gpu_process = (unsigned int)options->gpu_process;

     int min_score    =  (int)options->min_score;
     float match      =  (float)options->match;
     float mismatch   =  (float)options->mismatch;
     float gap_open   =  (float)options->gap_open;
     float gap_extend =  (float)options->gap_extend;

     printf("\n");
     printf("+===============================================================+\n");
     printf("|             HPG ALIGNER PARAMETERS CONFIGURATION              |\n");
     printf("+===============================================================+\n");
     //     printf("Num gpu threads %d\n", num_gpu_threads);
     //     printf("GPU Process: %s\n",  gpu_process == 0 ? "Disable":"Enable");
     printf("HPG Aligner version: %s\n", HPG_ALIGNER_VERSION);
     printf("Command line: %s\n", options->cmdline);
     printf("\n");
     printf("General options\n");
     printf("\tMode: %s\n", options->str_mode);
     if (in_filename2) {
       printf("\tInput FastQ filename, pair #1: %s\n", in_filename);
       printf("\tInput FastQ filename, pair #2: %s\n", in_filename2);
     } else {
       printf("\tInput %s filename: %s\n", FORMAT_STR(input_format), in_filename);
     }
     if (input_format == FASTQ_FORMAT) {
       printf("\tFastQ gzip mode: %s\n", options->gzip == 1 ? "Enable" : "Disable");
     }
     printf("\tIndex directory name: %s\n", bwt_dirname);
     printf("\tOutput directory name: %s\n", output_name);
          
     printf("\tOutput file format: %s\n", 
	    (options->bam_format || options->realignment || options->recalibration) ? "BAM" : "SAM");
     printf("\tAdapter: %s\n", (adapter ? adapter : "Not present"));
     printf("\n");

     printf("Architecture options:\n");
     printf("\tNumber of cpu threads: %d\n",  num_cpu_threads);
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
     if (options->mode == DNA_MODE) {
       printf("\tNum. seeds: %d\n",  num_seeds);
     }
     if (options->mode == RNA_MODE && fast_mode) {
       printf("\tMin CAL coverage form a read (0 to 100): %d \n",  min_cal_size);
     } else {
       printf("\tMin CAL size: %d\n",  min_cal_size);
     }
     printf("\n");

     //printf("Mapping filters\n");
     //printf("\tFor reads: %d mappings maximum, otherwise discarded\n", options->filter_read_mappings);
     //printf("\tFor seeds: %d mappings maximum, otherwise discarded\n", options->filter_seed_mappings);
     //printf("\n");

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

     if (options->mode == RNA_MODE) {
       printf("RNA parameters\n");
       printf("\tMode: %s\n", fast_mode ? "SA":"BWT");
       if (!fast_mode) {
	 printf("\tSeed size: %d\n",  seed_size);
       }
       printf("\tMax intron length: %d\n", max_intron_length);
       printf("\tMin intron length: %d\n", min_intron_length);
       printf("\tMin score        : %d\n", min_score);
     }

     if (options->realignment || options->recalibration) {
       printf("Post-processing\n");
       if (options->realignment) {
	 printf("\tIndel realignment\n");
       }
       if (options->recalibration) {
	 printf("\tRecalibration\n");
       }
     }
     printf("+===============================================================+\n");
     
     free(in_filename);
     if (in_filename2 != NULL) free(in_filename2);
     free(bwt_dirname);
     free(genome_filename);
     free(output_name);
     if (adapter) free(adapter);
}

//--------------------------------------------------------------------

void options_to_file(options_t *options, FILE *fd) {
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

  unsigned int num_cpu_threads =  (unsigned int)options->num_cpu_threads;

  unsigned int min_cal_size =  (unsigned int)options->min_cal_size; 

  unsigned int batch_size =  (unsigned int)options->batch_size; 


  unsigned int seed_size =  (unsigned int)options->seed_size;
  unsigned int num_seeds =  (unsigned int)options->num_seeds;

  unsigned int max_intron_length =  (unsigned int)options->max_intron_length;

  unsigned int pair_mode =  (unsigned int)options->pair_mode;
  unsigned int pair_min_distance =  (unsigned int)options->pair_min_distance;
  unsigned int pair_max_distance =  (unsigned int)options->pair_max_distance;
  unsigned int min_intron_length =  (unsigned int)options->min_intron_length;
  unsigned int fast_mode =   (unsigned int)options->fast_mode;
  char* adapter =  NULL;
  if (options->adapter) {
    adapter = strdup(options->adapter);
  }

  int input_format = options->input_format;

  //unsigned int gpu_process = (unsigned int)options->gpu_process;

  int min_score    =  (int)options->min_score;
  float match      =  (float)options->match;
  float mismatch   =  (float)options->mismatch;
  float gap_open   =  (float)options->gap_open;
  float gap_extend =  (float)options->gap_extend;
   
  fprintf(fd, "= HPG Aligner version: %s\n", HPG_ALIGNER_VERSION);
  fprintf(fd, "= Command line: %s\n", options->cmdline);
  fprintf(fd, "=------------------------------------=\n");
  fprintf(fd, "= G E N E R A L    P A R A M E T E R S \n");
  fprintf(fd, "=------------------------------------=\n");
  fprintf(fd, "= Mode: %s\n", options->str_mode);
  if (in_filename2) {
    fprintf(fd, "= Input FastQ filename, pair #1: %s\n", in_filename);
    fprintf(fd, "= Input FastQ filename, pair #2: %s\n", in_filename2);
  } else {
    fprintf(fd, "= Input %s filename: %s\n", FORMAT_STR(input_format), in_filename);
  }
  if (input_format == FASTQ_FORMAT) {
    fprintf(fd, "= FastQ gzip mode: %s\n", options->gzip == 1 ? "Enable" : "Disable");
  }
  fprintf(fd, "= Index directory name: %s\n", bwt_dirname);
  fprintf(fd, "= Output directory name: %s\n", output_name);

  fprintf(fd, "= Output file format: %s\n", 
	 (options->bam_format || options->realignment || options->recalibration) ? "SAM" : "BAM");
  fprintf(fd, "= Adapter: %s\n", (adapter ? adapter : "Not present"));
  fprintf(fd, "\n\n");


  fprintf(fd, "=  A R C H I T E C T U R E    P A R A M E T E R S\n"); 
  fprintf(fd, "=-----------------------------------------------=\n");
  fprintf(fd, "= Number of cpu threads %d\n",  num_cpu_threads);
  fprintf(fd, "= Batch size: %d bytes\n",  batch_size);
  fprintf(fd, "\n\n");


  fprintf(fd, "= R E P O R T    P A R A M E T E R S\n");
  fprintf(fd, "=-----------------------------------------=\n");
  fprintf(fd, "= Report all hits: %s\n",  report_all == 0 ? "Disable":"Enable");
  fprintf(fd, "= Report n best hits: %d\n",  report_n_best);
  fprintf(fd, "= Report n hits: %d\n",  report_n_hits);
  fprintf(fd, "= Report best hits: %s\n",  report_best == 0 ? "Disable":"Enable");
  fprintf(fd, "= Report unpaired reads: %s\n",  report_only_paired == 0 ? "Enable":"Disable");
  fprintf(fd, "\n");


  fprintf(fd, "= S E E D I N G    A N D    C A L    P A R A M E T E R S \n");
  fprintf(fd, "=------------------------------------------------------=\n");
  if (options->mode == DNA_MODE) {
    fprintf(fd, "= Num. seeds: %d\n",  num_seeds);
  }

  if (options->mode == RNA_MODE && fast_mode) {
    fprintf(fd, "= Min CAL coverage form a read (0 to 100): %d \n",  min_cal_size);
  } else {
    fprintf(fd, "= Min CAL size: %d\n",  min_cal_size);
  }

  fprintf(fd, "\n\n");


  fprintf(fd, "= M A P P I N G    F I L T E R S \n");
  fprintf(fd, "=------------------------------=\n");
  fprintf(fd, "= For reads: %d mappings maximum, otherwise discarded\n", options->filter_read_mappings);
  fprintf(fd, "= For seeds: %d mappings maximum, otherwise discarded\n", options->filter_seed_mappings);
  fprintf(fd, "\n\n");


  fprintf(fd, "= P A I R - M O D E    P A R A M E T E R S \n");
  fprintf(fd, "=----------------------------------------=\n");
  fprintf(fd, "= Pair mode: %d\n", pair_mode);
  fprintf(fd, "= Min. distance: %d\n", pair_min_distance);
  fprintf(fd, "= Max. distance: %d\n", pair_max_distance);
  fprintf(fd, "\n\n");


  fprintf(fd, "= S M I T H - W A T E R M A N    P A R A M E T E R S \n");
  fprintf(fd, "=--------------------------------------------------=\n");
  fprintf(fd, "= Match      : %0.4f\n", match);
  fprintf(fd, "= Mismatch   : %0.4f\n", mismatch);
  fprintf(fd, "= Gap open   : %0.4f\n", gap_open);
  fprintf(fd, "= Gap extend : %0.4f\n", gap_extend);
  fprintf(fd, "\n\n");


  if (options->mode == RNA_MODE) {
    fprintf(fd, "= R N A    P A R A M E T E R S \n");
    fprintf(fd, "=----------------------------=\n");
    fprintf(fd, "= Mode: %s\n", fast_mode ? "Fast":"Slow");
    if (!fast_mode) { fprintf(fd, "= Seed size: %d\n",  seed_size); }
    fprintf(fd, "= Max intron length: %d\n", max_intron_length);
    fprintf(fd, "= Min intron length: %d\n", min_intron_length);
    fprintf(fd, "= Min score        : %d\n", min_score);
  }
  fprintf(fd, "\n");

  free(in_filename);
  free(bwt_dirname);

  if (options->in_filename2 != NULL) {
    free(in_filename2);
  }

  if (options->genome_filename != NULL) {
    free(genome_filename);
  }

  free(output_name);
  if (adapter) free(adapter);
}

//--------------------------------------------------------------------

void** argtable_options_new(int mode) {

  int num_options = NUM_OPTIONS;

  if      (mode == DNA_MODE) num_options += NUM_DNA_OPTIONS;
  else if (mode == RNA_MODE) num_options += NUM_RNA_OPTIONS;

  // NUM_OPTIONS +1 to allocate end structure
  void **argtable = (void**)malloc((num_options + 1) * sizeof(void*));	
  
  // NOTICE that order cannot be changed as is accessed by index in other functions
  int count = 0;
  argtable[count++] = arg_file1("f", "fq,fastq", NULL, "Reads file input. For more than one file: f1.fq,f2.fq,...");
  argtable[count++] = arg_file0("j", "fq2,fastq2", NULL, "Reads file input #2 (for paired mode)");
  argtable[count++] = arg_lit0("z", "gzip", "FastQ input files are gzip");
  argtable[count++] = arg_file1("i", "index", NULL, "Index directory name");
  argtable[count++] = arg_file0("o", "outdir", NULL, "Output directory");
  argtable[count++] = arg_int0(NULL, "filter-read-mappings", NULL, "Reads that map in more than <n> locations are discarded");
  argtable[count++] = arg_int0(NULL, "filter-seed-mappings", NULL, "Seeds that map in more than <n> locations are discarded");
  argtable[count++] = arg_int0(NULL, "min-cal-size", NULL, "Minimum CAL size");
  argtable[count++] = arg_int0("t", "cpu-threads", NULL, "Number of CPU Threads");
  argtable[count++] = arg_int0(NULL, "read-batch-size", NULL, "Batch Size");
  argtable[count++] = arg_dbl0(NULL, "sw-match", NULL, "Match value for Smith-Waterman algorithm");
  argtable[count++] = arg_dbl0(NULL, "sw-mismatch", NULL, "Mismatch value for Smith-Waterman algorithm");
  argtable[count++] = arg_dbl0(NULL, "sw-gap-open", NULL, "Gap open penalty for Smith-Waterman algorithm");
  argtable[count++] = arg_dbl0(NULL, "sw-gap-extend", NULL, "Gap extend penalty for Smith-Waterman algorithm");
  argtable[count++] = arg_int0(NULL, "min-score", NULL, "Minimum score for valid mappings (0 to 100 for RNA)"); //TODO: and DNA?
  //  argtable[count++] = arg_int0(NULL, "paired-mode", NULL, "Pair mode: 0 = single-end, 1 = paired-end, 2 = mate-pair [Default 0]");
  argtable[count++] = arg_int0(NULL, "paired-min-distance", NULL, "Minimum distance between pairs");
  argtable[count++] = arg_int0(NULL, "paired-max-distance", NULL, "Maximum distance between pairs");
  argtable[count++] = arg_lit0(NULL, "report-best", "Report all alignments with best score");
  argtable[count++] = arg_lit0(NULL, "report-all", "Report all alignments");
  argtable[count++] = arg_int0(NULL, "report-n-best", NULL, "Report the <n> best alignments");
  argtable[count++] = arg_int0(NULL, "report-n-hits", NULL, "Report <n> hits");
  argtable[count++] = arg_lit0(NULL, "report-only-paired", "Report only the paired reads");
  argtable[count++] = arg_str0(NULL, "prefix", NULL, "File prefix name");
  argtable[count++] = arg_int0("l", "log-level", NULL, "Log debug level");
  argtable[count++] = arg_lit0("h", "help", "Help option");
  argtable[count++] = arg_str0(NULL, "output-format", NULL, "BAM output format (otherwise, SAM format. This option is only available for SA mode, BWT mode always report in BAM format), this option turn the process slow");
  argtable[count++] = arg_lit0(NULL, "indel-realignment", "Indel-based realignment");
  argtable[count++] = arg_lit0(NULL, "recalibration", "Base quality score recalibration");
  argtable[count++] = arg_str0("a", "adapter", NULL, "Adapter sequence in the read");
  argtable[count++] = arg_str0(NULL, "input-format", NULL, "Input file format: fastq or bam. Default: fastq");
  argtable[count++] = arg_lit0("v", "version", "Display the HPG Aligner version");

  if (mode == DNA_MODE) {
    argtable[count++] = arg_int0(NULL, "num-seeds", NULL, "Number of seeds");
  } else if (mode == RNA_MODE) {
    argtable[count++] = arg_int0(NULL, "max-distance-seeds", NULL, "Maximum distance between seeds");
    argtable[count++] = arg_file0(NULL, "transcriptome-file", NULL, "Transcriptome file to help search splice junctions");
    argtable[count++] = arg_int0(NULL, "seed-size", NULL, "Number of nucleotides in a seed");
    argtable[count++] = arg_int0(NULL, "max-intron-size", NULL, "Maximum intron size");
    argtable[count++] = arg_int0(NULL, "min-intron-size", NULL, "Minimum intron size");
  }

  argtable[num_options] = arg_end(count);
     
  return argtable;
}


void argtable_options_free(void **argtable, int num_options) {
  if(argtable != NULL) {
    arg_freetable(argtable, num_options + 1);	// struct end must also be freed
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
  int count = -1;
  if (((struct arg_file*)argtable[++count])->count) { options->in_filename = strdup(*(((struct arg_file*)argtable[count])->filename)); }
  if (((struct arg_file*)argtable[++count])->count) { 
    options->pair_mode = 1;
    options->in_filename2 = strdup(*(((struct arg_file*)argtable[count])->filename)); 
  }
  if (((struct arg_int*)argtable[++count])->count) { options->gzip = ((struct arg_int*)argtable[count])->count; }
  if (((struct arg_file*)argtable[++count])->count) { options->bwt_dirname = strdup(*(((struct arg_file*)argtable[count])->filename)); }
  if (((struct arg_file*)argtable[++count])->count) { free(options->output_name); options->output_name = strdup(*(((struct arg_file*)argtable[count])->filename)); }  
  if (((struct arg_int*)argtable[++count])->count) { options->filter_read_mappings = *(((struct arg_int*)argtable[count])->ival); }
  if (((struct arg_int*)argtable[++count])->count) { options->filter_seed_mappings = *(((struct arg_int*)argtable[count])->ival); }
  if (((struct arg_int*)argtable[++count])->count) { 
    options->min_cal_size = *(((struct arg_int*)argtable[count])->ival); 
    options->set_cal = 1;
  }
  if (((struct arg_int*)argtable[++count])->count) { options->num_cpu_threads = *(((struct arg_int*)argtable[count])->ival); }
  if (((struct arg_int*)argtable[++count])->count) { options->batch_size = *(((struct arg_int*)argtable[count])->ival); }
  if (((struct arg_dbl*)argtable[++count])->count) { options->match = *((struct arg_dbl*)argtable[count])->dval; }
  if (((struct arg_dbl*)argtable[++count])->count) { options->mismatch = *(((struct arg_dbl*)argtable[count])->dval); }
  if (((struct arg_dbl*)argtable[++count])->count) { options->gap_open = *(((struct arg_dbl*)argtable[count])->dval); }
  if (((struct arg_dbl*)argtable[++count])->count) { options->gap_extend = *(((struct arg_dbl*)argtable[count])->dval); }
  if (((struct arg_int*)argtable[++count])->count) { options->min_score = *(((struct arg_int*)argtable[count])->ival); }
  //  if (((struct arg_int*)argtable[++count])->count) { options->pair_mode = *(((struct arg_int*)argtable[count])->ival); }
  if (((struct arg_int*)argtable[++count])->count) { options->pair_min_distance = *(((struct arg_int*)argtable[count])->ival); }
  if (((struct arg_int*)argtable[++count])->count) { options->pair_max_distance = *(((struct arg_int*)argtable[count])->ival); }
  if (((struct arg_int*)argtable[++count])->count) { options->report_best = (((struct arg_int*)argtable[count])->count); }
  if (((struct arg_file*)argtable[++count])->count) { options->report_all = (((struct arg_int *)argtable[count])->count); }
  if (((struct arg_int*)argtable[++count])->count) { options->report_n_best = *(((struct arg_int*)argtable[count])->ival); }
  if (((struct arg_int*)argtable[++count])->count) { options->report_n_hits = *(((struct arg_int*)argtable[count])->ival); }
  if (((struct arg_int*)argtable[++count])->count) { options->report_only_paired = (((struct arg_int*)argtable[count])->count); }
  if (((struct arg_str*)argtable[++count])->count) { options->prefix_name = strdup(*(((struct arg_str*)argtable[count])->sval)); }
  if (((struct arg_file*)argtable[++count])->count) { options->log_level = *(((struct arg_int*)argtable[count])->ival); }
  if (((struct arg_int*)argtable[++count])->count) { options->help = ((struct arg_int*)argtable[count])->count; }

  options->bam_format = 0; 
  options->set_bam_format = 0;
  options->output_format = SAM_FORMAT;
  if (((struct arg_int*)argtable[++count])->count) {
    char *format = (char *) (*((struct arg_str*)argtable[count])->sval);
    if (!strcmp(format, "bam") || !strcmp(format, "BAM")) {
      options->output_format = BAM_FORMAT;
      options->bam_format = 1; 
      options->set_bam_format = 1;
    }
  }

  if (((struct arg_int*)argtable[++count])->count) { options->realignment = ((struct arg_int*)argtable[count])->count; }
  if (((struct arg_int*)argtable[++count])->count) { options->recalibration = ((struct arg_int*)argtable[count])->count; }
  if (((struct arg_str*)argtable[++count])->count) { options->adapter = strdup(*(((struct arg_str*)argtable[count])->sval)); }

  if (options->adapter) options->adapter_length = strlen(options->adapter);

  if (((struct arg_int*)argtable[++count])->count) {
    char *format = (char *) (*((struct arg_str*)argtable[count])->sval);
    if (!strcmp(format, "sam") || !strcmp(format, "SAM")) {
      options->input_format = SAM_FORMAT;
    } else if (!strcmp(format, "bam") || !strcmp(format, "BAM")) {
      options->input_format = BAM_FORMAT;
    } else {
      options->input_format = FASTQ_FORMAT;
    }
  }

  if (((struct arg_int*)argtable[++count])->count) { options->version = ((struct arg_int*)argtable[count])->count; }

  if (options->mode == DNA_MODE) {
    if (((struct arg_int*)argtable[++count])->count) { options->num_seeds = *(((struct arg_int*)argtable[count])->ival); }
  } else if (options->mode == RNA_MODE) {
    if (((struct arg_int*)argtable[++count])->count) { options->seeds_max_distance = *(((struct arg_int*)argtable[count])->ival); }
    if (((struct arg_file*)argtable[++count])->count) { options->transcriptome_filename = strdup(*(((struct arg_file*)argtable[count])->filename)); }
    if (((struct arg_int*)argtable[++count])->count) { options->seed_size = *(((struct arg_int*)argtable[count])->ival); }
    if (((struct arg_int*)argtable[++count])->count) { options->max_intron_length = *(((struct arg_int*)argtable[count])->ival); }
    if (((struct arg_int*)argtable[++count])->count) { options->min_intron_length = *(((struct arg_int*)argtable[count])->ival); }

    if (((struct arg_file*)argtable[++count])->count) { options->fast_mode = (((struct arg_int *)argtable[count])->count); }

    /*
    if (!set_bam_format) {
      if (options->fast_mode) {
	options->bam_format = 0;
      } else {
	options->bam_format = 1;
      }
    }
    */
  }

  return options;
}

//--------------------------------------------------------------------

options_t *parse_options(int argc, char **argv) {
  //	struct arg_end *end = arg_end(10);
  //	void **argtable = argtable_options_get(argtable_options, end);

  int mode, num_options = NUM_OPTIONS;
  if (strcmp(argv[0], "dna") == 0) {
    mode = DNA_MODE;
    num_options += NUM_DNA_OPTIONS;
  } else if (strcmp(argv[0], "rna") == 0) {
    mode = RNA_MODE;
    num_options += NUM_RNA_OPTIONS;
  } else {
    LOG_FATAL("Command unknown.\nValid commands are:\n\tdna: to map DNA sequences\n\trna: to map RNA sequences\n\tbs: to map BS sequences\n\tbuild-index: to create the genome index.\nUse -h or --help to display options.\n");
  }

  void **argtable = argtable_options_new(mode);
  options_t *options = options_new();
  options->mode = mode;

  if (argc < 2) {
    display_help(argtable, mode);
    exit(EXIT_FAILURE);
  } else {
    int num_errors = arg_parse(argc, argv, argtable);   
    if (num_errors > 0) {
      fprintf(stdout, "\nError:\n");
      // struct end is always allocated in the last position
      arg_print_errors(stdout, argtable[num_options], "\t");
      display_help(argtable, mode);
      exit(EXIT_FAILURE);
    } else {
      options = read_CLI_options(argtable, options);
      if (options->help) {
	display_help(argtable, mode);
	argtable_options_free(argtable, num_options);
	options_free(options);
	exit(EXIT_SUCCESS);
      }
      // Check if 'help' option has been provided.
    }
    
  }
  argtable_options_free(argtable, num_options);

  options->mode = mode;
  return options;
}

//--------------------------------------------------------------------

void usage(void **argtable, int mode) {
  printf("\nUsage:\n\t%s %s <options>\n", HPG_ALIGNER_BIN, (mode == RNA_MODE ? "rna" : "dna"));
  //  arg_print_syntaxv(stdout, argtable, "\n");
  printf("\nOptions:\n");
  arg_print_glossary(stdout, argtable, "\t%-50s\t%s\n");
}

//--------------------------------------------------------------------

void usage_cli(int mode) {
  void **argtable = argtable_options_new(mode);
  display_help(argtable, mode);
}

//--------------------------------------------------------------------

void display_version() {
  printf("%s version %s\n", HPG_ALIGNER_NAME, HPG_ALIGNER_VERSION);
  printf("\n");
  printf("Source code at https://github.com/opencb/hpg-aligner\n");
  printf("Documentation at https://github.com/opencb/hpg-aligner/wiki/\n");
}

//--------------------------------------------------------------------

void display_help(void **argtable, int mode) {
  if (mode == DNA_MODE) {
    display_dna_help();
    exit(EXIT_SUCCESS);
  } else if (mode == RNA_MODE) {
    //    usage(argtable, mode);
    display_rna_help();
    exit(EXIT_SUCCESS);
  } else {
    printf("Error: unknown command\n");
    exit(EXIT_FAILURE);
  }
}

//--------------------------------------------------------------------

void display_options(options_t *options, FILE *file) {
  if (options->mode == DNA_MODE) {
    display_dna_options(options, file);
  } else if (options->mode == RNA_MODE) {
    display_rna_options(options, file);
    //    if (!file) {
    //      options_display(options);
    //    } else {
    //      options_to_file(options, file);
    //    }
  } else {
    printf("Error: unknown command\n");
    exit(EXIT_FAILURE);
  }
}

//--------------------------------------------------------------------

char* create_cmdline(int argc, char **argv) {
  char *cmdline = NULL;

  int size = 0;
  for (int i = 0; i < argc; i++) {
    size += strlen(argv[i]);
  }

  if (size) {
    size += 100;
    
    cmdline = (char *) calloc(size, sizeof(char));
    sprintf(cmdline, "%s ", HPG_ALIGNER_BIN);
    for (int i = 0; i < argc; i++) {
      sprintf(cmdline, "%s %s", cmdline, argv[i]);
    }
  }
  
  return cmdline;
}

//--------------------------------------------------------------------
// DNA
//--------------------------------------------------------------------

void display_dna_help() {
  printf("\n");
  printf("+===============================================================+\n");
  printf("|                HPG-Aligner help for DNA mapping               |\n");
  printf("+===============================================================+\n");
  printf("Usage:\n");
  printf("\t%s dna -f=<input filename> -i=<index dirname> [options]\n", HPG_ALIGNER_BIN);
  printf("\n");
  printf("Mandatory parameters:\n");
  printf("\t-f, --fq, --fastq=<file>         Input filename. For more than one file: f1.fq,f2.fq,...\n");
  printf("\t-i, --index=<file>               Index directory name\n");
  printf("\n");

  printf("General options:\n");
  printf("\t-j, --fq2, --fastq2=<file>        Input filename for mate #2 in paired-end mode\n");
  printf("\t-o, --outdir=<file>               Output directory name\n");
  printf("\t--prefix=<string>                 Prefix for the output filename\n");
  printf("\t-z, --gzip                        FastQ input files are gzipped\n");
  printf("\t--input-format=<string>           Input file format (fastq or bam) [fastq]\n");
  printf("\t--output-format=<string>          Output file format (sam or bam) [sam]\n");
  printf("\t-v,--version                      Display the current HPG-Aligner version\n");
  printf("\t-h,--help                         Display this help\n");
  printf("\n");

  printf("Architecture options:\n");
  printf("\t-t, --cpu-threads=<int>            Number of CPU threads [%lu]\n", get_optimal_cpu_num_threads());
  printf("\t--read-batch-size=<int>            Batch size in bytes [%i]\n", DEFAULT_DNA_READ_BATCH_SIZE);
  printf("\n");

  printf("Paired-end:\n");
  printf("\t--paired-min-distance=<int>        Minimum distance between pairs [%i]\n", DEFAULT_PAIR_MIN_DISTANCE);
  printf("\t--paired-max-distance=<int>        Maximum distance between pairs [%i]\n", DEFAULT_PAIR_MAX_DISTANCE);
  printf("\n");

  printf("Seeding and CAL options:\n");
  printf("\t--num-seeds=<int>                  Number of seeds [%i]\n", DEFAULT_DNA_NUM_SEEDS);
  printf("\t--min-cal-size=<int>               Minimum CAL size [%i]\n", DEFAULT_DNA_MIN_CAL_SIZE);
  //    --filter-read-mappings=<int>                      Reads that map in more than <n> locations are discarded
  //    --filter-seed-mappings=<int>                      Seeds that map in more than <n> locations are discarded
  printf("\n");

  printf("Smith-Waterman options:\n");
  printf("\t--sw-match=<double>                Match score for Smith-Waterman algorithm [%0.2f]\n", DEFAULT_SW_MATCH);
  printf("\t--sw-mismatch=<double>             Mismatch penalty for Smith-Waterman algorithm [%0.2f]\n", DEFAULT_SW_MISMATCH);
  printf("\t--sw-gap-open=<double>             Gap open penalty for Smith-Waterman algorithm [%0.2f]\n", DEFAULT_SW_GAP_OPEN);
  printf("\t--sw-gap-extend=<double>           Gap extend penalty for Smith-Waterman algorithm [%0.2f]\n", DEFAULT_SW_GAP_EXTEND);
  printf("\t--min-score=<int>                  Minimum Smith-Waterman score for valid mappings (normalized from 0 to 100) [%i]\n", DEFAULT_MIN_SCORE);
  printf("\n");

  printf("Report options:\n");
  printf("\t--report-best                      Report all alignments with best score\n");
  printf("\t--report-all                       Report all alignments\n");
  printf("\t--report-n-best=<int>              Report the <n> best alignments\n");
  printf("\t--report-n-hits=<int>              Report <n> hits\n");
  printf("\t--report-only-paired               Report only the paired reads\n");
  printf("\n");

  printf("Pre-processing options:\n");
  printf("\t-a,--adapter=<string>              Adapter sequence in the read to remove before mapping\n");
  printf("\n");

  printf("Post-processing options:\n");
  printf("\t--indel-realignment                Perform the indel-based realignment after mapping\n");
  printf("\t--recalibration                    Perform the base-quality recalibration after mapping\n");  
  printf("+===============================================================+\n");
}

//--------------------------------------------------------------------

void display_dna_options(options_t *options, FILE *f) {
  FILE *file = (f == NULL ? stdout : f);

  fprintf(file, "\n");
  fprintf(file, "+===============================================================+\n");
  fprintf(file, "|             HPG ALIGNER PARAMETERS CONFIGURATION              |\n");
  fprintf(file, "+===============================================================+\n");
  fprintf(file, "HPG Aligner version: %s\n", HPG_ALIGNER_VERSION);
  fprintf(file, "Command line       : %s\n", options->cmdline);
  fprintf(file, "\n");
  fprintf(file, "General parameters:\n");
  fprintf(file, "\tType of sequence: %s\n", options->str_mode);
  if (options->in_filename2) {
    fprintf(file, "\tInput FastQ filename #1: %s\n", options->in_filename);
    fprintf(file, "\tInput FastQ filename #2: %s\n", options->in_filename2);
    fprintf(file, "\tPaired-end:\n");
    fprintf(file, "\t\tMin. distance: %d\n", options->pair_min_distance);
    fprintf(file, "\t\tMax. distance: %d\n", options->pair_max_distance);
  } else {
    fprintf(file, "\tInput %s filename: %s\n", FORMAT_STR(options->input_format), options->in_filename);
    if (options->input_format == FASTQ_FORMAT) {
      fprintf(file, "\tSingle-end\n");
    }
  }
  printf("\n");
  if (options->input_format == FASTQ_FORMAT) {
    fprintf(file, "\tFastQ gzip mode    : %s\n", options->gzip == 1 ? "Enabled" : "Disabled");
  }
  fprintf(file, "\tIndex directory name : %s\n", options->bwt_dirname);
  fprintf(file, "\tOutput directory name: %s\n", options->output_name);
  if (options->prefix_name) {
    fprintf(file, "\tOutput prefix      : %s\n", options->prefix_name);
  }
  
  fprintf(file, "\tOutput file format: %s\n", 
	 (options->bam_format || options->realignment || options->recalibration) ? "BAM" : "SAM");
  fprintf(file, "\n");
  
  fprintf(file, "Architecture parameters:\n");
  fprintf(file, "\tNumber of cpu threads: %d\n",  options->num_cpu_threads);
  fprintf(file, "\tBatch size (in bytes): %d\n",  options->batch_size);
  fprintf(file, "\n");

  fprintf(file, "Seeding and CAL parameters:\n");
  fprintf(file, "\tNum. seeds  : %d\n",  options->num_seeds);
  fprintf(file, "\tMin CAL size: %d\n",  options->min_cal_size);
  fprintf(file, "\n");
  
  //fprintf(file, "Mapping filters\n");
  //fprintf(file, "\tFor reads: %d mappings maximum, otherwise discarded\n", options->filter_read_mappings);
  //fprintf(file, "\tFor seeds: %d mappings maximum, otherwise discarded\n", options->filter_seed_mappings);
  //fprintf(file, "\n");
    
  fprintf(file, "Smith-Waterman parameters:\n");
  fprintf(file, "\tMatch score       : %0.4f\n", options->match);
  fprintf(file, "\tMismatch penalty  : %0.4f\n", options->mismatch);
  fprintf(file, "\tGap open penalty  : %0.4f\n", options->gap_open);
  fprintf(file, "\tGap extend penalty: %0.4f\n", options->gap_extend);
  fprintf(file, "\n");
  
  fprintf(file, "Report parameters:\n");
  fprintf(file, "\tReport all hits      : %s\n",  options->report_all == 0 ? "Disabled" : "Enabled");
  fprintf(file, "\tReport n best hits   : %d\n",  options->report_n_best);
  fprintf(file, "\tReport n hits        : %d\n",  options->report_n_hits);
  fprintf(file, "\tReport best hits     : %s\n",  options->report_best == 0 ? "Disabled" : "Enabled");
  fprintf(file, "\tReport unpaired reads: %s\n",  options->report_only_paired == 0 ? "Enabled" : "Disabled");
  fprintf(file, "\n");
  
  fprintf(file, "Pre-processing:\n");
  fprintf(file, "\tAdapter sequence: %s\n",
	  (options->adapter ? options->adapter : "Not present"));
  fprintf(file, "\n");
  
  fprintf(file, "Post-processing:\n");
  fprintf(file, "\tIndel realignment         : %s\n",
	  (options->realignment ? "Enabled" : "Disabled"));
  fprintf(file, "\tBase-quality recalibration: %s\n",
	  (options->recalibration ? "Enabled" : "Disabled"));
  fprintf(file, "+===============================================================+\n");
}

//--------------------------------------------------------------------
// RNA
//--------------------------------------------------------------------

void display_rna_help() {
  printf("\n");
  printf("+===============================================================+\n");
  printf("|                HPG-Aligner help for RNA mapping               |\n");
  printf("+===============================================================+\n");
  printf("Usage:\n");
  printf("\t%s rna -f=<input filename> -i=<index dirname> [options]\n", HPG_ALIGNER_BIN);
  printf("\n");
  printf("Mandatory parameters:\n");
  printf("\t-f, --fq, --fastq=<file>         Input filename. For more than one file: f1.fq,f2.fq,...\n");
  printf("\t-i, --index=<file>               Index directory name\n");
  printf("\n");

  printf("General options:\n");
  printf("\t-j, --fq2, --fastq2=<file>        Input filename for mate #2 in paired-end mode\n");
  printf("\t-o, --outdir=<file>               Output directory name\n");
  printf("\t--prefix=<string>                 Prefix for the output filename\n");
  printf("\t-z, --gzip                        FastQ input files are gzipped\n");
  printf("\t--output-format=<string>          Output file format (sam or bam) [sam]. Only available for SA index (for BWT index, output format is always bam)\n");
  printf("\t-l,--log-level                    Set log debug level\n");
  printf("\t-v,--version                      Display the current HPG-Aligner version\n");
  printf("\t-h,--help                         Display this help\n");
  printf("\n");

  printf("Architecture options:\n");
  printf("\t-t, --cpu-threads=<int>            Number of CPU threads [%lu]\n", get_optimal_cpu_num_threads());
  printf("\t--read-batch-size=<int>            Batch size in bytes [%i]\n", DEFAULT_DNA_READ_BATCH_SIZE);
  printf("\n");

  printf("Paired-end:\n");
  printf("\t--paired-min-distance=<int>        Minimum distance between pairs [%i]\n", DEFAULT_PAIR_MIN_DISTANCE);
  printf("\t--paired-max-distance=<int>        Maximum distance between pairs [%i]\n", DEFAULT_PAIR_MAX_DISTANCE);
  printf("\n");

  printf("Seeding and CAL options:\n");
  printf("\t--num-seeds=<int>                  Number of seeds [%i]\n", DEFAULT_DNA_NUM_SEEDS);
  printf("\t--min-cal-size=<int>               Minimum CAL size [%i]\n", DEFAULT_DNA_MIN_CAL_SIZE);
  printf("\t--seed-size=<int>                  Seed size [%i]. Only available for BWT index (for SA index, seed size is always 18)\n", DEFAULT_RNA_SEED_SIZE);
  printf("\n");

  printf("Filtering options:\n");
  printf("\t--filter-read-mappings=<int>       Reads that map in more than <n> locations are discarded [%i]\n", DEFAULT_FILTER_READ_MAPPINGS);
  printf("\t--filter-seed-mappings=<int>       Seeds that map in more than <n> locations are discarded [%i]\n", DEFAULT_FILTER_SEED_MAPPINGS);
  printf("\n");

  printf("Intron and transcriptome options:\n");
  printf("\t--transcriptome-file=<file>        Transcriptome file to help searching splice junctions\n");
  printf("\t--max-intron-size=<int>            Maximum intron size [%i]\n", DEFAULT_MAX_INTRON_LENGTH);
  printf("\t--min-intron-size=<int>            Maximum intron size [%i]\n", DEFAULT_MIN_INTRON_LENGTH);
  printf("\n");

  printf("Smith-Waterman options:\n");
  printf("\t--sw-match=<double>                Match score for Smith-Waterman algorithm [%0.2f]\n", DEFAULT_SW_MATCH);
  printf("\t--sw-mismatch=<double>             Mismatch penalty for Smith-Waterman algorithm [%0.2f]\n", DEFAULT_SW_MISMATCH);
  printf("\t--sw-gap-open=<double>             Gap open penalty for Smith-Waterman algorithm [%0.2f]\n", DEFAULT_SW_GAP_OPEN);
  printf("\t--sw-gap-extend=<double>           Gap extend penalty for Smith-Waterman algorithm [%0.2f]\n", DEFAULT_SW_GAP_EXTEND);
  printf("\t--min-score=<int>                  Minimum Smith-Waterman score for valid mappings (normalized from 0 to 100) [%i]\n", DEFAULT_MIN_SCORE);
  printf("\n");

  printf("Report options:\n");
  printf("\t--report-best                      Report all alignments with best score\n");
  printf("\t--report-all                       Report all alignments\n");
  printf("\t--report-n-best=<int>              Report the <n> best alignments\n");
  printf("\t--report-n-hits=<int>              Report <n> hits\n");
  printf("\t--report-only-paired               Report only the paired reads\n");
  printf("+===============================================================+\n");
}

//--------------------------------------------------------------------

void display_rna_options(options_t *options, FILE *f) {
  FILE *file = (f == NULL ? stdout : f);

  fprintf(file, "\n");
  fprintf(file, "+===============================================================+\n");
  fprintf(file, "|             HPG ALIGNER PARAMETERS CONFIGURATION              |\n");
  fprintf(file, "+===============================================================+\n");
  fprintf(file, "HPG Aligner version: %s\n", HPG_ALIGNER_VERSION);
  fprintf(file, "Command line       : %s\n", options->cmdline);
  fprintf(file, "\n");
  fprintf(file, "General parameters:\n");
  fprintf(file, "\tType of sequence: %s\n", options->str_mode);
  if (options->in_filename2) {
    fprintf(file, "\tInput FastQ filename #1: %s\n", options->in_filename);
    fprintf(file, "\tInput FastQ filename #2: %s\n", options->in_filename2);
    fprintf(file, "\tPaired-end:\n");
    fprintf(file, "\t\tMin. distance: %d\n", options->pair_min_distance);
    fprintf(file, "\t\tMax. distance: %d\n", options->pair_max_distance);
  } else {
    fprintf(file, "\tInput %s filename: %s\n", FORMAT_STR(options->input_format), options->in_filename);
    fprintf(file, "\tSingle-end\n");
  }
  printf("\n");
  if (options->input_format == FASTQ_FORMAT) {
    fprintf(file, "\tFastQ gzip mode    : %s\n", options->gzip == 1 ? "Enabled" : "Disabled");
  }
  fprintf(file, "\tIndex directory name : %s\n", options->bwt_dirname);
  fprintf(file, "\tOutput directory name: %s\n", options->output_name);
  if (options->prefix_name) {
    fprintf(file, "\tOutput prefix      : %s\n", options->prefix_name);
  }
  
  fprintf(file, "\tOutput file format: %s\n", 
	 (options->bam_format || options->realignment || options->recalibration) ? "BAM" : "SAM");
  fprintf(file, "\n");
  
  fprintf(file, "Architecture parameters:\n");
  fprintf(file, "\tNumber of cpu threads: %d\n",  options->num_cpu_threads);
  fprintf(file, "\tBatch size (in bytes): %d\n",  options->batch_size);
  fprintf(file, "\n");

  fprintf(file, "Seeding and CAL parameters:\n");
  fprintf(file, "\tNum. seeds  : %d\n",  options->num_seeds);
  fprintf(file, "\tMin CAL size: %d\n",  options->min_cal_size);
  fprintf(file, "\n");
  
  fprintf(file, "Filtering options:\n");
  fprintf(file, "\tFor reads: %d mappings maximum, otherwise discarded\n", options->filter_read_mappings);
  fprintf(file, "\tFor seeds: %d mappings maximum, otherwise discarded\n", options->filter_seed_mappings);
  fprintf(file, "\n");
    
  fprintf(file, "Intron and transcriptome options:\n");
  fprintf(file, "\tTranscriptome filename: %s\n", 
	  (options->transcriptome_filename ? options->transcriptome_filename : "Not present"));
  fprintf(file, "\tMinimum intron size   : %i\n", options->min_intron_length);
  fprintf(file, "\tMaximum intron size   : %i\n", options->max_intron_length);
  fprintf(file, "\n");

  fprintf(file, "Smith-Waterman parameters:\n");
  fprintf(file, "\tMatch score       : %0.2f\n", options->match);
  fprintf(file, "\tMismatch penalty  : %0.2f\n", options->mismatch);
  fprintf(file, "\tGap open penalty  : %0.2f\n", options->gap_open);
  fprintf(file, "\tGap extend penalty: %0.2f\n", options->gap_extend);
  fprintf(file, "\n");
  
  fprintf(file, "Report parameters:\n");
  fprintf(file, "\tReport all hits      : %s\n",  options->report_all == 0 ? "Disabled" : "Enabled");
  fprintf(file, "\tReport n best hits   : %d\n",  options->report_n_best);
  fprintf(file, "\tReport n hits        : %d\n",  options->report_n_hits);
  fprintf(file, "\tReport best hits     : %s\n",  options->report_best == 0 ? "Disabled" : "Enabled");
  fprintf(file, "\tReport unpaired reads: %s\n",  options->report_only_paired == 0 ? "Enabled" : "Disabled");
  fprintf(file, "+===============================================================+\n");
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
