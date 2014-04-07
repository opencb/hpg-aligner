#include "rna_aligner.h"

#define NUM_SECTIONS_TIME 		8

//--------------------------------------------------------------------
// workflow input                                                                                                                                                    
//--------------------------------------------------------------------  
extern int splice_junction_type(char nt_start_1, char nt_start_2, char nt_end_1, char nt_end_2);



int w1_function(void *data) {
  batch_t *batch = (batch_t *) data;
  apply_bwt_rna(batch->bwt_input, batch);
  apply_caling_rna(batch->cal_input, batch);
  apply_sw_rna(batch->sw_input, batch);

  return CONSUMER_STAGE;
}

int w2_function(void *data) {
  batch_t *batch = (batch_t *) data;
  apply_rna_last(batch->sw_input, batch);
  
  return CONSUMER_STAGE;
}
 
int w3_function(void *data) {
  batch_t *batch = (batch_t *) data;
  apply_rna_last_hc(batch->sw_input, batch);

  return CONSUMER_STAGE;
}




void rna_aligner(options_t *options) {
  int path_length = strlen(options->output_name);
  int prefix_length = 0;
  
  if (options->prefix_name) {
    prefix_length = strlen(options->prefix_name);
  }

  char *reads_results = (char *)calloc((60 + prefix_length), sizeof(char));
  char *log_results = (char *)calloc((60 + prefix_length), sizeof(char));
  char *extend_junctions = (char *)calloc((60 + prefix_length), sizeof(char));
  char *exact_junctions = (char *)calloc((60 + prefix_length), sizeof(char));

  char *output_filename = (char *)calloc((path_length + prefix_length + 60), sizeof(char));
  char *log_filename = (char *)calloc((path_length + prefix_length + 60), sizeof(char));
  char *extend_filename = (char *)calloc((path_length + prefix_length + 60), sizeof(char));
  char *exact_filename = (char *)calloc((path_length + prefix_length + 60), sizeof(char));

  if (options->prefix_name) {
    strcat(reads_results, "/");
    strcat(reads_results, options->prefix_name);
    if (!options->fast_mode) {
      strcat(reads_results, "_alignments.bam");  
    } else {
      strcat(reads_results, "_alignments.sam");  
    }
    strcat(log_results, "/");
    strcat(log_results, options->prefix_name);
    strcat(log_results, "_hpg-aligner.log");  

    strcat(extend_junctions, "/");
    strcat(extend_junctions, options->prefix_name);
    strcat(extend_junctions, "_extend_junctions.bed");

    strcat(exact_junctions, "/");
    strcat(exact_junctions, options->prefix_name);
    strcat(exact_junctions, "_exact_junctions.bed");
 
  } else {
    if (!options->fast_mode) {
      strcat(reads_results, "/alignments.bam");
    } else {
      strcat(reads_results, "/alignments.sam");
    }
    strcat(log_results, "/hpg-aligner.log");
    strcat(extend_junctions, "/extend_junctions.bed");
    strcat(exact_junctions, "/exact_junctions.bed");
  } 

  strcat(output_filename, options->output_name);
  strcat(output_filename, reads_results);
  free(reads_results);

  strcat(log_filename, options->output_name);
  strcat(log_filename, log_results);
  free(log_results);

  strcat(extend_filename, options->output_name);
  strcat(extend_filename, extend_junctions);
  free(extend_junctions);

  strcat(exact_filename, options->output_name);
  strcat(exact_filename, exact_junctions);
  free(exact_junctions);

  FILE *fd_log;
  fd_log = fopen(log_filename, "w");

  LOG_DEBUG("Auto Thread Configuration Done !");

  // timing
  //if (time_on) { 

  char* labels_time[NUM_SECTIONS_TIME] = {"FASTQ Reader               ", 
					  "BWT Server                 ", 
					  "REGION Seeker              ", 
					  "CAL Seeker                 ", 
					  "RNA Preprocess             ", 
					  "RNA Server                 ",
					  "BAM Writer                 ", 
					  "TOTAL Time                 "};    

  //timing = timing_new((char**) labels_time, NUM_SECTIONS_TIME);
  //}

  //validate_options(options);
  // display selected options
  LOG_DEBUG("Displaying options...\n");
  options_display(options);

  //time_on =  (unsigned int) options->timming;
  //statistics_on =  (unsigned int) options->statistics;

  struct timeval time_genome_s, time_genome_e;
  double time_genome;

  metaexons_t *metaexons;
  genome_t *genome;
  
  bwt_index_t *bwt_index;
  sa_index3_t *sa_index;
  int num_chromosomes;
  // load index for dna/rna or for bisulfite case

  start_timer(time_genome_s);

  // genome parameters 
  if (!options->fast_mode) {
    //////////////// LOAD BWT INDEX //////////////////////    
    // BWT index
    LOG_DEBUG("Reading bwt index...");
    bwt_index = bwt_index_new(options->bwt_dirname, false);
    LOG_DEBUG("Reading bwt index done !!");
    
    LOG_DEBUG("Reading genome...");
    genome = genome_new("dna_compression.bin", options->bwt_dirname, BWT_MODE);  
    LOG_DEBUG("Done !!");
    //////////////////////////////////////////////////////
  } else {    
    ///////////////// LOAD SA INDEX ////////////////////// 
    LOG_DEBUG("Loading SA tables...");
    sa_index = sa_index3_new(options->bwt_dirname);
    sa_index3_display(sa_index);

    LOG_DEBUG("Reading genome...");
    genome = genome_new("dna_compression.bin", options->bwt_dirname, SA_MODE);

    genome->num_chromosomes = sa_index->genome->num_chroms;
    genome->chr_name = (char **) calloc(genome->num_chromosomes, sizeof(char *));
    genome->chr_size = (size_t *) calloc(genome->num_chromosomes, sizeof(size_t));
    genome->chr_offset = (size_t *) calloc(genome->num_chromosomes, sizeof(size_t));
    size_t offset = 0;

    for (int c = 0; c < genome->num_chromosomes; c++) {
      genome->chr_size[c] = sa_index->genome->chrom_lengths[c];
      genome->chr_name[c] = strdup(sa_index->genome->chrom_names[c]);
      genome->chr_offset[c] = offset;
      offset += genome->chr_size[c];
    }

    LOG_DEBUG("Done !!");
    //metaexons = metaexons_new(sa_index->genome->num_chroms, 
    //			      sa_index->genome->chrom_lengths);    
    //////////////////////////////////////////////////////
  }

  num_chromosomes = genome->num_chromosomes;

  // Metaexons structure
  metaexons = metaexons_new(genome->num_chromosomes, 
			    genome->chr_size);
  
  stop_timer(time_genome_s, time_genome_e, time_genome);

  //============================= INPUT INITIALIZATIONS =========================//  
  //BWT parameters
  
  
  bwt_optarg_t *bwt_optarg = bwt_optarg_new(1, 0,
					    //500,
					    options->filter_read_mappings, 
					    options->filter_seed_mappings);
  
  // CAL parameters
  //printf("%i\n", options->min_cal_size);
  cal_optarg_t *cal_optarg = cal_optarg_new(options->min_cal_size, 
					    options->seeds_max_distance, 
					    options->num_seeds, 
					    options->min_num_seeds_in_cal,
					    options->seed_size, 
					    options->min_seed_size, 
					    options->cal_seeker_errors, 
					    options->max_intron_length, 
					    options->min_intron_length);
  
  // paired mode parameters
  pair_mng_t *pair_mng = pair_mng_new(options->pair_mode, options->pair_min_distance, 
				      options->pair_max_distance, options->report_only_paired);
  
  // report parameters
  report_optarg_t *report_optarg = report_optarg_new(options->report_all,
						     options->report_n_best,
						     options->report_n_hits, 
						     options->report_only_paired,
						     options->report_best);  

  avls_list_t* avls_list = avls_list_new(num_chromosomes);

  if (options->transcriptome_filename != NULL) {
    printf("Loading transcriptome...\n");
    load_transcriptome(options->transcriptome_filename, genome, avls_list, metaexons);
    printf("Load done!\n");
  }

  linked_list_t *buffer    = linked_list_new(COLLECTION_MODE_SYNCHRONIZED);
  linked_list_t *buffer_hc = linked_list_new(COLLECTION_MODE_SYNCHRONIZED);

  FILE *f_sa, *f_hc;

  f_sa = fopen("buffer_sa.tmp", "w+b");
  if (f_sa == NULL) {
    LOG_FATAL("Error opening file 'buffer_sa.tmp' \n");
  }

  f_hc = fopen("buffer_hc.tmp", "w+b");
  if (f_hc == NULL) {
    LOG_FATAL("Error opening file 'buffer_hc.tmp' \n");
  }

  fastq_batch_reader_input_t reader_input;
  fastq_batch_reader_input_init(options->in_filename, options->in_filename2, 
				options->pair_mode, options->batch_size, 
				NULL, options->gzip, &reader_input);  

  //buffer_reader_input_t buffer_reader_input;
  //buffer_reader_input_init(&reader_input,
  //buffer,
  //buffer_reader_input);

  //if (options->pair_mode == SINGLE_END_MODE) {
  //reader_input.fq_file1 = fastq_fopen(options->in_filename);
  //} else {
  //reader_input.fq_file1 = fastq_fopen(options->in_filename);
  //reader_input.fq_file2 = fastq_fopen(options->in_filename2);
  //}

  linked_list_t *alignments_list = linked_list_new(COLLECTION_MODE_SYNCHRONIZED);
  linked_list_set_flag(options->pair_mode, alignments_list);

  sw_optarg_t sw_optarg;
  sw_optarg_init(options->gap_open, options->gap_extend, 
		 options->match, options->mismatch, &sw_optarg);

  bwt_server_input_t bwt_input;
  bwt_server_input_init(NULL,  0,  bwt_optarg, 
			bwt_index, NULL,  0, 
			NULL, metaexons, &sw_optarg, 
			genome, &bwt_input);
      
  region_seeker_input_t region_input;
  region_seeker_input_init(NULL, cal_optarg, 
			   bwt_optarg, bwt_index, NULL, 
			   0, 0, options->min_seed_padding_left, 
			   options->min_seed_padding_right, genome, metaexons, &region_input);
  
  cal_seeker_input_t cal_input;
  cal_seeker_input_init(NULL, cal_optarg, NULL, 0, NULL, NULL, genome, 
			bwt_optarg, bwt_index, metaexons, &cal_input);
  
  preprocess_rna_input_t preprocess_rna;  
  preprocess_rna_input_init(options->max_intron_length, options->flank_length, 
			    options->seeds_max_distance, options->seed_size, genome, 
			    &preprocess_rna);

  int pair_mode = pair_mng->pair_mode;
  sw_server_input_t sw_input;
  sw_server_input_init(NULL, NULL, 0,  options->match,  
		       options->mismatch,  options->gap_open, options->gap_extend,  
		       options->min_score,  options->flank_length, genome,  
		       options->max_intron_length, options->min_intron_length,  
		       options->seeds_max_distance,  bwt_optarg, avls_list, 
		       cal_optarg, bwt_index, metaexons, buffer, buffer_hc, 
		       f_sa, f_hc, pair_mode, &sw_input);
  

  pair_server_input_t pair_input;
  pair_server_input_init(pair_mng, report_optarg, NULL, NULL, NULL, &pair_input);

  batch_writer_input_t writer_input;
  batch_writer_input_init( output_filename,
			   exact_filename, 
			   extend_filename, 
			   alignments_list, 
			   genome, &writer_input);

  if (!options->fast_mode) {
    bam_header_t *bam_header = create_bam_header_by_genome(genome);
    //printf("%s\n", output_filename);
    writer_input.bam_file = bam_fopen_mode(output_filename, bam_header, "w");
    bam_fwrite_header(bam_header, writer_input.bam_file);
  } else {
    writer_input.bam_file = (bam_file_t *) fopen(output_filename, "w"); 
    write_sam_header(sa_index->genome, (FILE *) writer_input.bam_file);
  }

  extra_stage_t extra_stage_input;
  
  int workflow_enable = 1;

  batch_t *batch = batch_new(&bwt_input, &region_input, &cal_input, 
			     &pair_input, &preprocess_rna, &sw_input, &writer_input, RNA_MODE, NULL);

  //  fastq_batch_reader_input_t reader_input;

  //Parse input for more than one file
  char *fq_list1 = options->in_filename, *fq_list2 = options->in_filename2;
  char token[2] = ",";
  char *ptr;

  array_list_t *files_fq1 = array_list_new(50,
					   1.25f,
					   COLLECTION_MODE_ASYNCHRONIZED);
  array_list_t *files_fq2 = array_list_new(50,
					   1.25f,
					   COLLECTION_MODE_ASYNCHRONIZED);
  int num_files1 = 0;
  int num_files2 = 0;

  if (fq_list1) {
    ptr = strtok(fq_list1, token);    // Primera llamada => Primer token
    array_list_insert(strdup(ptr), files_fq1);
    while( (ptr = strtok( NULL, token)) != NULL ) {    // Posteriores llamadas
      array_list_insert(strdup(ptr), files_fq1);
    }
    num_files1 = array_list_size(files_fq1);
  }

  if (fq_list2) {
    ptr = strtok(fq_list2, token);    // Primera llamada => Primer token
    array_list_insert(strdup(ptr), files_fq2);
    while( (ptr = strtok( NULL, token)) != NULL ) {    // Posteriores llamadas
      array_list_insert(strdup(ptr), files_fq2);
    }    
    num_files2 = array_list_size(files_fq2);
  }

  
  if (fq_list2 && (num_files1 != num_files2)) {
    LOG_FATAL("Diferent number of files in paired-end/mate-pair mode");
  }
  
  extern double main_time;
  double time_alig;
  struct timeval time_start_alig, time_end_alig;  
  char *file1, *file2;

  printf("\n\n\t\t==================================================\n");
  printf("\t\t|................. M A P P I N G ................|\n");
  printf("\t\t==================================================\n");
  printf("\t\t |                ___  ___  ___                 |\n");
  printf("\t\t |              \\/ H \\/ P \\/ G \\/               |\n");
  printf("\t\t |              /\\___/\\___/\\___/\\               |\n");
  printf("\t\t |      ___  ___  ___  ___  ___  ___  ___       |\n");
  printf("\t\t |    \\/ A \\/ L \\/ I \\/ G \\/ N \\/ E \\/ R \\/     |\n");
  printf("\t\t |    /\\___/\\___/\\___/\\___/\\___/\\___/\\___/\\     |\n");
  printf("\t\t |                                              |\n");
  printf("\t\t==================================================\n");

  for (int f = 0; f < num_files1; f++) {
    file1 = array_list_get(f, files_fq1);

    if (num_files2) {
      file2 = array_list_get(f, files_fq2);
    } else {
      file2 = NULL;
    }

    fastq_batch_reader_input_init(file1, file2, 
				  options->pair_mode, options->batch_size, 
				  NULL, options->gzip, &reader_input);  
    
    if (options->pair_mode == SINGLE_END_MODE) {
      if (options->gzip) {
	reader_input.fq_gzip_file1 = fastq_gzopen(file1);
      } else {
	reader_input.fq_file1 = fastq_fopen(file1);
      }
    } else {
      if (options->gzip) {
	reader_input.fq_gzip_file1 = fastq_gzopen(file1);
	reader_input.fq_gzip_file2 = fastq_gzopen(file2);
      } else {
	reader_input.fq_file1 = fastq_fopen(file1);
	reader_input.fq_file2 = fastq_fopen(file2);
      }
    }


    f_sa = fopen("buffer_sa.tmp", "w+b");
    if (f_sa == NULL) {
      LOG_FATAL("Error opening file 'buffer_sa.tmp' \n");
    }
    
    f_hc = fopen("buffer_hc.tmp", "w+b");
    if (f_hc == NULL) {
      LOG_FATAL("Error opening file 'buffer_hc.tmp' \n");
    }

    sw_input.f_sa = f_sa;
    sw_input.f_hc = f_hc;

    if (!options->fast_mode) {
      ///////////////// BWT INDEX WORKFLOW //////////////////////
      //===================================================================================
      //-----------------------------------------------------------------------------------
      // workflow management
      //
      //
      // timing
      //struct timeval start, end;

      wf_input_t *wf_input = wf_input_new(&reader_input, batch);
      wf_input_file_t *wf_input_file = wf_input_file_new(f_sa, batch);   
      wf_input_file_t *wf_input_file_hc = wf_input_file_new(f_hc, batch);  

      //create and initialize workflow
      workflow_t *wf = workflow_new();
      workflow_stage_function_t stage_functions[] = {bwt_stage, cal_stage, 
						     sw_stage, post_pair_stage};
      char *stage_labels[] = {"BWT", "CAL", "SW", "POST PAIR"};
      workflow_set_stages(4, (workflow_stage_function_t *)&stage_functions, stage_labels, wf);
      // optional producer and consumer functions
      workflow_set_producer((workflow_producer_function_t *)fastq_reader, "FastQ reader", wf);
      workflow_set_consumer((workflow_consumer_function_t *)bam_writer, "BAM writer", wf);
  
      workflow_t *wf_last = workflow_new();
      workflow_stage_function_t stage_functions_last[] = {rna_last_stage, post_pair_stage};
      char *stage_labels_last[] = {"RNA LAST STAGE", "POST PAIR"};
      workflow_set_stages(2, (workflow_stage_function_t *)&stage_functions_last, stage_labels_last, wf_last);
      workflow_set_producer((workflow_producer_function_t *)file_reader, "Buffer reader", wf_last);
      workflow_set_consumer((workflow_consumer_function_t *)bam_writer, "BAM writer", wf_last);

      workflow_t *wf_hc = workflow_new();
      workflow_stage_function_t stage_functions_hc[] = {rna_last_hc_stage, post_pair_stage};
      char *stage_labels_hc[] = {"RNA HARD CLIPPINGS", "POST PAIR"};
      workflow_set_stages(2, (workflow_stage_function_t *)&stage_functions_hc, stage_labels_hc, wf_hc);
      workflow_set_producer((workflow_producer_function_t *)file_reader_2, "Buffer reader", wf_hc);
      workflow_set_consumer((workflow_consumer_function_t *)bam_writer, "BAM writer", wf_hc);
     

      // Create new thread POSIX for search extra Splice Junctions
      //============================================================
      pthread_attr_t attr;
      pthread_t thread;
      void *status;
      int ret;

      //Run workflow
      size_t tot_reads_in;
      size_t tot_reads_out;

      tot_reads_in = 0;
      tot_reads_out = 0;


      start_timer(time_start_alig);
      printf("\t\t|................................................|\n");
      printf("\t\t|            W O R K W F L O W    1              |\n");
      workflow_run_with(options->num_cpu_threads, wf_input, wf);

      printf("\t\t|................................................|\n");
      printf("\t\t|            W O R K W F L O W    2              |\n");
      rewind(f_sa);
      workflow_run_with(options->num_cpu_threads, wf_input_file, wf_last);

      printf("\t\t|................................................|\n");
      printf("\t\t|            W O R K W F L O W    3              |\n");
      rewind(f_hc);
      workflow_run_with(options->num_cpu_threads, wf_input_file_hc, wf_hc);
      printf("\t\t==================================================\n\n");
      stop_timer(time_start_alig, time_end_alig, time_alig);
      //start_timer(time_start_alig);

      basic_statistics_display(basic_st, 1, time_alig / 1000000, time_genome / 1000000);  

      // free memory
      workflow_free(wf);
      workflow_free(wf_last);
      workflow_free(wf_hc);

      wf_input_free(wf_input);
      wf_input_file_free(wf_input_file);
      wf_input_file_free(wf_input_file_hc);
    
    } else {
      ///////////////// SA INDEX WORKFLOW //////////////////////
      //--------------------------------------------------------------------------------------
      // workflow management
      //
      sw_optarg_t sw_optarg;
      sw_optarg_init(options->gap_open, options->gap_extend, 
		     options->match, options->mismatch, &sw_optarg);

      sa_rna_input_t sa_rna;
      sa_rna.genome    = genome;
      sa_rna.avls_list = avls_list;
      sa_rna.metaexons = metaexons;
      sa_rna.sw_optarg = &sw_optarg;
      sa_rna.file1 = f_sa;
      sa_rna.file2 = f_hc;
      //printf("FILEEEEEEEEEEEEEEEE %p\n", f_sa);
      sa_wf_batch_t *wf_batch = sa_wf_batch_new(NULL, (void *)sa_index, &writer_input, NULL, &sa_rna);
      sa_wf_input_t *wf_input = sa_wf_input_new(0, &reader_input, wf_batch);
      
      // create and initialize workflow
      workflow_t *wf = workflow_new();      
      workflow_stage_function_t stage_functions[] = {sa_rna_mapper};
      char *stage_labels[] = {"SA mapper"};
      workflow_set_stages(1, (workflow_stage_function_t *)&stage_functions, stage_labels, wf);      
      // optional producer and consumer functions
      workflow_set_producer(sa_fq_reader_rna, "FastQ reader", wf);
      workflow_set_consumer(sa_sam_writer_rna, "SAM writer", wf);
      
      //Create and initialize second workflow
      workflow_t *wf_last = workflow_new();
      workflow_stage_function_t stage_functions_last[] = {sa_rna_mapper_last};
      char *stage_labels_last[] = {"SA mapper last stage"};
      workflow_set_stages(1, (workflow_stage_function_t *)&stage_functions_last, stage_labels_last, wf_last);
      workflow_set_producer(sa_alignments_reader_rna, "FastQ reader", wf_last);
      workflow_set_consumer(sa_sam_writer_rna, "SAM writer", wf_last);

      double time_total;
      struct timeval time_s1, time_e1;

      printf("Run workflow with %i threads\n", options->num_cpu_threads);

      start_timer(time_s1);
      workflow_run_with(options->num_cpu_threads, wf_input, wf);
      stop_timer(time_s1, time_e1, time_total);

      printf("Time W1: %f(s)\n", time_total / 1000000);
      
      printf("=============== SECOND WORKFLOW ================\n");
      start_timer(time_s1);
      rewind(f_sa);
      workflow_run_with(options->num_cpu_threads, wf_input, wf_last);      
      stop_timer(time_s1, time_e1, time_total);

      printf("Time W2: %f(s)\n", time_total / 1000000);

      //printf("TOTAL TIME : %f\n", time_total / 1000000);
      //printf("----------------------------------------------\n");
      //workflow_display_timing(wf);
      //printf("----------------------------------------------\n");   
      
      // free memory
      sa_wf_input_free(wf_input);
      sa_wf_batch_free(wf_batch);
      workflow_free(wf);      
      workflow_free(wf_last);      

      //printf("\n");
      //for (int x = 0; x <= 10; x++) {
      //printf("%i CALs: %i reads (%f)\n", x, tot_cals[x], ((float)tot_cals[x]*100)/(float)total_reads );
      //}
      
    }

    if (file1) { free(file1); }
    if (file2) { free(file2); }

    fclose(f_sa);
    fclose(f_hc);

    //closing files
    if (options->gzip) {
      if (options->pair_mode == SINGLE_END_MODE) {
	fastq_gzclose(reader_input.fq_gzip_file1);
      } else {
	fastq_gzclose(reader_input.fq_gzip_file1);
	fastq_gzclose(reader_input.fq_gzip_file2);
      }
    } else {
      if (options->pair_mode == SINGLE_END_MODE) {
	fastq_fclose(reader_input.fq_file1);
      } else {
	fastq_fclose(reader_input.fq_file1);
	fastq_fclose(reader_input.fq_file2);
      }
    }    
  }

  array_list_free(files_fq1, NULL);
  array_list_free(files_fq2, NULL);

  //Write chromosome avls
  write_chromosome_avls(extend_filename,
			exact_filename, num_chromosomes, avls_list);
  
  


  /*  
  printf("= = = = T I M I N G    W O R K F L O W    '1' = = = =\n");
  workflow_display_timing(wf);
  printf("= = = = - - - - - - - - - - - - - - - - - - - = = = =\n\n");
  
  printf("= = = = T I M I N G    W O R K F L O W    '2' = = = =\n");
  workflow_display_timing(wf_last);
  printf("= = = = - - - - - - - - - - - - - - - - - - - = = = =\n\n");
  
  printf("= = = = T I M I N G    W O R K F L O W    '3' = = = =\n");
  workflow_display_timing(wf_hc); 
  printf("= = = = - - - - - - - - - - - - - - - - - - - = = = =\n\n");
  */
  

  
  //Write statistics to file 

  //free(basic_st);
  
  //if (time_on){ timing_free(timing); }
  metaexons_free(metaexons);
  /*  
  //========================= O P E N M P    P I P E L I N E =============================//
  
  batch_t *batch = batch_new(&bwt_input, &region_input, &cal_input, 
  &pair_input, &preprocess_rna, &sw_input, &writer_input, RNA_MODE, NULL);

  
  int num_threads = options->num_cpu_threads - 2;
  
  list_t fastq_list;
  list_init("fastq-list", 1, options->num_cpu_threads * 3, &fastq_list);
  
  list_t bam_list;
  list_init("bam-list", num_threads, options->num_cpu_threads * 3, &bam_list);
  //  list_init("bam-list", 1, num_threads * 3, &bam_list);
  
  //omp_set_num_threads(num_threads);
  
  int num_cpus = 64;
  int cpuArray[num_cpus];    
  for (int i = 0; i < num_cpus; i++) {
  cpuArray[i] = i;
  }
  omp_set_nested(1);
    
  double time_alig;
  struct timeval time_start_alig, time_end_alig;
  start_timer(time_start_alig);
  
  #pragma omp parallel sections num_threads (3)
  {      
  #pragma omp section
  {
  printf("OMP_THREAD (%i): START READER\n", omp_get_thread_num());
  //testing_reader(fastq_filename, &reader_input, in);	  
  // FastQ batch reader
  //struct timeval start_time, end_time;
  //double reading_time = 0, total_reading_time = 0;	  
  //int id = omp_get_thread_num();
  cpu_set_t cpu_set;
  CPU_ZERO(&cpu_set);
  CPU_SET(0, &cpu_set);
  sched_setaffinity(syscall(SYS_gettid), sizeof(cpu_set), &cpu_set);
  
  list_item_t *item;
  void *data;
  size_t num_batches = 0;
  while (1) {
  //reading_time = 0;
  data = fastq_reader(wf_input);
  
  if ((data) == NULL) break;
  
  item = list_item_new(num_batches, 0, data);
  list_insert_item(item, &fastq_list);
  num_batches++;
  }
  //printf("OMP_THREAD: END READER\n");
  list_decr_writers(&fastq_list);
  //printf("Reading time: %0.4f sec\n", total_reading_time / 1000000.0f);
  }
  #pragma omp section
  {
  // batch mapper
  #pragma omp parallel num_threads(num_threads)
  {
  cpu_set_t cpu_set;
  CPU_ZERO( &cpu_set);
  int id = omp_get_thread_num() + 2;
  CPU_SET( cpuArray[ id % num_cpus], &cpu_set);
  sched_setaffinity(syscall(SYS_gettid), sizeof(cpu_set), &cpu_set);
  
  printf("START WORKER: %i\n", omp_get_thread_num());
  list_item_t *item;
  while ((item = list_remove_item(&fastq_list)) != NULL) {
  bwt_stage(item->data_p);
  cal_stage(item->data_p);
  sw_stage(item->data_p);
  post_pair_stage(item->data_p);
  list_insert_item(item, &bam_list);
  }
  //printf("END WORKER %i\n", omp_get_thread_num());
  list_decr_writers(&bam_list);
  }
  }
  #pragma omp section
  {
  cpu_set_t cpu_set;
  //int id = omp_get_thread_num();
  CPU_ZERO(&cpu_set);
  CPU_SET(1, &cpu_set);
  sched_setaffinity(syscall(SYS_gettid), sizeof(cpu_set), &cpu_set);
  //printf("OMP_THREAD (%i): START WRITER\n", omp_get_thread_num());
  list_item_t *item;
  //wf_batch_t *wf_batch;
  size_t num_batches = 0;
  while ((item = list_remove_item(&bam_list)) != NULL) {
  //writing_time = 0;
  //start_timer(start_time);
  //bam_writer1(item->data_p);
  //stop_timer(start_time, end_time, writing_time);
  //total_writing_time += writing_time;
  bam_writer(item->data_p);
  list_item_free(item);
  }
  //printf("OMP_THREAD: END WRITER\n");
  //printf("Writing time: %0.4f sec\n", total_writing_time / 1000000.0f);
  }
  }
  
  rewind(f_sa);
  fastq_list.writers = 1;
  bam_list.writers  = num_threads;
  
  #pragma omp parallel sections num_threads (3)
  {
  #pragma omp section
  {
  printf("OMP_THREAD: START READER\n");
  //testing_reader(fastq_filename, &reader_input, in);	  
  // FastQ batch reader
  //struct timeval start_time, end_time;
  //double reading_time = 0, total_reading_time = 0;	  
  cpu_set_t cpu_set;
  CPU_ZERO(&cpu_set);
  CPU_SET(0, &cpu_set);
  sched_setaffinity(syscall(SYS_gettid), sizeof(cpu_set), &cpu_set);
  
  list_item_t *item;
  void *data;
  size_t num_batches = 0;
  while (1) {
  //reading_time = 0;
  data = file_reader(wf_input_file);
  if ((data) == NULL) break;
  
  item = list_item_new(num_batches, 0, data);
  list_insert_item(item, &fastq_list);
  num_batches++;
  }
  printf("OMP_THREAD: END READER\n");
  list_decr_writers(&fastq_list);
  //printf("Reading time: %0.4f sec\n", total_reading_time / 1000000.0f);
  }
  #pragma omp section
  {
  // batch mapper
  #pragma omp parallel num_threads(num_threads)
  {
  cpu_set_t cpu_set;
  CPU_ZERO( &cpu_set);
  int id = omp_get_thread_num() + 2;
  CPU_SET( cpuArray[ id % num_cpus], &cpu_set);
  sched_setaffinity(syscall(SYS_gettid), sizeof(cpu_set), &cpu_set);
  
  printf("START WORKER: %i\n", omp_get_thread_num());
  list_item_t *item;
  while ((item = list_remove_item(&fastq_list)) != NULL) {
  rna_last_stage(item->data_p);
  post_pair_stage(item->data_p);
  list_insert_item(item, &bam_list);
  }
  printf("END WORKER %i\n", omp_get_thread_num());
  list_decr_writers(&bam_list);
  }
  }
  #pragma omp section
  {
  // BAM batch writer
  //struct timeval start_time, end_time;
  //double writing_time = 0, total_writing_time = 0;
  cpu_set_t cpu_set;
  //int id = omp_get_thread_num();
  CPU_ZERO(&cpu_set);
  CPU_SET(1, &cpu_set);
  sched_setaffinity(syscall(SYS_gettid), sizeof(cpu_set), &cpu_set);
  
  printf("OMP_THREAD: START WRITER\n");
  list_item_t *item;
  //wf_batch_t *wf_batch;
  size_t num_batches = 0;
  while ((item = list_remove_item(&bam_list)) != NULL) {
  //writing_time = 0;
  //start_timer(start_time);
  //bam_writer1(item->data_p);
  //stop_timer(start_time, end_time, writing_time);
  //total_writing_time += writing_time;
  bam_writer(item->data_p);
  list_item_free(item);
  }
  printf("OMP_THREAD: END WRITER\n");
  //printf("Writing time: %0.4f sec\n", total_writing_time / 1000000.0f);
  }
  }
  
  rewind(f_hc);
  fastq_list.writers = 1;
  bam_list.writers  = num_threads;
  
  #pragma omp parallel sections num_threads (3)
  {
  #pragma omp section
  {
  cpu_set_t cpu_set;
  CPU_ZERO(&cpu_set);
  CPU_SET(0, &cpu_set);
  sched_setaffinity(syscall(SYS_gettid), sizeof(cpu_set), &cpu_set);

  printf("OMP_THREAD: START READER\n");
  //testing_reader(fastq_filename, &reader_input, in);	  
  // FastQ batch reader
  //struct timeval start_time, end_time;
  //double reading_time = 0, total_reading_time = 0;	  
  list_item_t *item;
  void *data;
  size_t num_batches = 0;
  while (1) {
  //reading_time = 0;
  data = file_reader_2(wf_input_file_hc);
  
  if ((data) == NULL) break;
  
  item = list_item_new(num_batches, 0, data);
  list_insert_item(item, &fastq_list);
  num_batches++;
  }
  printf("OMP_THREAD: END READER\n");
  list_decr_writers(&fastq_list);
  //printf("Reading time: %0.4f sec\n", total_reading_time / 1000000.0f);
  }
  #pragma omp section
  {
  // batch mapper
  #pragma omp parallel num_threads(num_threads)
  {
  cpu_set_t cpu_set;
  CPU_ZERO( &cpu_set);
  int id = omp_get_thread_num() + 2;
  CPU_SET( cpuArray[ id % num_cpus], &cpu_set);
  sched_setaffinity(syscall(SYS_gettid), sizeof(cpu_set), &cpu_set);
  
  printf("START WORKER: %i\n", omp_get_thread_num());
  list_item_t *item;
  while ((item = list_remove_item(&fastq_list)) != NULL) {
  rna_last_hc_stage(item->data_p);
  post_pair_stage(item->data_p);
  list_insert_item(item, &bam_list);
  }
  printf("END WORKER %i\n", omp_get_thread_num());
  list_decr_writers(&bam_list);
  }
  }
  #pragma omp section
  {
  // BAM batch writer
  //struct timeval start_time, end_time;
  //double writing_time = 0, total_writing_time = 0;
  cpu_set_t cpu_set;
  //int id = omp_get_thread_num();
  CPU_ZERO(&cpu_set);
  CPU_SET(1, &cpu_set);
  sched_setaffinity(syscall(SYS_gettid), sizeof(cpu_set), &cpu_set);
  
  printf("OMP_THREAD: START WRITER\n");
  list_item_t *item;
  //wf_batch_t *wf_batch;
  size_t num_batches = 0;
  while ((item = list_remove_item(&bam_list)) != NULL) {
  //writing_time = 0;
  //start_timer(start_time);
  //bam_writer1(item->data_p);
  //stop_timer(start_time, end_time, writing_time);
  //total_writing_time += writing_time;
  bam_writer(item->data_p);
  list_item_free(item);
  }
  printf("OMP_THREAD: END WRITER\n");
  //printf("Writing time: %0.4f sec\n", total_writing_time / 1000000.0f);
  }
  }
  //Write chromosome avls
  write_chromosome_avls(extend_filename,
  exact_filename, genome->num_chromosomes, avls_list);
  stop_timer(time_start_alig, time_end_alig, time_alig);
  
  //}
  */

  if (!options->fast_mode) {
    // free memory
    if (bwt_index) { bwt_index_free(bwt_index); }
    bam_fclose(writer_input.bam_file);
  } else {
    // free memory
    if (sa_index)  { sa_index3_free(sa_index);  }
    fclose((FILE *) writer_input.bam_file);
  }
 
  batch_free(batch);

  //
  //
  // end of workflow management
  //--------------------------------------------------------------------------------------

  if (genome) {
    genome_free(genome);
  }

  bwt_optarg_free(bwt_optarg);
  cal_optarg_free(cal_optarg);
  pair_mng_free(pair_mng);
  report_optarg_free(report_optarg);


  free(log_filename);
  free(output_filename);
  free(exact_filename);
  free(extend_filename);
  linked_list_free(alignments_list, (void *)NULL);
  linked_list_free(buffer, (void *)NULL);
  linked_list_free(buffer_hc, (void *)NULL);

}

//--------------------------------------------------------------------------------------
//        T R A N S C R I P T O M E     M A N A G E M E N T
//--------------------------------------------------------------------------------------

static inline char *parse_attribute(char *name, char *attrs) {
  int name_len = strlen(name);
  char *p1 = strstr(attrs, name);
  char * p2 = strstr(p1 + name_len, "\";");
  int len = p2 - p1 - name_len;
  char *res = (char *) malloc(len + 2);
  strncpy(res, p1 + name_len, len);
  res[len] = 0;
  return res;
}

//--------------------------------------------------------------------------------------

static inline exon_t *parse_exon_line(FILE *f, genome_t *genome) {
  const int MAX_LENGTH = 8192;

  int field, found;
  char line[MAX_LENGTH], *token, *str, *p1, *p2;
  
  int chr, strand, start, end, exon_number;
  char *chr_name, *gene_id, *transcript_id, *exon_id, *exon_number_str;

  while (fgets(line, MAX_LENGTH, f) != NULL) {
    //    LOG_DEBUG_F("%s", line);
    str = strdup(line);
    token = strtok(str, "\t");
    field = 0;
    while (token != NULL) {
      if (field == 0) { // seqname          

	found = 0;
	for (chr = 0; chr < genome->num_chromosomes; chr++) {
	  if (strcmp(token, genome->chr_name[chr]) == 0) { found = 1; break; }
	}
	if (!found) break;

	chr_name = strdup(token);
      } else if (field == 2) { // feature type name: exon
	if (strcmp(token, "exon") != 0) {
	  found = 0;
	  free(chr_name);
	  break;
	}
      } else if (field == 3) { // start position of the feature (starting at 1)
	start = atoi(token);
      } else if (field == 4) { // end position of the feature (starting at 1)  
        end = atoi(token);
      } else if (field == 6) { // strand + (forward) or - (reverse)
        strand = ((strcmp(token, "-") == 0) ? 1 : 0);
	//	printf("---> strand = %i, token = %s\n", strand, token);
      } else if (field == 8) { // attributes, a semicolon-separated list of tag-value pairs
	// gene_id, transcript_id, exon_id
	gene_id = parse_attribute("gene_id \"", token);
	transcript_id = parse_attribute("transcript_id \"", token);
	exon_id = parse_attribute("exon_id \"", token);

	// exon_number
	exon_number_str = parse_attribute("exon_number \"", token);
	if (exon_number_str) {
	  exon_number = atoi(exon_number_str);
	  free(exon_number_str);
	}
      }
      
      token = strtok(NULL, "\t");
      field++;
    }
    free(str);
    if (found) {
      return exon_new(chr, strand, start, end, exon_number,
		      chr_name, gene_id, transcript_id, exon_id);
    }
  }
  return NULL;
}

//--------------------------------------------------------------------------------------

void load_transcriptome(char *filename, genome_t *genome, 
			avls_list_t *avls_list, metaexons_t *metaexons) {

  FILE *f = fopen(filename, "r");

  int pos, direction, count = 0;
  exon_t *exon = NULL, *exon1 = NULL, *exon2 = NULL;

  size_t g_start, g_end;
  int type, splice_strand, strand;
  char nt_start[2], nt_end[2];
  char *transcript_id;

  array_list_t *list = array_list_new(100, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  
  while (1) {

    // read the first exon
    if (!exon) {
      exon = parse_exon_line(f, genome);  
      if (!exon) {
	break;
      }
    }

    array_list_insert(exon, list);
    strand = exon->strand;
    transcript_id = exon->transcript_id;

    while ((exon = parse_exon_line(f, genome)) && strcmp(exon->transcript_id, transcript_id) == 0) {
      array_list_insert(exon, list);
    }
    
    // process the whole list
    if (strand) {
      // - strand
      pos = list->size - 1;
      direction = -1;
    } else {
      // + strand
      pos = 0;
      direction = 1;
    }

    exon1 = NULL;
    exon2 = NULL;

    for (int i = 0; i < list->size; i++, pos += direction) {

      //      printf("---> pos = %i\n", pos);

      // get exon1
      if (!exon1) {
	exon1 = array_list_get(pos, list);
	count++;
	continue;
      }

      // get exon2
      if (!exon2) {
	exon2 = array_list_get(pos, list);
	if (exon2) {
	  //	exon_display(exon2);
	  count++;
	}
      }

      if (exon1 && exon2) {
	//	printf("Process Exon1[%lu-%lu] and Exon2[%lu-%lu]\n", exon1->start, exon1->end, exon2->start, exon2->end);
	// exons belonging to the same transcript
	// process exon1 and exon2, to init the splice junction
	g_start = exon1->end + 1;
	g_end = g_start + 1;
	//	printf("chr = %s (%i), strand %i, g_start = %i, g_end = %i\n", exon1->chr_name, exon1->chr, exon1->strand, g_start, g_end);
	genome_read_sequence_by_chr_index(nt_start, exon1->strand, exon1->chr, &g_start, &g_end, genome);
	//	printf("\tnt_start = %s\n", nt_start);
	
	g_end = exon2->start - 1;
	g_start = g_end - 1;
	//	printf("chr = %s (%i), strand %i, g_start = %i, g_end = %i\n", exon2->chr_name, exon2->chr, exon2->strand, g_start, g_end);
	genome_read_sequence_by_chr_index(nt_end, exon2->strand, exon2->chr, &g_start, &g_end, genome);
	//	printf("\tnt_end = %s\n", nt_end);
	
	type = splice_junction_type(nt_start[0], nt_start[1], nt_end[0], nt_end[1]);
	int splice_strand;
	avl_node_t *avl_node_start, *avl_node_end;
	
	splice_strand = 0;
	if (type == CT_AC_SPLICE || type == GT_AT_SPLICE || type == CT_GC_SPLICE ) {
	  splice_strand = 1;
	}
	
	//LOG_FATAL_F("nt_start = %s, nt_end = %s, type = %i\n", nt_start, nt_end, type);
	
	//LOG_DEBUG_F("start_splice = %lu - end_splice = %lu (%s, %s, %s)\n", exon1->end, exon2->start, exon1->transcript_id, exon2->transcript_id, transcript_id);
	if (exon2->start - 1 < exon1->end + 1) {
	  LOG_FATAL_F("start_splice = %lu - end_splice = %lu (%s, %s, %s)\n", exon1->end, exon2->start, exon1->transcript_id, exon2->transcript_id, transcript_id);
	}
	
	allocate_start_node(exon1->chr, // startint at 0
			    splice_strand,
			    exon1->end + 1,   // splice start
			    exon2->start - 1, // splice_end,
			    exon1->end + 1,   // splice start
			    exon2->start - 1, // splice_end,
			    FROM_FILE,
			    type,
			    NULL, 
			    &avl_node_start,
			    &avl_node_end, 
			    avls_list);
	
	metaexon_insert(0/*exon1->strand*/, exon1->chr, exon1->start, exon1->end, 40,
			METAEXON_RIGHT_END, avl_node_start, metaexons);
	
	metaexon_insert(0/*exon2->strand*/, exon2->chr, exon2->start, exon2->end, 40,
			METAEXON_LEFT_END, avl_node_end, metaexons);
	
	// ...and then, free exon1 and update it to exon2
	exon_free(exon1);
	exon1 = exon2;
	exon2 = NULL;
      } else if (exon1) {
	// process the last exon (exon1) but only if it's the first (no splice)
	// (otherwise, it was already processed, there was a splice)
	if (exon1->exon_number == 1) {
	  //	  printf("Process Exon1[%lu-%lu]\n", exon1->start, exon1->end);
	  metaexon_insert(0/*exon1->strand*/, exon1->chr, exon1->start, exon1->end, 40,
			  METAEXON_NORMAL, NULL, metaexons);
	}
	// and free and exit
	exon_free(exon1);
	exon1 = NULL;
	exon2 = NULL;
	break;
      }
    } // end for

    if (exon1) {
      // process the last exon (exon1) but only if it's the first (no splice)
      // (otherwise, it was already processed, there was a splice)
      if (exon1->exon_number == 1) {
	//	printf("Process Exon1[%lu-%lu]\n", exon1->start, exon1->end);
	metaexon_insert(0/*exon1->strand*/, exon1->chr, exon1->start, exon1->end, 40,
			METAEXON_NORMAL, NULL, metaexons);
      }
      // and free and exit
      exon_free(exon1);
    }
    
    array_list_clear(list, NULL);
  } // end while

  fclose(f);
  array_list_free(list, NULL);

  LOG_DEBUG_F("Number of processed exons: %i", count);

}

//--------------------------------------------------------------------
//--------------------------------------------------------------------










