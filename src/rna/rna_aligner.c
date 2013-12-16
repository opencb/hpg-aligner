#include "rna_aligner.h"

#define NUM_SECTIONS_TIME 		8

//--------------------------------------------------------------------
// workflow input                                                                                                                                                     
//--------------------------------------------------------------------  
extern int splice_junction_type(char nt_start_1, char nt_start_2, char nt_end_1, char nt_end_2);


//--------------------------------------------------------------------


/*
  buffer_item_t *buffer_item_new(fastq_read_t *fq_read, array_list_t *array_list) {
  buffer_item_t *item = (buffer_item_t *)malloc(sizeof(buffer_item_t));
  item->read = fq_read;
  item->alignments_list = array_list;

  return item;
}

buffer_pair_item_t *buffer_pair_item_new(fastq_read_t *fq_read_1, array_list_t *array_list_1, 
					 fastq_read_t *fq_read_2, array_list_t *array_list_2) {
  buffer_pair_item_t *item = (buffer_pair_item_t *)malloc(sizeof(buffer_pair_item_t));
  item->read_1 = fq_read_1;
  item->alignments_list_1 = array_list_1;

  item->read_2 = fq_read_2;
  item->alignments_list_2 = array_list_2;

  return item;
}

int convert_alignments_to_CAL(array_list_t *alignments_list, array_list_t *cal_list) {
  alignment_t *alginment;
  size_t num_items = array_list_size(alignments_list);
  int found = 0;

  for (int i = num_items - 1; i >= 0; i--) {
    alignment = array_list_get_at(i, alignments_list);
    if (alignment->large_hard_clipping) {
      found = 1;
      array_list_remove_at(i, alignments_list);
      //Extract strand, chr, start and search in AVL. Nex, we will form CALs and insert to cal_list
    } else {
      alignment->primary_alignment = 0;
    }
  }
  return found;
}

void thread_function(extra_stage_t *extra_stage_input) {
  linked_list_t *list = extra_stage_input->align_list;
  linked_list_iterator_init(align_list, itr);
  workflow_t *wf = extra_stage_input->workflow;
  pair_mng_t *pair_mng = extra_stage_input->pair_mng;
  batch_t *batch;
  size_t num_alignments;
  alignment_t *align;
  void *buffer_item;
  buffer_pair_item_t *item_pair;
  buffer_item_t *item;
  mapping_batch_t *mapping_batch = mapping_batch_new_by_num(MAX_READS_RNA + 1, pair_mng);
  array_list_t *fq_batch = array_list_new(MAX_READS_RNA, 
					  1.25f, 
					  COLLECTION_MODE_ASYNCHRONIZED);
  size_t num_reads = 0;
  size_t num_targets = 0;
  int found;
  global_status = WORKFLOW_STATUS_RUNNING;
  wf->complete_extra_stage = 0;

  printf("Extrem search STARTs...\n");
  
  while(workflow_get_simple_status(wf) == WORKFLOW_STATUS_RUNNING) {
    pthread_cond_wait(&cond_sp, &mutex_sp);
    while(buffer_item = linked_list_remove_last(list) {
	if (linked_list_get_flag(list) != SINGLE_END_MODE) {
	  buffer_item_t *item = (buffer_item_t *)buffer_item;
	  array_list_insert(item->read, fq_batch);
	  convert_alignments_to_CAL(item->alignment_list, mapping_batch->mapping_lists[num_reads]);
	  mapping_batch->old_mapping_list[num_reads] = item->alignment_list;
	  mapping_batch->targets[num_targets++] = 1;
	} else {
	  buffer_pair_item_t *item = (buffer_pair_item_t *)buffer_item;
	  array_list_insert(item->read_1, fq_batch); 
	  if (convert_alignments_to_CAL(item->alignment_list_1, mapping_batch->mapping_lists[num_reads])) {
	    mapping_batch->targets[num_targets++] = 1;
	  }
	  mapping_batch->old_mapping_list[num_reads++] = item->alignment_list_1;

	  array_list_insert(item->read_2, fq_batch);
	  if (convert_alignments_to_CAL(item->alignment_list_2, mapping_batch->mapping_lists[num_reads])) {
	    mapping_batch->targets[num_targets++] = 1;
	  }
	  mapping_batch->old_mapping_list[num_reads++] = item->alignment_list_2;
	}
	if (array_list_size(fq_batch) >= MAX_READS_RNA ) {
	  mapping_batch->num_targets = num_targets;
	  //Insert in the corrected site
	  mapping_batch = mapping_batch_new_by_num(MAX_READS_RNA + 1, pair_mng);
	}
      }
  }

  wf->complete_extra_stage = 1;
  
  printf("Finish search!\n");
}
*/

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



void run_rna_aligner(genome_t *genome, bwt_index_t *bwt_index, pair_mng_t *pair_mng,
		     bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
		     report_optarg_t *report_optarg, metaexons_t *metaexons, options_t *options) {
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
    strcat(reads_results, "_alignments.bam");  

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
    strcat(reads_results, "/alignments.bam");
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


  extern FILE *fd_log;
  fd_log = fopen(log_filename, "w");

  LOG_DEBUG("Auto Thread Configuration Done !");

  // timing
  if (time_on) { 
    char* labels_time[NUM_SECTIONS_TIME] = {"FASTQ Reader               ", 
                                            "BWT Server                 ", 
					    "REGION Seeker              ", 
					    "CAL Seeker                 ", 
					    "RNA Preprocess             ", 
					    "RNA Server                 ",
					    "BAM Writer                 ", 
					    "TOTAL Time                 "};    
    timing = timing_new((char**) labels_time, NUM_SECTIONS_TIME);
  }

  // display selected options
  LOG_DEBUG("Displaying options...\n");
  options_display(options);

  //============================= INPUT INITIALIZATIONS =========================
  //allocate_splice_elements_t chromosome_avls[genome->num_chromosomes];
  //init_allocate_splice_elements(chromosome_avls, genome->num_chromosomes);
  avls_list_t* avls_list = avls_list_new(genome->num_chromosomes);
  printf("Loading transcriptome...\n");
  if (options->transcriptome_filename != NULL) {
    load_transcriptome(options->transcriptome_filename, genome, avls_list, metaexons);
    //    LOG_FATAL_F("transcriptome filename = %s\n", options->transcriptome_filename);
  }
  printf("Load done!\n");
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
				NULL, &reader_input);  

  //buffer_reader_input_t buffer_reader_input;
  //buffer_reader_input_init(&reader_input,
  //buffer,
  //buffer_reader_input);

  if (options->pair_mode == SINGLE_END_MODE) {
    reader_input.fq_file1 = fastq_fopen(options->in_filename);
  } else {
    reader_input.fq_file1 = fastq_fopen(options->in_filename);
    reader_input.fq_file2 = fastq_fopen(options->in_filename2);
  }

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
			   alignments_list, genome, &writer_input);

  bam_header_t *bam_header = create_bam_header_by_genome(genome);
  writer_input.bam_file = bam_fopen_mode(output_filename, bam_header, "w");
  bam_fwrite_header(bam_header, writer_input.bam_file);
  
  extra_stage_t extra_stage_input;
  
  int workflow_enable = options->workflow_enable;

  batch_t *batch = batch_new(&bwt_input, &region_input, &cal_input, 
			     &pair_input, &preprocess_rna, &sw_input, &writer_input, RNA_MODE, NULL);

  wf_input_t *wf_input = wf_input_new(&reader_input, batch);
  wf_input_file_t *wf_input_file = wf_input_file_new(f_sa, batch);   
  wf_input_file_t *wf_input_file_hc = wf_input_file_new(f_hc, batch);  
  
  if (workflow_enable) {
    //===================================================================================
    //-----------------------------------------------------------------------------------
    // workflow management
    //
    //
    // timing
    //struct timeval start, end;
    extern double main_time;

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
     
    /*
    workflow_t *wf = workflow_new();
    workflow_stage_function_t stage_functions[] = { w1_function };
    char *stage_labels[] = {"W1"};
    workflow_set_stages(1, (workflow_stage_function_t *)&stage_functions, stage_labels, wf);
    // optional producer and consumer functions
    workflow_set_producer((workflow_producer_function_t *)fastq_reader, "FastQ reader", wf);
    workflow_set_consumer((workflow_consumer_function_t *)bam_writer, "BAM writer", wf);
  
    workflow_t *wf_last = workflow_new();
    workflow_stage_function_t stage_functions_last[] = { w2_function };
    char *stage_labels_last[] = { "W2" };
    workflow_set_stages(1, (workflow_stage_function_t *)&stage_functions_last, stage_labels_last, wf_last);
    workflow_set_producer((workflow_producer_function_t *)file_reader, "Buffer reader", wf_last);
    workflow_set_consumer((workflow_consumer_function_t *)bam_writer, "BAM writer", wf_last);

    workflow_t *wf_hc = workflow_new();
    workflow_stage_function_t stage_functions_hc[] = { w3_function };
    char *stage_labels_hc[] = { "W3" };
    workflow_set_stages(1, (workflow_stage_function_t *)&stage_functions_hc, stage_labels_hc, wf_hc);
    workflow_set_producer((workflow_producer_function_t *)file_reader_2, "Buffer reader", wf_hc);
    workflow_set_consumer((workflow_consumer_function_t *)bam_writer, "BAM writer", wf_hc);
    */

    // Create new thread POSIX for search extra Splice Junctions
    //============================================================
    pthread_attr_t attr;
    pthread_t thread;
    void *status;
    int ret;

    //Run workflow
    //extern size_t num_reads_map;
    //extern size_t num_reads;
  
    //num_reads_map = 0;
    //num_reads     = 0;
    extern size_t tot_reads_in;
    extern size_t tot_reads_out;

    tot_reads_in = 0;
    tot_reads_out = 0;

    extern double time_alig;
    extern struct timeval time_start_alig, time_end_alig;

    start_timer(time_start_alig);
    fprintf(stderr, "START WORKWFLOW '1ph'\n");
    workflow_run_with(options->num_cpu_threads, wf_input, wf);
    //fprintf(stderr, "TOTAL READS MAP %lu / %lu\n", num_reads_map, num_reads);
    fprintf(stderr, "END WORKWFLOW '1ph'\n\n");
  
    //fprintf(stderr, "TOTAL READS PROCESS IN: %lu\n", tot_reads_in);
    //fprintf(stderr, "TOTAL READS PROCESS OUT: %lu\n", tot_reads_out);

    tot_reads_in = 0;
    tot_reads_out = 0;  
    //num_reads_map = 0;
    //num_reads     = 0;
    fprintf(stderr, "START WORKWFLOW '2ph'\n");
    rewind(f_sa);
    workflow_run_with(options->num_cpu_threads, wf_input_file, wf_last);
    //fprintf(stderr, "TOTAL READS MAP %lu / %lu\n", num_reads_map, num_reads);
    fprintf(stderr, "END WORKWFLOW '2ph'\n\n");

    //fprintf(stderr, "TOTAL READS PROCESS IN: %lu\n", tot_reads_in);
    //fprintf(stderr, "TOTAL READS PROCESS OUT: %lu\n", tot_reads_out);
  
    //num_reads_map = 0;
    //num_reads     = 0;
    tot_reads_in = 0;
    tot_reads_out = 0;  
    fprintf(stderr, "START WORKWFLOW '3ph'\n");
    rewind(f_hc);
    workflow_run_with(options->num_cpu_threads, wf_input_file_hc, wf_hc);
    //fprintf(stderr, "TOTAL READS MAP %lu / %lu\n", num_reads_map, num_reads);
    fprintf(stderr, "END WORKWFLOW '3ph'\n\n");
  
    //fprintf(stderr, "TOTAL READS PROCESS IN: %lu\n", tot_reads_in);
    //fprintf(stderr, "TOTAL READS PROCESS OUT: %lu\n", tot_reads_out);

    //extern size_t w2_3_r;    
    //extern size_t w2_r;
    //extern size_t w3_r;
    //fprintf(stderr, "w2_r = %lu, w3_r = %lu, w2_3_r = %lu\n", w2_r, w3_r, w2_3_r);

    /*options->num_cpu_threads*/

    //Write chromosome avls
    write_chromosome_avls(extend_filename,
			  exact_filename, genome->num_chromosomes, avls_list);

    /*if (time_on) { 
      stop_timer(start, end, time);
      timing_add(time, TOTAL_TIME, timing);
      }*/
    stop_timer(time_start_alig, time_end_alig, time_alig);

    printf("= = = = T I M I N G    W O R K F L O W    '1' = = = =\n");
    workflow_display_timing(wf);
    printf("= = = = - - - - - - - - - - - - - - - - - - - = = = =\n\n");

    printf("= = = = T I M I N G    W O R K F L O W    '2' = = = =\n");
    workflow_display_timing(wf_last);
    printf("= = = = - - - - - - - - - - - - - - - - - - - = = = =\n\n");

    printf("= = = = T I M I N G    W O R K F L O W    '3' = = = =\n");
    workflow_display_timing(wf_hc); 
    printf("= = = = - - - - - - - - - - - - - - - - - - - = = = =\n\n");

    extern size_t TOTAL_SW,
      TOTAL_READS_PROCESS,
      TOTAL_READS_SEEDING,
      TOTAL_READS_SEEDING2,
      TOTAL_READS_SA;

    printf("TOTAL READS PROCESS = %lu,\n TOTAL READS SEEDING x1 = %lu,\n TOTAL READS SEEDING x2 = %lu,\n TOTAL SW = %lu,\n TOTAL READS SINGLE ANCHOR FINAL = %lu\n\n",
	   TOTAL_READS_PROCESS, TOTAL_READS_SEEDING, TOTAL_READS_SEEDING2, TOTAL_SW, TOTAL_READS_SA);

    // free memory
    workflow_free(wf);
    workflow_free(wf_last);
    workflow_free(wf_hc);
    
  } else {
    //**************************************************************************************//
    //========================= O P E N M P    P I P E L I N E =============================//
    //**************************************************************************************//
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
    /*
    #pragma omp parallel 
    {
      cpu_set_t cpu_set;
      int id = omp_get_thread_num();
      printf("Thread %i to CPU %i\n", id, cpuArray[ id % num_cpus]);
      CPU_ZERO( &cpu_set);
      CPU_SET( cpuArray[ id % num_cpus], &cpu_set);
      sched_setaffinity(syscall(SYS_gettid), sizeof(cpu_set), &cpu_set);
    }
*/
    omp_set_nested(1);
    extern double time_alig;
    extern struct timeval time_start_alig, time_end_alig;
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
	 /*cpu_set_t cpu_set;
	 CPU_ZERO( &cpu_set);
	 int id = omp_get_thread_num();
	 CPU_SET( cpuArray[ id % num_cpus], &cpu_set);
	 sched_setaffinity(syscall(SYS_gettid), sizeof(cpu_set), &cpu_set);	 
	 // BAM batch writer
	 //struct timeval start_time, end_time;
	 //double writing_time = 0, total_writing_time = 0;*/
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

  }

  wf_input_free(wf_input);
  wf_input_file_free(wf_input_file);
  wf_input_file_free(wf_input_file_hc);

  //wf_in_free(in);
  //wf_batch_free(batch);
  //
  // end of workflow management
  //--------------------------------------------------------------------------------------
  
  //     printf("***** (num_sws, num_ext_sws, num_gaps) = (%i, %i, %i)\n", num_sws, num_ext_sws, num_gaps);
  
  //closing files
  //fastq_fclose(reader_input.fq_file1);
  //bam_fclose(writer_input.bam_file);



  
  //closing files
  if (options->pair_mode == SINGLE_END_MODE) {
    fastq_fclose(reader_input.fq_file1);
  } else {
    fastq_fclose(reader_input.fq_file1);
    fastq_fclose(reader_input.fq_file2);
  }
  
  bam_fclose(writer_input.bam_file);
    

  
  //avls_list_free();
  batch_free(batch);

  //
  //
  // end of workflow management
  //--------------------------------------------------------------------------------------

  //printf("========== FINAL BUFFER %i =========\n", linked_list_size(buffer));

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










