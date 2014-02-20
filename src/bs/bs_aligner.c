#include "bs_aligner.h"

//--------------------------------------------------------------------
// run dna aligner
//--------------------------------------------------------------------

//void run_bs_aligner(genome_t *genome, genome_t *genome1, genome_t *genome2,
void run_bs_aligner(options_t *options) {
  genome_t *genome2, *genome1, *genome;
  bwt_index_t *bwt_index2, *bwt_index1;
  //bwt_optarg_t *bwt_optarg, cal_optarg_t *cal_optarg, 
  //pair_mng_t *pair_mng, report_optarg_t *report_optarg, 
  //options_t *options
  //printf("index1 %s\n", bwt_index1->nucleotides);
  //printf("index2 %s\n", bwt_index2->nucleotides);

  int path_length = strlen(options->output_name);
  int prefix_length = 0;
  if (options->prefix_name) {
    prefix_length = strlen(options->prefix_name);
  }
  
  char *reads_results = (char *) calloc((60 + prefix_length), sizeof(char));
  char *output_filename = (char *) calloc((path_length + prefix_length + 60), sizeof(char));
  
  if (options->prefix_name) {
    strcat(reads_results, "/");
    strcat(reads_results, options->prefix_name);
    strcat(reads_results, "_alignments.bam");  
  } else {
    strcat(reads_results, "/alignments.bam");
  }
  
  strcat(output_filename, options->output_name);
  strcat(output_filename, reads_results);
  free(reads_results);
  
  // display selected options
  LOG_DEBUG("Displaying options...\n");
  options_display(options);

    char bs_dir1[256];
    sprintf(bs_dir1, "%s/AGT_index", options->bwt_dirname);
    char bs_dir2[256];
    sprintf(bs_dir2, "%s/ACT_index", options->bwt_dirname);
    // genome parameters
    LOG_DEBUG("Reading genomes...");

    ////////////////////////
    //descomentar la siguiente linea para probar con el alfabeto de 4 letras
    ////////////////////////
    //printf("dir %s\n", options->bwt_dirname);
    genome  = genome_new("dna_compression.bin", options->bwt_dirname, BWT_MODE);
    //genome = genome_new("dna_compression.bin", bs_dir1);
    genome1 = genome_new("dna_compression.bin", bs_dir1, BWT_MODE);
    genome2 = genome_new("dna_compression.bin", bs_dir2, BWT_MODE);

    LOG_DEBUG("Done !!");
    
    // BWT index
    //if (time_on) { timing_start(INIT_BWT_INDEX, 0, timing_p); }

    //printf("Load index1 (%s)\n", bs_dir1);
    LOG_DEBUG("Loading AGT index...");
    bwt_index1 = bwt_index_new(bs_dir1, false);
    //printf("Load index1 done\n");
    /*
    printf("+++bwt_index1 -> %s\n" ,bwt_index1->nucleotides);
    bwt_index1->nucleotides = strdup(readNucleotide(bs_dir1, "Nucleotide"));
    bwt_init_replace_table(bwt_index1->nucleotides, bwt_index1->table, bwt_index1->rev_table);
    printf("***bwt_index1 -> %s\n" ,bwt_index1->nucleotides);
    */
    LOG_DEBUG("Loading AGT index done !!");

    //printf("Load index2 (%s)\n", bs_dir2);
    LOG_DEBUG("Loading ACT index...");
    bwt_index2 = bwt_index_new(bs_dir2, false);
    //printf("Load index2 done\n");
    /*
    printf("+++bwt_index2 -> %s\n", bwt_index2->nucleotides);
    bwt_index2->nucleotides = strdup(readNucleotide(bs_dir2, "Nucleotide"));
    bwt_init_replace_table(bwt_index2->nucleotides, bwt_index2->table, bwt_index2->rev_table);
    printf("***bwt_index2 -> %s\n", bwt_index2->nucleotides);
    */
    LOG_DEBUG("Loading ACT index done !!");
    //exit(-1);
  
  //BWT parameters
  bwt_optarg_t *bwt_optarg = bwt_optarg_new(1, 0,
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

  
  // preparing input FastQ file
  fastq_batch_reader_input_t reader_input;
  fastq_batch_reader_input_init(options->in_filename, options->in_filename2, 
				options->pair_mode, options->batch_size, 
				NULL, options->gzip, &reader_input);
  
  if (options->pair_mode == SINGLE_END_MODE) {
    reader_input.fq_file1 = fastq_fopen(options->in_filename);
  } else {
    reader_input.fq_file1 = fastq_fopen(options->in_filename);
    reader_input.fq_file2 = fastq_fopen(options->in_filename2);
  }
  
  // preparing output BAM file
  batch_writer_input_t writer_input;
  batch_writer_input_init(output_filename, NULL, NULL, NULL, genome, &writer_input);

  metil_file_t *metil_file = calloc(1, sizeof(metil_file_t));
  //printf("out = %s\n", options->output_name);
  //metil_file_init(metil_file, "/home/pascual/tmp/", genome);
  metil_file_init(metil_file, options->output_name, genome);
  writer_input.metil_file = metil_file;
  
  bam_header_t *bam_header = create_bam_header_by_genome(genome);
  writer_input.bam_file = bam_fopen_mode(output_filename, bam_header, "w");
  bam_fwrite_header(bam_header, writer_input.bam_file);
  
  // preparing workflow stage input
  bwt_server_input_t bwt_input;
  bwt_server_input_init(NULL, 0, bwt_optarg, bwt_index1, 
			NULL, 0, NULL, NULL, NULL, genome, &bwt_input);
  bwt_input.bwt_index2_p = bwt_index2;
  
  region_seeker_input_t region_input;
  region_seeker_input_init(NULL, cal_optarg, bwt_optarg, 
			   bwt_index1, NULL, 0, 0, 0, 0, 
			   genome, NULL, &region_input);
  region_input.bwt_index2_p = bwt_index2;
  
  cal_seeker_input_t cal_input;
  cal_seeker_input_init(NULL, cal_optarg, NULL, 0, 
			NULL, NULL, genome1, bwt_optarg, bwt_index1, NULL, &cal_input);
  cal_input.index2 = bwt_index2;
  cal_input.genome2 = genome2;
  
  pair_server_input_t pair_input;
  pair_server_input_init(pair_mng, report_optarg, NULL, NULL, NULL, &pair_input);
  
  sw_server_input_t sw_input;
  sw_server_input_init(NULL, NULL, 0, options->match, options->mismatch, 
		       options->gap_open, options->gap_extend, options->min_score, 
		       options->flank_length, genome, 0, 0, 0,  bwt_optarg, NULL, 
		       NULL, NULL, NULL, NULL, NULL, NULL, NULL, options->pair_mode, &sw_input);
  sw_input.genome1_p = genome1;
  sw_input.genome2_p = genome2;

  sw_input.valuesCT = (unsigned long long **)malloc(50 * sizeof(unsigned long long *));
  sw_input.valuesGA = (unsigned long long **)malloc(50 * sizeof(unsigned long long *));
  load_encode_context(options->bwt_dirname, sw_input.valuesCT, sw_input.valuesGA);
  //printf("-------- CT = %llu -------- GA = %llu\n", sw_input.valuesCT[0][100000], sw_input.valuesGA[0][100000]);
  //return;

  //--------------------------------------------------------------------------------------
  // workflow management
  //
  //
  // timing
  struct timeval start, end;
  double main_time;

  batch_t *batch = batch_new(&bwt_input, &region_input, &cal_input, 
			     &pair_input, NULL, &sw_input, &writer_input, BS_MODE, NULL);
  
  wf_input_t *wf_input = wf_input_new(&reader_input, batch);
  
  // create and initialize workflow
  workflow_t *wf = workflow_new();
  
  // workflow definition for the analysis without the postprocess

  workflow_stage_function_t stage_functions[] = {bwt_stage_bs, cal_stage_bs, 
						 pre_pair_stage, sw_stage_bs, post_pair_stage_bs};
  char *stage_labels[] = {"BWT", "CAL", "PRE PAIR", "SW", "POST PAIR"};
  workflow_set_stages(5, stage_functions, stage_labels, wf);

  // workflow definition for the analysis with the postprocess
  /*  
  workflow_stage_function_t stage_functions[] = {bwt_stage_bs, seeding_stage_bs, cal_stage_bs, 
  						 pre_pair_stage, sw_stage_bs, post_pair_stage_bs, bs_status_stage};
  char *stage_labels[] = {"BWT", "SEEDING", "CAL", "PRE PAIR", "SW", "POST PAIR", "BS STATUS"};
  workflow_set_stages(7, &stage_functions, stage_labels, wf);
  */
  // optional producer and consumer functions
  workflow_set_producer(fastq_reader, "FastQ reader", wf);
  workflow_set_consumer(bs_writer, "BAM BS writer", wf);
  
  //if (time_on) {
  start_timer(start);
  //}
  
  workflow_run_with(options->num_cpu_threads, wf_input, wf);
  
  //if (time_on) {
  stop_timer(start, end, main_time);
  //printf("Total Time: %4.04f sec\n", time / 1000000);
  //}
  
  // free memory
  workflow_free(wf);
  wf_input_free(wf_input);
  batch_free(batch);


  // show statistics for cytosines methylated/unmethylated

  FILE * STAT = metil_file->STAT;
  if (STAT == NULL) {
    printf("reopen Statistics file\n");
    STAT = fopen(metil_file->filenameSTAT, "a");
  }

  size_t c_analysed = 
    metil_file->MUT_methyl +
    metil_file->CpG_methyl + metil_file->CpG_unmethyl +
    metil_file->CHG_methyl + metil_file->CHG_unmethyl +
    metil_file->CHH_methyl + metil_file->CHH_unmethyl;

  //--------------------------------------------------------------------------------------
  // Same data as Bismark
  printf("\n\n");
  printf("======================================\n");
  printf("Final Cytosine Methylation Report\n");
  printf("======================================\n");
  printf("Total number of C's analysed: %lu\n", c_analysed);
  printf("\n");
  printf("Total methylated C's in CpG context: %lu\n", metil_file->CpG_methyl);
  printf("Total methylated C's in CHG context: %lu\n", metil_file->CHG_methyl);
  printf("Total methylated C's in CHH context: %lu\n", metil_file->CHH_methyl);
  printf("\n");
  printf("Total C to T conversions in CpG context: %lu\n", metil_file->CpG_unmethyl);
  printf("Total C to T conversions in CHG context: %lu\n", metil_file->CHG_unmethyl);
  printf("Total C to T conversions in CHH context: %lu\n", metil_file->CHH_unmethyl);

  printf("\n");
  printf("C methylated in CpG context: %5.2f%c\n", (metil_file->CpG_methyl + metil_file->CpG_unmethyl == 0) ? 0.0 :
	 (float) metil_file->CpG_methyl / (metil_file->CpG_methyl + metil_file->CpG_unmethyl)  * 100, '%');
  printf("C methylated in CHG context: %5.2f%c\n", (metil_file->CHG_methyl + metil_file->CHG_unmethyl == 0) ? 0.0 :
	 (float) metil_file->CHG_methyl / (metil_file->CHG_methyl + metil_file->CHG_unmethyl)  * 100, '%');
  printf("C methylated in CHH context: %5.2f%c\n", (metil_file->CHH_methyl + metil_file->CHH_unmethyl == 0) ? 0.0 :
	 (float) metil_file->CHH_methyl / (metil_file->CHH_methyl + metil_file->CHH_unmethyl)  * 100, '%');

  fprintf(STAT, "\n\n");
  fprintf(STAT, "======================================\n");
  fprintf(STAT, "Final Cytosine Methylation Report\n");
  fprintf(STAT, "======================================\n");
  fprintf(STAT, "Total number of C's analysed: %lu\n", c_analysed);
  fprintf(STAT, "\n");
  fprintf(STAT, "Total methylated C's in CpG context: %lu\n", metil_file->CpG_methyl);
  fprintf(STAT, "Total methylated C's in CHG context: %lu\n", metil_file->CHG_methyl);
  fprintf(STAT, "Total methylated C's in CHH context: %lu\n", metil_file->CHH_methyl);
  fprintf(STAT, "\n");
  fprintf(STAT, "Total C to T conversions in CpG context: %lu\n", metil_file->CpG_unmethyl);
  fprintf(STAT, "Total C to T conversions in CHG context: %lu\n", metil_file->CHG_unmethyl);
  fprintf(STAT, "Total C to T conversions in CHH context: %lu\n", metil_file->CHH_unmethyl);

  fprintf(STAT, "\n");
  fprintf(STAT, "C methylated in CpG context: %5.2f%c\n", (metil_file->CpG_methyl + metil_file->CpG_unmethyl == 0) ? 0.0 :
	  (float) metil_file->CpG_methyl / (metil_file->CpG_methyl + metil_file->CpG_unmethyl)  * 100, '%');
  fprintf(STAT, "C methylated in CHG context: %5.2f%c\n", (metil_file->CHG_methyl + metil_file->CHG_unmethyl == 0) ? 0.0 :
	  (float) metil_file->CHG_methyl / (metil_file->CHG_methyl + metil_file->CHG_unmethyl)  * 100, '%');
  fprintf(STAT, "C methylated in CHH context: %5.2f%c\n", (metil_file->CHH_methyl + metil_file->CHH_unmethyl == 0) ? 0.0 :
	  (float) metil_file->CHH_methyl / (metil_file->CHH_methyl + metil_file->CHH_unmethyl)  * 100, '%');
  //
  //--------------------------------------------------------------------------------------
  


  metil_file_free(metil_file);
  //
  //
  // end of workflow management
  //--------------------------------------------------------------------------------------
  
  
  //closing files
  if (options->pair_mode == SINGLE_END_MODE) {
    fastq_fclose(reader_input.fq_file1);
  } else {
    fastq_fclose(reader_input.fq_file1);
    fastq_fclose(reader_input.fq_file2);
  }
  bam_fclose(writer_input.bam_file);
  
  free(output_filename);
    genome_free(genome);
    bwt_index_free(bwt_index1);
    genome_free(genome1);
    bwt_index_free(bwt_index2);
    genome_free(genome2);
  bwt_optarg_free(bwt_optarg);
  cal_optarg_free(cal_optarg);
  pair_mng_free(pair_mng);
  report_optarg_free(report_optarg);

  
  if (0) {
    size_t total_item = 0;
    double max_time = 0, total_throughput = 0;
    printf("\nBWT time:\n");
    for (int i = 0; i < options->num_cpu_threads; i++) {
      printf("\tThread %d: %0.4f s (%d batches, %d items) -> %0.2f BWT/s (reads)\n", 
	     i, bwt_time[i] / 1e6, thr_batches[i], thr_bwt_items[i], 1e6 * thr_bwt_items[i] / bwt_time[i]);
      total_item += thr_bwt_items[i];
      total_throughput += (1e6 * thr_bwt_items[i] / bwt_time[i]);
      if (max_time < bwt_time[i]) max_time = bwt_time[i];
    }
    printf("\n\tTotal BWTs: %lu, Max time = %0.4f, Throughput = %0.2f BWT/s\n", 
	   total_item, max_time / 1e6, total_throughput);
    
    total_item = 0; max_time = 0; total_throughput = 0;
    printf("\nSeeding time:\n");
    for (int i = 0; i < options->num_cpu_threads; i++) {
      printf("\tThread %d: %0.4f s (%d batches, %d items) -> %0.2f BWT/s (seeds)\n", 
	     i, seeding_time[i] / 1e6, thr_batches[i], thr_seeding_items[i], 
	     1e6 * thr_seeding_items[i] / seeding_time[i]);
      total_item += thr_seeding_items[i];
      total_throughput += (1e6 * thr_seeding_items[i] / seeding_time[i]);
      if (max_time < seeding_time[i]) max_time = seeding_time[i];
    }
    printf("\n\tTotal BWTs: %lu, Max time = %0.4f, Throughput = %0.2f BWT/s\n", 
	   total_item, max_time / 1e6, total_throughput);
    
    total_item = 0; max_time = 0; total_throughput = 0;
    printf("\nCAL time:\n");
    for (int i = 0; i < options->num_cpu_threads; i++) {
      printf("\tThread %d: %0.4f s (%d batches, %d items) -> %0.2f CAL/s)\n", 
	     i, cal_time[i] / 1e6, thr_batches[i], thr_cal_items[i], 1e6 * thr_cal_items[i] / cal_time[i]);
      total_item += thr_cal_items[i];
      total_throughput += (1e6 * thr_cal_items[i] / cal_time[i]);
      if (max_time < cal_time[i]) max_time = cal_time[i];
    }
    printf("\n\tTotal CALs: %lu, Max time = %0.4f, Throughput = %0.2f CAL/s\n", 
	   total_item, max_time / 1e6, total_throughput);
    
    total_item = 0; max_time = 0; total_throughput = 0;
    printf("\nSW time:\n");
    for (int i = 0; i < options->num_cpu_threads; i++) {
      printf("\tThread %d: %0.4f s (%d batches, %d items) -> %0.2f SW/s)\n", 
	     i, sw_time[i] / 1e6, thr_batches[i], thr_sw_items[i], 1e6 * thr_sw_items[i] / sw_time[i]);
      total_item += thr_sw_items[i];
      total_throughput += (1e6 * thr_sw_items[i] / sw_time[i]);
      if (max_time < sw_time[i]) max_time = sw_time[i];
    }
    printf("\n\tTotal SWs: %lu, Max time = %0.4f, Throughput = %0.2f SW/s\n", 
	   total_item, max_time / 1e6, total_throughput);
  }
}
