#include "rna_aligner.h"

//#include "extrae_user_events.h" 

#define NUM_SECTIONS_TIME 		8

//--------------------------------------------------------------------
// workflow input                                                                                        //--------------------------------------------------------------------  




void *write_sam_header_BWT(genome_t *genome, FILE *f) {
  fprintf(f, "@PG\tID:HPG-Aligner\tVN:2.0\n");
  for (int i = 0; i < genome->num_chromosomes; i++) {
    fprintf(f, "@SQ\tSN:%s\tLN:%lu\n", genome->chr_name[i], genome->chr_size[i] + 1);
    //printf("%iName %s\n", i, genome->chr_name[i]);
  }
  //fclose(f);
  //exit(-1);
}


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



//--------------------------------------------------------------------------------------

void sa_index3_parallel_genome_new(char *sa_index_dirname, int num_threads,
				   sa_index3_t **sa_index_out, genome_t **genome_out) {  

  FILE *f_tab;
  char line[1024], filename_tab[strlen(sa_index_dirname) + 1024];
  char *prefix;
  uint k_value, pre_length, A_items, IA_items, num_suffixes, genome_len, num_chroms, num_items;

  struct timeval stop, start, end;

  PREFIX_TABLE_NT_VALUE['A'] = 0;
  PREFIX_TABLE_NT_VALUE['N'] = 0;
  PREFIX_TABLE_NT_VALUE['C'] = 1;
  PREFIX_TABLE_NT_VALUE['G'] = 2;
  PREFIX_TABLE_NT_VALUE['T'] = 3;

  sprintf(filename_tab, "%s/params.txt", sa_index_dirname, prefix);
  //printf("reading %s\n", filename_tab);

  f_tab = fopen(filename_tab, "r");
  if (!f_tab) {
    fprintf(stderr, "Error opening file %s!\n", filename_tab);
    exit(-1);
  }

  // prefix
  fgets(line, 1024, f_tab);
  line[strlen(line) - 1] = 0;
  prefix = strdup(line);
  // k_value
  fgets(line, 1024, f_tab);
  k_value = atoi(line);
  // pre_length
  fgets(line, 1024, f_tab);
  pre_length = atoi(line);
  // A_items
  fgets(line, 1024, f_tab);
  A_items = atoi(line);
  // IA_items
  fgets(line, 1024, f_tab);
  IA_items = atol(line);
  // num_suffixes
  fgets(line, 1024, f_tab);
  num_suffixes = atoi(line);
  // genome_length
  fgets(line, 1024, f_tab);
  genome_len = atoi(line);
  // num_chroms
  fgets(line, 1024, f_tab);
  num_chroms = atoi(line);

  size_t *chrom_lengths = (size_t *) malloc(num_chroms * sizeof(size_t));
  char **chrom_names = (char **) malloc(num_chroms * sizeof(char *));
  char chrom_name[1024];
  size_t chrom_len;
		  
  for (int i = 0; i < num_chroms; i++) {
    fgets(line, 1024, f_tab);
    sscanf(line, "%s %lu\n", chrom_name, &chrom_len);
    //printf("chrom_name: %s, chrom_len: %lu\n", chrom_name, chrom_len);
    chrom_names[i] = strdup(chrom_name);
    chrom_lengths[i] = chrom_len;
  }

  fclose(f_tab);

  char *S, *CHROM;
  unsigned char *JA;
  sa_genome3_t *genome;
  uint *SA, *PRE, *A, *IA;
  genome_t *genome_;

  //printf("Parametro: %i num threads \n", num_threads);

  // Compressed Row Storage (IA table)
  //struct timeval start_f, end_f;
  //double time_f = 0;

  int final_threads = num_threads >= 2 ? 2 : num_threads;

  {
    #pragma omp parallel sections num_threads(final_threads)
    {
      #pragma omp section
      {
	struct timeval stop, start;
	FILE *f_tab;
	char filename_tab[strlen(sa_index_dirname) + 1024];

	// read S from file
	sprintf(filename_tab, "%s/%s.S", sa_index_dirname, prefix);
	f_tab = fopen(filename_tab, "rb");
	if (f_tab == NULL) {
	  printf("Error: could not open %s to write\n", filename_tab);
	  exit(-1);
	}
	//  printf("genome: filename %s, length = %lu\n", filename_tab, genome_len);
	S = (char *) malloc(genome_len);
  
	//printf("\nreading S from file %s...\n", filename_tab);
	gettimeofday(&start, NULL);
	num_items = fread(S, sizeof(char), genome_len, f_tab);
	if (num_items != genome_len) {
	  printf("Error: (%s) mismatch num_items = %lu vs length = %lu\n", 
		 filename_tab, num_items, genome_len);
	  exit(-1);
	}
	gettimeofday(&stop, NULL);
	//printf("end of reading S (%lu len) from file %s in %0.2f s\n", 
	//	     genome_len, filename_tab,
	//	     (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);      
	fclose(f_tab);
      
	genome = sa_genome3_new(genome_len, num_chroms, 
				chrom_lengths, chrom_names, S);
      
	for (size_t i = 0; i < genome->length; i++) {
	  if (genome->S[i] == 'N' || genome->S[i] == 'n') {
	    genome->S[i] = 'A';
	  }
	}

	// Compressed Row Storage (A table)
	A = NULL;
	sprintf(filename_tab, "%s/%s.A", sa_index_dirname, prefix);
	f_tab = fopen(filename_tab, "rb");
	if (f_tab) {
	  A = (uint *) malloc(A_items * sizeof(uint));
	
	  //	printf("\nreading A table (Compression Row Storage) from file %s...\n", filename_tab);
	  gettimeofday(&start, NULL);
	  if ((num_items = fread(A, sizeof(uint), A_items, f_tab)) != A_items) {
	    printf("Error: (%s) mismatch read num_items = %lu (it must be %lu)\n", 
		   filename_tab, num_items, A_items);
	    exit(-1);
	  }
	  gettimeofday(&stop, NULL);
	  //	printf("end of reading A table (%lu num_items) from file %s in %0.2f s\n", 
	  //	       num_items, filename_tab,
	  //	       (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
	  fclose(f_tab);
	}

	// Compressed Row Storage (IA table)
	IA = NULL;
	sprintf(filename_tab, "%s/%s.IA", sa_index_dirname, prefix);
	f_tab = fopen(filename_tab, "rb");
	if (f_tab) {
	  IA = (uint *) malloc(IA_items * sizeof(uint));
	
	  //	printf("\nreading IA table (Compression Row Storage) from file %s...\n", filename_tab);
	  gettimeofday(&start, NULL);
	  if ((num_items = fread(IA, sizeof(uint), IA_items, f_tab)) != IA_items) {
	    printf("Error: (%s) mismatch read num_items = %lu (it must be %lu)\n", 
		   filename_tab, num_items, IA_items);
	    exit(-1);
	  }
	  gettimeofday(&stop, NULL);
	  //	printf("end of reading IA table (%lu num_items) from file %s in %0.2f s\n", 
	  //	       num_items, filename_tab,
	  //	       (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
	  fclose(f_tab);
	}
      }

      #pragma omp section
      {
	struct timeval stop, start;
	FILE *f_tab;
	char filename_tab[strlen(sa_index_dirname) + 1024];

	// read SA table from file
	sprintf(filename_tab, "%s/%s.SA", sa_index_dirname, prefix);
	f_tab = fopen(filename_tab, "rb");
	if (f_tab == NULL) {
	  printf("Error: could not open %s to write\n", filename_tab);
	  exit(-1);
	}
	//      printf("SA: filename %s, num_suffixes = %lu\n", filename_tab, num_suffixes);
	SA = (uint *) malloc(num_suffixes * sizeof(uint));
  
	//      printf("\nreading SA table from file %s...\n", filename_tab);
	gettimeofday(&start, NULL);
	num_items = fread(SA, sizeof(uint), num_suffixes, f_tab);
	if (num_items != num_suffixes) {
	  printf("Error: (%s) mismatch num_items = %lu vs num_suffixes = %lu\n", 
		 filename_tab, num_items, num_suffixes);
	  exit(-1);
	}
	gettimeofday(&stop, NULL);
	//      printf("end of reading SA table (%lu num_suffixes) from file %s in %0.2f s\n", 
	//	     num_suffixes, filename_tab,
	//	     (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);      
	fclose(f_tab);

	// read CHROM table from file
	sprintf(filename_tab, "%s/%s.CHROM", sa_index_dirname, prefix);
	f_tab = fopen(filename_tab, "rb");
	if (f_tab == NULL) {
	  printf("Error: could not open %s to write\n", filename_tab);
	  exit(-1);
	}
	//      printf("CHROM: filename %s, num_suffixes = %lu\n", filename_tab, num_suffixes);
	CHROM = (char *) malloc(num_suffixes * sizeof(char));
      
	//      printf("\nreading CHROM table from file %s...\n", filename_tab);
	gettimeofday(&start, NULL);
	num_items = fread(CHROM, sizeof(char), num_suffixes, f_tab);
	if (num_items != num_suffixes) {
	  printf("Error: (%s) mismatch num_items = %lu vs num_suffixes = %lu\n", 
		 filename_tab, num_items, num_suffixes);
	  exit(-1);
	}
	gettimeofday(&stop, NULL);
	//      printf("end of reading CHROM table (%lu num_suffixes) from file %s in %0.2f s\n", 
	//	     num_suffixes, filename_tab,
	//	     (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);      
	fclose(f_tab);

	// read PRE table from file
	PRE = NULL;
	sprintf(filename_tab, "%s/%s.PRE", sa_index_dirname, prefix);
	f_tab = fopen(filename_tab, "rb");
	if (f_tab) {
	  num_items = 1LLU << (2 * k_value);
	  assert(num_items == pre_length);
	  PRE = (uint *) malloc(num_items * sizeof(uint));
	
	  //	printf("\nreading PRE table from file %s...\n", filename_tab);
	  gettimeofday(&start, NULL);
	  if (num_items != fread(PRE, sizeof(uint), num_items, f_tab)) {
	    printf("Error: (%s) mismatch num_items = %lu\n", 
		   filename_tab, num_items);
	    exit(-1);
	  }
	  gettimeofday(&stop, NULL);
	  //	printf("end of reading PRE table (%lu num_items) from file %s in %0.2f s\n", 
	  //	       num_items, filename_tab,
	  //	       (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
	  fclose(f_tab);
	}
   
	// Compressed Row Storage (JA table)
	JA = NULL;
	sprintf(filename_tab, "%s/%s.JA", sa_index_dirname, prefix);
	f_tab = fopen(filename_tab, "rb");
	if (f_tab) {
	  JA = (unsigned char *) malloc(A_items * sizeof(unsigned char));
	
	  //	printf("\nreading JA table (Compression Row Storage) from file %s...\n", filename_tab);
	  gettimeofday(&start, NULL);
	  if ((num_items = fread(JA, sizeof(unsigned char), A_items, f_tab)) != A_items) {
	    printf("Error: (%s) mismatch read num_items = %lu (it must be %lu)\n", 
		   filename_tab, num_items, A_items);
	    exit(-1);
	  }
	  gettimeofday(&stop, NULL);
	  //	printf("end of reading JA table (%lu num_items) from file %s in %0.2f s\n", 
	  //	       num_items, filename_tab,
	  //	       (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
	  fclose(f_tab);
	}

	genome_ = genome_new("dna_compression.bin", sa_index_dirname, SA_MODE);	

	genome_->num_chromosomes = num_chroms;
	genome_->chr_name = (char **) calloc(genome_->num_chromosomes, sizeof(char *));
	genome_->chr_size = (size_t *) calloc(genome_->num_chromosomes, sizeof(size_t));
	genome_->chr_offset = (size_t *) calloc(genome_->num_chromosomes, sizeof(size_t));
	size_t offset = 0;
	
	for (int c = 0; c < num_chroms; c++) {
	  genome_->chr_size[c] = chrom_lengths[c];
	  genome_->chr_name[c] = strdup(chrom_names[c]);
	  genome_->chr_offset[c] = offset;
	  offset += genome_->chr_size[c];
	}
	
      }
    }  
  }

  goto final;

  //======================================================================
  {
    #pragma omp parallel sections num_threads(final_threads)
    {
      #pragma omp section
      {
	struct timeval stop, start;
	FILE *f_tab;
	char filename_tab[strlen(sa_index_dirname) + 1024];

	// read S from file
	sprintf(filename_tab, "%s/%s.S", sa_index_dirname, prefix);
	f_tab = fopen(filename_tab, "rb");
	if (f_tab == NULL) {
	  printf("Error: could not open %s to write\n", filename_tab);
	  exit(-1);
	}
	// printf("genome: filename %s, length = %lu\n", filename_tab, genome_len);
	S = (char *) malloc(genome_len);
  
	//printf("\nreading S from file %s...\n", filename_tab);
	gettimeofday(&start, NULL);
	num_items = fread(S, sizeof(char), genome_len, f_tab);
	if (num_items != genome_len) {
	  printf("Error: (%s) mismatch num_items = %lu vs length = %lu\n",
		 filename_tab, num_items, genome_len);
	  exit(-1);
	}
	gettimeofday(&stop, NULL);
	//printf("end of reading S (%lu len) from file %s in %0.2f s\n",
	// genome_len, filename_tab,
	// (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);
	fclose(f_tab);
      
	genome = sa_genome3_new(genome_len, num_chroms,
				chrom_lengths, chrom_names, S);
      
	for (size_t i = 0; i < genome->length; i++) {
	  if (genome->S[i] == 'N' || genome->S[i] == 'n') {
	    genome->S[i] = 'A';
	  }
	}


	// Compressed Row Storage (IA table)
	IA = NULL;
	sprintf(filename_tab, "%s/%s.IA", sa_index_dirname, prefix);
	f_tab = fopen(filename_tab, "rb");
	if (f_tab) {
	  IA = (uint *) malloc(IA_items * sizeof(uint));

	  // printf("\nreading IA table (Compression Row Storage) from file %s...\n", filename_tab);
	  gettimeofday(&start, NULL);
	  if ((num_items = fread(IA, sizeof(uint), IA_items, f_tab)) != IA_items) {
	    printf("Error: (%s) mismatch read num_items = %lu (it must be %lu)\n",
		   filename_tab, num_items, IA_items);
	    exit(-1);
	  }
	  gettimeofday(&stop, NULL);
	  // printf("end of reading IA table (%lu num_items) from file %s in %0.2f s\n",
	  // num_items, filename_tab,
	  // (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);
	  fclose(f_tab);
	}
      }
      #pragma omp section
      {
	struct timeval stop, start;
	FILE *f_tab;
	char filename_tab[strlen(sa_index_dirname) + 1024];

	// Compressed Row Storage (A table)
	A = NULL;
	sprintf(filename_tab, "%s/%s.A", sa_index_dirname, prefix);
	f_tab = fopen(filename_tab, "rb");
	if (f_tab) {
	  A = (uint *) malloc(A_items * sizeof(uint));

	  // printf("\nreading A table (Compression Row Storage) from file %s...\n", filename_tab);
	  gettimeofday(&start, NULL);
	  if ((num_items = fread(A, sizeof(uint), A_items, f_tab)) != A_items) {
	    printf("Error: (%s) mismatch read num_items = %lu (it must be %lu)\n",
		   filename_tab, num_items, A_items);
	    exit(-1);
	  }
	  gettimeofday(&stop, NULL);
	  // printf("end of reading A table (%lu num_items) from file %s in %0.2f s\n",
	  // num_items, filename_tab,
	  // (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);
	  fclose(f_tab);
	}

      }
      #pragma omp section
      {
	struct timeval stop, start;
	FILE *f_tab;
	char filename_tab[strlen(sa_index_dirname) + 1024];

	// read SA table from file
	sprintf(filename_tab, "%s/%s.SA", sa_index_dirname, prefix);
	f_tab = fopen(filename_tab, "rb");
	if (f_tab == NULL) {
	  printf("Error: could not open %s to write\n", filename_tab);
	  exit(-1);
	}
	// printf("SA: filename %s, num_suffixes = %lu\n", filename_tab, num_suffixes);
	SA = (uint *) malloc(num_suffixes * sizeof(uint));
  
	// printf("\nreading SA table from file %s...\n", filename_tab);
	gettimeofday(&start, NULL);
	num_items = fread(SA, sizeof(uint), num_suffixes, f_tab);
	if (num_items != num_suffixes) {
	  printf("Error: (%s) mismatch num_items = %lu vs num_suffixes = %lu\n",
		 filename_tab, num_items, num_suffixes);
	  exit(-1);
	}
	gettimeofday(&stop, NULL);
	// printf("end of reading SA table (%lu num_suffixes) from file %s in %0.2f s\n",
	// num_suffixes, filename_tab,
	// (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);
	fclose(f_tab);

	genome_ = genome_new("dna_compression.bin", sa_index_dirname, SA_MODE);	
	genome_->num_chromosomes = num_chroms;
	genome_->chr_name = (char **) calloc(genome_->num_chromosomes, sizeof(char *));
	genome_->chr_size = (size_t *) calloc(genome_->num_chromosomes, sizeof(size_t));
	genome_->chr_offset = (size_t *) calloc(genome_->num_chromosomes, sizeof(size_t));
	size_t offset = 0;
	
	for (int c = 0; c < num_chroms; c++) {
	  genome_->chr_size[c] = chrom_lengths[c];
	  genome_->chr_name[c] = chrom_names[c];
	  genome_->chr_offset[c] = offset;
	  offset += genome_->chr_size[c];
	}
	
	//stop_timer(start, end, time);
	//printf("************* Tiempo GENOME : %f\n", time / 1000000);
	
      }

      #pragma omp section 
      {
	struct timeval stop, start;
	FILE *f_tab;
	char filename_tab[strlen(sa_index_dirname) + 1024];

	// read CHROM table from file
	sprintf(filename_tab, "%s/%s.CHROM", sa_index_dirname, prefix);
	f_tab = fopen(filename_tab, "rb");
	if (f_tab == NULL) {
	  printf("Error: could not open %s to write\n", filename_tab);
	  exit(-1);
	}
	// printf("CHROM: filename %s, num_suffixes = %lu\n", filename_tab, num_suffixes);
	CHROM = (char *) malloc(num_suffixes * sizeof(char));
      
	// printf("\nreading CHROM table from file %s...\n", filename_tab);
	gettimeofday(&start, NULL);
	num_items = fread(CHROM, sizeof(char), num_suffixes, f_tab);
	if (num_items != num_suffixes) {
	  printf("Error: (%s) mismatch num_items = %lu vs num_suffixes = %lu\n",
		 filename_tab, num_items, num_suffixes);
	  exit(-1);
	}
	gettimeofday(&stop, NULL);
	// printf("end of reading CHROM table (%lu num_suffixes) from file %s in %0.2f s\n",
	// num_suffixes, filename_tab,
	// (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);
	fclose(f_tab);

	// read PRE table from file
	PRE = NULL;
	sprintf(filename_tab, "%s/%s.PRE", sa_index_dirname, prefix);
	f_tab = fopen(filename_tab, "rb");
	if (f_tab) {
	  num_items = 1LLU << (2 * k_value);
	  assert(num_items == pre_length);
	  PRE = (uint *) malloc(num_items * sizeof(uint));

	  // printf("\nreading PRE table from file %s...\n", filename_tab);
	  gettimeofday(&start, NULL);
	  if (num_items != fread(PRE, sizeof(uint), num_items, f_tab)) {
	    printf("Error: (%s) mismatch num_items = %lu\n",
		   filename_tab, num_items);
	    exit(-1);
	  }
	  gettimeofday(&stop, NULL);
	  // printf("end of reading PRE table (%lu num_items) from file %s in %0.2f s\n",
	  // num_items, filename_tab,
	  // (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);
	  fclose(f_tab);
	}
   
	// Compressed Row Storage (JA table)
	JA = NULL;
	sprintf(filename_tab, "%s/%s.JA", sa_index_dirname, prefix);
	f_tab = fopen(filename_tab, "rb");
	if (f_tab) {
	  JA = (unsigned char *) malloc(A_items * sizeof(unsigned char));

	  // printf("\nreading JA table (Compression Row Storage) from file %s...\n", filename_tab);
	  gettimeofday(&start, NULL);
	  if ((num_items = fread(JA, sizeof(unsigned char), A_items, f_tab)) != A_items) {
	    printf("Error: (%s) mismatch read num_items = %lu (it must be %lu)\n",
		   filename_tab, num_items, A_items);
	    exit(-1);
	  }
	  gettimeofday(&stop, NULL);
	  // printf("end of reading JA table (%lu num_items) from file %s in %0.2f s\n",
	  // num_items, filename_tab,
	  // (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);
	  fclose(f_tab);
	}

      }
    }
  }
  
  goto final;

  {
    
    struct timeval start_f, end_f;
    double time_f;

    start_timer(start_f); 
    
    
    // read SA table from file
    FILE *f_tab, *f_tab2;
    char filename_tab[strlen(sa_index_dirname) + 1024];
    struct timeval end, start;
    double time;

    uint num_items0, num_items1, num_items2, num_items3;
    
    start_timer(start); 
    sprintf(filename_tab, "%s/%s.SA", sa_index_dirname, prefix);
    printf("SA: filename %s, num_suffixes = %lu\n", filename_tab, num_suffixes);
    SA = (uint *) malloc(num_suffixes * sizeof(uint));      
    #pragma omp parallel sections num_threads(2)
    {
        #pragma omp section
        {
	  f_tab = fopen(filename_tab, "rb");
	  size_t offset = sizeof(uint)*(num_suffixes/2);
	  uint extra = num_suffixes % 2;
	  fseek(f_tab, offset, SEEK_SET);

	  num_items0 = fread(&SA[num_suffixes/2], sizeof(uint), num_suffixes/2 + extra, f_tab);
	}
        #pragma omp section
        {
	  f_tab2 = fopen(filename_tab, "rb");
	  num_items1 = fread(SA, sizeof(uint), num_suffixes/2, f_tab2);
	}
    }
    fclose(f_tab);
    fclose(f_tab2);
    stop_timer(start, end, time);
    printf("****************Tiempo SA (%lu =? %lu): %f\n", num_items0 + num_items1, num_suffixes, time / 1000000);            


    start_timer(start); 
    sprintf(filename_tab, "%s/%s.A", sa_index_dirname, prefix);
    A = (uint *) malloc(A_items * sizeof(uint));
    #pragma omp parallel sections num_threads(2)
    {
        #pragma omp section
        {
	  f_tab = fopen(filename_tab, "rb");
	  size_t offset = sizeof(uint)*(A_items/2);
	  uint extra = A_items % 2;
	  fseek(f_tab, offset, SEEK_CUR);
	  num_items0 = fread(&A[A_items/2], sizeof(uint), A_items/2 + extra, f_tab);
	}
        #pragma omp section
        {
	  f_tab2 = fopen(filename_tab, "rb");
	  num_items1 = fread(A, sizeof(uint), A_items/2, f_tab2);
	}
    }	
    fclose(f_tab);    
    stop_timer(start, end, time);
    printf("**************Tiempo A (%lu =? %lu): %f\n", num_items0 + num_items1, A_items, time / 1000000);


    
    //exit(-1);
    
    char filename_tab_2[2048];
    sprintf(filename_tab, "%s/%s.S", sa_index_dirname, prefix);
    sprintf(filename_tab_2, "%s/%s.CHROM", sa_index_dirname, prefix);

    S = (char *) malloc(genome_len);
    CHROM = (char *) malloc(num_suffixes);
    #pragma omp parallel sections num_threads(4)
    {
      #pragma omp section
      {
	  f_tab = fopen(filename_tab, "rb");
	  size_t offset = sizeof(char)*(genome_len/2);
	  uint extra = genome_len % 2;
	  fseek(f_tab, offset, SEEK_CUR);
	  num_items0 = fread(&S[genome_len/2], sizeof(char), genome_len/2 + extra, f_tab);
      }

      #pragma omp section
      {
	f_tab2 = fopen(filename_tab, "rb");
	num_items1 = fread(S, sizeof(char), genome_len/2, f_tab2);
      }
      
      #pragma omp section
      {
	FILE *f_tab3 = fopen(filename_tab_2, "rb");
	size_t offset = sizeof(char)*(num_suffixes/2);
	uint extra = num_suffixes % 2;
	fseek(f_tab, offset, SEEK_CUR);
	num_items2 = fread(&CHROM[num_suffixes/2], sizeof(char), num_suffixes/2 + extra, f_tab3);
      }
      
      #pragma omp section
      {
	FILE *f_tab4 = fopen(filename_tab_2, "rb");
	num_items3 = fread(CHROM, sizeof(char), num_suffixes/2, f_tab4);
      }	
    }
    fclose(f_tab);
    stop_timer(start, end, time);
    printf("**************** Tiempo S (%lu =? %lu) | CHROM (%lu =? %lu): %f\n",  num_items0 + num_items1, genome_len,  num_items2 + num_items3, num_suffixes, time / 1000000);



    

    start_timer(start); 
    sprintf(filename_tab, "%s/%s.JA", sa_index_dirname, prefix);
    JA = (unsigned char *) malloc(A_items * sizeof(unsigned char));
    #pragma omp parallel sections num_threads(2)
    {
        #pragma omp section
        {
	  f_tab = fopen(filename_tab, "rb");
	  size_t offset = sizeof(unsigned char)*(A_items/2);
	  uint extra = A_items % 2;
	  fseek(f_tab, offset, SEEK_CUR);
	  num_items0 = fread(&JA[A_items/2], sizeof(unsigned char), A_items/2 + extra, f_tab);
	}
        #pragma omp section
        {
	  f_tab2 = fopen(filename_tab, "rb");
	  num_items1 = fread(JA, sizeof(unsigned char), A_items/2, f_tab2);
	}
    }	
    fclose(f_tab);    
    stop_timer(start, end, time);
    printf("**************Tiempo JA (%lu =? %lu): %f\n", num_items0 + num_items1, A_items, time / 1000000);




    
    start_timer(start); 
    IA = NULL;
    sprintf(filename_tab, "%s/%s.IA", sa_index_dirname, prefix);
    f_tab = fopen(filename_tab, "rb");
    IA = (uint *) malloc(IA_items * sizeof(uint));
    #pragma omp parallel sections num_threads(2)
    {
        #pragma omp section
        {
	  f_tab = fopen(filename_tab, "rb");
	  size_t offset = sizeof(uint)*(IA_items/2);
	  uint extra = IA_items % 2;
	  fseek(f_tab, offset, SEEK_CUR);
	  num_items0 = fread(&IA[IA_items/2], sizeof(uint), IA_items/2 + extra, f_tab);
	}
        #pragma omp section
        {
	  f_tab2 = fopen(filename_tab, "rb");
	  num_items1 = fread(IA, sizeof(uint), IA_items/2, f_tab2);
	}
    }
    fclose(f_tab);    
    stop_timer(start, end, time);
    printf("**************Tiempo IA (%lu =? %lu): %f\n", num_items0 + num_items1, IA_items, time / 1000000);
	



    #pragma omp parallel sections num_threads(3)
    {
        #pragma omp section
        {
	  PRE = NULL;
	  sprintf(filename_tab, "%s/%s.PRE", sa_index_dirname, prefix);
	  f_tab = fopen(filename_tab, "rb");
	  if (f_tab) {
	    num_items = 1LLU << (2 * k_value);
	    assert(num_items == pre_length);
	    PRE = (uint *) malloc(num_items * sizeof(uint));
	    if (num_items != fread(PRE, sizeof(uint), num_items, f_tab)) {
	      printf("Error: (%s) mismatch num_items = %lu\n", 
		     filename_tab, num_items);
	      exit(-1);
	    }
	    fclose(f_tab);
	  }
	}
        #pragma omp section
        {
	  genome = sa_genome3_new(genome_len, num_chroms, 
				  chrom_lengths, chrom_names, S);
	  
	  for (size_t i = 0; i < genome->length; i++) {
	    if (genome->S[i] == 'N' || genome->S[i] == 'n') {
	      genome->S[i] = 'A';
	    }
	  }
	}	
        #pragma omp section
        {
	  genome_ = genome_new("dna_compression.bin", sa_index_dirname, SA_MODE);
      
	  genome_->num_chromosomes = num_chroms;
	  genome_->chr_name = (char **) calloc(genome_->num_chromosomes, sizeof(char *));
	  genome_->chr_size = (size_t *) calloc(genome_->num_chromosomes, sizeof(size_t));
	  genome_->chr_offset = (size_t *) calloc(genome_->num_chromosomes, sizeof(size_t));
	  size_t offset = 0;
      
	  for (int c = 0; c < num_chroms; c++) {
	    genome_->chr_size[c] = chrom_lengths[c];
	    genome_->chr_name[c] = chrom_names[c];
	    genome_->chr_offset[c] = offset;
	    offset += genome_->chr_size[c];
	  }

	  stop_timer(start, end, time);
	  //printf("************* Tiempo GENOME : %f\n", time / 1000000);

	}
    }
    
    stop_timer(start_f, end_f, time_f);
    printf("FINAL TIME : %f\n", time_f / 1000000);
    
  }
  //======================================================================

 final:
  free(prefix);

  {
    // creating the sa_index_t structure
    sa_index3_t *p = (sa_index3_t *) malloc(sizeof(sa_index3_t));

    p->num_suffixes = num_suffixes;
    p->prefix_length = pre_length;
    p->A_items = A_items;
    p->IA_items = IA_items;
    p->k_value = k_value;
    p->SA = SA;
    p->CHROM = CHROM;
    p->PRE = PRE;
    p->A = A;
    p->IA = IA;
    p->JA = JA;
    p->genome = genome;
    
    *sa_index_out = p;
    *genome_out = genome_;
    
    return;
    
  }

  exit(-1);

    /*    
    start_timer(start); 
    sprintf(filename_tab, "%s/%s.CHROM", sa_index_dirname, prefix);
    CHROM = (char *) malloc(num_suffixes * sizeof(char));
    #pragma omp parallel sections num_threads(2)
    {
        #pragma omp section
        {
	  f_tab = fopen(filename_tab, "rb");
	  int offset = sizeof(char)*(num_suffixes/2);
	  int extra = num_suffixes % 2;
	  fseek(f_tab, offset, SEEK_CUR);
	  fread(&CHROM[num_suffixes/2], sizeof(char), num_suffixes/2 + extra, f_tab);
	}
        #pragma omp section
        {
	  f_tab2 = fopen(filename_tab, "rb");
	  fread(CHROM, sizeof(char), num_suffixes/2, f_tab2);
	}
    }    
    fclose(f_tab);
    stop_timer(start, end, time);
    printf("**************Tiempo CHROM : %f\n", time / 1000000);        
    */

  //start_timer(start_f); 
  #pragma omp parallel sections num_threads(final_threads)
  {
    //printf("Parallel section %i\n", omp_get_num_threads());
    #pragma omp section
    {
      FILE *f_tab;
      char filename_tab[strlen(sa_index_dirname) + 1024];
      struct timeval end, start;      
      double time;  
      time = 0;
      start_timer(start); 

      IA = NULL;
      sprintf(filename_tab, "%s/%s.IA", sa_index_dirname, prefix);
      f_tab = fopen(filename_tab, "rb");
      if (f_tab) {
	IA = (uint *) malloc(IA_items * sizeof(uint));
	
	//	printf("\nreading IA table (Compression Row Storage) from file %s...\n", filename_tab);
	//gettimeofday(&start, NULL);
	if ((num_items = fread(IA, sizeof(uint), IA_items, f_tab)) != IA_items) {
	  printf("Error: (%s) mismatch read num_items = %lu (it must be %lu)\n", 
		 filename_tab, num_items, IA_items);
	  exit(-1);
	}
	//gettimeofday(&stop, NULL);
	//	printf("end of reading IA table (%lu num_items) from file %s in %0.2f s\n", 
	//	       num_items, filename_tab,
	//	       (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
	fclose(f_tab);
      }
      //stop_timer(start, end, time);
      //printf("Tiempo IA : %f\n", time / 1000000);


      // read CHROM table from file
      //FILE *f_tab;
      //char filename_tab[strlen(sa_index_dirname) + 1024];
      //struct timeval end, start;
      //double time;
      //time = 0;
      //start_timer(start); 

      sprintf(filename_tab, "%s/%s.CHROM", sa_index_dirname, prefix);
      f_tab = fopen(filename_tab, "rb");
      if (f_tab == NULL) {
	printf("Error: could not open %s to write\n", filename_tab);
	exit(-1);
      }
      //      printf("CHROM: filename %s, num_suffixes = %lu\n", filename_tab, num_suffixes);
      CHROM = (char *) malloc(num_suffixes * sizeof(char));
      
      //      printf("\nreading CHROM table from file %s...\n", filename_tab);
      //gettimeofday(&start, NULL);
      num_items = fread(CHROM, sizeof(char), num_suffixes, f_tab);
      if (num_items != num_suffixes) {
	printf("Error: (%s) mismatch num_items = %lu vs num_suffixes = %lu\n", 
	       filename_tab, num_items, num_suffixes);
	exit(-1);
      }
      //gettimeofday(&stop, NULL);
      //      printf("end of reading CHROM table (%lu num_suffixes) from file %s in %0.2f s\n", 
      //	     num_suffixes, filename_tab,
      //	     (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);      
      fclose(f_tab);
      //stop_timer(start, end,  time);
      //printf("Tiempo CHROM : %f\n", time / 1000000);


      // read PRE table from file
      //FILE *f_tab;
      //char filename_tab[strlen(sa_index_dirname) + 1024];
      //struct timeval end, start;
      //double time;
      //time = 0;
      //start_timer(start); 

      PRE = NULL;
      sprintf(filename_tab, "%s/%s.PRE", sa_index_dirname, prefix);
      f_tab = fopen(filename_tab, "rb");
      if (f_tab) {
	num_items = 1LLU << (2 * k_value);
	assert(num_items == pre_length);
	PRE = (uint *) malloc(num_items * sizeof(uint));
	
	//	printf("\nreading PRE table from file %s...\n", filename_tab);
	//gettimeofday(&start, NULL);
	if (num_items != fread(PRE, sizeof(uint), num_items, f_tab)) {
	  printf("Error: (%s) mismatch num_items = %lu\n", 
		 filename_tab, num_items);
	  exit(-1);
	}
	//gettimeofday(&stop, NULL);
	//	printf("end of reading PRE table (%lu num_items) from file %s in %0.2f s\n", 
	//	       num_items, filename_tab,
	//	       (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
	fclose(f_tab);
      }
      //stop_timer(start, end, time);
      //printf("Tiempo PRE : %f\n", time / 1000000);


      // Compressed Row Storage (JA table)
      //FILE *f_tab;
      //char filename_tab[strlen(sa_index_dirname) + 1024];
      //struct timeval end, start;
      //double time;
      //time = 0;
      //start_timer(start); 

      JA = NULL;
      sprintf(filename_tab, "%s/%s.JA", sa_index_dirname, prefix);
      f_tab = fopen(filename_tab, "rb");
      if (f_tab) {
	JA = (unsigned char *) malloc(A_items * sizeof(unsigned char));
	
	//	printf("\nreading JA table (Compression Row Storage) from file %s...\n", filename_tab);
	//gettimeofday(&start, NULL);
	if ((num_items = fread(JA, sizeof(unsigned char), A_items, f_tab)) != A_items) {
	  printf("Error: (%s) mismatch read num_items = %lu (it must be %lu)\n", 
		 filename_tab, num_items, A_items);
	  exit(-1);
	}
	//gettimeofday(&stop, NULL);
	//	printf("end of reading JA table (%lu num_items) from file %s in %0.2f s\n", 
	//	       num_items, filename_tab,
	//	       (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
	fclose(f_tab);
      }
      //stop_timer(start, end, time);
      //printf("Tiempo JA : %f\n", time / 1000000);


      //FILE *f_tab;
      //char filename_tab[strlen(sa_index_dirname) + 1024];
      //struct timeval end, start;
      //double time;
      //time = 0;
      //start_timer(start); 

      genome_ = genome_new("dna_compression.bin", sa_index_dirname, SA_MODE);
      
      genome_->num_chromosomes = num_chroms;
      genome_->chr_name = (char **) calloc(genome_->num_chromosomes, sizeof(char *));
      genome_->chr_size = (size_t *) calloc(genome_->num_chromosomes, sizeof(size_t));
      genome_->chr_offset = (size_t *) calloc(genome_->num_chromosomes, sizeof(size_t));
      size_t offset = 0;
      
      for (int c = 0; c < num_chroms; c++) {
	genome_->chr_size[c] = chrom_lengths[c];
	genome_->chr_name[c] = chrom_names[c];
	genome_->chr_offset[c] = offset;
	offset += genome_->chr_size[c];
      }
      stop_timer(start, end, time);
      printf("************* Tiempo GENOME : %f\n", time / 1000000);
    }

    #pragma omp section
    {
      //struct timeval stop, start;
      //FILE *f_tab;
      FILE *f_tab;
      char filename_tab[strlen(sa_index_dirname) + 1024];
      struct timeval end, start;
      double time;
      time = 0;
      start_timer(start); 
      
      // read S from file
      sprintf(filename_tab, "%s/%s.S", sa_index_dirname, prefix);
      f_tab = fopen(filename_tab, "rb");
      if (f_tab == NULL) {
	printf("Error: could not open %s to write\n", filename_tab);
	exit(-1);
      }
      //  printf("genome: filename %s, length = %lu\n", filename_tab, genome_len);
      S = (char *) malloc(genome_len);
  
      //printf("\nreading S from file %s...\n", filename_tab);
      //gettimeofday(&start, NULL);
      num_items = fread(S, sizeof(char), genome_len, f_tab);
      if (num_items != genome_len) {
	printf("Error: (%s) mismatch num_items = %lu vs length = %lu\n", 
	       filename_tab, num_items, genome_len);
	exit(-1);
      }
      //gettimeofday(&stop, NULL);
      //printf("end of reading S (%lu len) from file %s in %0.2f s\n", 
      //	     genome_len, filename_tab,
      //	     (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);      
      fclose(f_tab);
      
      genome = sa_genome3_new(genome_len, num_chroms, 
			      chrom_lengths, chrom_names, S);
      
      for (size_t i = 0; i < genome->length; i++) {
	if (genome->S[i] == 'N' || genome->S[i] == 'n') {
	  genome->S[i] = 'A';
	}
      }

      stop_timer(start, end, time);
      printf("**************** Tiempo S : %f\n", time / 1000000);

    }

    #pragma omp section
    {
      // Compressed Row Storage (A table)
      FILE *f_tab;
      char filename_tab[strlen(sa_index_dirname) + 1024];
      struct timeval end, start;
      double time;
      start_timer(start); 

      A = NULL;
      sprintf(filename_tab, "%s/%s.A", sa_index_dirname, prefix);
      f_tab = fopen(filename_tab, "rb");
      if (f_tab) {
	A = (uint *) malloc(A_items * sizeof(uint));
	
	//	printf("\nreading A table (Compression Row Storage) from file %s...\n", filename_tab);
	//gettimeofday(&start, NULL);
	if ((num_items = fread(A, sizeof(uint), A_items, f_tab)) != A_items) {
	  printf("Error: (%s) mismatch read num_items = %lu (it must be %lu)\n", 
		 filename_tab, num_items, A_items);
	  exit(-1);
	}
	//gettimeofday(&stop, NULL);
	//	printf("end of reading A table (%lu num_items) from file %s in %0.2f s\n", 
	//	       num_items, filename_tab,
	//	       (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
	fclose(f_tab);
      }
      stop_timer(start, end, time);
      printf("**************Tiempo A : %f\n", time / 1000000);
    }

    #pragma omp section
    {
      // read SA table from file
      FILE *f_tab;
      char filename_tab[strlen(sa_index_dirname) + 1024];
      struct timeval end, start;
      double time;
      start_timer(start); 
      //gettimeofday(&start, NULL);

      sprintf(filename_tab, "%s/%s.SA", sa_index_dirname, prefix);
      f_tab = fopen(filename_tab, "rb");
      if (f_tab == NULL) {
	printf("Error: could not open %s to write\n", filename_tab);
	exit(-1);
      }
      //      printf("SA: filename %s, num_suffixes = %lu\n", filename_tab, num_suffixes);
      SA = (uint *) malloc(num_suffixes * sizeof(uint));  
      //      printf("\nreading SA table from file %s...\n", filename_tab);
      num_items = fread(SA, sizeof(uint), num_suffixes, f_tab);
      if (num_items != num_suffixes) {
	printf("Error: (%s) mismatch num_items = %lu vs num_suffixes = %lu\n", 
	       filename_tab, num_items, num_suffixes);
	exit(-1);
      }
      stop_timer(start, end, time);
      printf("****************Tiempo SA : %f\n", time / 1000000);
      //      printf("end of reading SA table (%lu num_suffixes) from file %s in %0.2f s\n", 
      //	     num_suffixes, filename_tab,
      //	     (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);      
      fclose(f_tab);
    }
  }
  //stop_timer(start_f, end_f, time_f);
  //printf("tiempo final %f :)\n", time_f / 1000000);

  free(prefix);

  // creating the sa_index_t structure
  sa_index3_t *p = (sa_index3_t *) malloc(sizeof(sa_index3_t));

  p->num_suffixes = num_suffixes;
  p->prefix_length = pre_length;
  p->A_items = A_items;
  p->IA_items = IA_items;
  p->k_value = k_value;
  p->SA = SA;
  p->CHROM = CHROM;
  p->PRE = PRE;
  p->A = A;
  p->IA = IA;
  p->JA = JA;
  p->genome = genome;

  *sa_index_out = p;
  *genome_out = genome_;

  //exit(-1);

  return p;
}





void rna_aligner(options_t *options) {
  //End fill 

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
    if (options->bam_format) {
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
    //sa_index = sa_index3_new(options->bwt_dirname);
    start_timer(time_genome_s);
    sa_index3_parallel_genome_new(options->bwt_dirname, options->num_cpu_threads, &sa_index, &genome);
    
    /*
    sa_index3_display(sa_index);
    stop_timer(time_genome_s, time_genome_e, time_genome);
    printf("1.TOTAL Carga genoma: %f(s)\n", time_genome / 1000000);
    
    start_timer(time_genome_s);
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
    stop_timer(time_genome_s, time_genome_e, time_genome);
    printf("2.TOTAL Carga genoma: %f(s)\n", time_genome / 1000000);

    LOG_DEBUG("Done !!");
    */
  }

  num_chromosomes = genome->num_chromosomes;

  //start_timer(time_genome_s);
  // Metaexons structure
  metaexons = metaexons_new(genome->num_chromosomes, 
			    genome->chr_size);
  
  stop_timer(time_genome_s, time_genome_e, time_genome);
  
  printf("TOTAL Carga genoma: %f(s)\n", time_genome / 1000000);
  //exit(-1);

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
  batch_writer_input_init(output_filename,
			  exact_filename, 
			  extend_filename, 
			  alignments_list, 
			  genome, 
			  &writer_input);

  writer_input.bam_format = options->bam_format;
  if (options->bam_format) {
    bam_header_t *bam_header;
    if (options->fast_mode) {
      bam_header = create_bam_header(sa_index->genome);
    } else {
      bam_header = create_bam_header_by_genome(genome);
    }
    writer_input.bam_file = bam_fopen_mode(output_filename, bam_header, "w");
    bam_fwrite_header(bam_header, writer_input.bam_file);
  } else {
    writer_input.bam_file = (bam_file_t *) fopen(output_filename, "w"); 
    if (options->fast_mode) {
      write_sam_header(sa_index->genome, (FILE *) writer_input.bam_file);
    } else {
      write_sam_header_BWT(genome, (FILE *) writer_input.bam_file);
    }
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

  double time_total_1, time_total_2;
  struct timeval time_s1, time_e1, time_s2, time_e2;

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

      if (options->bam_format) {
	workflow_set_consumer((workflow_consumer_function_t *)bam_writer, "BAM writer", wf);
      } else {
	workflow_set_consumer((workflow_consumer_function_t *)sam_writer, "SAM writer", wf);
      }

      workflow_t *wf_last = workflow_new();
      workflow_stage_function_t stage_functions_last[] = {rna_last_stage, post_pair_stage};
      char *stage_labels_last[] = {"RNA LAST STAGE", "POST PAIR"};
      workflow_set_stages(2, (workflow_stage_function_t *)&stage_functions_last, stage_labels_last, wf_last);
      workflow_set_producer((workflow_producer_function_t *)file_reader, "Buffer reader", wf_last);

      if (options->bam_format) {
	workflow_set_consumer((workflow_consumer_function_t *)bam_writer, "BAM writer", wf_last);
      } else {
	workflow_set_consumer((workflow_consumer_function_t *)sam_writer, "SAM writer", wf_last);
      }


      workflow_t *wf_hc = workflow_new();
      workflow_stage_function_t stage_functions_hc[] = {rna_last_hc_stage, post_pair_stage};
      char *stage_labels_hc[] = {"RNA HARD CLIPPINGS", "POST PAIR"};
      workflow_set_stages(2, (workflow_stage_function_t *)&stage_functions_hc, stage_labels_hc, wf_hc);
      workflow_set_producer((workflow_producer_function_t *)file_reader_2, "Buffer reader", wf_hc);

      if (options->bam_format) {
	workflow_set_consumer((workflow_consumer_function_t *)bam_writer, "BAM writer", wf_hc);
      } else {
	workflow_set_consumer((workflow_consumer_function_t *)sam_writer, "SAM writer", wf_hc);
      }
     

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
      sa_wf_input_t *wf_input = sa_wf_input_new(options->bam_format, &reader_input, wf_batch);
      
      // create and initialize workflow
      workflow_t *wf = workflow_new();      
      workflow_stage_function_t stage_functions[] = {sa_rna_mapper};
      char *stage_labels[] = {"SA mapper"};
      workflow_set_stages(1, (workflow_stage_function_t *)&stage_functions, stage_labels, wf);      
      // optional producer and consumer functions
      workflow_set_producer(sa_fq_reader_rna, "FastQ reader", wf);

      if (options->bam_format) {
	//workflow_set_consumer(sa_bam_writer_rna, "BAM writer", wf);
	workflow_set_consumer(write_to_file, "SAM writer", wf);
      } else {
	//workflow_set_consumer(sa_sam_writer_rna, "SAM writer", wf);
	workflow_set_consumer(write_to_file, "SAM writer", wf);
      }

      //Create and initialize second workflow
      workflow_t *wf_last = workflow_new();
      workflow_stage_function_t stage_functions_last[] = {sa_rna_mapper_last};
      char *stage_labels_last[] = {"SA mapper last stage"};
      workflow_set_stages(1, (workflow_stage_function_t *)&stage_functions_last, stage_labels_last, wf_last);
      workflow_set_producer(sa_alignments_reader_rna, "FastQ reader", wf_last);
      
      if (options->bam_format) {
	//workflow_set_consumer(sa_bam_writer_rna, "BAM writer", wf_last);
	workflow_set_consumer(write_to_file, "SAM writer", wf_last);
      } else {
	//workflow_set_consumer(sa_sam_writer_rna, "SAM writer", wf_last);
	workflow_set_consumer(write_to_file, "SAM writer", wf_last);
      }

      printf("Run workflow with %i threads\n", options->num_cpu_threads);

      //Extrae_init(); 

      start_timer(time_s1);
      workflow_run_with(options->num_cpu_threads, wf_input, wf);
      stop_timer(time_s1, time_e1, time_total_1);


      //extern size_t total_reads;
      //extern size_t reads_no_map;
      //printf("===== W1 %f(s) =====\n", time_total_1 / (double) 1000000);

      //printf("Total Reads           : %lu\n", total_reads);
      //printf("Total Reads Mapped    : %lu (%f)\n", total_reads - reads_no_map, (double)((total_reads - reads_no_map) * 100) / (double)(total_reads) );
      //printf("Total Reads No Mapped : %lu (%f)\n", reads_no_map, (double)((reads_no_map) * 100) / (double)(total_reads) );
      //printf("====================\n");

      printf("= = = = T I M I N G    W O R K F L O W    '1' = = = =\n");
      workflow_display_timing(wf);
      printf("= = = = - - - - - - - - - - - - - - - - - - - = = = =\n\n");

      //Extrae_event(6000019, 9); 
      //sleep(2);
      //Extrae_event(6000019, 0); 
      //extern double time_free_alig, time_free_list, time_free_batch;
      //printf("Time free Writer Alig-free=%f(s), list-free=%f(s), batch-free=%f\n", time_free_alig / (double)1000000, time_free_list / (double)1000000, time_free_batch / (double)1000000);

      //extern double time_timer0;
      //extern double time_timer1;
      //extern double time_timer2;
      //extern double time_timer3;

      //printf("Metaexon times: Timer0 = %f, Timer1 = %f, Timer2 = %f, Timer3 = %f\n", time_timer0 / (double)1000000, 
      //     time_timer1 / (double)1000000, time_timer2 / (double)1000000, time_timer3 / (double)1000000);



      start_timer(time_s2);
      rewind(f_sa);
      workflow_run_with(options->num_cpu_threads, wf_input, wf_last);      
      stop_timer(time_s2, time_e2, time_total_2);


      //Extrae_fini();
      //printf("===== W2 %f(s) =====\n", time_total_2 / (double)1000000);
      //printf("Total Reads           : %i\n", total_reads);
      //printf("Total Reads Mapped    : %i (%f)\n", total_reads - reads_no_map, (double)((total_reads - reads_no_map) * 100) / (double)(total_reads) );
      //printf("Total Reads No Mapped : %i (%f)\n", reads_no_map, (double)((reads_no_map) * 100) / (double)(total_reads) );
      //printf("====================\n");
      
      printf("= = = = T I M I N G    W O R K F L O W    '2' = = = =\n");
      workflow_display_timing(wf_last);
      printf("= = = = - - - - - - - - - - - - - - - - - - - = = = =\n\n");

      //extern double time_read_fq  ;
      //extern double time_read_fq_process ;
      //extern double time_read_alig ;
      //extern double time_read_proc ;
      
      //printf("Tiempo Lectura fq = %f, Tiempo procesamiento lectura = %f, Tiempo lectura Alig & buffer = %f, Tiempo procesamiento = %f\n", time_read_fq / (double)1000000, time_read_fq_process / (double)1000000, time_read_alig / (double)1000000, time_read_proc / (double)1000000);
      
      //printf("Time W2: %f(s)\n", time_total / 1000000);
      //fprintf(stderr, "========== FINISHED MAPPING in %f(s) =========\n", (time_total_1 + time_total_2) / 1000000);
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
      extern size_t search_calls;
      extern size_t insert_calls;
      extern double time_search;
      extern double time_insert;

     
      printf("Load Genome Time : %f(s)\n", time_genome / 1000000);
      printf("         W1 Time : %f(s)\n", time_total_1 / 1000000);
      printf("         W2 Time : %f(s)\n", time_total_2 / 1000000);
      printf("         W TOTAL : %f(s)\n", (time_total_1 + time_total_2) / 1000000);
      printf("         TOTAL   : %f(s)\n", (time_genome + time_total_1 + time_total_2) / 1000000);
     

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
  //size_t total_meta = 0;
  //for (unsigned int i = 0; i < num_chromosomes; i++) {
  //printf("Chr %i : %lu\n", i + 1, linked_list_size(metaexons->metaexons_list[i]));
  //total_meta += linked_list_size(metaexons->metaexons_list[i]);
  //}
  //printf("Total Metaexons = %lu\n", total_meta);
  //metaexons_free(metaexons);
  
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
  } else {
    // free memory
    if (sa_index)  { sa_index3_free(sa_index);  }
  }

  if (options->bam_format) {
    bam_fclose(writer_input.bam_file);
  } else {
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










