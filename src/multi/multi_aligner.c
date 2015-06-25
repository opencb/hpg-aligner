#include <multi/multi_aligner.h>

#define TMP_PATH     "TMP_FILES"
#define OUTPUT_PATH  "MPI_output"

#define SAM_FILE 1
#define BAM_FILE 2

#define MAX_HEAD_BUFFER 10485760
#define MAX_BUFFER      104857600

pthread_mutex_t mutex_sp;
int min_intron, max_intron;
char convert_ASCII[128];


int total_writes;

int mapper_mode;

    /*
    char* optional_fields = (char*) bam1_aux(bam_line);
    int i = 0;
    printf("Size %i: ", bam_line->l_aux);
    while (i < bam_line->l_aux) {
        printf("%c", optional_fields[i++]);
	printf("%c:%c: ", optional_fields[i++], optional_fields[i++]);
    }
    printf("\n");
    */
    
    /*
    char rnext[4] = "*\0";    
    int pnext = 0, tlen = 0;
    
    convert_to_quality_string_length(qual_tmp, bam1_qual(bam_line), bam_line->core.l_qseq, 33);
    BAM_convert_to_sequence(seq_tmp, bam1_seq(bam_line), bam_line->core.l_qseq);

    
    sprintf(SAM_line_out, "%s\t%i\t%s\t%i\t%i\t%s\t%s\t%i\t%i\t%s\t%s\t%s\n", 
	    bam1_qname(bam_line), 
	    bam_line->core.flag,
	    bam_header->target_name[bam_line->core.tid],
	    bam_line->core.pos + 1,
	    bam_line->core.qual,
	    cigar,
	    rnext,
	    pnext,
	    tlen,
	    seq_tmp,
	    qual_tmp,
	    op_tmp);
    */
    
    //printf("%s\n", SAM_line_out);
    //exit(-1);


void sa_index3_parallel_genome_new(char *sa_index_dirname, int num_threads,
				   sa_index3_t **sa_index_out, genome_t **genome_out) {  

  float load_progress = 0.0;

  printf("\nLoad Genome Status\n");

  extern pthread_mutex_t mutex_sp;
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

  sprintf(filename_tab, "%s/params.txt", sa_index_dirname);
  //printf("reading %s\n", filename_tab);

  f_tab = fopen(filename_tab, "r");
  if (!f_tab) {
    fprintf(stderr, "Error opening file %s!\n", filename_tab);
    exit(-1);
  }

  // prefix
  char *res = fgets(line, 1024, f_tab);
  line[strlen(line) - 1] = 0;
  prefix = strdup(line);
  // k_value
  res = fgets(line, 1024, f_tab);
  k_value = atoi(line);
  // pre_length
  res = fgets(line, 1024, f_tab);
  pre_length = atoi(line);
  // A_items
  res = fgets(line, 1024, f_tab);
  A_items = atoi(line);
  // IA_items
  res = fgets(line, 1024, f_tab);
  IA_items = atol(line);
  // num_suffixes
  res = fgets(line, 1024, f_tab);
  num_suffixes = atoi(line);
  // genome_length
  res = fgets(line, 1024, f_tab);
  genome_len = atoi(line);
  // num_chroms
  res = fgets(line, 1024, f_tab);
  num_chroms = atoi(line);
  
  size_t *chrom_lengths = (size_t *) malloc(num_chroms * sizeof(size_t));
  char **chrom_names = (char **) malloc(num_chroms * sizeof(char *));
  char chrom_name[1024];
  size_t chrom_len;
		  
  for (int i = 0; i < num_chroms; i++) {
    res = fgets(line, 1024, f_tab);
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
	  printf("Error: (%s) mismatch num_items = %i vs length = %i\n", 
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

	pthread_mutex_lock(&mutex_sp);
	load_progress += 10;

	pthread_mutex_unlock(&mutex_sp);

	// Compressed Row Storage (A table)
	A = NULL;
	sprintf(filename_tab, "%s/%s.A", sa_index_dirname, prefix);
	f_tab = fopen(filename_tab, "rb");
	if (f_tab) {
	  A = (uint *) malloc(A_items * sizeof(uint));
	
	  //	printf("\nreading A table (Compression Row Storage) from file %s...\n", filename_tab);
	  gettimeofday(&start, NULL);
	  if ((num_items = fread(A, sizeof(uint), A_items, f_tab)) != A_items) {
	    printf("Error: (%s) mismatch read num_items = %i (it must be %i)\n", 
		   filename_tab, num_items, A_items);
	    exit(-1);
	  }
	  gettimeofday(&stop, NULL);
	  //	printf("end of reading A table (%lu num_items) from file %s in %0.2f s\n", 
	  //	       num_items, filename_tab,
	  //	       (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
	  fclose(f_tab);
	}

	pthread_mutex_lock(&mutex_sp);
	load_progress += 27.1;

	pthread_mutex_unlock(&mutex_sp);

	// Compressed Row Storage (IA table)
	IA = NULL;
	sprintf(filename_tab, "%s/%s.IA", sa_index_dirname, prefix);
	f_tab = fopen(filename_tab, "rb");
	if (f_tab) {
	  IA = (uint *) malloc(IA_items * sizeof(uint));
	
	  //	printf("\nreading IA table (Compression Row Storage) from file %s...\n", filename_tab);
	  gettimeofday(&start, NULL);
	  if ((num_items = fread(IA, sizeof(uint), IA_items, f_tab)) != IA_items) {
	    printf("Error: (%s) mismatch read num_items = %i (it must be %i)\n", 
		   filename_tab, num_items, IA_items);
	    exit(-1);
	  }
	  gettimeofday(&stop, NULL);
	  //	printf("end of reading IA table (%lu num_items) from file %s in %0.2f s\n", 
	  //	       num_items, filename_tab,
	  //	       (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
	  fclose(f_tab);
	}

	pthread_mutex_lock(&mutex_sp);
	load_progress += 3.5;

	pthread_mutex_unlock(&mutex_sp);

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
	  printf("Error: (%s) mismatch num_items = %i vs num_suffixes = %i\n", 
		 filename_tab, num_items, num_suffixes);
	  exit(-1);
	}
	gettimeofday(&stop, NULL);
	//      printf("end of reading SA table (%lu num_suffixes) from file %s in %0.2f s\n", 
	//	     num_suffixes, filename_tab,
	//	     (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);      
	fclose(f_tab);

	pthread_mutex_lock(&mutex_sp);
	load_progress += 40;

	pthread_mutex_unlock(&mutex_sp);

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
	  printf("Error: (%s) mismatch num_items = %i vs num_suffixes = %i\n", 
		 filename_tab, num_items, num_suffixes);
	  exit(-1);
	}
	gettimeofday(&stop, NULL);
	//      printf("end of reading CHROM table (%lu num_suffixes) from file %s in %0.2f s\n", 
	//	     num_suffixes, filename_tab,
	//	     (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);      
	fclose(f_tab);


	pthread_mutex_lock(&mutex_sp);
	load_progress += 9.6;

	pthread_mutex_unlock(&mutex_sp);

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
	    printf("Error: (%s) mismatch num_items = %i\n", 
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
	    printf("Error: (%s) mismatch read num_items = %i (it must be %i)\n", 
		   filename_tab, num_items, A_items);
	    exit(-1);
	  }
	  gettimeofday(&stop, NULL);
	  //	printf("end of reading JA table (%lu num_items) from file %s in %0.2f s\n", 
	  //	       num_items, filename_tab,
	  //	       (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
	  fclose(f_tab);
	}


	pthread_mutex_lock(&mutex_sp);
	load_progress += 6.8;

	pthread_mutex_unlock(&mutex_sp);

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
	
	pthread_mutex_lock(&mutex_sp);
	load_progress += 3.5;

	pthread_mutex_unlock(&mutex_sp);

      }
    }
  }


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
}

//--------------------------------------------------------------------------------------

/*
//Return TMP_PATH name
int split_input_file(char *fq_str, char *tmp_path, int numprocs) {

  FILE *tmp_fd;
  FILE *fq_file = fopen(fq_str, "r");

  //Prepare PATH tmp files
  create_directory(tmp_path);
  
  int count;
  char header1[MAX_READ_ID_LENGTH];
  char sequence[MAX_READ_SEQUENCE_LENGTH];
  char header2[MAX_READ_ID_LENGTH];
  char qualities[MAX_READ_SEQUENCE_LENGTH];

  char *res;

  char tmp_file[1024];
  int file_id = 0;

  int reads_split = ;
  
  res = fgets(header1, MAX_READ_ID_LENGTH, fq_file);
  
  while (res) {    
    sprintf(tmp_file, "%s/%i.tmp", tmp_path, file_id);
    tmp_fd = fopen(tmp_file, "w");
    count = 0;
    
    while (count < reads_split) {      
      res = fgets(sequence, MAX_READ_SEQUENCE_LENGTH, fq_file);
      res = fgets(header2, MAX_READ_ID_LENGTH, fq_file);
      res = fgets(qualities, MAX_READ_SEQUENCE_LENGTH, fq_file);

      fputs(header1, tmp_fd);
      fputs(sequence, tmp_fd);
      fputs(header2, tmp_fd);
      fputs(qualities, tmp_fd);
      
      count++;

      res = fgets(header1, MAX_READ_ID_LENGTH, fq_file);
      if (!res) break;
    }
    
    fclose(tmp_fd);   
    //if (!res) break;    
    file_id++;    
    
  }
  
  
  
  fclose(fq_file);

  //file_id++;
  //printf("Split file in %i\n", file_id);
  
  return file_id;
  
}
*/


int split_input_file(char *fq_str, char *tmp_path, int numprocs) {
  FILE *tmp_fd;
  FILE *fq_file = fopen(fq_str, "r");
  if (!fq_file) { printf("ERROR opening file\n"); exit(-1); }

  //Prepare PATH tmp files
  char cmd[2048];
  sprintf(cmd, "rm -rf %s", tmp_path);
  system(cmd);  
  create_directory(tmp_path);

  //Calculate the number of reads per file
  int reads_split;
  size_t reads_positions[numprocs + 1];

  fseek(fq_file, 0L, SEEK_END);
  size_t fd_total_bytes = ftell(fq_file);
  fseek(fq_file, 0L, SEEK_SET);
    
  int i;
  size_t seek_bytes = 0;
  size_t inc_bytes = fd_total_bytes / numprocs;
  size_t offset;
  char line[2048];
    
  for (i = 0; i < numprocs; i++) {
    offset = 0;
    fseek(fq_file, seek_bytes, SEEK_SET);
    fgets(line, 2048, fq_file);
    while (line[0] != '@') {
      offset += strlen(line);
      fgets(line, 2048, fq_file);
    }      
    seek_bytes += offset;
    reads_positions[i] = seek_bytes;
    seek_bytes += inc_bytes;
  }
  
  reads_positions[i] = fd_total_bytes;
  fclose(fq_file);

  fq_file = fopen(fq_str, "r");
  char sequence[MAX_READ_ID_LENGTH];
  size_t total_size = 0;
  char tmp_file[1024];
  char *res;

  //printf("SPLIT FILES\n");
  for (int rank = 0; rank < numprocs; rank++) {
    sprintf(tmp_file, "%s/%i.tmp", tmp_path, rank);
    tmp_fd = fopen(tmp_file, "w");
    //printf("[%i] %s => %lu < %lu\n", rank, tmp_file, total_size, reads_positions[rank + 1]);
    
    while (total_size < reads_positions[rank + 1]) {
      fgets(sequence, MAX_READ_ID_LENGTH, fq_file);
      fputs(sequence, tmp_fd);
      total_size += strlen(sequence);
    }
    
    fclose(tmp_fd);
    
  }

  fclose(fq_file);

}


char* BAM_convert_to_sequence(char *sequence_string, uint8_t* sequence_p, int sequence_length) {

  //char* sequence_string = (char*) calloc(1, sequence_length + 1); //each byte codes two nts ( 1 nt = 4 bits)

  for (int i = 0; i < sequence_length; i++) {
    switch (bam1_seqi(sequence_p, i)) {
    case 1:
      sequence_string[i] = 'A';
      break;
    case 2:
      sequence_string[i] = 'C';
      break;
    case 4:
      sequence_string[i] = 'G';
      break;
    case 8:
      sequence_string[i] = 'T';
      break;
    case 15:
      sequence_string[i] = 'N';
      break;
    }
    
  }
  
  sequence_string[sequence_length] = '\0';
  
  return sequence_string;

}


int align_to_file(int line_format, char *align_str, char *align_tmp, FILE *fd_buffer,
		  int second_phase, avls_list_t *avls_list, metaexons_t *metaexons,
		  genome_t *genome, bam1_t* bam_line, bam_header_t* bam_header,
		  char *SAM_line_out) {

  int to_file = 1;
  int update_meta = 0;
  int cigar_null = 0;
  int no_map = 0;
  
  char *pch;
  int itr = 0;
  char *cigar = NULL;
  size_t start_map;
  int chromosome;
  int strand;
  int flag;
  char *id, *sequence, *quality;
  cigar_code_t *cigar_code = NULL;
  int hard_clipping = 0;
  
  //printf("%s", align_str);

  if (mapper_mode == RNA_MODE) {
    if (line_format == SAM_FILE) {
      pch = strtok (align_str, "\t");
      while (pch != NULL) {
	itr++;
	if (itr == 2) {
	  flag = atoi(pch);
	  strand = (((uint32_t) flag) & BAM_FREVERSE) ? 1 : 0;
	} else if (itr == 3) {
	  //TODO: When chromosome is not a number or X,Y and MT, we have a problem... Solve this ;-)
	  if (pch[0] == 'c' && pch[1] == 'h' && pch[2] == 'r') {
	    pch += 3;
	  }      
	  if (strcmp(pch,"X") == 0) {
	    chromosome = 23;
	  } else if (strcmp(pch,"Y") == 0) {
	    chromosome = 24;
	  } else if (strcmp(pch,"MT") == 0) {
	    chromosome = 25;
	  } else {
	    chromosome = atoi(pch);
	  }
	} else if (itr == 4) {
	  start_map = atol(pch);
	  start_map--;
	} else if (itr == 6) {
	  cigar = strdup(pch);
	  break;
	}
	pch = strtok(NULL, "\t");	  
      }
    } else {
      flag = bam_line->core.flag;
      strand = (((uint32_t) flag) & BAM_FREVERSE) ? 1 : 0;

      if (flag != 4) {
	char *chr_aux = bam_header->target_name[bam_line->core.tid];
        
	if (chr_aux[0] == 'c' && chr_aux[1] == 'h' && chr_aux[2] == 'r') {
	  chr_aux += 3;
	}
    
	if (strcmp(chr_aux,"X") == 0) {
	  chromosome = 23;
	} else if (strcmp(chr_aux,"Y") == 0) {
	  chromosome = 24;
	} else if (strcmp(chr_aux,"MT") == 0) {
	  chromosome = 25;
	} else {
	  chromosome = atoi(chr_aux);
	}

	start_map = bam_line->core.pos;
	uint32_t *cigar_bam = bam1_cigar(bam_line);
	cigar = convert_to_cigar_string(cigar_bam, bam_line->core.n_cigar);
      } else {
	cigar = strdup("*");
      }    
      //printf("BAM: %i %i %i\n", flag, strand, bam_line->core.tid);
    }
  } else {

    if (line_format == SAM_FILE) {
      pch = strtok (align_str, "\t");
      while (pch != NULL) {
	itr++;
	if (itr == 2) {
	  flag = atoi(pch);
	  strand = (((uint32_t) flag) & BAM_FREVERSE) ? 1 : 0;
	  break;
	}
	pch = strtok(NULL, "\t");	  
      }      
    } else {
      flag = bam_line->core.flag;
    }

    if (flag == 4 && second_phase) {
      to_file = 0;
    }
    
  }


  if (mapper_mode == RNA_MODE) {
    if (cigar) {
      //printf("\tcigar: %s\n", cigar);
      int len = strlen(cigar);
      int br = 0;
    
      for (int i = 0; i < len; i++) {
	if (cigar[i] == '*') {
	  no_map = 1;
	  if (second_phase) { to_file = 0; }
	  br = 1;
	  break;
	} else if (cigar[i] == 'S') {
	  if (second_phase) { to_file = 0; }
	  br = 1;
	  break;
	} else if (cigar[i] == 'H') {
	  br = 1;
	  break;
	}
      }
    
      //if (!br) { update_meta = 1; }    
    
    } else {
      cigar_null = 1;
      to_file = 0;
      no_map = 1;
    }
  
  
    if (!no_map) {
      cigar_code = cigar_code_new_by_string(cigar);
    
      size_t sj_start, sj_end;
      size_t offset = start_map;
      size_t exon_start = offset;
      size_t exon_end;
    
      const int MAX_SJ = 100;
      int exons_array[MAX_SJ*4];
      int num_sj = 0;
      int pos = 0;
    
      avl_node_t *node_avl_start, *node_avl_end;

      for (int i = 0; i < array_list_size(cigar_code->ops); i++) {
	cigar_op_t *op = array_list_get(i, cigar_code->ops);
      
	if (op->name == 'D' || op->name == 'M') {
	  offset += op->number;
	}
    
	if (op->name == 'N') {
	  exon_end = offset;
	
	  exons_array[pos++] = exon_start;
	  exons_array[pos++] = exon_end;
	  num_sj++;
	
	  offset += op->number;
	  exon_start = offset;
	}    
      }

      if (num_sj > 0) {
	exon_end = offset;
	exons_array[pos++] = exon_start;
	exons_array[pos++] = exon_end;
      
	pos = 0;
      
	int sp_type = 1;
	//char sj_start_ref[10];
	//char sj_end_ref[10];
	char sj_ref[10];
	size_t genome_start, genome_end;
	int sj_strand;
	int report;

	for (int i = 0; i < num_sj; i++) {
	  sj_start = exons_array[pos + 1] + 1;
	  sj_end   = exons_array[pos + 2];
	
	  //printf("%lu : %i\n", start_map, chromosome);	
	  genome_start = sj_start;
	  genome_end   = genome_start + 2;      
	  genome_read_sequence_by_chr_index(sj_ref, 0, chromosome - 1,
					    &genome_start, &genome_end, genome);
	
	  //printf("SP_START: %c%c\n", sj_ref[0], sj_ref[1]);

	  sj_ref[2] = '-';

	  genome_start = sj_end - 1;
	  genome_end   = genome_start + 1;
	  genome_read_sequence_by_chr_index(&sj_ref[3], 0, chromosome - 1,
					    &genome_start, &genome_end, genome);
		
	  sj_ref[5] = '\0';

	  //printf("SP_END: %c%c\n", sj_ref[3], sj_ref[4]);
	
	  //sj_strand = strand;
	  //report = 0;
	
	  //report = 0;
	  if ((sj_ref[0] == 'C' && sj_ref[1] == 'T' && sj_ref[3] == 'A' && sj_ref[4] == 'C') ||
	      (sj_ref[0] == 'G' && sj_ref[1] == 'T' && sj_ref[3] == 'A' && sj_ref[4] == 'T') || 
	      (sj_ref[0] == 'C' && sj_ref[1] == 'T' && sj_ref[3] == 'G' && sj_ref[4] == 'C')) {
	    //report = 1;
	    sj_strand = 1;
	  } else if ((sj_ref[0] == 'G' && sj_ref[1] == 'T' && sj_ref[3] == 'A' && sj_ref[4] == 'G') ||
		     (sj_ref[0] == 'A' && sj_ref[1] == 'T' && sj_ref[3] == 'A' && sj_ref[4] == 'C') || 
		     (sj_ref[0] == 'G' && sj_ref[1] == 'C' && sj_ref[3] == 'A' && sj_ref[4] == 'G')) {
	    sj_strand = 0;
	    //report = 1;
	  } else {
	    sj_strand = strand;
	  }

	  //if (report) {
	  //printf("Insert SJ\n");
	
	  allocate_start_node(chromosome - 1,
			      sj_strand,
			      sj_start,
			      sj_end,
			      sj_start,
			      sj_end,
			      FROM_READ,
			      sp_type,
			      sj_ref,
			      &node_avl_start,
			      &node_avl_end,
			      avls_list);
	
	  exon_start = exons_array[pos] + 1;
	  exon_end = sj_start - 1;
	  //printf("Meta-L: %lu - %lu\n", exon_start, exon_end);
	
	  metaexon_insert(strand, chromosome - 1,
			  exon_start, exon_end, 40,
			  METAEXON_RIGHT_END, node_avl_start,
			  metaexons);
	  
	
	  exon_start = sj_end + 1;
	  exon_end   = exons_array[pos + 3];
	  //printf("Meta-R: %lu - %lu\n", exon_start, exon_end);
	
	  
	  metaexon_insert(strand, chromosome - 1,
			  exon_start, exon_end, 40,
			  METAEXON_LEFT_END, node_avl_end,
			  metaexons);
	  
	  //size_t num_sj = get_total_items(25, avls_list);
	  //printf("Num SJ = %i\n", num_sj);
      
	  pos += 2;
	}
      }
    }
  }

  if (!to_file) {    
    
    total_writes++;
    //======= Select buffer to write =========//
    char read_tmp[4098*3];
    int len_tmp = 1;
    
    read_tmp[0] = '@';
    
    if (line_format == SAM_FILE) {
      itr = 0;
      pch = strtok(align_tmp, "\t");
      while (pch != NULL) {
	itr++;
	if (itr == 1) {
	  strcpy(&read_tmp[1], pch);
	  len_tmp += strlen(pch);
	  read_tmp[len_tmp] = '\n';
	  len_tmp++;
	} else if (itr == 10) {
	  strcpy(&read_tmp[len_tmp], pch);
	  len_tmp += strlen(pch);
	  strcpy(&read_tmp[len_tmp], "\n+\n\0");
	  len_tmp += 3;
	} else if (itr == 11) {
	  strcpy(&read_tmp[len_tmp], pch);
	  len_tmp += strlen(pch);
	  //read_tmp[strlen(read_tmp)] = '\n';
	  break;
	} 
	pch = strtok(NULL, "\t");	
      }

      if (read_tmp[len_tmp - 1] != '\n') {
	strcpy(&read_tmp[len_tmp], "\n");
	len_tmp++;
      }
    
      read_tmp[len_tmp] = '\0';
      
    } else {
      char *id_aux = bam1_qname(bam_line);
      //printf("id: %s, (%lu)\n", id_aux, strlen(id_aux));
      
      strcpy(&read_tmp[1], id_aux);
      len_tmp += strlen(id_aux);
      read_tmp[len_tmp] = '\n';
      len_tmp++;
      
      BAM_convert_to_sequence(&read_tmp[len_tmp], bam1_seq(bam_line), bam_line->core.l_qseq);
      
      len_tmp += bam_line->core.l_qseq;
      strcpy(&read_tmp[len_tmp], "\n+\n\0");
      len_tmp += 3;
      
      convert_to_quality_string_length(&read_tmp[len_tmp], bam1_qual(bam_line), bam_line->core.l_qseq, 33);
      len_tmp += bam_line->core.l_qseq;

      if (read_tmp[len_tmp - 1] != '\n') {
	strcpy(&read_tmp[len_tmp], "\n");
	len_tmp++;
      }
      
      read_tmp[len_tmp] = '\0';

    }

    //printf("%s\n", read_tmp);
    //exit(-1);
    
    fwrite(read_tmp, strlen(read_tmp), sizeof(char), fd_buffer);
    
  } else if (line_format == BAM_FILE) {
    char qual_tmp[2048];
    char seq_tmp[2048];
    uint8_t op_tmp[bam_line->l_aux];
    char fop_tmp[bam_line->l_aux];

    char *txt = bam_format1_core(bam_header, bam_line, BAM_OFDEC);
    strcpy(SAM_line_out, txt);
    SAM_line_out[strlen(txt)] = '\n';
    SAM_line_out[strlen(txt) + 1] = '\0';
    
  }
  
  return to_file;
  
}





/*
void *multialigner_second_phase_reader(void *input) {
  wf_input_t *wf_input = (wf_input_t *) input;
  batch_t *new_batch = NULL;
  batch_t *batch = wf_input->batch;
  fastq_batch_reader_input_t *fq_reader_input = wf_input->fq_reader_input;

  const int MAX_READS = 200;
  array_list_t *reads = array_list_new(MAX_READS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  int n_reads = 0;
  
  FILE *fd = fq_reader_input.file1;

  while (fgets(line, 4096, fd)) {
    if (line[0] == '@') { continue; }

    n_reads++;
    if (n_reads >= MAX_READS) { break; }
  }
  
  strcpy(line_strtok, line);
  strcpy(line_tmp, line);
  if (is_correct_or_to_buffer(line_strtok, line_tmp, second_phase, avls_list, metaexons, genome)) {
    //Update metaexon and AVL
	strcpy(&buffer[len_buffer], line);
	len_buffer += strlen(line);	  
  }
  
  while (fgets(line, 4096, fd)) {
    strcpy(line_strtok, line);
    strcpy(line_tmp, line);
    if (is_correct_or_to_buffer(line_strtok, line_tmp, second_phase, avls_list, metaexons, genome)) {	    
      if (len_buffer + strlen(line) >= MAX_BUFFER) {
	MPI_Send(buffer, MAX_BUFFER, MPI_CHAR, w_rank, 1, MPI_COMM_WORLD);
	len_buffer = 0;
      }
      strcpy(&buffer[len_buffer], line);
      len_buffer += strlen(line);
    }
  }
  

     size_t num_reads = array_list_size(reads);

     if (num_reads == 0) {
	  array_list_free(reads, (void *)fastq_read_free);
     } else {
	  mapping_batch_t *mapping_batch = mapping_batch_new(reads, 
							     batch->pair_input->pair_mng);

	  new_batch = batch_new(batch->bwt_input, batch->region_input, batch->cal_input, 
				batch->pair_input, batch->preprocess_rna, batch->sw_input, batch->writer_input, 
				batch->mapping_mode, mapping_batch, batch->data_out);
     }

     //if (time_on) { stop_timer(start, end, time); timing_add(time, FASTQ_READER, timing); }
     //printf("Read batch %i\n", num_reads);
     
     return new_batch;

  for (int i = 0; i < array_list_size(files_found); i++) {
    char *file = array_list_get(i, files_found);
    int type   = array_list_get(i, files_type);
      
    char file_path[strlen(file) + strlen(node_out) + 1024];
    sprintf(file_path, "%s/%s", node_out, file);
      
    FILE *fd = fopen(file_path, "r");
    char line[4096], line_strtok[4096], line_tmp[4096];
    char *point;
    size_t nlines = 0;
      
    if (type == 1) {
	
    } else if (type == 2) { 
      //TODO: Read BAM
      continue;
    }
      
  } //end for
  
}
*/

/*
void genome_fasta(char *path) {

  DIR *dir;
  struct dirent *ent;
  int n_files = 0;
  char *file_name;
  
  if ((dir = opendir (path)) != NULL) {
    while ((ent = readdir (dir)) != NULL) {
      //printf ("\t%s\n", ent->d_name);
      if (strstr(ent->d_name, "fa") != NULL) {
	if (n_files == 1) {
	  printf("Error more than one fasta files found\n");
	  exit(-1);
	}
	file_name = strdup(ent->d_name);
	n_files++;
      }
    }
    closedir (dir);
  } else {
    printf("Directory not found\n");
    exit(-1);
  }

  char path_fa[strlen(path) + strlen(file_name) + 1024];
  sprintf(path_fa, "%s/%s", path, file_name);


  FILE *fd = fopen(path_fa, "r");
  const int MAXLINE = 4096;
  size_t max_len_genome = 3200000000; //3GB aprox.
  size_t max_chr = 30;
  
  char *genome = (char *)calloc(max_len_genome, sizeof(char));
  size_t *chr_offset = (size_t *)malloc(sizeof(size_t)*max_chr);
  char   **chr_name   = (char **)malloc(sizeof(char *)*max_chr);

  
  char line[MAXLINE];
  char chr_name[256];
  int i, j;
  
  while (fgets(line, MAXLINE, fd)) {
    if (line[0] == '>') {
      //Line start with >
      i = 1;
      j = 0;
      while (line[i] != ' ') {
	chr_name[j++] = line[i];
      }
      chr_name[j] = '\0';
    } else {
      strcpy(genome, line);
    }
    
  }
    fclose(fd);

  
}
*/


int hpg_multialigner_main(int argc, char *argv[]) {

  int numprocs;
  int rank;  
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  int w_rank;
  MPI_Status status;
  int num_chromosomes_header = 0;
  int second_phase;
  
  struct timeval time_start, time_stop;
  struct timeval total_start, total_stop;
  
  double time_split   = 0.0;
  double time_mapping = 0.0;
  double time_merge   = 0.0;
  double time_total   = 0.0;  

  MPI_Init(&argc, &argv);
  
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Get_processor_name(processor_name, &name_len);  
  
  //printf("[%i/%i] %s\n", rank, numprocs, processor_name);

  ////////////////////////////////// MULTI ALIGNER VALIDATE //////////////////////////////////////////
  
  if (argc <= 1) {
    if (rank == 0) {
      printf("Missing command.\nValid commands are:\n\tdna: to map DNA sequences\n\trna: to map RNA sequences\nUse -h or --help to display hpg-aligner options.\n");
    }
    MPI_Finalize();
    exit(-1);
  }
  
  if (argc <= 1 || strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
    if (rank == 0) {
      usage_cli(mapper_mode);
    }
    MPI_Finalize();
    exit(-1);
  }
  
  char *command = argv[1];
  
  argc -= 1;
  argv += 1;
  
  if(strcmp(command, "dna") != 0 && 
     strcmp(command, "rna") != 0) {
    if (rank == 0) {
      printf("Missing command.\nValid commands are:\n\tdna: to map DNA sequences\n\trna: to map RNA sequences\nUse -h or --help to display hpg-aligner options.\n");
    }
    MPI_Finalize();
    exit(-1);
  }

  if (strcmp(command, "dna") == 0) {
    mapper_mode = DNA_MODE;
  } else if (strcmp(command, "rna") == 0) {
    mapper_mode = RNA_MODE;
  }

  int mode, num_options = NUM_OPTIONS;

  mode = RNA_MODE;

  void **argtable = argtable_options_new(mode);
  options_t *options = options_new();
  options->mode = mode;

  if (argc < 2) {
    if (rank == 0) {
      usage(argtable);
    }
    MPI_Finalize();
    exit(-1);
  } else {
    int num_errors = arg_parse(argc, argv, argtable);   
    //if (((struct arg_int*)argtable[24])->count) {
    //if (rank == 0) {
    //	usage(argtable);
    //}
    //argtable_options_free(argtable, num_options);
    //options_free(options);
    //MPI_Finalize();
    //exit(0);
    //}

    
    if (num_errors > 0) {
      if (rank == 0) {
	fprintf(stdout, "Errors:\n");
	// struct end is always allocated in the last position
	arg_print_errors(stdout, argtable[num_options], "hpg-multialigner");      
	usage(argtable);
      }
      MPI_Finalize();
      exit(-1);
    }else {
      options = read_CLI_options(argtable, options);
      if(options->help) {
	if (rank == 0) {
	  usage(argtable);
	  argtable_options_free(argtable, num_options);
	  options_free(options);
	}
	MPI_Finalize();
	exit(0);
      }
    }    
  }

  argtable_options_free(argtable, num_options);
  options->mode = mode;

  validate_options(options);
  
  second_phase = options->second_phase;
  int hpg_enable = 0;
  
  char *second_cli = options->second_command;
  
  if (second_phase) {
    if (!second_cli) {
      hpg_enable = 1;
    } else {
      hpg_enable = 0;
    }
  }


  if (!options->in_filename) {
    if (rank == 0) {
      printf("Filename input is missing. Please, insert it with option '-f FILENAME'.\n");
      usage_cli(mapper_mode);
    }
    MPI_Finalize();
    exit(-1);
  }

  if (!options->command) {
    if (rank == 0) {
      printf("Command line of mapper is missing. Please, insert it with option '-c \"CLI MAPPER\"'.\n");
      usage_cli(mapper_mode);
    }
    MPI_Finalize();
    exit(-1);
  }

    
  if (mapper_mode == RNA_MODE) {
    if (!options->bwt_dirname) {
      if (rank == 0) {
	printf("Index directory is missing. Please, insert it with option '-i DIRNAME'.");
	if (!second_phase || (second_phase && second_cli)) {
	  printf("Second phase is disable or enable with a second command line, therefore you only need dna_compression.bin file.\n");
	}
	usage_cli(mapper_mode);
      }
      MPI_Finalize();
      exit(-1);
    }
  } else {
    if (second_phase && !second_cli) {
      if (rank == 0) {
	printf("DNA mode is enable but second CLI is missing. Please, insert it with option '--second-command \"CLI MAPPER 2\"' or with option '-i DIRNAME'.\n");
      }
      MPI_Finalize();
      exit(-1);
    }
  } 
  /*
  if (!options->tmp_path && !options->tmp_file) {
    printf("Not temporal file or not temporal path found! Please, insert it with '--tmp-file' or '--tmp-path' option\n");
    MPI_Finalize();
    exit(-1);
  }

  if (options->tmp_path && options->tmp_file) {
    printf("Only one of '--tmp-file' and '--tmp-path' options may be set\n");
    MPI_Finalize();
    exit(-1);
  }

  if (second_phase && (!options->second_tmp_path && !options->second_tmp_file)) {
    printf("Not temporal file or not temporal path found! Please, insert it with '--tmp-file' or '--tmp-path' option\n");
    MPI_Finalize();
    exit(-1);
  }
  */
  if (options->tmp_file && (strstr(options->tmp_file, "sam") == NULL) && (strstr(options->tmp_file, "bam") == NULL)) {
    printf("Extension of temporal file name must be '.sam' or '.bam'");
    MPI_Finalize();
    exit(-1);    
  }
  /*
  if (second_phase && (options->second_tmp_path && options->second_tmp_file)) {
    printf("Only one of '--tmp-file' and '--tmp-path' options may be set\n");
    MPI_Finalize();
    exit(-1);
  }
  */
  
  ///////////////////////////////////////////////////////////////////////////////////////////

  extern int min_intron, max_intron;
  min_intron = options->min_intron_length;
  max_intron = options->max_intron_length;
  
  numprocs--;
  w_rank = numprocs;
  
  MPI_Barrier(MPI_COMM_WORLD);

  array_list_t *files_found;
  array_list_t *files_type;
  array_list_t *files_splice;
  
  char *head_buffer = (char *)malloc(sizeof(char)*MAX_HEAD_BUFFER);
  size_t len_head_buffer = 0;
  char *mapper_cli = strdup(options->command);
  char *file_input = strdup(options->in_filename);
  char *path_tmp;
  char *second_path_tmp = NULL;
  
  char pwd[2048];
  if (getcwd(pwd, sizeof(pwd)) == NULL) {
    perror("getcwd() error");
  }
  
  char *path_output;
  char *file_output = (char *)malloc(sizeof(char)*(strlen(OUTPUT_PATH) + strlen(pwd) + 1024));
  char *exact_filename = (char *)malloc(sizeof(char)*(strlen(OUTPUT_PATH) + strlen(pwd) + 1024));
  
  if (!options->output_name) {
    path_output = (char *)malloc(sizeof(char)*(strlen(OUTPUT_PATH) + strlen(pwd) + 1024));
    sprintf(path_output, "%s/%s/", pwd, OUTPUT_PATH);
  } else {
    path_output = strdup(options->output_name);
  }
  
  if (options->tmp_path) {
    path_tmp = strdup(options->tmp_path);
  } else {
    path_tmp = strdup(path_output);
  }
  
  if (second_phase) {
    if (options->second_tmp_path) {
      second_path_tmp = strdup(options->second_tmp_path);
    } else {
      second_path_tmp = strdup(path_output);
    }
  }
  
  sprintf(file_output, "%s/%s", path_output, "alignments.sam");
  sprintf(exact_filename, "%s/%s", path_output, "SJ.bed");
  

  FILE *fd_out;
  
  if (rank == 0) {
    printf("=====================================================\n"); 
    printf("START Multi Aligner in %s mode\n", mapper_mode == RNA_MODE ? "RNA" : "DNA"); 
    printf("=====================================================\n"); 
    printf("MAPPER COMMAND LINE :  '%s'\n", mapper_cli);
    printf("INPUT  FILE         :  '%s'\n", file_input);

    if (options->tmp_file) {
      printf("TMP OUTPUT PATH     :  '%s/%s'\n", path_tmp, options->tmp_file);
    } else {
      printf("TMP OUTPUT FILE     :  '%s'\n", path_tmp);
    }

    if (second_phase) {
      if (options->second_tmp_file) {
	printf("SECOND PHASE TMP OUTPUT PATH     :  '%s/%s'\n", second_path_tmp, options->second_tmp_file);
      } else {
	printf("SECOND PHASE TMP OUTPUT FILE     :  '%s'\n", path_tmp);
      }
    }    

    printf("=====================================================\n");
    printf("OUTPUT FILE         :  '%s' of node '%s'\n", file_output, processor_name);
    printf("=====================================================\n"); 
  }

  //exit(-1);
  
  if (rank == 0) {
    start_timer(total_start);
  }
  
  char node_out[strlen(path_output) + 1024];
  
  size_t max_len_cli;
  char *cli_in;
  char *cli_out;
  char *cli_end;
  unsigned char first_in = 0;

  char *tmp_input_path;

  if (options->tmp_input == NULL) {
    tmp_input_path = strdup(TMP_PATH);
  } else {
    tmp_input_path = strdup(options->tmp_input);
  }
  
  //if (rank != numprocs) {    
  //================ First split input File ====================//
  int nfiles;  
  int ntasks, extra_tasks;
  
  if (rank == 0) {
    start_timer(time_start);
    printf("[%i] Spliting input file...\n", rank);
    nfiles = split_input_file(file_input, tmp_input_path, numprocs);
    stop_timer(time_start, time_stop, time_split);
    printf("[%i] Spliting END! (%0.2f s)\n", rank, time_split);
  }
  
  MPI_Bcast(&nfiles, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  ntasks = nfiles / numprocs;
  extra_tasks = nfiles % numprocs;
  
  if (rank < extra_tasks) {
    ntasks++;
  }
  //============================================================//


  //================= SPLIT INPUT CLI AND TRANSFORM IT ==========================//
  max_len_cli = strlen(mapper_cli) + strlen(path_output) + strlen(tmp_input_path) + 1024;
  char *cli_tmp = (char *)malloc(sizeof(char)*max_len_cli);
  
  cli_in  = (char *)malloc(sizeof(char)*max_len_cli);
  cli_out = (char *)malloc(sizeof(char)*max_len_cli);
  cli_end = (char *)malloc(sizeof(char)*max_len_cli);
  
  size_t len_cli = 0;
  unsigned char find_mark = 0, first_out = 0;
    
  //if (rank == 0) {
  char *str = mapper_cli;
  char *pch;
  
  pch = strtok (str," ");
  while (pch != NULL) {
    if (strcmp(pch, "%I") == 0 || strcmp(pch, "%i") == 0) {	
      strcpy(cli_in, cli_tmp);
      len_cli = 0;
      cli_tmp[0] = '\0';
      find_mark++;
      first_in = !first_out;
    } else if (strcmp(pch, "%O") == 0 || strcmp(pch, "%o") == 0) {
      strcpy(cli_out, cli_tmp);
      len_cli = 0;
      cli_tmp[0] = '\0';
      find_mark++;
      first_out = !first_in;
    } else {
      strcpy(cli_tmp + len_cli, pch);
      len_cli = strlen(cli_tmp);
      strcpy(cli_tmp + len_cli, " \0");
      len_cli++;
    }
    
    pch = strtok(NULL, " ");
    
  }
  
  strcpy(cli_end, cli_tmp);
  
  if (find_mark != 2) {
    if (rank == 0) {
      printf("ERROR in CLI Mapper\n");
      exit(-1);
    }  
  }
  //}
  
  //=============================================================================//
  
  
  //================================ MPI PROCESS =======================================//
  //COMMON OUTPUT PATH 

  
  if (rank == 0) {
    char cmd[strlen(path_output) + 1024];
    sprintf(cmd, "rm -rf %s", path_output);
    system(cmd);
    create_directory(path_output);  

  }
    
  if (rank != numprocs) {
    char path_out_tmp[strlen(path_tmp) + 1024];
    if (options->tmp_path) {
      sprintf(path_out_tmp, "%s/%i.out/", path_tmp, rank);

      char cmd[strlen(path_out_tmp) + 1024];
      sprintf(cmd, "rm -rf %s", path_out_tmp);
      system(cmd);
      
      create_directory(path_out_tmp);      

    }     
  }  
  
  
  //PROCESS FILES
  MPI_Barrier(MPI_COMM_WORLD);
  
  
  if (rank == numprocs) {
    fd_out = fopen(file_output, "w");
  }

  MPI_Barrier(MPI_COMM_WORLD);

  char buffer_node_out[strlen(path_tmp) + 1024];
  
  if (rank != numprocs) {
    char *mapper_run = (char *)malloc(sizeof(char)*max_len_cli);

    if (options->tmp_file) {
      sprintf(node_out, "%s/%i.out/%s", path_tmp, rank, options->tmp_file);
    } else {
      sprintf(node_out, "%s/%i.out/", path_tmp, rank);
    }

    sprintf(buffer_node_out, "%s/%i.out/", path_tmp, rank);
    
    if (first_in) {
      sprintf(mapper_run, "%s %s/%i.tmp %s %s %s", cli_in, tmp_input_path, rank, cli_out, node_out, cli_end);
    } else {
      sprintf(mapper_run, "%s %s %s %s/%i.tmp %s", cli_out, node_out, cli_in, tmp_input_path, rank, cli_end);
    }

    
    start_timer(time_start);
    system(mapper_run);
    stop_timer(time_start, time_stop, time_mapping);

    printf("@@@@ [%i] (%0.2f): %s\n", rank, time_mapping / 1000000, mapper_run);
    //====================================================================================//
  }

  int num_chromosomes;
  metaexons_t *metaexons = NULL;
  avls_list_t* avls_list = NULL;
  int fast_mode;
  bwt_index_t *bwt_index = NULL;
  sa_index3_t *sa_index = NULL;
  genome_t *genome = NULL;
  sa_wf_batch_t *sa_wf_batch;      
  workflow_t *wf;          
  workflow_SA_t *sa_wf;
  wf_input_t *wf_input;
  sa_wf_input_t *sa_wf_input;      
  FILE *f_sa, *f_hc;
  batch_t *batch;
  int tag;
  //batch_buffer_t *data_out;      
  sa_rna_input_t sa_rna;
  
  if ((rank != numprocs) && (mapper_mode == RNA_MODE)) {
    printf("LOAD GENOME...\n");
    //First, search the head files
    //*************************** S E C O N D    P H A S E *******************************//
    //convert ASCII fill
    convert_ASCII['a'] = 'T';
    convert_ASCII['A'] = 'T';
      
    convert_ASCII['c'] = 'G';
    convert_ASCII['C'] = 'G';
      
    convert_ASCII['g'] = 'C';
    convert_ASCII['G'] = 'C';
      
    convert_ASCII['t'] = 'a';
    convert_ASCII['T'] = 'A';
      
    convert_ASCII['n'] = 'N';
    convert_ASCII['N'] = 'N';
    
    if (second_phase && !second_cli) {
      //printf("Second\n");
      
      char filename_tab[strlen(options->bwt_dirname) + 1024];  
      sprintf(filename_tab, "%s/params.info", options->bwt_dirname);
      FILE *fd = fopen(filename_tab, "r");      

      if (fd) { 
	fast_mode = 1; 
	fclose(fd);
      } else { 
	fast_mode = 0; 
      }

      if (hpg_enable) {
	// genome parameters 
	printf("HPG ENABLE LOADING GENOME...\n");
	if (!fast_mode) {
	  //////////////// LOAD BWT INDEX //////////////////////    
	  bwt_index = bwt_index_new(options->bwt_dirname, false);
	  genome = genome_new("dna_compression.bin", options->bwt_dirname, BWT_MODE);  
	  //////////////////////////////////////////////////////
	} else {    
	  ///////////////// LOAD SA INDEX ////////////////////// 
	  sa_index3_parallel_genome_new(options->bwt_dirname, options->num_cpu_threads, &sa_index, &genome);	
	}
      }
    } else {
      printf("Load dna...\n");
      genome = genome_new("dna_compression.bin", options->bwt_dirname, BWT_MODE);  
    }
    
    num_chromosomes = genome->num_chromosomes;
      
    //************************* M E T A E X O N    B U I L D *******************************//
    avls_list = avls_list_new(num_chromosomes);
    //if (options->transcriptome_filename != NULL) {
    //load_transcriptome(options->transcriptome_filename, genome, avls_list, metaexons);
    //}

    //You can obtain the chr_size array reading fastq file header
    metaexons = metaexons_new(num_chromosomes,
			      genome->chr_size);
    //**************************************************************************************//
  }
    
  //====================================================================================//
  // READ OUTPUT SAM
  //printf("%s:\n", node_out);
  if (rank != numprocs) {
    DIR *dir;
    struct dirent *ent;
    files_found  = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    files_type   = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  
    if ((dir = opendir (buffer_node_out)) != NULL) {
      while ((ent = readdir (dir)) != NULL) {
	//printf ("\t%s\n", ent->d_name);
	if (strstr(ent->d_name, "sam") != NULL) {
	  array_list_insert(strdup(ent->d_name), files_found);
	  array_list_insert(SAM_FILE, files_type);
	} else if (strstr(ent->d_name, "bam") != NULL) {
	  array_list_insert(strdup(ent->d_name), files_found);
	  array_list_insert(BAM_FILE, files_type);
	}
      }
      closedir (dir);
    } else {
      printf("Directory not found\n");
      exit(-1);
    }
  }
  
    /*
      } else {
      printf("Aqui file\n");
      if (strstr(path_tmp, "sam") != NULL) {
      printf("SAM file\n");
      array_list_insert(strdup(path_tmp), files_found);
      array_list_insert(SAM_FILE, files_type);
      } else if (strstr(path_tmp, "bam") != NULL) {
      array_list_insert(strdup(path_tmp), files_found);
      array_list_insert(BAM_FILE, files_type);
      }
      }
    */
  //}

  int head_found = 0;
  printf("PROCESS FILE %s\n", path_tmp);
  if (rank == 0) {
    //First, search the head files
    for (int i = 0; i < array_list_size(files_found); i++) {
      char *file = array_list_get(i, files_found);
      int type   = array_list_get(i, files_type);
      
      char file_path[strlen(file) + strlen(node_out) + 1024];
      sprintf(file_path, "%s/%s", buffer_node_out, file);
      
      char line[4096];
      char *point;
      size_t nlines = 0;
      
      if (type == SAM_FILE) {
	FILE *fd = fopen(file_path, "r");
	while (fgets(line, 4096, fd)) {
	  if (line[0] == '@') {
	    if (line[1] == 'S' && line[2]  == 'Q') {
	      num_chromosomes_header++;
	    }
	    head_found = 1;
	    strcpy(&head_buffer[len_head_buffer], line);
	    len_head_buffer += strlen(line);
	    if (len_head_buffer >= MAX_HEAD_BUFFER) {
	      printf("ERROR MAX HEAD BUFFER OVERFLOW\n");
	      exit(-1);
	    }
	  } else {
	    break;
	  }
	  ///free(str);
	} // End while
	fclose(fd);
      } else 	if (type == BAM_FILE) { 
	//Read BAM
	bam_file_t* bam_file_p =  bam_fopen(file_path);
	
	if (bam_file_p->bam_header_p->text[0] == '@') {
	  head_found = 1;
	  if (bam_file_p->bam_header_p->l_text >= MAX_HEAD_BUFFER) {
	    printf("ERROR MAX HEAD BUFFER OVERFLOW\n");
	    exit(-1);
	  }
	  strcpy(head_buffer, bam_file_p->bam_header_p->text);	    
	}
	
	bam_fclose(bam_file_p);
	
      } //End if
      
      if (head_found) { break; }
      
    }
    
    if (head_found) {
      MPI_Send(head_buffer, MAX_HEAD_BUFFER, MPI_CHAR, w_rank, 1, MPI_COMM_WORLD);
    } else  {
      printf("ERROR Not head found\n");
      exit(-1);
    }
    
  } else if (rank == numprocs) { //rank == numprocs is the writer
    MPI_Recv(head_buffer, MAX_HEAD_BUFFER, MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
    fwrite(head_buffer, strlen(head_buffer), sizeof(char), fd_out);
  }
  
  free(head_buffer);
  
  printf("MERGE OUTPUT DONE! : %i\n", second_phase);
  
  //MPI_Barrier(MPI_COMM_WORLD);
  
  if (rank == 0) {
    start_timer(time_start);
  }

  char *buffer = (char *)malloc(sizeof(char)*MAX_BUFFER);
  size_t len_buffer = 0;
  FILE *fd_buffer;
  char buffer_path[strlen(buffer_node_out) + 1024];
  sprintf(buffer_path, "%s/reads_buffer.fq", buffer_node_out);    
  char line[4096], line_strtok[4096], line_tmp[4096];
  linked_list_t *buffer_no_map = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
  
  printf("[%i]MERGE DATA, tmp %s\n", rank, buffer_path);
  
  size_t nlines = 0;  
  total_writes = 0;

  MPI_Barrier(MPI_COMM_WORLD);
  
  //if (!second_phase) {
  if (rank != numprocs) {
    //First, search the head files
    fd_buffer = fopen(buffer_path, "w");
    
    for (int i = 0; i < array_list_size(files_found); i++) {
      char *file = array_list_get(i, files_found);
      int type   = array_list_get(i, files_type);
      
      char file_path[strlen(file) + strlen(node_out) + 1024];
      sprintf(file_path, "%s/%s", buffer_node_out, file);
      
      FILE *fd = fopen(file_path, "r");
      char *point;      
      
      printf("process file %s and type = %i, second_phase = %i\n", file_path, type, second_phase);
      
      if (type == SAM_FILE) {
	int fend = 1;
	while (fgets(line, 4096, fd)) {
	  if (line[0] != '@') { fend = 0; break; }
	}

	if (!fend) {	
	  nlines++;
	  strcpy(line_strtok, line);
	  strcpy(line_tmp, line);
	  
	  if (align_to_file(SAM_FILE, line_strtok, line_tmp, fd_buffer, second_phase,
			    avls_list, metaexons, genome, NULL, NULL, NULL)) {
	    //Update metaexon and AVL
	    strcpy(&buffer[len_buffer], line);
	    len_buffer += strlen(line);
	  }
	
	  while (fgets(line, 4096, fd)) {
	    nlines++;
	    strcpy(line_strtok, line);
	    strcpy(line_tmp, line);

	    if (align_to_file(SAM_FILE, line_strtok, line_tmp, fd_buffer, second_phase,
			      avls_list, metaexons, genome, NULL, NULL, NULL)) {
	      if (len_buffer + strlen(line) >= MAX_BUFFER) {
		MPI_Send(buffer, MAX_BUFFER, MPI_CHAR, w_rank, 1, MPI_COMM_WORLD);
		len_buffer = 0;
	      }
	      strcpy(&buffer[len_buffer], line);
	      len_buffer += strlen(line);	  
	    }
	  }
	} 
      } else if (type == BAM_FILE) { 
	//TODO: Read BAM	  
	bam_file_t* bam_file_p =  bam_fopen(file_path);
	int read_bytes;
	char *id_aux;
	bam1_t* bam_line = bam_init1();
	char SAM_line[4096];
	bam_header_t* bam_header = bam_file_p->bam_header_p;
	size_t nlines = 0;
	while ((read_bytes = bam_read1(bam_file_p->bam_fd, bam_line)) > 0) {
	  nlines++;
	  
	  if (align_to_file(BAM_FILE, line_strtok, line_tmp, fd_buffer, second_phase,
			    avls_list, metaexons, genome, bam_line, bam_header, SAM_line)) {
	    
	    if (len_buffer + strlen(line) >= MAX_BUFFER) {
	      MPI_Send(buffer, MAX_BUFFER, MPI_CHAR, w_rank, 1, MPI_COMM_WORLD);
	      len_buffer = 0;
	    } 
	    strcpy(&buffer[len_buffer], SAM_line);
	    len_buffer += strlen(SAM_line);	  
	  }
	  
	}
      }
    }

    printf("Finished MERGE 1\n");
    
    if (len_buffer) {
      MPI_Send(buffer, MAX_BUFFER, MPI_CHAR, w_rank, 1, MPI_COMM_WORLD);
    }
    
    strcpy(buffer, "END");
    MPI_Send(buffer, MAX_BUFFER, MPI_CHAR, w_rank, 1, MPI_COMM_WORLD);
    
    printf("Finished MERGE\n");
    
    fclose(fd_buffer);    

  } else { //rank == numprocs is the writer
    int numprocs_end = 0;
    while (1) {
      MPI_Recv(buffer, MAX_BUFFER, MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
      if (strcmp("END", buffer) == 0) {
	numprocs_end++;
	if (numprocs_end == numprocs) {
	  break;
	}
      } else {
	//printf("Write %lu\n", strlen(buffer));
	fwrite(buffer, strlen(buffer), sizeof(char), fd_out);
      }
    }  
  }
  
  
  if (rank == 0) {
    time_merge = 0;
    stop_timer(time_start, time_stop, time_merge);
  }

  
  //fclose(fd_buffer);
  MPI_Barrier(MPI_COMM_WORLD);
  //exit(-1);
  
  //printf("CLOSE with %i lines\n", total_writes);
  printf("END MERGE DATA OUTPUT : %i\n", second_phase);

  //MPI_Barrier(MPI_COMM_WORLD);
  
  //exit(-1);

  //MERGE METAEXON AND AVL STRUCTURE FOR SECOND PHASE      
  //create new communicator
  //--------------------------------------------------------------  
  MPI_Comm new_comm; 
  MPI_Group orig_group, new_group;
  MPI_Comm_group(MPI_COMM_WORLD, &orig_group); 
  
  /* Divide tasks into two distinct groups based upon rank */ 
  //if (rank != numprocs -1) {
  int ranks_group[numprocs];
  for (int i = 0; i < numprocs; i++) {
    ranks_group[i] = i;
  }
  MPI_Group_incl(orig_group, numprocs, ranks_group, &new_group);
  
  /* Create new communicator and then perform collective communications */
  MPI_Comm_create(MPI_COMM_WORLD, new_group, &new_comm);   
  //--------------------------------------------------------------

  size_t len_s = 1024;
  if (second_path_tmp != NULL) {
    len_s += strlen(second_path_tmp);
  }
  
  char second_path_out_tmp[len_s];
  char second_path_out_search[len_s];

  if (second_phase) {
    if (rank != numprocs) {
      //if (options->second_tmp_path) {
      sprintf(second_path_out_tmp, "%s/%i.second.out/", path_tmp, rank);
      sprintf(second_path_out_search, "%s/%i.second.out/", path_tmp, rank);

      printf("::::::::::::::::::::::::SECOND PHASE ENABLE:::::::::::::::: %s\n", second_path_out_tmp);

      char cmd[strlen(second_path_out_tmp) + 1024];
      sprintf(cmd, "rm -rf %s", second_path_out_tmp);
      system(cmd);

      create_directory(second_path_out_tmp);

      //} else {
      //sprintf(second_path_out_tmp, "%s", second_path_tmp);
      //}
      //printf("Create path %s\n", second_path_out_tmp);

    }
  }

  MPI_Barrier(MPI_COMM_WORLD);


  if (mapper_mode == RNA_MODE) {
    /////////////////////// FOR RNA SECOND PHASE //////////////////////////////////  
    if (second_phase) {
      
      char alig_second[strlen(second_path_out_tmp) + 1024];
      sprintf(alig_second, "%s/alig_second.sam", second_path_out_tmp);
      
      if ((rank != numprocs)) {
	if (hpg_enable) {
	  printf("MERGE METAEXON\n");
	  if (numprocs >= 2) {
	    //================================= MERGE METAEXON AND AVL ============================//             
	    unsigned long num_sj;        
	    int numprocs_merge = numprocs - 1;
	    int dec_s = 1;
	    int dec_r = 2;
	    int send, recv;
	    int dec = 2;
    
	    while (1) {
	      send = numprocs_merge - dec_s;
	      recv = numprocs_merge - dec_r;
	      if (send > 0 && recv < 0)
		recv = 0;
	      if (send <= 0 && recv <= 0)
		break;
	
	      for (int r = numprocs_merge - 1; r >= 0; r--) {
		if (rank == send) {
		  //SEND AVL
		  num_sj = get_total_items(num_chromosomes, avls_list);
		  unsigned long  *list_sj = (unsigned long *)malloc(sizeof(unsigned long)*num_sj*5); 	  
		  MPI_avl_package(num_chromosomes, avls_list, num_sj, list_sj);	  
	    
	  
		  //SEND METAEXON
		  unsigned long num_meta = 0;
		  unsigned long num_left_breaks = 0, num_right_breaks = 0;
		  for (int c = 0; c < metaexons->num_chromosomes; c++) {
		    linked_list_t *metaexon_list = metaexons->metaexons_list[c];
		    linked_list_item_t *item;
		    for (item = metaexon_list->first; item != NULL; item = item->next) {
		      metaexon_t *metaexon = item->item;
		      num_left_breaks  += array_list_size(metaexon->left_breaks);
		      num_right_breaks += array_list_size(metaexon->right_breaks);
		      num_meta++;	      
		    }
		  }
	  
		  //printf("Pack meta %i, %i, %i...\n", num_meta, num_left_breaks, num_right_breaks);
		  unsigned long *list_meta = (unsigned long *)malloc(sizeof(unsigned long)*num_meta*5);
		  unsigned long *list_left_breaks  = (unsigned long *)malloc(sizeof(unsigned long)*num_left_breaks*2);
		  unsigned long *list_right_breaks = (unsigned long *)malloc(sizeof(unsigned long)*num_right_breaks*2);
	  
		  MPI_metaexon_package(metaexons, num_meta, list_meta, 
				       num_left_breaks, list_left_breaks, 
				       num_right_breaks, list_right_breaks);
	  
		  //printf("RECV-RANK0: %i, %i, %i, %i, %i\n", list_meta[num_meta], list_meta[num_meta + 1], list_meta[num_meta + 2], list_meta[num_meta + 3], list_meta[num_meta + 4]);
	  
		  MPI_Send(&num_sj, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);	  
		  MPI_Send(&num_meta, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
		  MPI_Send(&num_left_breaks, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
		  MPI_Send(&num_right_breaks, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
		  //printf("[%i]SEND VALUES : %lu, %lu, %lu, %lu\n", rank, num_sj, num_meta, num_left_breaks, num_right_breaks);	  
	  
		  MPI_Send(list_sj, num_sj*5, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
		  MPI_Send(list_meta, num_meta*5, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
		  MPI_Send(list_left_breaks, num_left_breaks*2, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
		  MPI_Send(list_right_breaks, num_right_breaks*2, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  
		  //printf("[%i]SEND LISTS : %lu, %lu, %lu, %lu\n", rank, list_sj[num_sj], list_meta[num_meta], list_sj[num_sj + 1], list_meta[num_meta + 1]);
		  //printf("SEND: %i, %i, %i, %i, %i\n", list_right_breaks[0], list_right_breaks[1], list_right_breaks[2], list_right_breaks[3], list_right_breaks[4]);
	  
		  free(list_sj);
		  free(list_meta);
		  free(list_left_breaks);
		  free(list_right_breaks);
	  
		}
	
		if (rank == recv) {
		  //Recv AVLs
		  unsigned long recv_num_sj;	    
		  unsigned long recv_num_meta;
		  unsigned long recv_num_left_breaks;
		  unsigned long recv_num_right_breaks;
	  
		  //printf("RECV: %i, %i, %i, %i, %i\n", recv_list_meta[recv_num_meta], recv_list_meta[recv_num_meta + 1], recv_list_meta[recv_num_meta + 2], recv_list_meta[recv_num_meta + 3], recv_list_meta[recv_num_meta + 4]);	  
	    
		  MPI_Recv(&recv_num_sj, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);	
		  MPI_Recv(&recv_num_meta, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
		  MPI_Recv(&recv_num_left_breaks, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
		  MPI_Recv(&recv_num_right_breaks, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
		  //printf("[%i]RECV VALUES : %lu, %lu, %lu, %lu\n", rank, recv_num_sj, recv_num_meta, recv_num_left_breaks,recv_num_right_breaks);	  
	  
		  unsigned long *recv_list_sj = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_sj*5);
		  unsigned long *recv_list_meta = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_meta*5);
		  unsigned long *recv_list_left_breaks = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_left_breaks*2);
		  unsigned long *recv_list_right_breaks = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_right_breaks*2);
	  
		  MPI_Recv(recv_list_sj, recv_num_sj*5, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);	  	
		  MPI_Recv(recv_list_meta, recv_num_meta*5, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
		  MPI_Recv(recv_list_left_breaks, recv_num_left_breaks*2, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
		  MPI_Recv(recv_list_right_breaks, recv_num_right_breaks*2, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);
	  
	  
		  //Merge AVLs
		  avl_node_t *node_avl_start, *node_avl_end;
		  MPI_splice_t MPI_sj;
		  for (size_t sj = 0; sj < recv_num_sj; sj++) {
		    MPI_sj.start          = recv_list_sj[sj];
		    MPI_sj.end            = recv_list_sj[recv_num_sj + sj];
		    MPI_sj.reads_number   = recv_list_sj[recv_num_sj*2 + sj];
		    MPI_sj.strand         = recv_list_sj[recv_num_sj*3 + sj];
		    MPI_sj.chr            = recv_list_sj[recv_num_sj*4 + sj];
	    
		    //printf("\t[%i] %i/%i Merge (%i) %i: %lu - %lu\n", rank, sj, recv_num_sj, MPI_sj.strand, MPI_sj.chr, MPI_sj.start, MPI_sj.end);
		    allocate_start_node(MPI_sj.chr - 1,
					MPI_sj.strand,
					MPI_sj.start,
					MPI_sj.end,
					MPI_sj.start,
					MPI_sj.end,
					FROM_READ,
					1,
					"",
					&node_avl_start,
					&node_avl_end,
					avls_list);
		    //Refresh MPI_sj.reads_number in nodes
		    //end_data_t *end_data = node_avl_start->data;	    
		  }
		  //Merge Metaexon
		  //sleep(1);
	  
		  //printf("[%i]Merge %i META:\n", rank, recv_num_meta);
		  size_t l_pos = 0, r_pos = 0;
		  MPI_metaexon_t mpi_meta;
		  for (size_t m = 0; m < recv_num_meta; m++) {
		    mpi_meta.chromosome = recv_list_meta[m];
		    mpi_meta.start      = recv_list_meta[recv_num_meta   + m];
		    mpi_meta.end        = recv_list_meta[recv_num_meta*2 + m];
		    mpi_meta.n_starts   = recv_list_meta[recv_num_meta*3 + m];
		    mpi_meta.n_ends     = recv_list_meta[recv_num_meta*4 + m];
	    
		    //printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends);
		    MPI_breaks_t bk;
		    for (size_t s = 0; s < mpi_meta.n_starts; s++) {
		      bk.pos    = recv_list_right_breaks[r_pos];
		      bk.strand = recv_list_right_breaks[recv_num_right_breaks + r_pos];
		      r_pos++;
		      //MPI_breaks_t bk = MPI_breaks[se_pos++];
		      //printf("\tM-STARTS-(%i)%lu\n", bk.strand, bk.pos);
		      cp_avltree *avl = avls_list->avls[bk.strand][mpi_meta.chromosome].avl;
		      avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
		      assert(avl_node);
		      //if (!avl_node) {
		      //printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends);
		      //printf("%i:%i:%i\n", m, r_pos - 1, bk.pos);
		      //exit(-1);
		      //}
		      metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
				      METAEXON_RIGHT_END, avl_node, metaexons);
		    }
		
		    for (size_t s = 0; s < mpi_meta.n_ends; s++) {
		      bk.pos    = recv_list_left_breaks[l_pos];
		      bk.strand = recv_list_left_breaks[recv_num_left_breaks + l_pos];
		      l_pos++;
		      //MPI_breaks_t bk = MPI_breaks[se_pos++];
		      //printf("\t[%i]M-ENDS-(%i)%lu\n", rank, bk.strand, bk.pos);
		      cp_avltree *avl = avls_list->ends_avls[bk.strand][mpi_meta.chromosome].avl;
		      avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
		      assert(avl_node);
		      metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
				      METAEXON_LEFT_END, avl_node, metaexons);
		    }
		  }
	  
		  free(recv_list_sj);
		  free(recv_list_meta);
		  free(recv_list_left_breaks);
		  free(recv_list_right_breaks);
	  
		}
		send -= dec;
		recv -= dec;
		if (send <= 0) { send = -1; }
	      }
	      dec   *= 2;
	      dec_s *= 2;
	      dec_r *= 2;
	    }
	
	    printf("[%i]MERGE NODO 0 WORKFLOW 1 END...BROADCAST TO OTHER NODES\n", rank);
	
	    //int num_sj;
	    unsigned long recv_num_sj, recv_n_metaexons, recv_n_starts_ends;
    
	    unsigned long *list_sj;
	    unsigned long *list_meta;
	    unsigned long *list_left_breaks;
	    unsigned long *list_right_breaks;
    
	    unsigned long num_meta = 0;
	    unsigned long num_left_breaks = 0, num_right_breaks = 0;    
    
	    if (rank == 0) {
	      //PACK AVL
	      num_sj = get_total_items(num_chromosomes, avls_list);
	      list_sj = (unsigned long *)malloc(sizeof(unsigned long)*num_sj*5); 	  
	      MPI_avl_package(num_chromosomes, avls_list, num_sj, list_sj);	  
      
	      //PACK METAEXON
	      for (int c = 0; c < metaexons->num_chromosomes; c++) {
		linked_list_t *metaexon_list = metaexons->metaexons_list[c];
		linked_list_item_t *item;
		for (item = metaexon_list->first; item != NULL; item = item->next) {
		  metaexon_t *metaexon = item->item;
		  num_left_breaks  += array_list_size(metaexon->left_breaks);
		  num_right_breaks += array_list_size(metaexon->right_breaks);
		  num_meta++;	      
		}
	      }
      
	      list_meta = (unsigned long *)malloc(sizeof(unsigned long)*num_meta*5);
	      list_left_breaks  = (unsigned long *)malloc(sizeof(unsigned long)*num_left_breaks*2);
	      list_right_breaks = (unsigned long *)malloc(sizeof(unsigned long)*num_right_breaks*2);
      
	      MPI_metaexon_package(metaexons, num_meta, list_meta, 
				   num_left_breaks, list_left_breaks, 
				   num_right_breaks, list_right_breaks);
	    }
    
	    MPI_Bcast(&num_sj, 1, MPI_UNSIGNED_LONG, 0, new_comm);
	    MPI_Bcast(&num_meta, 1, MPI_UNSIGNED_LONG, 0, new_comm);
	    MPI_Bcast(&num_left_breaks, 1, MPI_UNSIGNED_LONG, 0, new_comm);
	    MPI_Bcast(&num_right_breaks, 1, MPI_UNSIGNED_LONG, 0, new_comm);
    
	    size_t xn = num_sj*5 + num_meta*5 + num_left_breaks*2 + num_right_breaks*2;
	    unsigned long *pack_send = (unsigned long *)malloc(sizeof(unsigned long)*xn);
    
	    if (rank == 0) {
	      memcpy(pack_send, list_sj, sizeof(unsigned long)*num_sj*5);
	      memcpy(&pack_send[num_sj*5], list_meta, sizeof(unsigned long)*num_meta*5);
	      memcpy(&pack_send[num_sj*5 + num_meta*5], list_left_breaks, sizeof(unsigned long)*num_left_breaks*2);
	      memcpy(&pack_send[num_sj*5 + num_meta*5 + num_left_breaks*2], list_right_breaks, sizeof(unsigned long)*num_right_breaks*2);
      
	      free(list_sj);
	      free(list_meta);
	      free(list_left_breaks);
	      free(list_right_breaks);
	    }
    
	    MPI_Bcast(pack_send, xn, MPI_UNSIGNED_LONG, 0, new_comm);
    
	    list_sj           = pack_send;
	    list_meta         = &pack_send[num_sj*5];
	    list_left_breaks  = &pack_send[num_sj*5 + num_meta*5];
	    list_right_breaks = &pack_send[num_sj*5 + num_meta*5 + num_left_breaks*2];
    
	    //Merge AVLs
	    avl_node_t *node_avl_start, *node_avl_end;
	    MPI_splice_t MPI_sj;
    
	    for (size_t sj = 0; sj < num_sj; sj++) {
	      MPI_sj.start          = list_sj[sj];
	      MPI_sj.end            = list_sj[num_sj + sj];
	      MPI_sj.reads_number   = list_sj[num_sj*2 + sj];
	      MPI_sj.strand         = list_sj[num_sj*3 + sj];
	      MPI_sj.chr            = list_sj[num_sj*4 + sj];
      
	      //printf("[%i]Merge 2 (%i) %i: %lu - %lu\n", rank, MPI_sj.strand, MPI_sj.chr, MPI_sj.start, MPI_sj.end);
	      allocate_start_node(MPI_sj.chr - 1,
				  MPI_sj.strand,
				  MPI_sj.start,
				  MPI_sj.end,
				  MPI_sj.start,
				  MPI_sj.end,
				  FROM_READ,
				  1,
				  "",
				  &node_avl_start,
				  &node_avl_end,
				  avls_list);
	      //Refresh MPI_sj.reads_number in nodes
	      end_data_t *end_data = node_avl_start->data;	    
	    }
      
	    //printf("[%i]Merge Meta %i\n", rank, num_meta);
      
	    //Merge Metaexon
	    size_t se_pos = 0;
	    size_t r_pos = 0, l_pos = 0;
	    MPI_metaexon_t mpi_meta;
    
	    for (size_t m = 0; m < num_meta; m++) {
	      mpi_meta.chromosome = list_meta[m];
	      mpi_meta.start      = list_meta[num_meta   + m];
	      mpi_meta.end        = list_meta[num_meta*2 + m];
	      mpi_meta.n_starts   = list_meta[num_meta*3 + m];
	      mpi_meta.n_ends     = list_meta[num_meta*4 + m];
	
	      //printf("[%i]M-%i:%lu-%lu  (%i|%i):\n", rank, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, mpi_meta.n_starts, mpi_meta.n_ends);
	      MPI_breaks_t bk;
	      for (int s = 0; s < mpi_meta.n_starts; s++) {
		bk.pos    = list_right_breaks[r_pos];
		bk.strand = list_right_breaks[num_right_breaks + r_pos];
		r_pos++;
	  
		//printf("\tM-STARTS-(%i)%lu\n", bk.strand, bk.pos);
		cp_avltree *avl = avls_list->avls[bk.strand][mpi_meta.chromosome].avl;
		avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
		//assert(avl_node);
		if (!avl_node) {
		  printf("ERROR: (%i|%i)%lu\n", m, s, bk.pos);
		  exit(-1);
		}
		metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
				METAEXON_RIGHT_END, avl_node, metaexons);
	      }
	      for (int s = 0; s < mpi_meta.n_ends; s++) {
		bk.pos    = list_left_breaks[l_pos];
		bk.strand = list_left_breaks[num_left_breaks + l_pos];
		l_pos++;
	  
		//printf("\t[%i]M-ENDS-(%i)%lu\n", rank, bk.strand, bk.pos);
		cp_avltree *avl = avls_list->ends_avls[bk.strand][mpi_meta.chromosome].avl;
		avl_node_t *avl_node = (avl_node_t *)cp_avltree_get(avl, (void *)bk.pos);
		assert(avl_node);
		metaexon_insert(0, mpi_meta.chromosome, mpi_meta.start, mpi_meta.end, 40,
				METAEXON_LEFT_END, avl_node, metaexons);
	      }
	    }
      
	    printf("Frees merge\n");
	    free(pack_send);
      
	    //====                E N D    M E R G E                ====//
	    //==========================================================//
	  }	

	
	  //================================= INPUT INITIALIZATIONS =============================//  
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
	  fastq_batch_reader_input_t reader_input;
	  fastq_batch_reader_input_init(options->in_filename, options->in_filename2, 
					options->pair_mode, options->batch_size, 
					NULL, options->gzip, &reader_input);  
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
			       cal_optarg, bwt_index, metaexons, NULL, NULL, 
			       NULL, NULL, pair_mode, &sw_input);
	  pair_server_input_t pair_input;
	  pair_server_input_init(pair_mng, report_optarg, NULL, NULL, NULL, &pair_input);
	  batch_writer_input_t writer_input;
	  batch_writer_input_init(NULL,
				  exact_filename, 
				  NULL, 
				  alignments_list, 
				  genome, 
				  &writer_input);
	  writer_input.bam_format = 0;      
	  batch_t *batch = batch_new(&bwt_input, &region_input, &cal_input, 
				     &pair_input, &preprocess_rna, &sw_input, &writer_input, RNA_MODE,
				     NULL, NULL);
	  fastq_batch_reader_input_init(NULL, NULL, 
					options->pair_mode, options->batch_size, 
					NULL, options->gzip, &reader_input);
	
	  writer_input.bam_file = (bam_file_t *) fopen(alig_second, "w"); 
	  reader_input.fq_file1 = fastq_fopen(buffer_path);

	  char buffer_sa_path[strlen(second_path_out_tmp) + 1024];
	  sprintf(buffer_sa_path, "%s/buffer_sa.tmp", second_path_out_tmp);    
	  f_sa = fopen(buffer_sa_path, "w+b");
	  if (f_sa == NULL) {
	    LOG_FATAL("Error opening file 'buffer_sa.tmp' \n");
	  }

	  char buffer_hc_path[strlen(second_path_out_tmp) + 1024];
	  sprintf(buffer_hc_path, "%s/buffer_hc.tmp", second_path_out_tmp);    
	  f_hc = fopen(buffer_hc_path, "w+b");
	  if (f_hc == NULL) {
	    LOG_FATAL("Error opening file 'buffer_hc.tmp' \n");
	  }
      
	  sw_input.f_sa = f_sa;
	  sw_input.f_hc = f_hc;
      
	  if (!fast_mode) {
	    ///////////////// BWT INDEX WORKFLOW //////////////////////
	    //===================================================================================
	    //-----------------------------------------------------------------------------------
	    wf_input_t *wf_input = wf_input_new(&reader_input, batch);
	    wf_input_file_t *wf_input_file = wf_input_file_new(f_sa, batch);   
	    wf_input_file_t *wf_input_file_hc = wf_input_file_new(f_hc, batch);  

	    //create and initialize workflow
	    workflow_t *wf = workflow_new();
	    workflow_stage_function_t stage_functions[] = {rna_wq};
	    char *stage_labels[] = {"RNA"};
	    workflow_set_stages(1, (workflow_stage_function_t *)&stage_functions, stage_labels, wf);

	    // optional producer and consumer functions
	    workflow_set_producer((workflow_producer_function_t *)fastq_reader, "FastQ reader", wf);
	    workflow_set_consumer((workflow_consumer_function_t *)sam_writer, "SAM writer", wf);
      
	    workflow_t *wf_last = workflow_new();
	    workflow_stage_function_t stage_functions_last[] = {rna_wq_w2};
	    char *stage_labels_last[] = {"RNA LAST STAGE"};
	    workflow_set_stages(1, (workflow_stage_function_t *)&stage_functions_last, stage_labels_last, wf_last);      
	    workflow_set_producer((workflow_producer_function_t *)file_reader, "Buffer reader", wf_last);
	    workflow_set_consumer((workflow_consumer_function_t *)sam_writer, "SAM writer", wf_last);

	    workflow_t *wf_hc = workflow_new();
	    workflow_stage_function_t stage_functions_hc[] = {rna_wq_w3};
	    char *stage_labels_hc[] = {"RNA HARD CLIPPINGS"};
	    workflow_set_stages(1, (workflow_stage_function_t *)&stage_functions_hc, stage_labels_hc, wf_hc);
	    workflow_set_producer((workflow_producer_function_t *)file_reader_2, "Buffer reader", wf_hc);
	    workflow_set_consumer((workflow_consumer_function_t *)sam_writer, "SAM writer", wf_hc);

	    printf("[%i]BWT run workflow %s...\n", rank, buffer_path);
	    workflow_run_with(options->num_cpu_threads, wf_input, wf);
	    printf("W1 END \n");	  
	    rewind(f_sa);
	    workflow_run_with(options->num_cpu_threads, wf_input_file, wf_last);
	    rewind(f_hc);
	    workflow_run_with(options->num_cpu_threads, wf_input_file_hc, wf_hc);
	    printf("[%i]BWT run workflow end\n", rank);

	    // free memory
	    workflow_free(wf);
	    //workflow_free(wf_last);
	    //workflow_free(wf_hc);	
	    wf_input_free(wf_input);
	    //wf_input_file_free(wf_input_file);
	    //wf_input_file_free(wf_input_file_hc);
	
	  } else {
	    ///////////////// SA INDEX WORKFLOW //////////////////////
	    //--------------------------------------------------------------------------------------
	    // workflow management
	    //
	    sw_optarg_t sw_optarg;
	    sw_optarg_init(options->gap_open, options->gap_extend, 
			   options->match, options->mismatch, &sw_optarg);
      
	    sa_rna_input_t sa_rna;
	    sa_rna.cal_optarg   = cal_optarg;
	    sa_rna.genome    = genome;
	    sa_rna.avls_list = avls_list;
	    sa_rna.metaexons = metaexons;
	    sa_rna.sw_optarg = &sw_optarg;
	    sa_rna.file1 = f_sa;
	    sa_rna.file2 = f_hc;
	    sa_rna.pair_input = &pair_input;
	    sa_rna.min_score = options->min_score;
	    sa_rna.max_alig = options->filter_read_mappings;

	    sa_wf_batch_t *wf_batch = sa_wf_batch_new(NULL, (void *)sa_index, &writer_input, NULL, &sa_rna);
	    sa_wf_input_t *wf_input = sa_wf_input_new(options->bam_format, &reader_input, wf_batch);
      
	    // create and initialize workflow
	    workflow_SA_t *wf = workflow_SA_new();      
	    workflow_stage_function_SA_t stage_functions[] = {sa_rna_mapper};
	    char *stage_labels[] = {"SA mapper"};
	    workflow_set_stages_SA(1, stage_functions, stage_labels, wf);      
	    // optional producer and consumer functions

	    workflow_set_producer_SA((workflow_producer_function_SA_t *)sa_fq_reader_rna, "FastQ reader", wf);
	    workflow_set_consumer_SA((workflow_consumer_function_SA_t *)write_to_file, "SAM writer", wf);      

	    //Create and initialize second workflow
	    workflow_SA_t *wf_last = workflow_SA_new();
	    workflow_stage_function_SA_t stage_functions_last[] = {sa_rna_mapper_last};
	    char *stage_labels_last[] = {"SA mapper last stage"};
	    workflow_set_stages_SA(1, (workflow_stage_function_SA_t *)&stage_functions_last, stage_labels_last, wf_last);
	    workflow_set_producer_SA((workflow_producer_function_SA_t *)sa_alignments_reader_rna, "FastQ reader", wf_last);      
	    workflow_set_consumer_SA((workflow_consumer_function_SA_t *)write_to_file, "SAM writer", wf_last);

	    printf("[%i]SA run workflow...\n", rank);
	    workflow_run_with_SA(options->num_cpu_threads, wf_input, wf);      
	    rewind(f_sa);
	    workflow_run_with_SA(options->num_cpu_threads, wf_input, wf_last);
	    printf("[%i]SA run workflow end\n", rank);
	  
	    // free memory
	    sa_wf_input_free(wf_input);
	    sa_wf_batch_free(wf_batch);
	    workflow_SA_free(wf);      
	    workflow_SA_free(wf_last);            
	  }
      
      
	  //FREE MEMORY
	  //metaexons_free(metaexons);      
	  bwt_optarg_free(bwt_optarg);
	  cal_optarg_free(cal_optarg);
	  pair_mng_free(pair_mng);
	  report_optarg_free(report_optarg);
	  batch_free(batch);
	  fclose(writer_input.bam_file);

	  if (!fast_mode) {
	    // free memory
	    if (bwt_index) { bwt_index_free(bwt_index); }
	  } else {
	    // free memory
	    if (sa_index)  { sa_index3_free(sa_index);  }
	  }   
    
	  //if (genome) {
	  //genome_free(genome);
	  //}    

	} else { //hpg_enable ???
	  //system("CLI second aligner");

	  //================= SPLIT INPUT CLI AND TRANSFORM IT ==========================//
	  max_len_cli = strlen(second_cli) + strlen(path_output) + strlen(tmp_input_path) + 1024;
	  char *cli_tmp = (char *)malloc(sizeof(char)*max_len_cli);
    
	  cli_in  = (char *)malloc(sizeof(char)*max_len_cli);
	  cli_out = (char *)malloc(sizeof(char)*max_len_cli);
	  cli_end = (char *)malloc(sizeof(char)*max_len_cli);
    
	  size_t len_cli = 0;
	  unsigned char find_mark = 0, first_out = 0;
  
	  //if (rank == 0) {
	  char *str = second_cli;
	  char *pch;
    
	  pch = strtok (str," ");
	  while (pch != NULL) {
	    if (strcmp(pch, "%I") == 0 || strcmp(pch, "%i") == 0) {	
	      strcpy(cli_in, cli_tmp);
	      len_cli = 0;
	      cli_tmp[0] = '\0';
	      find_mark++;
	      first_in = !first_out;
	    } else if (strcmp(pch, "%O") == 0 || strcmp(pch, "%o") == 0) {
	      strcpy(cli_out, cli_tmp);
	      len_cli = 0;
	      cli_tmp[0] = '\0';
	      find_mark++;
	      first_out = !first_in;
	    } else {
	      strcpy(cli_tmp + len_cli, pch);
	      len_cli = strlen(cli_tmp);
	      strcpy(cli_tmp + len_cli, " \0");
	      len_cli++;
	    }
    
	    pch = strtok(NULL, " ");
    
	  }
    
	  strcpy(cli_end, cli_tmp);
    
	  if (find_mark != 2) {
	    if (rank == 0) {
	      printf("ERROR in CLI Mapper\n");
	      exit(-1);
	    }  
	  }
    
	  //================================ MPI PROCESS =======================================//
	  //COMMON OUTPUT PATH 

	  //PROCESS FILES
    
	  if (rank != numprocs) {
	    char *mapper_run = (char *)malloc(sizeof(char)*max_len_cli);
	    //sprintf(node_out, "%s/%i.out/", path_tmp, rank);      
	    
	    if (options->second_tmp_file) {
	      sprintf(second_path_out_tmp, "%s/%s", second_path_out_tmp, options->second_tmp_file);
	    } 
	    
	    if (first_in) {
	      sprintf(mapper_run, "%s %s %s %s %s", cli_in, buffer_path, cli_out, second_path_out_tmp, cli_end);
	    } else {
	      sprintf(mapper_run, "%s %s %s %s %s", cli_out, second_path_out_tmp, cli_in, buffer_path, cli_end);
	    }

	    printf("%s\n", mapper_run);
	    //exit(-1);
	  
	    system(mapper_run);
      
	  }
	  //===========================================================================//               		
	}
      }
       
      printf("[%i]MERGE NODO 0 WORKFLOW 1 END...BROADCAST TO OTHER NODES\n", rank);
        
      //
      // end of workflow management
      //--------------------------------------------------------------------------------------
    }
   } else { 

    /////////////////////// FOR DNA SECOND PHASE //////////////////////////////////
    
    if (second_phase) {
      if (second_cli) {
	printf("SECOND PHASE ENABLE\n");
	//================= SPLIT INPUT CLI AND TRANSFORM IT ==========================//
	max_len_cli = strlen(second_cli) + strlen(path_output) + strlen(tmp_input_path) + 1024;
	char *cli_tmp = (char *)malloc(sizeof(char)*max_len_cli);
      
	cli_in  = (char *)malloc(sizeof(char)*max_len_cli);
	cli_out = (char *)malloc(sizeof(char)*max_len_cli);
	cli_end = (char *)malloc(sizeof(char)*max_len_cli);
      
	size_t len_cli = 0;
	unsigned char find_mark = 0, first_out = 0;
      
	//if (rank == 0) {
	char *str = second_cli;
	char *pch;
      
	pch = strtok (str," ");
	while (pch != NULL) {
	  if (strcmp(pch, "%I") == 0 || strcmp(pch, "%i") == 0) {	
	    strcpy(cli_in, cli_tmp);
	    len_cli = 0;
	    cli_tmp[0] = '\0';
	    find_mark++;
	    first_in = !first_out;
	  } else if (strcmp(pch, "%O") == 0 || strcmp(pch, "%o") == 0) {
	    strcpy(cli_out, cli_tmp);
	    len_cli = 0;
	    cli_tmp[0] = '\0';
	    find_mark++;
	    first_out = !first_in;
	  } else {
	    strcpy(cli_tmp + len_cli, pch);
	    len_cli = strlen(cli_tmp);
	    strcpy(cli_tmp + len_cli, " \0");
	    len_cli++;
	  }
	
	  pch = strtok(NULL, " ");
	
	}
    
	strcpy(cli_end, cli_tmp);
      
	if (find_mark != 2) {
	  if (rank == 0) {
	    printf("ERROR in CLI Mapper\n");
	    exit(-1);
	  }  
	}
      
	//================================ MPI PROCESS =======================================//
	//COMMON OUTPUT PATH 
      
	//PROCESS FILES
      
	if (rank != numprocs) {
	  char *mapper_run = (char *)malloc(sizeof(char)*max_len_cli);
	  //sprintf(node_out, "%s/%i.out/", path_tmp, rank);      

	  if (options->second_tmp_file) {
	    sprintf(second_path_out_tmp, "%s/%s", second_path_out_tmp, options->second_tmp_file);
	  } 
	  
	  if (first_in) {
	    sprintf(mapper_run, "%s %s %s %s %s", cli_in, buffer_path, cli_out, second_path_out_tmp, cli_end);
	  } else {
	    sprintf(mapper_run, "%s %s %s %s %s", cli_out, second_path_out_tmp, cli_in, buffer_path, cli_end);
	  }
	
	  printf("%s\n", mapper_run);
	  //exit(-1);
	
	  system(mapper_run);
	
	}
	//===========================================================================//	
      } else { //DNA SA
	//NOTHING
      }
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  
  if (second_phase) {
    //MERGE OUTPUT
    //MPI_Barrier(MPI_COMM_WORLD);  
    //printf("[%i]Finish\n", rank);
    //exit(-1);    
    //printf("[%i]MERGE ALIGNMENTS MPI %s\n", rank, alig_second);

    //====================================================================================//
    // READ OUTPUT SAM
    //printf("%s:\n", node_out);
    DIR *dir;
    struct dirent *ent;
    array_list_t *files_found_2  = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    array_list_t *files_type_2   = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

    //printf("XXXXXX: %s\n", second_path_out_tmp);

    //if (options->tmp_path) {
    if (rank != numprocs) {
      if ((dir = opendir (second_path_out_search)) != NULL) {
	while ((ent = readdir (dir)) != NULL) {
	  //printf ("\t%s\n", ent->d_name);
	  if (strstr(ent->d_name, "sam") != NULL) {
	    array_list_insert(strdup(ent->d_name), files_found_2);
	    array_list_insert(SAM_FILE, files_type_2);
	  } else if (strstr(ent->d_name, "bam") != NULL) {
	    array_list_insert(strdup(ent->d_name), files_found_2);
	    array_list_insert(BAM_FILE, files_type_2);
	  }
	}
	closedir (dir);
      } else {
	printf("Directory not found\n");
	exit(-1);
      }
    }
    
    /*
    //} else {
      if (strstr(path_tmp, "sam") != NULL) {
	array_list_insert(strdup(path_tmp), files_found_2);
	array_list_insert(SAM_FILE, files_type_2);
      } else if (strstr(ent->d_name, "bam") != NULL) {
	array_list_insert(strdup(path_tmp), files_found_2);
	array_list_insert(BAM_FILE, files_type_2);
      }
    }
    */
    
    printf("************* [%i] After first loop\n", rank);
    
    second_phase = 0;
    if (rank != numprocs) {
      len_buffer= 0;
      for (int i = 0; i < array_list_size(files_found_2); i++) {
	char *file = array_list_get(i, files_found_2);
	int type   = array_list_get(i, files_type_2);

	char file_path[strlen(file) + strlen(second_path_out_tmp) + 1024];
	sprintf(file_path, "%s/%s", second_path_out_search, file);

	FILE *fd = fopen(file_path, "r");
	char *point;      
      
	printf("process file %s and type = %i, second_phase = %i\n", file_path, type, second_phase);
      
	if (type == SAM_FILE) {
	  int fend = 1;
	  while (fgets(line, 4096, fd)) {
	    if (line[0] != '@') { fend = 0; break; }
	  }

	  if (!fend) {	
	    nlines++;
	    strcpy(line_strtok, line);
	    strcpy(line_tmp, line);
	  
	    if (align_to_file(SAM_FILE, line_strtok, line_tmp, fd_buffer, second_phase,
			      avls_list, metaexons, genome, NULL, NULL, NULL)) {
	      //Update metaexon and AVL
	      strcpy(&buffer[len_buffer], line);
	      len_buffer += strlen(line);
	    }
	
	    while (fgets(line, 4096, fd)) {
	      nlines++;
	      strcpy(line_strtok, line);
	      strcpy(line_tmp, line);

	      if (align_to_file(SAM_FILE, line_strtok, line_tmp, fd_buffer, second_phase,
				avls_list, metaexons, genome, NULL, NULL, NULL)) {
		if (len_buffer + strlen(line) >= MAX_BUFFER) {
		  MPI_Send(buffer, MAX_BUFFER, MPI_CHAR, w_rank, 1, MPI_COMM_WORLD);
		  len_buffer = 0;
		} 
		strcpy(&buffer[len_buffer], line);
		len_buffer += strlen(line);	  
	      }
	    }
	  } 
	} else if (type == BAM_FILE) { 
	  //TODO: Read BAM	  
	  bam_file_t* bam_file_p =  bam_fopen(file_path);
	  int read_bytes;
	  char *id_aux;
	  bam1_t* bam_line = bam_init1();
	  char SAM_line[4096];
	  bam_header_t* bam_header = bam_file_p->bam_header_p;
	  size_t nlines = 0;
	  while ((read_bytes = bam_read1(bam_file_p->bam_fd, bam_line)) > 0) {
	    nlines++;
	  
	    if (align_to_file(BAM_FILE, line_strtok, line_tmp, fd_buffer, second_phase,
			      avls_list, metaexons, genome, bam_line, bam_header, SAM_line)) {
	    
	      if (len_buffer + strlen(line) >= MAX_BUFFER) {
		MPI_Send(buffer, MAX_BUFFER, MPI_CHAR, w_rank, 1, MPI_COMM_WORLD);
		len_buffer = 0;
	      } 
	      strcpy(&buffer[len_buffer], SAM_line);
	      len_buffer += strlen(SAM_line);	  
	    }
	  
	  }
	}
      }

      printf("Finished MERGE 1\n");
    
      if (len_buffer) {
	MPI_Send(buffer, MAX_BUFFER, MPI_CHAR, w_rank, 1, MPI_COMM_WORLD);
      }
    
      strcpy(buffer, "END");
      MPI_Send(buffer, MAX_BUFFER, MPI_CHAR, w_rank, 1, MPI_COMM_WORLD);
    
      printf("Finished MERGE\n");
    
      //fclose(fd_buffer);    
      
    } else { //rank == numprocs is the writer
      int numprocs_end = 0;
      while (1) {
	MPI_Recv(buffer, MAX_BUFFER, MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
	if (strcmp("END", buffer) == 0) {
	  numprocs_end++;
	  if (numprocs_end == numprocs) {
	    break;
	  }
	} else {
	  printf("Write %lu\n", strlen(buffer));
	  fwrite(buffer, strlen(buffer), sizeof(char), fd_out);
	}
      }
    }
  }
  
  if (rank != numprocs) {    
    if (mapper_mode == RNA_MODE) {
      printf("MERGE AVL\n");  
      //============================================================//
      //= M E R G E    M E T A E X O N - A V L    T O    W R I T E =//
      //============================================================//
      //== Package structures ==//    
      unsigned long num_sj;        
      int numprocs_merge = numprocs - 1;
      int dec_s = 1;
      int dec_r = 2;
      int send, recv;
      int dec = 2;
    
      while (1) {
	send = numprocs_merge - dec_s;
	recv = numprocs_merge - dec_r;
	if (send > 0 && recv < 0)
	  recv = 0;
	if (send <= 0 && recv <= 0)
	  break;
      
	for (int r = numprocs_merge - 1; r >= 0; r--) {
	  if (rank == send) {
	    //SEND AVL
	    num_sj = get_total_items(num_chromosomes, avls_list);
	    unsigned long  *list_sj = (unsigned long *)malloc(sizeof(unsigned long)*num_sj*5); 	  
	    MPI_avl_package(num_chromosomes, avls_list, num_sj, list_sj);	  
	  
	    MPI_Send(&num_sj, 1, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);	  	    
	    MPI_Send(list_sj, num_sj*5, MPI_UNSIGNED_LONG, recv, tag, MPI_COMM_WORLD);
	  
	    free(list_sj);
	  }
	
	  if (rank == recv) {
	    //Recv AVLs
	    unsigned long recv_num_sj; 	    
	    MPI_Recv(&recv_num_sj, 1, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);		    
	    unsigned long *recv_list_sj = (unsigned long *)malloc(sizeof(unsigned long)*recv_num_sj*5);	    
	    MPI_Recv(recv_list_sj, recv_num_sj*5, MPI_UNSIGNED_LONG, send, tag, MPI_COMM_WORLD, &status);	  	
	  
	  
	    //Merge AVLs
	    avl_node_t *node_avl_start, *node_avl_end;
	    MPI_splice_t MPI_sj;
	    for (size_t sj = 0; sj < recv_num_sj; sj++) {
	      MPI_sj.start          = recv_list_sj[sj];
	      MPI_sj.end            = recv_list_sj[recv_num_sj + sj];
	      MPI_sj.reads_number   = recv_list_sj[recv_num_sj*2 + sj];
	      MPI_sj.strand         = recv_list_sj[recv_num_sj*3 + sj];
	      MPI_sj.chr            = recv_list_sj[recv_num_sj*4 + sj];
	    
	      //printf("\t[%i] %i/%i Merge (%i) %i: %lu - %lu\n", rank, sj, recv_num_sj, MPI_sj.strand, MPI_sj.chr, MPI_sj.start, MPI_sj.end);
	      allocate_start_node(MPI_sj.chr - 1,
				  MPI_sj.strand,
				  MPI_sj.start,
				  MPI_sj.end,
				  MPI_sj.start,
				  MPI_sj.end,
				  FROM_READ,
				  1,
				  "",
				  &node_avl_start,
				  &node_avl_end,
				  avls_list);
	      //Refresh MPI_sj.reads_number in nodes
	      //end_data_t *end_data = node_avl_start->data;	    
	    }
	    //Merge Metaexon
	    //sleep(1);
	    free(recv_list_sj);
	  
	  }
	  send -= dec;
	  recv -= dec;
	  if (send <= 0) { send = -1; }
	}
	dec   *= 2;
	dec_s *= 2;
	dec_r *= 2;
      }
    
      if (rank == 0) {
	printf("WRITE AVL in %s\n", exact_filename);
	write_chromosome_avls(NULL, exact_filename, num_chromosomes, avls_list);
	printf("WRITE OK\n");
      }

    }    
  } else {    
    fclose(fd_out);    
  }

  
  MPI_Barrier(MPI_COMM_WORLD);
  printf("FINISHHHH!! :)\n");
  
  if (rank == 0) {
    stop_timer(total_start, total_stop, time_total);
    //stop_timer(time_start, time_stop, time_merge);    
    printf("============ T    I    M    I    N    G ===============\n");
    printf("Time Split : %0.2f\n", time_split / 1000000);
    printf("Time Merge : %0.2f\n", time_merge / 1000000);
    printf("Total Time : %0.2f\n", time_total / 1000000);
    printf("=======================================================\n");    
  }

  options_free(options);
  
  //free(buffer);
  //fwrite(recv_buff, strlen(recv_buff), sizeof(char), fd_out);
  
  //printf("[%i]%s\n", rank, buffer);

  //====================================================================================//
  
  
  //================= MPI for input aligner ====================//
  //printf("[%i] Run Tasks %i\n", rank, ntasks);
  //============================================================//  
  
  MPI_Finalize();
  
  //printf("[%i] Finished\n", rank);
  
  return 1;
  
}



  /*

    //First, search the head files
    for (int i = 0; i < array_list_size(files_found_2); i++) {
      char *file = array_list_get(i, files_found_2);
      int type   = array_list_get(i, files_type_2);
	
      char file_path[strlen(file) + strlen(second_path_out_tmp) + 1024];
      sprintf(file_path, "%s/%s", second_path_out_tmp, file);
	
      char line[4096];
      char *point;
      size_t nlines = 0;
      
      second_phase = 0;
      if (rank != numprocs) {
	//First, search the head files
	if (type == SAM_FILE) {
	  len_buffer = 0;
	  
	  FILE *fd_buffer = fopen(alig_second, "r");
	  if (!fd_buffer) { exit(-1); }
	  
	  while (fgets(line, 4096, fd_buffer)) {
	    strcpy(line_strtok, line);
	    strcpy(line_tmp, line);
	    //printf("%s\n", line);
	    if (align_to_file(SAM_FILE, line_strtok, line_tmp, fd_buffer, second_phase,
			      avls_list, metaexons, genome, NULL, NULL, NULL)) {
	      //Update metaexon and AVL
	      strcpy(&buffer[len_buffer], line);
	      len_buffer += strlen(line);
	    }

	    if (len_buffer + strlen(line) >= MAX_BUFFER) {
	      printf("[%i]SEND MERGE\n", rank);
	      MPI_Send(buffer, MAX_BUFFER, MPI_CHAR, w_rank, 1, MPI_COMM_WORLD);
	      len_buffer = 0;
	    }
	  } //end for
	} else {
	  bam_file_t* bam_file_p =  bam_fopen(file_path);
	  int read_bytes;
	  char *id_aux;
	  bam1_t* bam_line = bam_init1();
	  char SAM_line[4096];
	  bam_header_t* bam_header = bam_file_p->bam_header_p;
	  
	  while ((read_bytes = bam_read1(bam_file_p->bam_fd, bam_line)) > 0) {    	    
	    if (align_to_file(BAM_FILE, line_strtok, line_tmp, fd_buffer, second_phase,
			      avls_list, metaexons, genome, bam_line, bam_header, SAM_line)) {	      
	      if (len_buffer + strlen(line) >= MAX_BUFFER) {
		MPI_Send(buffer, MAX_BUFFER, MPI_CHAR, w_rank, 1, MPI_COMM_WORLD);
		len_buffer = 0;
	      } 
	      strcpy(&buffer[len_buffer], SAM_line);
	      len_buffer += strlen(SAM_line);	  
	    }	    
	  }		  
	}

	if (len_buffer) {
	  MPI_Send(buffer, MAX_BUFFER, MPI_CHAR, w_rank, 1, MPI_COMM_WORLD);
	}
      
	strcpy(buffer, "END");
	MPI_Send(buffer, MAX_BUFFER, MPI_CHAR, w_rank, 1, MPI_COMM_WORLD);
	
	fclose(fd_buffer);
      
      } else { //rank == numprocs is the writer
	int numprocs_end = 0;
	while (1) {
	  MPI_Recv(buffer, MAX_BUFFER, MPI_CHAR, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
	  if (strcmp("END", buffer) == 0) {
	    numprocs_end++;
	    if (numprocs_end == numprocs) {
	      break;
	    }
	  } else {
	    fwrite(buffer, strlen(buffer), sizeof(char), fd_out);
	    printf("[%i]Write Done!\n", rank);
	  } 
	}
	printf("END WHILE WRITER\n");
      }
    }
  */

