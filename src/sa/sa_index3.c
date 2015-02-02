#include "sa_index3.h"

#include "options.h"
 
#define PROGRESS 1000000
#define MAX_LINE_LENGTH 4096

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

alt_names_t *alt_names_new(char *alt_filename) {
  // open alternative names file
  FILE *f = fopen(alt_filename, "r");
  if (f == NULL) {
    printf("Error reading alternative names filename %s\n", alt_filename);
    exit(EXIT_FAILURE);
  }

  size_t count, num_alts = 0;
  char line[MAX_LINE_LENGTH], alt_name[MAX_LINE_LENGTH], chrom_name[MAX_LINE_LENGTH];
  while (fgets(line, MAX_LINE_LENGTH, f) != NULL) {
    num_alts++;
  }
  fseek(f, 0, SEEK_SET);	

  char **alt_names = (char **) calloc(num_alts, sizeof(char *));
  char **chrom_names = (char **) calloc(num_alts, sizeof(char *));
  count = 0;
  while (fgets(line, MAX_LINE_LENGTH, f) != NULL) {
    sscanf(line, "%s\t%s\n", alt_name, chrom_name);
    alt_names[count] = strdup(alt_name);
    chrom_names[count] = strdup(chrom_name);
    count++;
  }

  fclose(f);
  
  alt_names_t *p = (alt_names_t *) calloc(1, sizeof(alt_names_t));
  p->size = num_alts;
  p->alt_names = alt_names;
  p->chrom_names = chrom_names;
  return p;
}

//--------------------------------------------------------------------------------------

void alt_names_free(alt_names_t *p) {
  if (p) {
    for (size_t i = 0; i < p->size; i++) {
      if (p->alt_names && p->alt_names[i]) free(p->alt_names[i]);
      if (p->chrom_names && p->chrom_names[i]) free(p->chrom_names[i]);
    }
    if (p->alt_names) free(p->alt_names);
    if (p->chrom_names) free(p->chrom_names);
    free(p);
  }
}

//--------------------------------------------------------------------------------------

int alt_names_exists(char *alt_name, alt_names_t *p) {
  if (p) {
    for (size_t i = 0; i < p->size; i++) {
      if (strcmp(alt_name, p->alt_names[i]) == 0) {
	return 1;
      }
    }
  }
  return 0;
}

//--------------------------------------------------------------------------------------

char *alt_names_get_chrom_name(char *alt_name, alt_names_t *p) {
  if (p) {
    for (size_t i = 0; i < p->size; i++) {
      if (strcmp(alt_name, p->alt_names[i]) == 0) {
	return p->chrom_names[i];
      }
    }
  }
  return NULL;
}

//--------------------------------------------------------------------------------------

char *alt_names_display(alt_names_t *p) {
  if (p) {
    printf("%lu\n", p->size);
    for (size_t i = 0; i < p->size; i++) {
      printf("%lu\t%s\t%s\n", i, p->alt_names[i], p->chrom_names[i]);
    }
  }
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

int get_seq_flag(char *seq_name, alt_names_t *alt_names, array_list_t *decoy_names) {
  if (alt_names) {
    if (alt_names_exists(seq_name, alt_names)) {
      return ALT_FLAG;
    }
  }
  for (int i = 0; i < array_list_size(decoy_names); i++) {
    if (strcmp(seq_name, array_list_get(i, decoy_names)) == 0) {
      return DECOY_FLAG;
    }
  }
  return CHROM_FLAG;
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

sa_genome3_t *read_genome3(char *genome_filename) {
  return read_genome3_alt(genome_filename, NULL, NULL);
}

//--------------------------------------------------------------------------------------

sa_genome3_t *read_genome3_alt(char *genome_filename, char *alt_filename,
			       char *decoy_filename) {

  const int MAX_CHROM_NAME_LENGHT = 1024;
  uint reading_name, seq_name_count = 0;
  char seq_name[MAX_CHROM_NAME_LENGHT];

  alt_names_t *alt_names = NULL;
  if (alt_filename) {
    alt_names = alt_names_new(alt_filename);
    //alt_names_display(alt_names);
  }

  array_list_t *decoy_names = array_list_new(100, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  FILE *f, *genome_f, *decoy_f = NULL;

  // decoy
  size_t decoy_length = 0;
  if (decoy_filename) {
    decoy_f = fopen(decoy_filename, "r");
    if (decoy_f == NULL) {
      printf("Error reading decoy filename %s\n", decoy_filename);
      exit(EXIT_FAILURE);
    } 

    // get file length
    fseek(decoy_f, 0, SEEK_END);
    decoy_length = ftell(decoy_f);
    printf("file %s length = %lu\n", decoy_filename, decoy_length);
    fseek(decoy_f, 0, SEEK_SET);	

    get_seq_names(decoy_filename, decoy_names);
  }

  // open fasta genome file
  genome_f = fopen(genome_filename, "r");
  if (genome_f == NULL) {
    printf("Error reading genome filename %s\n", genome_filename);
    exit(EXIT_FAILURE);
  }

  // get file length
  fseek(genome_f, 0, SEEK_END);
  size_t file_length = ftell(genome_f);
  printf("file %s length = %lu\n", genome_filename, file_length);
  fseek(genome_f, 0, SEEK_SET);	

  // read genome
  size_t num_A = 0, num_C = 0, num_G = 0, num_N = 0, num_T = 0;
  size_t seq_length = 0;
  size_t num_seqs = 0, num_allocated_seqs = 100;
  int *seq_flags = (int *) calloc(num_allocated_seqs, sizeof(int));
  size_t *seq_chroms = (size_t *) calloc(num_allocated_seqs, sizeof(size_t));
  size_t *seq_starts = (size_t *) calloc(num_allocated_seqs, sizeof(size_t));
  size_t *seq_ends = (size_t *) calloc(num_allocated_seqs, sizeof(size_t));
  size_t *seq_lengths = (size_t *) calloc(num_allocated_seqs, sizeof(size_t));
  char **seq_names = (char **) calloc(num_allocated_seqs, sizeof(char *));
  char *S = (char *) calloc(file_length + decoy_length + 1, sizeof(char));

  int flag, skip = 0, hap = 0, initial_flank = 1;
  size_t l = 0, process = 0, hap_padding = 0;

  char c;

  seq_length = 0;

  for (int k = 0; k < 2; k++) {

    // init
    skip = 0;
    hap = 0;
    initial_flank = 1;
    process = 0;
    hap_padding = 0;

    f = (k == 0 ? genome_f : decoy_f);
    if (f == NULL) break;

    while ((c = getc(f)) != EOF) {
      process++;
      //    if (process % PROGRESS == 0) printf("reading %0.2f %c...\n", 100.0f * process / file_length, '%'); 
      //printf("len %lu: c = %c\n", len, c);
      if (c == '>') {
	printf("%c", c);
	if (skip == 0) {
	  num_seqs++;
	  reading_name = 1;
	  // to do: get and save sequence name
	  if (num_seqs > 1) {
	    if (hap && hap_padding) {
	      seq_ends[num_seqs - 2] = hap_padding;
	      num_N -= hap_padding; seq_length -= hap_padding; l -= hap_padding;
	    }
	    seq_lengths[num_seqs - 2] = seq_length;
	    
	    printf("setting %s seq %lu: length = %lu (but num_seqs = %lu)\n", seq_name, num_seqs - 2, seq_length, num_seqs);
	    hap = 0;
	    hap_padding = 0;
	    initial_flank = 1;
	    seq_length = 0;
	    if (num_seqs >= num_allocated_seqs) {
	    num_allocated_seqs += 50;
	    
	    seq_flags = (int *) realloc(seq_flags, num_allocated_seqs * sizeof(int));
	    seq_chroms = (size_t *) realloc(seq_chroms, num_allocated_seqs * sizeof(size_t));
	    seq_starts = (size_t *) realloc(seq_starts, num_allocated_seqs * sizeof(size_t));
	    seq_ends = (size_t *) realloc(seq_ends, num_allocated_seqs * sizeof(size_t));
	    
	    seq_lengths = (size_t *) realloc(seq_lengths, num_allocated_seqs * sizeof(size_t));
	    seq_names = (char **) realloc(seq_names, num_allocated_seqs * sizeof(char *));
	    }
	  }
	}
	skip = 1;
      } else if (c == '\n') {
	if (reading_name) {
	  seq_name[seq_name_count] = 0;
	  seq_names[num_seqs - 1] = strdup(seq_name);
	  
	  flag = get_seq_flag(seq_name, alt_names, decoy_names);
	  hap = (flag == ALT_FLAG ? 1 : 0);
	  seq_flags[num_seqs - 1] = flag;
	  
	  seq_name_count = 0;
	  reading_name = 0;
	}
	if (skip == 1) {
	  printf("%c", c);
	}
	skip = 0;
      } else if (c == ' ' || c == '\t') {
	printf("%c", c);
	if (reading_name) {
	  seq_name[seq_name_count] = 0;
	  seq_names[num_seqs - 1] = strdup(seq_name);
	  
	  flag = get_seq_flag(seq_name, alt_names, decoy_names);
	  hap = (flag == ALT_FLAG ? 1 : 0);
	  seq_flags[num_seqs - 1] = flag;
	  
	  seq_name_count = 0;
	  reading_name = 0;
	}
      } else {
	if (skip == 0) {
	  if (c == 'A' || c == 'a') { 
	    // A
	    S[l++] = 'A'; num_A++; seq_length++; 
	    if (hap && hap_padding && initial_flank) {
	      seq_starts[num_seqs - 1] = hap_padding;
	    }
	    hap_padding = 0;
	    initial_flank = 0;
	  } else if (c == 'C' || c == 'c') { 
	    // C
	    S[l++] = 'C'; num_C++; seq_length++; 
	    if (hap && hap_padding && initial_flank) {
	      seq_starts[num_seqs - 1] = hap_padding;
	    }
	    hap_padding = 0;
	    initial_flank = 0;
	  } else if (c == 'G' || c == 'g') { 
	    // G
	    S[l++] = 'G'; num_G++; seq_length++; 
	    if (hap && hap_padding && initial_flank) {
	      seq_starts[num_seqs - 1] = hap_padding;
	    }
	    hap_padding = 0;
	    initial_flank = 0;
	  } else if (c == 'T' || c == 't') { 
	    // T
	    S[l++] = 'T'; num_T++; seq_length++; 
	    if (hap && hap_padding && initial_flank) {
	      seq_starts[num_seqs - 1] = hap_padding;
	    }
	    hap_padding = 0;
	    initial_flank = 0;
	  } else if (c == 'N' || c == 'n') { 
	    // N
	    hap_padding++;
	    if (!hap || (hap && !initial_flank)) {
	      S[l++] = 'N'; num_N++; seq_length++; 
	    }
	  } else {
	    printf("Unknown character %c at %lu position\n", c, process);
	  }
	} else {
	  if (reading_name) {
	    seq_name[seq_name_count++] = c;
	    if (seq_name_count >= MAX_CHROM_NAME_LENGHT) {
	      seq_name[MAX_CHROM_NAME_LENGHT - 1] = 0;
	      printf("Sequence name (%s) exceeds max. length (%u)", seq_name, MAX_CHROM_NAME_LENGHT);
	      exit(EXIT_FAILURE);
	    }
	  }
	  printf("%c", c);
	}
      }
    }
    if (hap && hap_padding) {
      seq_ends[num_seqs - 1] = hap_padding;
      num_N -= hap_padding; seq_length -= hap_padding; l -= hap_padding;
    }
    seq_lengths[num_seqs - 1] = seq_length;
    printf("last: setting %s seq %lu: length = %lu (but num_seqs = %lu)\n", seq_name, num_seqs - 2, seq_length, num_seqs);

    fclose(f);
  }

  S[l++] = '$';

  printf("...genome reading done!\n");

  // update chromosomes for ALT sequences: chromosomes and flanks
  char *chrom_name;
  for (size_t i = 0; i < num_seqs; i++) {
    seq_chroms[i] = i;
    if (seq_flags[i] == ALT_FLAG) {
      chrom_name = alt_names_get_chrom_name(seq_names[i], alt_names);
      if (chrom_name) {
	for (size_t j = 0; j < num_seqs; j++) {
	  if (strcmp(chrom_name, seq_names[j]) == 0) {
	    // set chrom
	    seq_chroms[i] = j;
	    break;
	  }
	}
      }
    }
  }

  // free memory
  if (alt_names) alt_names_free(alt_names);

  // create sa_genome3_t and return it
  sa_genome3_t *genome = sa_genome3_new(l, num_seqs, seq_lengths, seq_flags, 
					seq_chroms, seq_starts, seq_ends, 
					seq_names, S);
  sa_genome3_set_nt_counters(num_A, num_C, num_G, num_N, num_T, genome);
  return genome;
}

//--------------------------------------------------------------------------------------

char *global_S;

//--------------------------------------------------------------------------------------

typedef struct suffix_tmp {
  uint value;
  unsigned short int seq;
} suffix_tmp_t;


int suffix_tmp_cmp3(void const *a, void const *b) { 
  suffix_tmp_t *item_a = (suffix_tmp_t *) a;
  suffix_tmp_t *item_b = (suffix_tmp_t *) b;

  return (strncmp(&global_S[item_a->value], &global_S[item_b->value], 1000));
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

//--------------------------------------------------------------------------------------
// create a SA index from genome
//--------------------------------------------------------------------------------------

void sa_index3_build(char *genome_filename, uint k_value, char *sa_index_dirname) {

  FILE *f_tab;
  char filename_tab[strlen(genome_filename) + strlen(sa_index_dirname) + 100];
  struct timeval stop, start;

  // getting prefix
  char *prefix = strrchr(genome_filename, '/');
  if (prefix == NULL) {
    prefix = genome_filename;
  } else {
    prefix++;
  }

  //-----------------------------------------
  // read genome S from file
  //-----------------------------------------

  printf("\nreading file genome %s...\n", genome_filename);
  gettimeofday(&start, NULL);
  sa_genome3_t *genome = read_genome3(genome_filename);
  gettimeofday(&stop, NULL);

  //  sa_genome3_display(genome);

  // write S to file
  sprintf(filename_tab, "%s/%s.S", sa_index_dirname, prefix);
  f_tab = fopen(filename_tab, "wb");
  if (f_tab == NULL) {
    printf("Error: could not open %s to write\n", filename_tab);
    exit(EXIT_FAILURE);
  }
  fwrite(genome->S, sizeof(char), genome->length, f_tab);
  fclose(f_tab);

  global_S = genome->S;

  //-----------------------------------------
  // compute SA table
  //-----------------------------------------
  uint num_suffixes = 0;
  printf("\ncomputing SA and CHROM table...\n");
  gettimeofday(&start, NULL);
  suffix_tmp_t *tmp[4];
  uint tmp_A = 0, tmp_C = 0, tmp_G = 0, tmp_T = 0, nts[4];

  nts[0] = genome->num_A;
  nts[1] = genome->num_C;
  nts[2] = genome->num_G;
  nts[3] = genome->num_T;

  for (int i = 0; i < 4; i++) {
    num_suffixes += nts[i];
    tmp[i] = (suffix_tmp_t *) malloc(nts[i] * sizeof(suffix_tmp_t));
  }

  char nt;
  uint c = 0;
  for (uint i = 0; i < genome->num_seqs; i++) {
    for (uint j = 0; j < genome->seq_lengths[i]; j++) {
      nt = genome->S[c];
      if (nt == 'A' || nt == 'a') {
	tmp[0][tmp_A].value = c;
	tmp[0][tmp_A].seq = (unsigned short int) i;
	tmp_A++;
      } else if (nt == 'C' || nt == 'c') {
	tmp[1][tmp_C].value = c;
	tmp[1][tmp_C].seq = (unsigned short int) i;
	tmp_C++;
      } else if (nt == 'G' || nt == 'g') {
	tmp[2][tmp_G].value = c;
	tmp[2][tmp_G].seq = (unsigned short int) i;
	tmp_G++;
      } else if (nt == 'T' || nt == 't') {
	tmp[3][tmp_T].value = c;
	tmp[3][tmp_T].seq = (unsigned short int) i;
	tmp_T++;
      }
      c++;
    }
  }

  // sanity check
  if (tmp_A != genome->num_A) {
    printf("mismatch: A's (genome, S) = (%lu, %lu)\n", genome->num_A, tmp_A);
    exit(EXIT_FAILURE);
  }
  if (tmp_T != genome->num_T) {
    printf("mismatch: T's (genome, S) = (%lu, %lu)\n", genome->num_T, tmp_T);
    exit(EXIT_FAILURE);
  }
  if (tmp_G != genome->num_G) {
    printf("mismatch: G's (genome, S) = (%lu, %lu)\n", genome->num_G, tmp_G);
    exit(EXIT_FAILURE);
  }
  if (tmp_C != genome->num_C) {
    printf("mismatch: C's (genome, S) = (%lu, %lu)\n", genome->num_C, tmp_C);
    exit(EXIT_FAILURE);
  }

  // sort
  #pragma omp parallel for num_threads(4)
  for (int i = 0; i < 4; i++) {
    qsort(tmp[i], nts[i], sizeof(suffix_tmp_t), suffix_tmp_cmp3);
  }

  gettimeofday(&stop, NULL);
  printf("end of computing SA and CHROM table in %0.2f s\n", 
	 (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  

  // write SA to file
  sprintf(filename_tab, "%s/%s.SA", sa_index_dirname, prefix);
  f_tab = fopen(filename_tab, "wb");
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < nts[i]; j++) {
      fwrite(&(tmp[i][j].value), sizeof(uint), 1, f_tab);
    }
  }
  fclose(f_tab);

  // write CHROM to file
  sprintf(filename_tab, "%s/%s.CHROM", sa_index_dirname, prefix);
  f_tab = fopen(filename_tab, "wb");
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < nts[i]; j++) {
      fwrite(&(tmp[i][j].seq), sizeof(unsigned short int), 1, f_tab);
    }
  }
  fclose(f_tab);

  //-----------------------------------------
  // compute PRE table
  //-----------------------------------------
  printf("\ncomputing PRE table...\n");
  gettimeofday(&start, NULL);

  PREFIX_TABLE_NT_VALUE['A'] = 0;
  PREFIX_TABLE_NT_VALUE['N'] = 0;
  PREFIX_TABLE_NT_VALUE['C'] = 1;
  PREFIX_TABLE_NT_VALUE['G'] = 2;
  PREFIX_TABLE_NT_VALUE['T'] = 3;

  uint pre_length = 1LLU << (2 * k_value);
  uint *PRE = (uint *) calloc(pre_length, sizeof(uint));

  size_t value, count = 0;
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < nts[i]; j++) {
      value = compute_prefix_value(&genome->S[tmp[i][j].value], k_value);
      if (PRE[value] == 0) {
	PRE[value] = count;
      }
      count++;
    }
  }
  gettimeofday(&stop, NULL);
  printf("end of computing PRE table in %0.2f s\n", 
	 (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  

  // write PRE to file for the next time
  sprintf(filename_tab, "%s/%s.PRE", sa_index_dirname, prefix);
  f_tab = fopen(filename_tab, "wb");
  if (f_tab == NULL) {
    printf("Error: could not open %s to write\n", filename_tab);
    exit(EXIT_FAILURE);
  } 
  fwrite(PRE, sizeof(uint), pre_length, f_tab);
  fclose(f_tab);
    
  //-----------------------------------------
  // save parameters
  //-----------------------------------------
  sprintf(filename_tab, "%s/params.txt", sa_index_dirname);
  f_tab = fopen(filename_tab, "w");
  fprintf(f_tab, "%s\n", prefix);
  fprintf(f_tab, "%i\n", k_value);
  fprintf(f_tab, "%i\n", pre_length);
  fprintf(f_tab, "%i\n", num_suffixes);
  fprintf(f_tab, "%lu\n", genome->length);
  fprintf(f_tab, "%lu\n", genome->num_seqs);
  for (size_t i = 0; i < genome->num_seqs; i++) {
    fprintf(f_tab, "%s\t%lu\n", 
	   (genome->seq_names ? genome->seq_names[i] : "no-name"), 
	    genome->seq_lengths[i]);
  }
  fclose(f_tab);
}

//--------------------------------------------------------------------------------------

void sa_index3_build_k18(char *genome_filename, uint k_value, char *sa_index_dirname) {
  sa_index3_build_k18_alt(genome_filename, NULL, NULL, k_value, sa_index_dirname);
}

//--------------------------------------------------------------------------------------

void sa_index3_build_k18_alt(char *genome_filename, char *alt_filename, char *decoy_filename,
			     uint k_value, char *sa_index_dirname) {

  //printf("\n***************** K value = 18 ***************************\n");
  k_value = 18;

  const size_t value4M = 4194304; // 4 M
  const size_t value16M = 16777216; // 16 M
  const size_t value64M = 67108864; // 64 M
  const size_t value256M = 268435456; // 256 M
  
  size_t M_items = 256LLU * value16M;
  size_t M_bytes = M_items * sizeof(uint);
  printf("allocating matrix M (%lu bytes: %lu items)\n", M_bytes, M_items);
  uint *M = (uint *) malloc(M_bytes);
  if (M == NULL) {
    printf("Error allocating memory for M matrix\n");
    exit(EXIT_FAILURE);
  }
  printf("end of allocating matrix M\n");

  PREFIX_TABLE_NT_VALUE['A'] = 0;
  PREFIX_TABLE_NT_VALUE['N'] = 0;
  PREFIX_TABLE_NT_VALUE['C'] = 1;
  PREFIX_TABLE_NT_VALUE['G'] = 2;
  PREFIX_TABLE_NT_VALUE['T'] = 3;

  FILE *f_tab;
  char filename_tab[strlen(genome_filename) + strlen(sa_index_dirname) + 100];
  struct timeval stop, start;

  // getting prefix
  char *prefix = strrchr(genome_filename, '/');
  if (prefix == NULL) {
    prefix = genome_filename;
  } else {
    prefix++;
  }

  //-----------------------------------------
  // read genome S from file
  //-----------------------------------------
  sprintf(filename_tab, "%s/%s.S", sa_index_dirname, prefix);
  printf("\nreading file genome %s...\n", genome_filename);
  gettimeofday(&start, NULL);
  sa_genome3_t *genome = read_genome3_alt(genome_filename, alt_filename, decoy_filename);
  gettimeofday(&stop, NULL);
  
  if (genome->length > MAX_GENOME_LENGTH || genome->num_seqs > MAX_NUM_SEQUENCES) {
    printf("Genome not supported: (%s)\n", genome_filename);
    printf("\tGenome length: %lu\n", genome->length);
    printf("\tNumber segments (chromosomes, scaffolds,...): %i\n", genome->num_seqs);
    printf("\n");
    printf("Genomes supported by HPG Aligner %s\n", HPG_ALIGNER_VERSION);
    printf("\tMax. genome length: %lu\n", MAX_GENOME_LENGTH);
    printf("\tMax. number segments (chromosomes, scaffolds,...): %i\n", MAX_NUM_SEQUENCES);
    exit(EXIT_FAILURE);
  }
  
  //  sa_genome3_display(genome);

  // write S to file
  f_tab = fopen(filename_tab, "wb");
  if (f_tab == NULL) {
    printf("Error: could not open %s to write\n", filename_tab);
    exit(EXIT_FAILURE);
  }
  fwrite(genome->S, sizeof(char), genome->length, f_tab);
  fclose(f_tab);

  //-----------------------------------------
  // compute SA table
  //-----------------------------------------

  char nt;
  uint A_counter = 0, IA_counter = 0;
  size_t num_suffixes = genome->num_A + genome->num_C + genome->num_G + genome->num_T;

  printf("\ncomputing SA and CHROM table...\n");
  gettimeofday(&start, NULL);
  suffix_tmp_t *tmp[4];
  size_t tmp_A = 0, tmp_C = 0, tmp_G = 0, tmp_T = 0, nts[4];
  
  nts[0] = genome->num_A;
  nts[1] = genome->num_C;
  nts[2] = genome->num_G;
  nts[3] = genome->num_T;
  
  for (int i = 0; i < 4; i++) {
    tmp[i] = (suffix_tmp_t *) malloc(nts[i] * sizeof(suffix_tmp_t));
  }
  
  size_t c = 0;
  for (size_t i = 0; i < genome->num_seqs; i++) {
    for (size_t j = 0; j < genome->seq_lengths[i]; j++) {
      nt = genome->S[c];
      if (nt == 'A' || nt == 'a') {
	tmp[0][tmp_A].value = c;
	tmp[0][tmp_A].seq = (unsigned short int) i;
	tmp_A++;
      } else if (nt == 'C' || nt == 'c') {
	tmp[1][tmp_C].value = c;
	tmp[1][tmp_C].seq = (unsigned short int) i;
	tmp_C++;
      } else if (nt == 'G' || nt == 'g') {
	tmp[2][tmp_G].value = c;
	tmp[2][tmp_G].seq = (unsigned short int) i;
	tmp_G++;
      } else if (nt == 'T' || nt == 't') {
	tmp[3][tmp_T].value = c;
	tmp[3][tmp_T].seq = (unsigned short int) i;
	tmp_T++;
      }
      c++;
    }
  }
  
  // sanity check
  if (tmp_A != genome->num_A) {
    printf("mismatch: A's (genome, S) = (%lu, %lu)\n", genome->num_A, tmp_A);
    exit(EXIT_FAILURE);
  }
  if (tmp_T != genome->num_T) {
    printf("mismatch: T's (genome, S) = (%lu, %lu)\n", genome->num_T, tmp_T);
    exit(EXIT_FAILURE);
  }
  if (tmp_G != genome->num_G) {
    printf("mismatch: G's (genome, S) = (%lu, %lu)\n", genome->num_G, tmp_G);
    exit(EXIT_FAILURE);
  }
  if (tmp_C != genome->num_C) {
    printf("mismatch: C's (genome, S) = (%lu, %lu)\n", genome->num_C, tmp_C);
    exit(EXIT_FAILURE);
  }
  
  for (size_t i = 0; i < genome->length; i++) {
    if (genome->S[i] == 'N' || genome->S[i] == 'n') {
      genome->S[i] = 'A';
    }
  }
  global_S = genome->S;
  
  // sort
  #pragma omp parallel for num_threads(4)
  for (size_t i = 0; i < 4; i++) {
    qsort(tmp[i], nts[i], sizeof(suffix_tmp_t), suffix_tmp_cmp3);
  }
  
  gettimeofday(&stop, NULL);
  printf("end of computing SA and CHROM table in %0.2f s\n", 
	 (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  
  // write SA to file
  sprintf(filename_tab, "%s/%s.SA", sa_index_dirname, prefix);
  f_tab = fopen(filename_tab, "wb");
  for (size_t i = 0; i < 4; i++) {
    for (size_t j = 0; j < nts[i]; j++) {
      fwrite(&(tmp[i][j].value), sizeof(uint), 1, f_tab);
    }
  }
  fclose(f_tab);
  
  // write CHROM to file
  sprintf(filename_tab, "%s/%s.CHROM", sa_index_dirname, prefix);
  f_tab = fopen(filename_tab, "wb");
  for (size_t i = 0; i < 4; i++) {
    for (size_t j = 0; j < nts[i]; j++) {
      fwrite(&(tmp[i][j].seq), sizeof(unsigned short int), 1, f_tab);
    }
  }
  fclose(f_tab);
  
  // free memory
  for (int i = 0; i < 4; i++) {
    free(tmp[i]);
  }

  sprintf(filename_tab, "%s/%s.SA", sa_index_dirname, prefix);
  f_tab = fopen(filename_tab, "rb");
  if (f_tab == NULL) {
    printf("Error: could not open %s to write\n", filename_tab);
    exit(EXIT_FAILURE);
  }
  printf("SA: filename %s, num_suffixes = %lu\n", filename_tab, num_suffixes);
  uint *SA = (uint *) malloc(num_suffixes * sizeof(uint));
  
  printf("\nreading SA table from file %s...\n", filename_tab);
  gettimeofday(&start, NULL);
  uint num_items = fread(SA, sizeof(uint), num_suffixes, f_tab);
  if (num_items != num_suffixes) {
    printf("Error: (%s) mismatch num_items = %i vs num_suffixes = %lu\n", 
	   filename_tab, num_items, num_suffixes);
    exit(EXIT_FAILURE);
  }
  gettimeofday(&stop, NULL);
  printf("end of reading SA table (%lu num_suffixes) from file %s in %0.2f s\n", 
  	 num_suffixes, filename_tab,
	 (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);      
  fclose(f_tab);

  //-----------------------------------------
  // compute Compressed Row Storage tables
  //-----------------------------------------
  printf("\ncomputing Compressed Row Storage tables...\n");
  gettimeofday(&start, NULL);

  size_t num_prefixes = 0;

  // A vector
  sprintf(filename_tab, "%s/%s.A", sa_index_dirname, prefix);
  FILE *f_A = fopen(filename_tab, "wb");

  // IA vector
  sprintf(filename_tab, "%s/%s.IA", sa_index_dirname, prefix);
  FILE *f_IA = fopen(filename_tab, "wb");
  uint first_in_row;

  // JA vectors, they'll be merged into a single JA vector
  sprintf(filename_tab, "%s/%s.JA", sa_index_dirname, prefix);
  FILE *f_JA = fopen(filename_tab, "wb");

  unsigned char ja;
  uint row, col, m_value;
  
  size_t matrix, matrix_items, value, min_value, max_value;
  // set matrix to 0
  matrix = 0;
  matrix_items = 0;
  memset(M, 255, M_bytes);

  // start filling matrix
  //printf("filling matrix (%lu)...\n", matrix);
  for (uint i = 0; i < num_suffixes; i++) {
    value = compute_prefix_value(&genome->S[SA[i]], k_value);
    //printf("value = %lu, (i = %i, j = %i)\n", value, i, j);
    if (value / M_items == matrix) { 
      if (M[value % M_items] == max_uint) {
	//printf("\t\tvalue = %lu -> (modulo: %lu) setting to %lu\n", value, value % M_items, count);
	M[value % M_items] = i;
	num_prefixes++;
	matrix_items++;
      }
    } else {
      // printf("end of filling sub-matrix (%lu): %lu items. Done !!\n", matrix, matrix_items);
      //printf("\t -----> new matrix due to the value SA[%lu] = %lu -> ", i, value);
      //display_prefix(&genome->S[SA[i]], k_value);
      //printf("\n");
      
      // read matrix and creating A, IA and JA vectors
      //printf("reading matrix %lu and updating CRS vectors...\n", matrix);
      for (row = 0; row < value16M; row++) {
	first_in_row = 1;
	for (col = 0; col < 256; col++) {
	  //printf("accesing to M[%lu]...\n", 256 * row + col);
	  m_value = M[256 * row + col];
	  //printf("accesing to M[%lu] = %lu\n", 256 * row + col, m_value);
	  if (m_value != max_uint) {
	    // handling A
	    //printf("A[%lu]: %lu\n", A_counter, m_value);
	    //fprintf(f_A, "%lu\n", m_value);
	    fwrite(&m_value, sizeof(uint), 1, f_A);
	    
	    // handling IA
	    if (first_in_row) {
	      //fprintf(f_IA, "%lu\n", A_counter);
	      fwrite(&A_counter, sizeof(uint), 1, f_IA);
	      IA_counter++;
	      first_in_row = 0;
	    }
	    
	    // handling JA
	    ja = col;
	    //	      fprintf(f_JA, "%lu\n", col);
	    fwrite(&ja, sizeof(unsigned char), 1, f_JA);
	    
	    A_counter++;
	  }
	}
	if (first_in_row == 1) {
	  uint max = max_uint;
	  fwrite(&max, sizeof(uint), 1, f_IA);
	  //fprintf(f_IA, "%lu\n", max_uint);
	  IA_counter++;
	  
	  //printf("**** matrix %u, row %u without elements !!!\n", matrix, row);
	  //	    exit(-1);
	}
      }

      //printf("end of reading matrix %lu and updating CRS vectors. Done !!\n", matrix);
      //	exit(-1);

      // initialize matrix to re-fill it
      matrix = value / M_items;
      matrix_items = 0;
      memset(M, 255, M_bytes);
      
      //printf("init and filling matrix (%lu)...\n", matrix);
      M[value % M_items] = i;
      matrix_items++;
    }
  }

  // read matrix and creating A, IA and JA vectors
  //printf("reading matrix %lu and updating CRS vectors...\n", matrix);
  for (row = 0; row < value16M; row++) {
    first_in_row = 1;
    for (col = 0; col < 256; col++) {
      m_value = M[256 * row + col];
      if (m_value != max_uint) {
	// handling A
	//fprintf(f_A, "%lu\n", m_value);
	fwrite(&m_value, sizeof(uint), 1, f_A);
	
	// handling IA
	if (first_in_row) {
	  //fprintf(f_IA, "%lu\n", A_counter);
	  fwrite(&A_counter, sizeof(uint), 1, f_IA);
	  IA_counter++;
	  first_in_row = 0;
	}
	
	// handling JA
	//fprintf(f_JA, "%lu\n", col);
	nt = col;
	fwrite(&nt, sizeof(char), 1, f_JA);
	
	A_counter++;
      }
    }
    if (first_in_row == 1) {
      uint max = max_uint;
      fwrite(&max, sizeof(uint), 1, f_IA);
      //fprintf(f_IA, "%lu\n", max_uint);
      IA_counter++;
      //printf("**** matrix %u, row %u without elements !!!\n", matrix, row);
      //	    exit(-1);
    }
  }

  //printf("end of reading matrix %lu and updating CRS vectors. Done !!\n", matrix);

  printf("A length = %u, IA length = %u (num. prefixes = %lu)\n", A_counter, IA_counter, num_prefixes);
  fclose(f_A);
  fclose(f_IA);
  fclose(f_JA);

  //-----------------------------------------
  // save parameters
  //-----------------------------------------
  sprintf(filename_tab, "%s/params.txt", sa_index_dirname);
  f_tab = fopen(filename_tab, "w");
  fprintf(f_tab, "%s\n", prefix);
  fprintf(f_tab, "%i\n", k_value);
  fprintf(f_tab, "0\n");
  fprintf(f_tab, "%i\n", A_counter);
  fprintf(f_tab, "%i\n", IA_counter);
  fprintf(f_tab, "%lu\n", num_suffixes);
  fprintf(f_tab, "%lu\n", genome->length);
  fprintf(f_tab, "%lu\n", genome->num_seqs);
  for (size_t i = 0; i < genome->num_seqs; i++) {
    fprintf(f_tab, "%s\t%lu\t%i\t%lu\t%lu\t%lu\t%lu\t%lu\n", 
	    (genome->seq_names ? genome->seq_names[i] : "no-name"), 
	    genome->seq_lengths[i],
	    genome->seq_flags[i],
	    genome->seq_chroms[i],
	    genome->seq_starts[i],
	    genome->seq_ends[i],
	    genome->left_flanks[i],
	    genome->right_flanks[i]);
  }
  fclose(f_tab);

  sprintf(filename_tab, "%s/params.info", sa_index_dirname);
  f_tab = fopen(filename_tab, "w");
  fprintf(f_tab, "1. filename prefix\n");
  fprintf(f_tab, "2. k value\n");
  fprintf(f_tab, "3. skip\n");
  fprintf(f_tab, "4. A table length (= JA table length). Be carefull, A 32-bit items, for JA 8-bit items\n");
  fprintf(f_tab, "5. IA table length\n");
  fprintf(f_tab, "6. Number of suffixes\n");
  fprintf(f_tab, "7. Genome length\n");
  fprintf(f_tab, "8. Number of sequencess\n");
  fprintf(f_tab, "9. One line per sequence: name, length, type, chrom, start, end, left and right flanks (the last five fields for ALT sequences)\n");
  fclose(f_tab);

  sprintf(filename_tab, "%s/index", sa_index_dirname);
  f_tab = fopen(filename_tab, "w");
  for (size_t i = 0; i < genome->num_seqs; i++) {
    fprintf(f_tab, ">%s %i %lu\n", 
	   (genome->seq_names ? genome->seq_names[i] : "no-name"), 
	    0,
	    genome->seq_lengths[i] - 1);
  }
  fclose(f_tab);
}

//--------------------------------------------------------------------------------------
// load a SA index in memory
//--------------------------------------------------------------------------------------

sa_index3_t *sa_index3_new(char *sa_index_dirname) {

  FILE *f_tab;
  char line[1024], filename_tab[strlen(sa_index_dirname) + 1024];
  char *prefix;
  size_t k_value, A_items, IA_items, num_suffixes, genome_len, num_seqs, num_items;

  struct timeval stop, start;

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
    exit(EXIT_FAILURE);
  }

  char *res;

  // prefix
  res = fgets(line, 1024, f_tab);
  line[strlen(line) - 1] = 0;
  prefix = strdup(line);

  // k_value
  res = fgets(line, 1024, f_tab);
  sscanf(line, "%lu\n", &k_value);

  // skip
  res = fgets(line, 1024, f_tab);
  sscanf(line, "%lu\n", &A_items);

  // A_items
  res = fgets(line, 1024, f_tab);
  sscanf(line, "%lu\n", &A_items);

  // IA_items
  res = fgets(line, 1024, f_tab);
  sscanf(line, "%lu\n", &IA_items);

  // num_suffixes
  res = fgets(line, 1024, f_tab);
  sscanf(line, "%lu\n", &num_suffixes);

  // genome_length
  res = fgets(line, 1024, f_tab);
  sscanf(line, "%lu\n", &genome_len);

  // num_seqs
  res = fgets(line, 1024, f_tab);
  sscanf(line, "%lu\n", &num_seqs);

  size_t *seq_lengths = (size_t *) malloc(num_seqs * sizeof(size_t));
  int *seq_flags = (int *) malloc(num_seqs * sizeof(int));
  size_t *seq_chroms = (size_t *) malloc(num_seqs * sizeof(size_t));
  size_t *seq_starts = (size_t *) malloc(num_seqs * sizeof(size_t));
  size_t *seq_ends = (size_t *) malloc(num_seqs * sizeof(size_t));
  size_t *left_flanks = (size_t *) calloc(num_seqs, sizeof(size_t));
  size_t *right_flanks = (size_t *) calloc(num_seqs, sizeof(size_t));

  char **seq_names = (char **) malloc(num_seqs * sizeof(char *));
  char seq_name[1024];
  char seq_flag;
  size_t seq_len, seq_chrom, seq_start, seq_end, left_flank, right_flank;
		  
  for (size_t i = 0; i < num_seqs; i++) {
    res = fgets(line, 1024, f_tab);
    sscanf(line, "%s\t%lu\t%i\t%lu\t%lu\t%lu\t%lu\t%lu\n", 
	   seq_name, &seq_len, &seq_flag, &seq_chrom, &seq_start, &seq_end, &left_flank, &right_flank);
    //printf("chrom_name: %s, chrom_len: %lu\n", chrom_name, chrom_len);
    seq_names[i] = strdup(seq_name);
    seq_lengths[i] = seq_len;
    seq_flags[i] = seq_flag;
    seq_chroms[i] = seq_chrom;
    seq_starts[i] = seq_start;
    seq_ends[i] = seq_end;
    left_flanks[i] = left_flank;
    right_flanks[i] = right_flank;
  }

  fclose(f_tab);

  char *S;
  unsigned short int *CHROM;
  unsigned char *JA;
  sa_genome3_t *genome;
  uint *SA, *PRE, *A, *IA;

  #pragma omp parallel sections num_threads(2)
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
	exit(EXIT_FAILURE);
      }
      //  printf("genome: filename %s, length = %lu\n", filename_tab, genome_len);
      S = (char *) malloc(genome_len);
  
      //printf("\nreading S from file %s...\n", filename_tab);
      gettimeofday(&start, NULL);
      num_items = fread(S, sizeof(char), genome_len, f_tab);
      if (num_items != genome_len) {
	printf("Error: (%s) mismatch num_items = %i vs length = %i\n", 
	       filename_tab, num_items, genome_len);
	exit(EXIT_FAILURE);
      }
      gettimeofday(&stop, NULL);
      //printf("end of reading S (%lu len) from file %s in %0.2f s\n", 
      //	     genome_len, filename_tab,
      //	     (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);      
      fclose(f_tab);
      
      genome = sa_genome3_new(genome_len, num_seqs, seq_lengths, seq_flags, seq_chroms, 
			      seq_starts, seq_ends, seq_names, S);
      
      for (size_t i = 0; i < genome->length; i++) {
	if (genome->S[i] == 'N' || genome->S[i] == 'n') {
	  genome->S[i] = 'A';
	}
      }

      // display genome
      //sa_genome3_display(genome);

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
	  exit(EXIT_FAILURE);
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
	  printf("Error: (%s) mismatch read num_items = %i (it must be %i)\n", 
		 filename_tab, num_items, IA_items);
	  exit(EXIT_FAILURE);
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
	exit(EXIT_FAILURE);
      }
      //      printf("SA: filename %s, num_suffixes = %lu\n", filename_tab, num_suffixes);
      SA = (uint *) malloc(num_suffixes * sizeof(uint));
  
      //      printf("\nreading SA table from file %s...\n", filename_tab);
      gettimeofday(&start, NULL);
      num_items = fread(SA, sizeof(uint), num_suffixes, f_tab);
      if (num_items != num_suffixes) {
	printf("Error: (%s) mismatch num_items = %i vs num_suffixes = %i\n", 
	       filename_tab, num_items, num_suffixes);
	exit(EXIT_FAILURE);
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
	exit(EXIT_FAILURE);
      }
      //      printf("CHROM: filename %s, num_suffixes = %lu\n", filename_tab, num_suffixes);
      CHROM = (unsigned short int *) malloc(num_suffixes * sizeof(unsigned short int));
      
      //      printf("\nreading CHROM table from file %s...\n", filename_tab);
      gettimeofday(&start, NULL);
      num_items = fread(CHROM, sizeof(unsigned short int), num_suffixes, f_tab);
      if (num_items != num_suffixes) {
	printf("Error: (%s) mismatch num_items = %i vs num_suffixes = %i\n", 
	       filename_tab, num_items, num_suffixes);
	exit(EXIT_FAILURE);
      }
      gettimeofday(&stop, NULL);
      //      printf("end of reading CHROM table (%lu num_suffixes) from file %s in %0.2f s\n", 
      //	     num_suffixes, filename_tab,
      //	     (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);      
      fclose(f_tab);

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
	  exit(EXIT_FAILURE);
	}
	gettimeofday(&stop, NULL);
	//	printf("end of reading JA table (%lu num_items) from file %s in %0.2f s\n", 
	//	       num_items, filename_tab,
	//	       (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
	fclose(f_tab);
      }
    }
  }
  free(prefix);

  // creating the sa_index_t structure
  sa_index3_t *p = (sa_index3_t *) malloc(sizeof(sa_index3_t));

  p->num_suffixes = num_suffixes;
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

  return p;
}

//--------------------------------------------------------------------------------------

void sa_index3_free(sa_index3_t *p) {
  if (p) {
    
    if (p->SA) free(p->SA);
    if (p->CHROM) free(p->CHROM);
    if (p->PRE) free(p->PRE);
    if (p->A) free(p->A);
    if (p->IA) free(p->IA);
    if (p->JA) free(p->JA);
    if (p->genome) sa_genome3_free(p->genome);

    free(p);
  }
}

//--------------------------------------------------------------------------------------

void get_seq_names(char *fasta_filename, array_list_t *seqs_names) {
  const int MAX_CHROM_NAME_LENGHT = 1024;
  uint reading_name = 0, chrom_name_count = 0;
  char chrom_name[MAX_CHROM_NAME_LENGHT];

  // open fasta file
  FILE *f = fopen(fasta_filename, "r");
  if (f == NULL) {
    printf("Error reading FASTA filename %s\n", fasta_filename);
    exit(EXIT_FAILURE);
  }

  char c;
  size_t num_chroms = 0, skip = 0, process = 0;

  while ((c = getc(f)) != EOF) {
    process++;
    //    if (process % PROGRESS == 0) printf("reading %0.2f %c...\n", 100.0f * process / file_length, '%'); 
   //printf("len %lu: c = %c\n", len, c);
    if (c == '>') {
      if (skip == 0) {
	num_chroms++;
	reading_name = 1;
      }
      skip = 1;
    } else if (c == '\n') {
      if (reading_name) {
	chrom_name[chrom_name_count] = 0;
	array_list_insert(strdup(chrom_name), seqs_names);
	//	printf("************ chrom. name = %s\n", chrom_name);
	chrom_name_count = 0;
	reading_name = 0;
      }
      if (skip == 1) printf("%c", c);
      skip = 0;
    } else if (c == ' ' || c == '\t') {
      if (reading_name) {
	chrom_name[chrom_name_count] = 0;
	array_list_insert(strdup(chrom_name), seqs_names);
	//	printf("******** chrom. name = %s\n", chrom_name);
	chrom_name_count = 0;
	reading_name = 0;
      }
    } else {
      if (skip) {
	if (reading_name) {
	  chrom_name[chrom_name_count++] = c;
	  if (chrom_name_count >= MAX_CHROM_NAME_LENGHT) {
	    chrom_name[MAX_CHROM_NAME_LENGHT - 1] = 0;
	    printf("Chromosome name (%s) exceeds max. length (%u)", chrom_name, MAX_CHROM_NAME_LENGHT);
	    exit(EXIT_FAILURE);
	  }
	}
      }
    }
  }
  fclose(f);
}

//--------------------------------------------------------------------------------------

void sa_index3_set_decoy_names(array_list_t *decoy_names, char *sa_index_dirname) {

  FILE *f_tab;
  char line[1024], filename_tab[strlen(sa_index_dirname) + 1024];
  char *prefix;
  size_t pre_length, IA_items, A_items, num_suffixes, genome_len, num_seqs;
  int k_value;

  // read meta info
  sprintf(filename_tab, "%s/params.txt", sa_index_dirname);
  //printf("reading %s\n", filename_tab);

  f_tab = fopen(filename_tab, "r");
  if (!f_tab) {
    fprintf(stderr, "Error opening file %s!\n", filename_tab);
    exit(EXIT_FAILURE);
  }

  char *res;

  // prefix
  res = fgets(line, 1024, f_tab);
  line[strlen(line) - 1] = 0;
  prefix = strdup(line);

  // k_value
  res = fgets(line, 1024, f_tab);
  sscanf(line, "%lu\n", &k_value);

  // skip
  res = fgets(line, 1024, f_tab);
  sscanf(line, "%lu\n", &A_items);

  // A_items
  res = fgets(line, 1024, f_tab);
  sscanf(line, "%lu\n", &A_items);

  // IA_items
  res = fgets(line, 1024, f_tab);
  sscanf(line, "%lu\n", &IA_items);

  // num_suffixes
  res = fgets(line, 1024, f_tab);
  sscanf(line, "%lu\n", &num_suffixes);

  // genome_length
  res = fgets(line, 1024, f_tab);
  sscanf(line, "%lu\n", &genome_len);

  // num_seqs
  res = fgets(line, 1024, f_tab);
  sscanf(line, "%lu\n", &num_seqs);

  size_t *seq_lengths = (size_t *) malloc(num_seqs * sizeof(size_t));
  int *seq_flags = (int *) malloc(num_seqs * sizeof(int));
  size_t *seq_chroms = (size_t *) malloc(num_seqs * sizeof(size_t));
  size_t *seq_starts = (size_t *) malloc(num_seqs * sizeof(size_t));
  size_t *seq_ends = (size_t *) malloc(num_seqs * sizeof(size_t));
  size_t *left_flanks = (size_t *) calloc(num_seqs, sizeof(size_t));
  size_t *right_flanks = (size_t *) calloc(num_seqs, sizeof(size_t));

  char **seq_names = (char **) malloc(num_seqs * sizeof(char *));
  char seq_name[1024];
  int seq_flag;
  size_t seq_len, seq_chrom, seq_start, seq_end, left_flank, right_flank;
		  
  int num_decoys = array_list_size(decoy_names);
  for (size_t i = 0; i < num_seqs; i++) {
    res = fgets(line, 1024, f_tab);
    sscanf(line, "%s\t%lu\t%i\t%lu\t%lu\t%lu\t%lu\t%lu\n", 
	   seq_name, &seq_len, &seq_flag, &seq_chrom, &seq_start, &seq_end, &left_flank, &right_flank);
    //printf("chrom_name: %s, chrom_len: %lu\n", chrom_name, chrom_len);
    seq_names[i] = strdup(seq_name);
    seq_lengths[i] = seq_len;
    seq_flags[i] = seq_flag;
    seq_chroms[i] = seq_chrom;
    seq_starts[i] = seq_start;
    seq_ends[i] = seq_end;
    left_flanks[i] = left_flank;
    right_flanks[i] = right_flank;
    for (int j = 0; j < num_decoys; j++) {
      if (strcmp(seq_name, array_list_get(j, decoy_names)) == 0) {
	seq_flags[i] = DECOY_FLAG;
	break;
      }
    }
  }

  fclose(f_tab);

  //  char filename_tab2[1024];
  //  sprintf(filename_tab2, "%s", filename
  strcat(filename_tab, "_2");

  // write meta info
  f_tab = fopen(filename_tab, "w");
  fprintf(f_tab, "%s\n", prefix);
  fprintf(f_tab, "%i\n", k_value);
  fprintf(f_tab, "%lu\n", pre_length);
  fprintf(f_tab, "%lu\n", A_items);
  fprintf(f_tab, "%lu\n", IA_items);
  fprintf(f_tab, "%lu\n", num_suffixes);
  fprintf(f_tab, "%lu\n", genome_len);
  fprintf(f_tab, "%lu\n", num_seqs);
  for (size_t i = 0; i < num_seqs; i++) {
    fprintf(f_tab, "%s\t%lu\t%i\t%lu\t%lu\t%lu\t%lu\t%lu\n", 
	    seq_names[i], seq_lengths[i], seq_flags[i], seq_chroms[i], 
	    seq_starts[i], seq_ends[i], left_flanks[i], right_flanks[i]);
  }
  fclose(f_tab);
}

//--------------------------------------------------------------------------------------
// merge genomes
//--------------------------------------------------------------------------------------

void sa_index3_set_decoy(char *decoy_genome, char *sa_index_dirname) {

  const int MAX_CHROM_NAME_LENGHT = 1024;
  uint reading_name = 0, chrom_name_count = 0;
  char chrom_name[MAX_CHROM_NAME_LENGHT];

  // open fasta file
  FILE *f = fopen(decoy_genome, "r");
  if (f == NULL) {
    printf("Error reading decoy filename %s\n", decoy_genome);
    exit(EXIT_FAILURE);
  }

  // read genome
  array_list_t *seqs_names = array_list_new(100, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  char c;
  size_t num_chroms = 0, skip = 0, process = 0;

  while ((c = getc(f)) != EOF) {
    process++;
    //    if (process % PROGRESS == 0) printf("reading %0.2f %c...\n", 100.0f * process / file_length, '%'); 
   //printf("len %lu: c = %c\n", len, c);
    if (c == '>') {
      if (skip == 0) {
	num_chroms++;
	reading_name = 1;
      }
      skip = 1;
    } else if (c == '\n') {
      if (reading_name) {
	chrom_name[chrom_name_count] = 0;
	array_list_insert(strdup(chrom_name), seqs_names);
	//	printf("************ chrom. name = %s\n", chrom_name);
	chrom_name_count = 0;
	reading_name = 0;
      }
      if (skip == 1) printf("%c", c);
      skip = 0;
    } else if (c == ' ' || c == '\t') {
      if (reading_name) {
	chrom_name[chrom_name_count] = 0;
	array_list_insert(strdup(chrom_name), seqs_names);
	//	printf("******** chrom. name = %s\n", chrom_name);
	chrom_name_count = 0;
	reading_name = 0;
      }
    } else {
      if (skip) {
	if (reading_name) {
	  chrom_name[chrom_name_count++] = c;
	  if (chrom_name_count >= MAX_CHROM_NAME_LENGHT) {
	    chrom_name[MAX_CHROM_NAME_LENGHT - 1] = 0;
	    printf("Chromosome name (%s) exceeds max. length (%u)", chrom_name, MAX_CHROM_NAME_LENGHT);
	    exit(EXIT_FAILURE);
	  }
	}
      }
    }
  }
  fclose(f);

  sa_index3_set_decoy_names(seqs_names, sa_index_dirname);
}


//--------------------------------------------------------------------------------------
// merge genomes
//--------------------------------------------------------------------------------------

void merge_genomes(char *genome1, char *genome2, char *out_genome) {
  FILE *f1 = fopen(genome1, "r");
  FILE *f2 = fopen(genome2, "r");
 
  if (f1 == NULL || f2 == NULL) { 
    printf("Error: could not open genomes files: %s, %s\n", genome1, genome2);
    exit(EXIT_FAILURE);
  }
 
  FILE *f3 = fopen(out_genome,"w");
 
  if (f3 == NULL ) {
    printf("Error: could not write file: %s\n", out_genome);
    exit(EXIT_FAILURE);
  }
 
  char c;
  while ((c = fgetc(f1)) != EOF) fputc(c, f3);
  while ((c = fgetc(f2)) != EOF) fputc(c, f3);
 
  fclose(f1);
  fclose(f2);
  fclose(f3); 
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
