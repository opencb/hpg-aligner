#include "sa_index3.h"

#define PROGRESS 1000000

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

sa_genome3_t *read_genome3(char *filename) {

  const int MAX_CHROM_NAME_LENGHT = 1024;
  uint reading_name, chrom_name_count = 0;
  char chrom_name[MAX_CHROM_NAME_LENGHT];

  // open fasta file
  FILE *f = fopen(filename, "r");
  if (f == NULL) {
    printf("Error reading filename %s\n", filename);
    exit(-1);
  }

  // get file length
  fseek(f, 0, SEEK_END);
  long int file_length = ftell(f);
  printf("file %s length = %lu\n", filename, file_length);
  fseek(f, 0, SEEK_SET);	

  // read genome
  size_t num_A = 0, num_C = 0, num_G = 0, num_N = 0, num_T = 0;
  size_t chrom_length = 0;
  size_t num_chroms = 0, num_allocated_chroms = 100;
  size_t *chrom_lengths = (size_t *) calloc(num_allocated_chroms, sizeof(size_t));
  char **chrom_names = (char **) calloc(num_allocated_chroms, sizeof(char *));
  char *S = (char *) calloc(file_length + 1, sizeof(char));

  uint skip = 0, l = 0, process = 0;
  char *c;
  while ((c = getc(f)) != EOF) {
    process++;
    //    if (process % PROGRESS == 0) printf("reading %0.2f %c...\n", 100.0f * process / file_length, '%'); 
   //printf("len %lu: c = %c\n", len, c);
    if (c == '>') {
      printf("%c", c);
      if (skip == 0) {
	num_chroms++;
	reading_name = 1;
	// to do: get and save chromosome name
	if (num_chroms > 1) {
	  chrom_lengths[num_chroms - 2] = chrom_length;
	  printf("setting chrom %lu: length = %lu (but num_chroms = %lu)\n", num_chroms - 2, chrom_length, num_chroms);
	  chrom_length = 0;
	  if (num_chroms >= num_allocated_chroms) {
	    num_allocated_chroms += 50;
	    chrom_lengths = (size_t *) realloc(chrom_lengths, num_allocated_chroms * sizeof(size_t));
	    chrom_names = (char **) realloc(chrom_lengths, num_allocated_chroms * sizeof(char *));
	  }
	}
      }
      skip = 1;
    } else if (c == '\n') {
      if (reading_name) {
	chrom_name[chrom_name_count] = 0;
	chrom_names[num_chroms - 1] = strdup(chrom_name);
	//	printf("************ chrom. name = %s (num. chrom = %u)\n", chrom_name, num_chroms);
	chrom_name_count = 0;
	reading_name = 0;
      }
      if (skip == 1) printf("%c", c);
      skip = 0;
    } else if (c == ' ' || c == '\t') {
      if (reading_name) {
	chrom_name[chrom_name_count] = 0;
	chrom_names[num_chroms - 1] = strdup(chrom_name);
	//	printf("******** chrom. name = %s (nu. chrom = %u)\n", chrom_name, num_chroms);
	chrom_name_count = 0;
	reading_name = 0;
      }
    } else {
      if (skip == 0) {
	if      (c == 'A' || c == 'a') { S[l++] = 'A'; num_A++; chrom_length++; }
	else if (c == 'C' || c == 'c') { S[l++] = 'C'; num_C++; chrom_length++; }
	else if (c == 'G' || c == 'g') { S[l++] = 'G'; num_G++; chrom_length++; }
	else if (c == 'T' || c == 't') { S[l++] = 'T'; num_T++; chrom_length++; }
	else if (c == 'N' || c == 'n') { S[l++] = 'N'; num_N++; chrom_length++; }
	else {
	  printf("Unknown character %c at %lu position\n", c, process);
	}
      } else {
	if (reading_name) {
	  chrom_name[chrom_name_count++] = c;
	  if (chrom_name_count >= MAX_CHROM_NAME_LENGHT) {
	    chrom_name[MAX_CHROM_NAME_LENGHT - 1] = 0;
	    printf("Chromosome name (%s) exceeds max. length (%u)", chrom_name, MAX_CHROM_NAME_LENGHT);
	    exit(-1);
	  }
	}
	printf("%c", c);
      }
    }
  }
  chrom_lengths[num_chroms - 1] = chrom_length;
  S[l++] = '$';

  fclose(f);
  printf("reading %0.2f %c...\n", 100.0f, '%');

  // create sa_genome3_t and return it
  sa_genome3_t *genome = sa_genome3_new(l, num_chroms, chrom_lengths, chrom_names, S);
  sa_genome3_set_nt_counters(num_A, num_C, num_G, num_N, num_T, genome);
  return genome;
}

//--------------------------------------------------------------------------------------

char *global_S;

//--------------------------------------------------------------------------------------

typedef struct suffix_tmp {
  uint value;
  char chrom;
} suffix_tmp_t;


int suffix_tmp_cmp3(void const *a, void const *b) { 
  suffix_tmp_t *item_a = (suffix_tmp_t *) a;
  suffix_tmp_t *item_b = (suffix_tmp_t *) b;

  return (strncmp(&global_S[item_a->value], &global_S[item_b->value], 1000));
}

//--------------------------------------------------------------------------------------

/*
void compute_LCP(uint num_suffixes, char *S, uint *SA, uint *LCP) {
  uint h, i, j, k;
  
  // Kasai et al linear-time construction
  // computes the inverse permutation of the suffix array
  uint *ISA = (uint *) malloc(num_suffixes * sizeof(uint));
  for (i = 0; i < num_suffixes; i++) {
    ISA[SA[i]] = i;
  }
  
  h=0;  // index to support result that lcp[rank[i]]>=lcp[rank[i-1]]-1
  for (i = 0; i < num_suffixes; i++) {
    k = ISA[i];
    if (k == 0) {
      LCP[k]= 0; //(-1);
    } else {
      j = SA[k-1];
      // attempt to extend lcp
      while (i + h < num_suffixes && 
	     j + h < num_suffixes && 
	     S[i+h] == S[j+h]) {
	h++;
      }
      LCP[k]=h;
    }
    if ( h > 0) {
      h--;  // decrease according to result
    }
  }

  // free memory
  free(ISA);
}
*/
//--------------------------------------------------------------------------------------
/*
void compute_PRE(uint num_suffixes, char *S, uint *SA, uint *PRE) {

  PREFIX_TABLE_NT_VALUE['A'] = 0;
  PREFIX_TABLE_NT_VALUE['N'] = 0;
  PREFIX_TABLE_NT_VALUE['C'] = 1;
  PREFIX_TABLE_NT_VALUE['G'] = 2;
  PREFIX_TABLE_NT_VALUE['T'] = 3;

  uint count = 0;

  uint value;
  for (uint i = 1; i < num_suffixes; i++) {
    value = compute_prefix_value(&S[SA[i]], PREFIX_TABLE_K_VALUE);
    if (PRE[value] == 0) {
      PRE[value] = i; //sa[i];
*/
      /*
      if (++count < 20) {
	printf("i = %i\tvalue = %lu\tPRE[value] = %lu\tSA[PRE[value]] = %lu\t", i, value, PRE[value], SA[PRE[value]]);
	char *p = &S[SA[PRE[value]]];
	for (size_t j = 0; j < PREFIX_TABLE_K_VALUE; j++) {
	  printf("%c", p[j]);
	}
	printf("\n");
      }
      */
/*
    }
  }
}
*/

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
  uint len = 0;
  printf("\nreading file genome %s...\n", genome_filename);
  gettimeofday(&start, NULL);
  sa_genome3_t *genome = read_genome3(genome_filename);
  gettimeofday(&stop, NULL);
  printf("end of reading file (%lu items) in %0.2f s\n", 
	 genome_filename,
	 (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);

  sa_genome3_display(genome);

  // write S to file
  sprintf(filename_tab, "%s/%s.S", sa_index_dirname, prefix);
  f_tab = fopen(filename_tab, "wb");
  if (f_tab == NULL) {
    printf("Error: could not open %s to write\n", filename_tab);
    exit(-1);
  }
  fwrite(genome->S, sizeof(char), genome->length, f_tab);
  fclose(f_tab);

  global_S = genome->S;

  //-----------------------------------------
  // compute SA table
  //-----------------------------------------
  uint num_suffixes;
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
  for (uint i = 0; i < genome->num_chroms; i++) {
    for (uint j = 0; j < genome->chrom_lengths[i]; j++) {
      nt = genome->S[c];
      if (nt == 'A' || nt == 'a') {
	tmp[0][tmp_A].value = c;
	tmp[0][tmp_A].chrom = i;
	tmp_A++;
      } else if (nt == 'C' || nt == 'c') {
	tmp[1][tmp_C].value = c;
	tmp[1][tmp_C].chrom = i;
	tmp_C++;
      } else if (nt == 'G' || nt == 'g') {
	tmp[2][tmp_G].value = c;
	tmp[2][tmp_G].chrom = i;
	tmp_G++;
      } else if (nt == 'T' || nt == 't') {
	tmp[3][tmp_T].value = c;
	tmp[3][tmp_T].chrom = i;
	tmp_T++;
      }
      c++;
    }
  }

  assert(tmp_A == genome->num_A);
  assert(tmp_C == genome->num_C);
  assert(tmp_G == genome->num_G);
  assert(tmp_T == genome->num_T);

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
      fwrite(&(tmp[i][j].chrom), sizeof(char), 1, f_tab);
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
    exit(-1);
  } 
  fwrite(PRE, sizeof(uint), pre_length, f_tab);
  fclose(f_tab);
    
  //-----------------------------------------
  // save parameters
  //-----------------------------------------
  sprintf(filename_tab, "%s/params.txt", sa_index_dirname);
  f_tab = fopen(filename_tab, "w");
  fprintf(f_tab, "%s\n", prefix);
  fprintf(f_tab, "%lu\n", k_value);
  fprintf(f_tab, "%lu\n", pre_length);
  fprintf(f_tab, "%lu\n", num_suffixes);
  fprintf(f_tab, "%lu\n", genome->length);
  fprintf(f_tab, "%lu\n", genome->num_chroms);
  for (size_t i = 0; i < genome->num_chroms; i++) {
    fprintf(f_tab, "%s\t%lu\n", 
	   (genome->chrom_names ? genome->chrom_names[i] : "no-name"), 
	    genome->chrom_lengths[i]);
  }
  fclose(f_tab);
}

//--------------------------------------------------------------------------------------

void sa_index3_build_k18(char *genome_filename, uint k_value, char *sa_index_dirname) {

  k_value = 18;
  printf("\n***************** K value = 18 ***************************\n");

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
    exit(-1);
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
  sa_genome3_t *genome = read_genome3(genome_filename);
  gettimeofday(&stop, NULL);
  printf("end of reading file (%lu items) in %0.2f s\n", 
	 genome_filename,
	 (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);
  
  sa_genome3_display(genome);

  // write S to file
  f_tab = fopen(filename_tab, "wb");
  if (f_tab == NULL) {
    printf("Error: could not open %s to write\n", filename_tab);
    exit(-1);
  }
  fwrite(genome->S, sizeof(char), genome->length, f_tab);
  fclose(f_tab);

  //-----------------------------------------
  // compute SA table
  //-----------------------------------------

  char nt;
  uint num_suffixes = genome->num_A + genome->num_C + genome->num_G + genome->num_T;

  sprintf(filename_tab, "%s/%s.SA", sa_index_dirname, prefix);
  f_tab = fopen(filename_tab, "rb");
  if (!f_tab) {
    printf("\ncomputing SA and CHROM table...\n");
    gettimeofday(&start, NULL);
    suffix_tmp_t *tmp[4];
    uint tmp_A = 0, tmp_C = 0, tmp_G = 0, tmp_T = 0, nts[4];
    
    nts[0] = genome->num_A;
    nts[1] = genome->num_C;
    nts[2] = genome->num_G;
    nts[3] = genome->num_T;
    
    for (int i = 0; i < 4; i++) {
      tmp[i] = (suffix_tmp_t *) malloc(nts[i] * sizeof(suffix_tmp_t));
    }
    
    uint c = 0;
    for (uint i = 0; i < genome->num_chroms; i++) {
      for (uint j = 0; j < genome->chrom_lengths[i]; j++) {
	nt = genome->S[c];
	if (nt == 'A' || nt == 'a') {
	  tmp[0][tmp_A].value = c;
	  tmp[0][tmp_A].chrom = i;
	  tmp_A++;
	} else if (nt == 'C' || nt == 'c') {
	  tmp[1][tmp_C].value = c;
	  tmp[1][tmp_C].chrom = i;
	  tmp_C++;
	} else if (nt == 'G' || nt == 'g') {
	  tmp[2][tmp_G].value = c;
	  tmp[2][tmp_G].chrom = i;
	  tmp_G++;
	} else if (nt == 'T' || nt == 't') {
	  tmp[3][tmp_T].value = c;
	  tmp[3][tmp_T].chrom = i;
	  tmp_T++;
	}
	c++;
      }
    }
    
    assert(tmp_A == genome->num_A);
    assert(tmp_C == genome->num_C);
    assert(tmp_G == genome->num_G);
    assert(tmp_T == genome->num_T);
    
    for (size_t i = 0; i < genome->length; i++) {
      if (genome->S[i] == 'N' || genome->S[i] == 'n') {
	genome->S[i] = 'A';
      }
    }
    global_S = genome->S;
    /*
      printf("before sorting...\n");
      for (int i = 0; i < nts[0]; i++) {
p      display_prefix(&genome->S[tmp[0][i].value], k_value);
      printf("\ttmp[0][%i] = %u -> prefix.value = %lu\n", 
      i, tmp[0][i].value, compute_prefix_value(&genome->S[tmp[0][i].value], k_value));
      }
    */
    #pragma omp parallel for num_threads(4)
    for (size_t i = 0; i < 4; i++) {
      qsort(tmp[i], nts[i], sizeof(suffix_tmp_t), suffix_tmp_cmp3);
    }
    /*
      printf("after sorting...\n");
      for (int i = 0; i < nts[0]; i++) {
      display_prefix(&genome->S[tmp[0][i].value], k_value);
      printf("\ttmp[0][%i] = %u -> prefix.value = %lu\n", 
      i, tmp[0][i].value, compute_prefix_value(&genome->S[tmp[0][i].value], k_value));
      }
    */
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
	fwrite(&(tmp[i][j].chrom), sizeof(char), 1, f_tab);
      }
    }
    fclose(f_tab);

    // free memory
    for (int i = 0; i < 4; i++) {
      free(tmp[i]);
    }
  } else {
    printf("updating S genome (N -> A)...\n");
    for (size_t i = 0; i < genome->length; i++) {
      if (genome->S[i] == 'N' || genome->S[i] == 'n') {
	genome->S[i] = 'A';
      }
    }
    printf("end of updating S genome. Done!\n");
  }

  sprintf(filename_tab, "%s/%s.SA", sa_index_dirname, prefix);
  f_tab = fopen(filename_tab, "rb");
  if (f_tab == NULL) {
    printf("Error: could not open %s to write\n", filename_tab);
    exit(-1);
  }
  printf("SA: filename %s, num_suffixes = %lu\n", filename_tab, num_suffixes);
  uint *SA = (uint *) malloc(num_suffixes * sizeof(uint));
  
  printf("\nreading SA table from file %s...\n", filename_tab);
  gettimeofday(&start, NULL);
  uint num_items = fread(SA, sizeof(uint), num_suffixes, f_tab);
  if (num_items != num_suffixes) {
    printf("Error: (%s) mismatch num_items = %lu vs num_suffixes = %lu\n", 
	   filename_tab, num_items, num_suffixes);
    exit(-1);
  }
  gettimeofday(&stop, NULL);
  printf("end of reading SA table (%lu num_suffixes) from file %s in %0.2f s\n", 
  	 num_suffixes, filename_tab,
	 (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);      
  fclose(f_tab);


  /*
  {
    char *seq;
    size_t num_prefixes, prev_value, value, pos;
    for (size_t k = 23; k <= 26; k++) {
      prev_value = 999999;
      num_prefixes = 0;
      for (size_t i = 0; i < num_suffixes; i++) {
	pos = SA[i];
	if (pos + k < genome->length) {
	  seq = &genome->S[pos];
	  value = compute_prefix_value(seq, k);
	  if (value != prev_value) {
	    prev_value = value;
	    num_prefixes++;
	  }
	}
      }
      printf("k = %lu -> num. prefixes = %lu\n", k, num_prefixes);
    }
    exit(-1);
  }
  */

  //-----------------------------------------
  // compute Compressed Row Storage tables
  //-----------------------------------------
  printf("\ncomputing Compressed Row Storage tables...\n");
  gettimeofday(&start, NULL);

  size_t num_prefixes = 0;

  // A vector
  sprintf(filename_tab, "%s/%s.A", sa_index_dirname, prefix);
  FILE *f_A = fopen(filename_tab, "wb");
  uint A_counter = 0;

  // IA vector
  sprintf(filename_tab, "%s/%s.IA", sa_index_dirname, prefix);
  FILE *f_IA = fopen(filename_tab, "wb");
  uint IA_counter = 0, first_in_row;

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
  printf("filling matrix (%lu)...\n", matrix);
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
      printf("end of filling sub-matrix (%lu): %lu items. Done !!\n", matrix, matrix_items);
      printf("\t -----> new matrix due to the value SA[%lu] = %lu -> ", i, value);
      display_prefix(&genome->S[SA[i]], k_value);
      printf("\n");
      
      // read matrix and creating A, IA and JA vectors
      printf("reading matrix %lu and updating CRS vectors...\n", matrix);
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

      printf("end of reading matrix %lu and updating CRS vectors. Done !!\n", matrix);
      //	exit(-1);

      // initialize matrix to re-fill it
      matrix = value / M_items;
      matrix_items = 0;
      memset(M, 255, M_bytes);
      
      printf("init and filling matrix (%lu)...\n", matrix);
      M[value % M_items] = i;
      matrix_items++;
    }
  }

  // read matrix and creating A, IA and JA vectors
  printf("reading matrix %lu and updating CRS vectors...\n", matrix);
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

  printf("end of reading matrix %lu and updating CRS vectors. Done !!\n", matrix);

  printf("A length = %u, IA length = %u (num. prefixes = %lu)\n", A_counter, IA_counter, num_prefixes);
  fclose(f_A);
  fclose(f_IA);
  fclose(f_JA);

  uint pre_length;// = 1LLU << (2 * k_value);
  uint *PRE;// = (uint *) calloc(pre_length, sizeof(uint));
/*

  gettimeofday(&stop, NULL);
  printf("end of computing PRE table in %0.2f s\n", 
	 (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  

  // write PRE to file for the next time
  sprintf(filename_tab, "%s/%s.PRE", sa_index_dirname, prefix);
  f_tab = fopen(filename_tab, "wb");
  if (f_tab == NULL) {
    printf("Error: could not open %s to write\n", filename_tab);
    exit(-1);
  } 
  fwrite(PRE, sizeof(uint), pre_length, f_tab);
  fclose(f_tab);
*/  
  //-----------------------------------------
  // save parameters
  //-----------------------------------------
  sprintf(filename_tab, "%s/params.txt", sa_index_dirname);
  f_tab = fopen(filename_tab, "w");
  fprintf(f_tab, "%s\n", prefix);
  fprintf(f_tab, "%lu\n", k_value);
  fprintf(f_tab, "%lu\n", pre_length);
  fprintf(f_tab, "%lu\n", A_counter);
  fprintf(f_tab, "%lu\n", IA_counter);
  fprintf(f_tab, "%lu\n", num_suffixes);
  fprintf(f_tab, "%lu\n", genome->length);
  fprintf(f_tab, "%lu\n", genome->num_chroms);
  for (size_t i = 0; i < genome->num_chroms; i++) {
    fprintf(f_tab, "%s\t%lu\n", 
	   (genome->chrom_names ? genome->chrom_names[i] : "no-name"), 
	    genome->chrom_lengths[i]);
  }
  fclose(f_tab);

  sprintf(filename_tab, "%s/params.info", sa_index_dirname);
  f_tab = fopen(filename_tab, "w");
  fprintf(f_tab, "1. filename prefix\n");
  fprintf(f_tab, "2. k value\n");
  fprintf(f_tab, "3. prefix table length: pow(2, k-value * 2)\n");
  fprintf(f_tab, "4. A table length (= JA table length). Be carefull, A 32-bit items, for JA 8-bit items\n");
  fprintf(f_tab, "5. IA table length\n");
  fprintf(f_tab, "6. Number of suffixes\n");
  fprintf(f_tab, "7. Genome length\n");
  fprintf(f_tab, "8. Number of chromosomes\n");
  fprintf(f_tab, "9. One line per chromsomome: name and length\n");
  fclose(f_tab);
}

//--------------------------------------------------------------------------------------
// load a SA index in memory
//--------------------------------------------------------------------------------------

sa_index3_t *sa_index3_new(char *sa_index_dirname) {

  FILE *f_tab;
  char line[1024], filename_tab[strlen(sa_index_dirname) + 1024];
  char *prefix;
  uint k_value, pre_length, A_items, IA_items, num_suffixes, genome_len, num_chroms, num_items;

  struct timeval stop, start;

  PREFIX_TABLE_NT_VALUE['A'] = 0;
  PREFIX_TABLE_NT_VALUE['N'] = 0;
  PREFIX_TABLE_NT_VALUE['C'] = 1;
  PREFIX_TABLE_NT_VALUE['G'] = 2;
  PREFIX_TABLE_NT_VALUE['T'] = 3;

  sprintf(filename_tab, "%s/params.txt", sa_index_dirname, prefix);

  printf("reading %s\n", filename_tab);

  f_tab = fopen(filename_tab, "r");
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

  #pragma omp parallel sections 
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
  
      printf("\nreading S from file %s...\n", filename_tab);
      gettimeofday(&start, NULL);
      num_items = fread(S, sizeof(char), genome_len, f_tab);
      if (num_items != genome_len) {
	printf("Error: (%s) mismatch num_items = %lu vs length = %lu\n", 
	       filename_tab, num_items, genome_len);
	exit(-1);
      }
      gettimeofday(&stop, NULL);
      printf("end of reading S (%lu len) from file %s in %0.2f s\n", 
	     genome_len, filename_tab,
	     (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);      
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
	
	printf("\nreading A table (Compression Row Storage) from file %s...\n", filename_tab);
	gettimeofday(&start, NULL);
	if ((num_items = fread(A, sizeof(uint), A_items, f_tab)) != A_items) {
	  printf("Error: (%s) mismatch read num_items = %lu (it must be %lu)\n", 
		 filename_tab, num_items, A_items);
	  exit(-1);
	}
	gettimeofday(&stop, NULL);
	printf("end of reading A table (%lu num_items) from file %s in %0.2f s\n", 
	       num_items, filename_tab,
	       (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
	fclose(f_tab);
      }

      // Compressed Row Storage (IA table)
      IA = NULL;
      sprintf(filename_tab, "%s/%s.IA", sa_index_dirname, prefix);
      f_tab = fopen(filename_tab, "rb");
      if (f_tab) {
	IA = (uint *) malloc(IA_items * sizeof(uint));
	
	printf("\nreading IA table (Compression Row Storage) from file %s...\n", filename_tab);
	gettimeofday(&start, NULL);
	if ((num_items = fread(IA, sizeof(uint), IA_items, f_tab)) != IA_items) {
	  printf("Error: (%s) mismatch read num_items = %lu (it must be %lu)\n", 
		 filename_tab, num_items, IA_items);
	  exit(-1);
	}
	gettimeofday(&stop, NULL);
	printf("end of reading IA table (%lu num_items) from file %s in %0.2f s\n", 
	       num_items, filename_tab,
	       (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
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
      printf("SA: filename %s, num_suffixes = %lu\n", filename_tab, num_suffixes);
      SA = (uint *) malloc(num_suffixes * sizeof(uint));
  
      printf("\nreading SA table from file %s...\n", filename_tab);
      gettimeofday(&start, NULL);
      num_items = fread(SA, sizeof(uint), num_suffixes, f_tab);
      if (num_items != num_suffixes) {
	printf("Error: (%s) mismatch num_items = %lu vs num_suffixes = %lu\n", 
	       filename_tab, num_items, num_suffixes);
	exit(-1);
      }
      gettimeofday(&stop, NULL);
      printf("end of reading SA table (%lu num_suffixes) from file %s in %0.2f s\n", 
	     num_suffixes, filename_tab,
	     (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);      
      fclose(f_tab);

      // read CHROM table from file
      sprintf(filename_tab, "%s/%s.CHROM", sa_index_dirname, prefix);
      f_tab = fopen(filename_tab, "rb");
      if (f_tab == NULL) {
	printf("Error: could not open %s to write\n", filename_tab);
	exit(-1);
      }
      printf("CHROM: filename %s, num_suffixes = %lu\n", filename_tab, num_suffixes);
      CHROM = (char *) malloc(num_suffixes * sizeof(char));
      
      printf("\nreading CHROM table from file %s...\n", filename_tab);
      gettimeofday(&start, NULL);
      num_items = fread(CHROM, sizeof(char), num_suffixes, f_tab);
      if (num_items != num_suffixes) {
	printf("Error: (%s) mismatch num_items = %lu vs num_suffixes = %lu\n", 
	       filename_tab, num_items, num_suffixes);
	exit(-1);
      }
      gettimeofday(&stop, NULL);
      printf("end of reading CHROM table (%lu num_suffixes) from file %s in %0.2f s\n", 
	     num_suffixes, filename_tab,
	     (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);      
      fclose(f_tab);

      // read PRE table from file
      PRE = NULL;
      sprintf(filename_tab, "%s/%s.PRE", sa_index_dirname, prefix);
      f_tab = fopen(filename_tab, "rb");
      if (f_tab) {
	num_items = 1LLU << (2 * k_value);
	assert(num_items == pre_length);
	PRE = (uint *) malloc(num_items * sizeof(uint));
	
	printf("\nreading PRE table from file %s...\n", filename_tab);
	gettimeofday(&start, NULL);
	if (num_items != fread(PRE, sizeof(uint), num_items, f_tab)) {
	  printf("Error: (%s) mismatch num_items = %lu\n", 
		 filename_tab, num_items);
	  exit(-1);
	}
	gettimeofday(&stop, NULL);
	printf("end of reading PRE table (%lu num_items) from file %s in %0.2f s\n", 
	       num_items, filename_tab,
	       (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
	fclose(f_tab);
      }
   
      // Compressed Row Storage (JA table)
      JA = NULL;
      sprintf(filename_tab, "%s/%s.JA", sa_index_dirname, prefix);
      f_tab = fopen(filename_tab, "rb");
      if (f_tab) {
	JA = (unsigned char *) malloc(A_items * sizeof(unsigned char));
	
	printf("\nreading JA table (Compression Row Storage) from file %s...\n", filename_tab);
	gettimeofday(&start, NULL);
	if ((num_items = fread(JA, sizeof(unsigned char), A_items, f_tab)) != A_items) {
	  printf("Error: (%s) mismatch read num_items = %lu (it must be %lu)\n", 
		 filename_tab, num_items, A_items);
	  exit(-1);
	}
	gettimeofday(&stop, NULL);
	printf("end of reading JA table (%lu num_items) from file %s in %0.2f s\n", 
	       num_items, filename_tab,
	       (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
	fclose(f_tab);
      }
    }
  }

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
//--------------------------------------------------------------------------------------




  /*
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
  {
    // start filling matrix
    size_t value, val;
    size_t num_rows = 16 * 1024 * 1024; // 16 M
    size_t num_cols = 256;

    size_t M_items = num_rows * num_cols;
    size_t M_bytes = M_items * sizeof(uint);
    printf("allocating matrix M (%lu bytes)\n", M_bytes);
    uint *M = (uint *) malloc(M_bytes);
    if (M == NULL) {
      printf("Error allocating memory for M matrix\n");
      exit(-1);
    }
    printf("end of allocating matrix M\n");
    
    const uint max_uint = 4294967295;
    memset(M, 0, M_bytes);

    size_t matrix_id = 0, num_prefixes = 0;
    size_t offset = matrix_id * M_items;
    printf("num. suffixes = %lu\n", num_suffixes);
    printf("processing sub-matrix %lu (offset %lu)\n", matrix_id, offset);
    for (uint i = 0; i < num_suffixes; i++) {
      value = compute_prefix_value(&genome->S[SA[i]], k_value);
      if (value / M_items == matrix_id) {
	if (M[value % M_items] == 0) { //max_uint) {
	  M[value % M_items] = 1;
	  num_prefixes++;
	  //	  printf("\t");
	  //	  display_prefix(&genome->S[SA[i]], k_value);
	  //	  printf("\tprefix: (num, value, val) = (%lu, %lu, %lu)\n", num_prefixes, value, val);
	}
      } else {
	printf("\tend of processing sub-matrix %lu: num. prefixes = %lu\n", matrix_id, num_prefixes);
	memset(M, 0, M_bytes);
	matrix_id = value / M_items;
	M[value % M_items] = 1;
	num_prefixes++;
	//	printf("\t");
	//	display_prefix(&genome->S[SA[i]], k_value);
	//	printf("\tprefix: (num, value, val) = (%lu, %lu, %lu)\n", num_prefixes, value, val);
	printf("processing sub-matrix %lu \n", matrix_id);
      }
    }
    printf("\tend of processing sub-matrix %lu: num. prefixes = %lu\n", matrix_id, num_prefixes);
    printf("---> num. prefixes = %lu\n", num_prefixes);
    exit(-1);
  }
  //////////////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////////////
*/

  /*  
  // read LCP table from file
  sprintf(filename_tab, "%s/%s.LCP", sa_index_dirname, prefix);
  f_tab = fopen(filename_tab, "rb");
  if (f_tab == NULL) {
    printf("Error: could not open %s to write\n", filename_tab);
    exit(-1);
  }
  uint *LCP = (uint *) malloc(num_suffixes * sizeof(uint));

  //  printf("\nreading LCP table from file %s...\n", filename_tab);
  gettimeofday(&start, NULL);
  num_items = fread(LCP, sizeof(uint), num_suffixes, f_tab);
  if (num_items != num_suffixes) {
    printf("Error: (%s) mismatch num_items = %lu vs num_suffixes = %lu\n", 
	   filename_tab, num_items, num_suffixes);
    exit(-1);
  }
  gettimeofday(&stop, NULL);
  //  printf("end of reading LCP table (%lu num_suffixes) from file %s in %0.2f s\n", 
  //	 num_suffixes, filename_tab,
  //	 (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  fclose(f_tab);
  */


  // change N to A
  //  for (uint i = 0; i < len; i++) {
  //    if (S[i] == 'N') S[i] = 'A';
  //  }
  
  //    size_t prev_value, curr_value, index = 0;
  //    curr_value = compute_prefix_value(&s[sa[1]], seed_size);
  //    prefetch_table[index] = curr_value;
  //    prev_value = curr_value;
  /*
  for (size_t i = 0; i < num_items; i++) {
    printf("i = %i\tPRE[i] = %lu\tSA[PRE[i]] = %lu\t", i, PRE[i], SA[PRE[i]]);
    char *p = &S[SA[PRE[i]]];
    for (size_t j = 0; j < PREFIX_TABLE_K_VALUE; j++) {
      printf("%c", p[j]);
    }
    printf("\n");
    if (i > 20) break;
  }
  */
  //	exit(-1);
