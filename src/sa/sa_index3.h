#ifndef SA_INDEX3_H
#define SA_INDEX3_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <assert.h>

#include "sa_tools.h"

//--------------------------------------------------------------------------------------

#define max_uint 4294967295

//--------------------------------------------------------------------------------------

typedef struct sa_genome3 {
  uint length;
  size_t num_chroms;
  size_t num_A;
  size_t num_C;
  size_t num_G;
  size_t num_N;
  size_t num_T;
  size_t *chrom_lengths;
  size_t *chrom_offsets;
  char **chrom_names;
  char *S;
} sa_genome3_t;

inline sa_genome3_t *sa_genome3_new(uint length, size_t num_chroms, 
				    size_t *chrom_lengths,
				    char **chrom_names, char *S) {
  sa_genome3_t *p = (sa_genome3_t *) calloc(1, sizeof(sa_genome3_t));
  p->length = length;
  p->num_chroms = num_chroms;
  p->chrom_lengths = chrom_lengths;
  if (num_chroms && chrom_lengths) {
    p->chrom_offsets = (size_t *) calloc(num_chroms, sizeof(size_t));
    size_t offset = 0;
    for (size_t i = 0; i < num_chroms; i++) {
      p->chrom_offsets[i] = offset;
      offset += chrom_lengths[i];
    }
  } else {
    p->chrom_offsets = NULL;
  }
  p->chrom_names = chrom_names;
  p->S = S;
  return p;
}

//--------------------------------------------------------------------------------------

inline void sa_genome3_free(sa_genome3_t *p) {
  if (p) {
    if (p->chrom_lengths) free(p->chrom_lengths);
    if (p->chrom_offsets) free(p->chrom_offsets);
    if (p->chrom_names) {
      for (int i = 0; i < p->num_chroms; i++) {
	free(p->chrom_names[i]);
      }
      free(p->chrom_names);
    }
    if (p->S) free(p->S);
    free(p);
  }
}

//--------------------------------------------------------------------------------------

inline char *sa_genome_get_sequence(int chrom, size_t start, size_t end, sa_genome3_t *p) {
  size_t pos, len = end - start + 1;
  char *seq = (char *) malloc((len + 1) * sizeof(char));
  for (size_t i = 0, pos = start + p->chrom_offsets[chrom]; i < len; i++, pos++) {
    seq[i] = p->S[pos];
  }
  seq[len] = 0;
  return seq;
}

//--------------------------------------------------------------------------------------

inline void sa_genome3_set_nt_counters(size_t num_A, size_t num_C, size_t num_G, 
				      size_t num_N, size_t num_T, sa_genome3_t *p) {
  if (p) {
    p->num_A = num_A;
    p->num_C = num_C;
    p->num_G = num_G;
    p->num_N = num_N;
    p->num_T = num_T;
  }
}

//--------------------------------------------------------------------------------------

inline void sa_genome3_display(sa_genome3_t *p) {
  if (!p) return;

  printf("Genome length: %u\n", p->length);
  printf("Number of chromosomes: %u\n", p->num_chroms);
  for (size_t i = 0; i < p->num_chroms; i++) {
    printf("\tChrom %u: (name, length, offset) = (%s, %u, %u)\n", 
	   i, (p->chrom_names ? p->chrom_names[i] : "no-name"), 
	   p->chrom_lengths[i], p->chrom_offsets[i]);
  }
  printf("Nucleotide counters:\n");
  printf("\tNumber of A: %u\n", p->num_A);
  printf("\tNumber of C: %u\n", p->num_C);
  printf("\tNumber of G: %u\n", p->num_G);
  printf("\tNumber of T: %u\n", p->num_T);
  printf("\tNumber of N: %u\n", p->num_N);
}

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------
/*
typedef struct sa3 {
  uint size; // in bytes (= num_suffixes * 5)
  uint num_suffixes;
  char *data; // 4 bytes (SA value) + 1 byte (chromosome)
} sa3_t;

//--------------------------------------------------------------------------------------

sa3_t *sa3_compute(sa_genome3_t *genome);

//--------------------------------------------------------------------------------------

inline sa3_t *sa3_new(uint num_suffixes) {
  sa3_t *p = (sa3_t *) malloc(sizeof(sa3_t));
  p->num_suffixes = num_suffixes;
  p->size = num_suffixes * (sizeof(uint) + 1);
  p->data = (char *) malloc(p->size);
  return p;
}

//--------------------------------------------------------------------------------------

inline void sa3_free(sa3_t *p) {
  if (p) {
    if (p->data) free(p->data);
    free(p);
  }
}

//--------------------------------------------------------------------------------------

inline void sa3_set(uint index, uint value, char chrom, sa3_t *SA) {
  uint pos = index * 5;
  *((uint *) &SA->data[pos]) = value;

  pos += 4;
  SA->data[pos] = chrom;
}

//--------------------------------------------------------------------------------------

inline void sa3_get(uint index, uint *value, char *chrom, sa3_t *SA) {
  uint pos = index * 5;
  *value = *((uint *) &SA->data[pos]);

  pos += 4;
  *chrom = SA->data[pos];
}

//--------------------------------------------------------------------------------------

inline uint sa3_get_value(uint index, sa3_t *SA) {
  return *((uint *) &SA->data[index * 5]);
}

//--------------------------------------------------------------------------------------

inline void sa3_save(char *filename, sa3_t *SA) {
  // write SA to file for the next time
  FILE *f = fopen(filename, "wb");
  if (f == NULL) {
    printf("Error: could not open %s to write\n", filename);
    exit(-1);
  }
  fwrite(SA->data, sizeof(char), SA->size, f);
  fclose(f);
}
*/
//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

typedef struct sa_index3 {
  uint k_value;
  uint num_suffixes;
  uint prefix_length;
  uint A_items; // JA_items = A_items
  uint IA_items;
  char *CHROM;
  uint *PRE;
  uint *SA;
  uint *A;
  uint *IA;
  unsigned char *JA;
  sa_genome3_t *genome;
} sa_index3_t;

//--------------------------------------------------------------------------------------

void sa_index3_build(char *genome_filename, uint k_value, char *sa_index_dirname);
void sa_index3_build_k18(char *genome_filename, uint k_value, char *sa_index_dirname);

//--------------------------------------------------------------------------------------

sa_index3_t *sa_index3_new(char *sa_index_dirname);
void sa_index3_free(sa_index3_t *sa_index);

//--------------------------------------------------------------------------------------

inline void sa_index3_display(sa_index3_t *p) {
  if (!p) return;

  printf("Num. suffixes          : %lu\n", p->num_suffixes);
  printf("Prefix table length    : %lu\n", p->prefix_length);
  printf("Prefix length (k-value): %lu\n", p->k_value);
  printf("A length (JA length)   : %lu\n", p->A_items);
  printf("AI length              : %lu\n", p->IA_items);

  if (p->genome) sa_genome3_display(p->genome);

}

//--------------------------------------------------------------------------------------

inline void display_prefix(char *S, int len) {
  for (int i = 0; i < len; i++) {
    printf("%c", S[i]);
  }
}



//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

#endif // SA_INDEX3_H
