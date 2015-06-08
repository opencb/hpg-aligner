#ifndef DOSCADFUN_H
#define DOSCADFUN_H

#include "dna/sa/sa_dna_commons.h"

//--------------------------------------------------------------------

#define DEL_FINAL       10
#define MAX_NUM_ERRORS   1

//--------------------------------------------------------------------
//--------------------------------------------------------------------

typedef struct alig_out {
  int map_len1;
  int map_len2;
  int match;
  int mismatch;
  int gap_open; 
  int gap_extend;
  int st_map_len;
  // char *st1_map;
  // char *st2_map;
  float score;
  cigar_t cigar;
} alig_out_t;

//--------------------------------------------------------------------

static inline void alig_out_init(alig_out_t *out) {
  out->map_len1 = 0;
  out->map_len2 = 0;
  out->match = 0;
  out->mismatch = 0;
  out->gap_open = 0;
  out->gap_extend = 0;
  out->st_map_len = 0;
  //  cigar_init(&out->cigar);
  cigar_set_zero(&out->cigar);
  // 	char *st1_map;	char *st2_map;
}

//--------------------------------------------------------------------

static inline void alig_out_set(int map_l1, int map_l2, int mm, int mism,
				int gap1, int gapmas, int st_mm_len, alig_out_t *out) {
  out->map_len1 = map_l1;
  out->map_len2 = map_l2;
  out->match = mm;
  out->mismatch = mism;
  out->gap_open = gap1;
  out->gap_extend = gapmas;
  out->st_map_len = st_mm_len;
  // 	char *st1_map;	char *st2_map;
}

//--------------------------------------------------------------------

static inline void alig_out_print(alig_out_t *out) {
  printf("matches = %i, mismatches = %i, open-gaps = %i, extended-gaps = %i, (query map len, ref. map.len) = (%i, %i)\n",
	 out->match, out->mismatch, out->gap_open, out->gap_extend,
	 out->map_len1, out->map_len2);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

float doscadfun(char* st1, int ll1, char* st2, int ll2, 
		float mismatch_perc, alig_out_t *out);

float doscadfun_inv(char* st1, int ll1, char* st2, int ll2, 
		    float mismatch_perc, alig_out_t *out);
//--------------------------------------------------------------------
//--------------------------------------------------------------------

#endif // DOSCADFUN_H
