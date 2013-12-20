#ifndef _SA_DNA_COMMONS_H
#define _SA_DNA_COMMONS_H

#include "bioformats/fastq/fastq_batch_reader.h"
#include "aligners/bwt/bwt.h"

#include "buffers.h"
#include "cal_seeker.h"

#include "sa/sa_index3.h"

//--------------------------------------------------------------------

#ifdef _TIMING

#define FUNC_SEARCH_SUFFIX             0
#define FUNC_SEARCH_PREFIX             1
#define FUNC_SEARCH_SA                 2
#define FUNC_CALS_FROM_EXACT_READ      3
#define FUNC_CALS_FROM_SUFFIXES        4
#define FUNC_INIT_CALS_FROM_SUFFIXES   5
#define FUNC_SET_POSITIONS             6
#define FUNC_SET_REF_SEQUENCE          7
#define FUNC_SKIP_SUFFIXES             8
#define FUNC_MINI_SW_RIGHT_SIDE        9
#define FUNC_MINI_SW_LEFT_SIDE        10
#define FUNC_SEED_NEW                 11
#define FUNC_SEED_LIST_INSERT         12
#define FUNC_CAL_NEW                  13
#define FUNC_CAL_MNG_INSERT           14
#define FUNC_CAL_MNG_TO_LIST          15
#define FUNC_FILTER_BY_READ_AREA      16
#define FUNC_FILTER_BY_NUM_MISMATCHES 17
#define FUNC_PRE_SW                   18
#define FUNC_SW                       19
#define FUNC_POST_SW                  20
#define FUNC_OTHER                    21
#define FUNC_CREATE_ALIGNMENTS        22


#define NUM_TIMING (FUNC_CREATE_ALIGNMENTS + 1)
double func_times[NUM_TIMING];
char func_names[NUM_TIMING][1024];

#endif

//--------------------------------------------------------------------

#define MISMATCH_PERC 0.10f

#define MAX_NUM_MISMATCHES   4
#define MAX_NUM_SUFFIXES   2000

//--------------------------------------------------------------------
// sa_mapping_batch
//--------------------------------------------------------------------

typedef struct sa_mapping_batch {
  size_t num_reads;
  #ifdef _TIMING
  double func_times[NUM_TIMING];
  #endif
  array_list_t *fq_reads;
  char **revcomp_seqs;
  array_list_t **mapping_lists;
} sa_mapping_batch_t;

//--------------------------------------------------------------------

inline sa_mapping_batch_t *sa_mapping_batch_new(array_list_t *fq_reads) {
  char *revcomp, c;
  fastq_read_t *read;
  size_t read_length, num_reads = array_list_size(fq_reads);

  sa_mapping_batch_t *p = (sa_mapping_batch_t *) malloc(sizeof(sa_mapping_batch_t));
  p->num_reads = num_reads;
  p->fq_reads = fq_reads;
  p->revcomp_seqs = (char **) malloc(num_reads * sizeof(char *));
  p->mapping_lists = (array_list_t **) malloc(num_reads * sizeof(array_list_t *));
  for (size_t i = 0; i < num_reads; i++) {
    read = array_list_get(i, fq_reads);
    read_length = read->length;
    p->mapping_lists[i] = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    // prepare reverse complementary sequence
    revcomp = (char *) malloc(read_length + 1);
    for (int src = 0, dest = read_length - 1; src < read_length; src++, dest--) {
      c = read->sequence[src];
      if      (c == 'A') { revcomp[dest] = 'T'; }
      else if (c == 'T') { revcomp[dest] = 'A'; }
      else if (c == 'G') { revcomp[dest] = 'C'; }
      else if (c == 'C') { revcomp[dest] = 'G'; }
      else               { revcomp[dest] = c; }
    }
    revcomp[read_length] = 0;
    p->revcomp_seqs[i] = revcomp;
  }

  #ifdef _TIMING
  for (int i = 0; i < NUM_TIMING; i++) {
    p->func_times[i] = 0;
  }
  #endif

  return p;
}  

//--------------------------------------------------------------------

inline void sa_mapping_batch_free(sa_mapping_batch_t *p) {
  if (p) {
    if (p->fq_reads) { array_list_free(p->fq_reads, (void *) fastq_read_free); }
    if (p->mapping_lists) { free(p->mapping_lists); }
    if (p->revcomp_seqs) {
      for (size_t i = 0; i < p->num_reads; i++) {
	free(p->revcomp_seqs[i]);
      }
      free(p->revcomp_seqs);
    }
    free(p);
  }
}  

//--------------------------------------------------------------------
// sa_wf_batch_t
//--------------------------------------------------------------------

typedef struct sa_wf_batch {
  sa_index3_t *sa_index;
  batch_writer_input_t *writer_input;
  sa_mapping_batch_t *mapping_batch;  
} sa_wf_batch_t;

//--------------------------------------------------------------------

inline sa_wf_batch_t *sa_wf_batch_new(sa_index3_t *sa_index,
				      batch_writer_input_t *writer_input,
				      sa_mapping_batch_t *mapping_batch) {
  
  sa_wf_batch_t *p = (sa_wf_batch_t *) malloc(sizeof(sa_wf_batch_t));
  p->sa_index = sa_index;
  p->writer_input = writer_input;
  p->mapping_batch = mapping_batch;
  return p;
}

//--------------------------------------------------------------------

inline void sa_wf_batch_free(sa_wf_batch_t *p) {
  if (p) free(p);
}

//--------------------------------------------------------------------
// sa_wf_input_t
//--------------------------------------------------------------------

typedef struct sa_wf_input {
  fastq_batch_reader_input_t *fq_reader_input;
  sa_wf_batch_t *wf_batch;
} sa_wf_input_t;

//--------------------------------------------------------------------

inline sa_wf_input_t *sa_wf_input_new(fastq_batch_reader_input_t *fq_reader_input,
				      sa_wf_batch_t *wf_batch) {
  sa_wf_input_t *p = (sa_wf_input_t *) malloc(sizeof(sa_wf_input_t));
  p->fq_reader_input = fq_reader_input;
  p->wf_batch = wf_batch;
  return p;
}

//--------------------------------------------------------------------

inline void sa_wf_input_free(sa_wf_input_t *p) {
  if (p) free(p);
}

//--------------------------------------------------------------------
// cigar_t
//--------------------------------------------------------------------

typedef struct cigar {
  uint32_t ops[100];
  int num_ops;
} cigar_t;

//--------------------------------------------------------------------

inline cigar_t *cigar_new(int value, int name) {
  cigar_t *p = (cigar_t *) malloc(sizeof(cigar_t));
  p->ops[0] = ((value << 8) | (name & 255));
  p->num_ops = 1;
  //  printf("+++++ cigar_new: p = %x\n", p);
  return p;
}

//--------------------------------------------------------------------

inline cigar_t *cigar_new_empty() {
  cigar_t *p = (cigar_t *) malloc(sizeof(cigar_t));
  p->num_ops = 0;
  //  printf("+++++ cigar_new_empty: p = %x\n", p);
  return p;
}

//--------------------------------------------------------------------

inline void cigar_init(cigar_t *p) {
  if (p) {
    p->num_ops = 0;
  }
}

//--------------------------------------------------------------------

inline void cigar_get_op(int index, int *value, int *name, cigar_t *p) {
  assert(index < p->num_ops);
  *name = (p->ops[index] & 255);
  *value = (p->ops[index] >> 8);
}

//--------------------------------------------------------------------

inline char *cigar_to_string(cigar_t *p) {
  char *str = (char *) malloc(p->num_ops * 10);
  int name, value;
  str[0] = 0;
  for (int i = 0; i < p->num_ops; i++) {
    cigar_get_op(i, &value, &name, p);
    sprintf(str, "%s%i%c", str, value, name);
  }
  return str;
}

//--------------------------------------------------------------------

inline void cigar_free(cigar_t *p) {
  //  printf("---------- cigar_free: p = %x (%s)\n", p, cigar_to_string(p));
  //  printf("---------- cigar_free: p = %x\n", p);
  if (p) free(p);
}

//--------------------------------------------------------------------

inline void cigar_set_op(int index, int value, int name, cigar_t *p) {
  p->ops[index] = ((value << 8) | (name & 255));
}

//--------------------------------------------------------------------

inline void _cigar_append_op(int value, int name, cigar_t *p) {
  cigar_set_op(p->num_ops, value, name, p);
  p->num_ops++;
}

//--------------------------------------------------------------------

inline void cigar_append_op(int value, int name, cigar_t *p) {
  if (p->num_ops == 0) {
    cigar_set_op(0, value, name, p);
    p->num_ops++;
  } else {
    int last_op_name, last_op_value;
    cigar_get_op(p->num_ops - 1, &last_op_value, &last_op_name, p);
    if (last_op_name == name) {
      cigar_set_op(p->num_ops - 1, last_op_value + value, name, p);      
    } else {
      cigar_set_op(p->num_ops, value, name, p);
      p->num_ops++;      
    }
  }
}

//--------------------------------------------------------------------

inline void cigar_concat(cigar_t *src, cigar_t *dst) {
  if (dst->num_ops == 0) {
    memcpy(dst->ops, src->ops, src->num_ops * sizeof(uint32_t));
    dst->num_ops = src->num_ops;
  } else {
    int first_op_name, first_op_value, last_op_name, last_op_value;
    cigar_get_op(0, &first_op_value, &first_op_name, src);
    cigar_get_op(dst->num_ops - 1, &last_op_value, &last_op_name, dst);
    if (first_op_name == last_op_name) {
      cigar_set_op(dst->num_ops - 1, last_op_value + first_op_value, last_op_name, dst);
      if (src->num_ops > 1) {
	memcpy(&dst->ops[dst->num_ops], &src->ops[1], (src->num_ops - 1) * sizeof(uint32_t));
	dst->num_ops += (src->num_ops - 1);
      }
    } else {
      memcpy(&dst->ops[dst->num_ops], src->ops, src->num_ops * sizeof(uint32_t));
      dst->num_ops += src->num_ops;
    }
  }
}

//--------------------------------------------------------------------

inline void cigar_copy(cigar_t *dst, cigar_t *src) {
  if (src->num_ops > 0) {
    dst->num_ops = src->num_ops;
    memcpy(dst->ops, src->ops, src->num_ops * sizeof(uint32_t));
  }
}

//--------------------------------------------------------------------

inline void cigar_revcopy(cigar_t *dst, cigar_t *src) {
  if (src->num_ops > 0) {
    dst->num_ops = src->num_ops;
    for (int i = 0, j = src->num_ops - 1; i < src->num_ops; i++, j--) {
      dst->ops[i] = src->ops[j];
    }
  }
}

//--------------------------------------------------------------------

inline void cigar_rev(cigar_t *p) {
  if (p->num_ops > 0) {
    cigar_t aux;
    cigar_copy(&aux, p);
    for (int i = 0, j = p->num_ops - 1; i < p->num_ops; i++, j--) {
      p->ops[i] = aux.ops[j];
    }
  }
}

//--------------------------------------------------------------------
// cigarset_t
//--------------------------------------------------------------------

typedef struct cigarset {
  int num_cigars;
  int *active;
  cigar_t **cigars;
} cigarset_t;

//--------------------------------------------------------------------

inline cigarset_t *cigarset_new(int num_cigars) {
  cigarset_t *p = (cigarset_t *) malloc(sizeof(cigarset_t));
  p->num_cigars = num_cigars;
  p->active = (int *) calloc(num_cigars, sizeof(int));
  p->cigars = (cigar_t **) malloc(num_cigars * sizeof(cigar_t*));
  return p;
}

//--------------------------------------------------------------------

inline void cigarset_free(cigarset_t *p) {
  if (p) {
    if (p->active) free(p->active);
    if (p->cigars) free(p->cigars);
    free(p);
  }
}

//--------------------------------------------------------------------
// seed_t
//--------------------------------------------------------------------

typedef struct seed {
  size_t read_start;
  size_t read_end;
  size_t genome_start;
  size_t genome_end;

  int strand;
  int chromosome_id;
  int num_mismatches;
  int num_open_gaps;
  int num_extend_gaps;

  cigar_t cigar;
} seed_t;

//--------------------------------------------------------------------

inline seed_t *seed_new(size_t read_start, size_t read_end, 
			size_t genome_start, size_t genome_end) {

  seed_t *p = (seed_t *) malloc(sizeof(seed_t));

  p->read_start = read_start;
  p->read_end = read_end;
  p->genome_start = genome_start;
  p->genome_end = genome_end;

  p->strand = 0;
  p->chromosome_id = 0;
  p->num_mismatches = 0;
  p->num_open_gaps = 0;
  p->num_extend_gaps = 0;


  cigar_init(&p->cigar);

  return p;
}

//--------------------------------------------------------------------

void seed_free(seed_t *p);
void cal_free_ex(cal_t *cal);

//--------------------------------------------------------------------
// utils
//--------------------------------------------------------------------

int get_min_num_mismatches(array_list_t *cal_list);
int get_max_read_area(array_list_t *cal_list);
void filter_cals_by_min_read_area(int read_area, array_list_t **list);
void filter_cals_by_max_read_area(int read_area, array_list_t **list);
void filter_cals_by_max_num_mismatches(int num_mismatches, array_list_t **list);

void create_alignments(array_list_t *cal_list, fastq_read_t *read, 
		       array_list_t *mapping_list);

void display_suffix_mappings(int strand, size_t r_start, size_t suffix_len, 
			     size_t low, size_t high, sa_index3_t *sa_index);
void print_seed(char *msg, seed_t *s);
void display_sequence(uint j, sa_index3_t *index, uint len);
char *get_subsequence(char *seq, size_t start, size_t len);
void display_cmp_sequences(fastq_read_t *read, char *revcomp_seq, sa_index3_t *sa_index);

//--------------------------------------------------------------------
//--------------------------------------------------------------------

#endif // _SA_DNA_COMMONS_H
