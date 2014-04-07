#ifndef _SA_DNA_COMMONS_H
#define _SA_DNA_COMMONS_H

#include "bioformats/fastq/fastq_batch_reader.h"
#include "aligners/bwt/bwt.h"

#include "options.h"
#include "buffers.h"
#include "cal_seeker.h"
#include "pair_server.h"

#include "sa/sa_index3.h"

//--------------------------------------------------------------------

#define NUM_COUNTERS 10
extern int counters[NUM_COUNTERS];

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

#define MAX_NUM_MISMATCHES    4
#define MAX_NUM_SUFFIXES   2000

//--------------------------------------------------------------------
// sa_mapping_batch
//--------------------------------------------------------------------

typedef struct sa_mapping_batch {

  int counters[NUM_COUNTERS];

  int bam_format;
  size_t num_reads;

  #ifdef _TIMING
  double func_times[NUM_TIMING];
  #endif

  options_t *options;

  array_list_t *fq_reads;
  array_list_t **mapping_lists;
} sa_mapping_batch_t;

//--------------------------------------------------------------------

static inline sa_mapping_batch_t *sa_mapping_batch_new(array_list_t *fq_reads) {
  fastq_read_t *read;
  size_t read_length, num_reads = array_list_size(fq_reads);

  sa_mapping_batch_t *p = (sa_mapping_batch_t *) malloc(sizeof(sa_mapping_batch_t));

  for (int i = 0; i < NUM_COUNTERS; i++) {
    p->counters[i] = 0;
  }

  p->bam_format = 0;
  p->num_reads = num_reads;

  p->fq_reads = fq_reads;
  p->mapping_lists = (array_list_t **) malloc(num_reads * sizeof(array_list_t *));
  for (size_t i = 0; i < num_reads; i++) {
    p->mapping_lists[i] = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  }

  #ifdef _TIMING
  for (int i = 0; i < NUM_TIMING; i++) {
    p->func_times[i] = 0;
  }
  #endif

  return p;
}  

//--------------------------------------------------------------------

static inline void sa_mapping_batch_free(sa_mapping_batch_t *p) {
  if (p) {
    if (p->fq_reads) { array_list_free(p->fq_reads, (void *) fastq_read_free); }
    if (p->mapping_lists) { free(p->mapping_lists); }
    free(p);
  }
}  

//--------------------------------------------------------------------
// sa_wf_batch_t
//--------------------------------------------------------------------

typedef struct sa_wf_batch {
  options_t *options;
  sa_index3_t *sa_index;
  batch_writer_input_t *writer_input;
  void *mapping_batch;  
  void *data;
} sa_wf_batch_t;

//--------------------------------------------------------------------

static inline sa_wf_batch_t *sa_wf_batch_new(options_t *options,
					     sa_index3_t *sa_index,
					     batch_writer_input_t *writer_input,
					     void *mapping_batch,
					     void *data) {  
  sa_wf_batch_t *p = (sa_wf_batch_t *) malloc(sizeof(sa_wf_batch_t));

  p->options = options;
  p->sa_index = sa_index;
  p->writer_input = writer_input;

  p->mapping_batch = mapping_batch;
  p->data          = data;

  return p;
}

//--------------------------------------------------------------------

static inline void sa_wf_batch_free(sa_wf_batch_t *p) {
  if (p) free(p);
}

//--------------------------------------------------------------------
// sa_wf_input_t
//--------------------------------------------------------------------

typedef struct sa_wf_input {
  int bam_format;
  fastq_batch_reader_input_t *fq_reader_input;
  sa_wf_batch_t *wf_batch;
} sa_wf_input_t;

//--------------------------------------------------------------------

static inline sa_wf_input_t *sa_wf_input_new(int bam_format,
					     fastq_batch_reader_input_t *fq_reader_input,
					     sa_wf_batch_t *wf_batch) {
  sa_wf_input_t *p = (sa_wf_input_t *) malloc(sizeof(sa_wf_input_t));
  p->bam_format = bam_format;
  p->fq_reader_input = fq_reader_input;
  p->wf_batch = wf_batch;
  return p;
}

//--------------------------------------------------------------------

static inline void sa_wf_input_free(sa_wf_input_t *p) {
  if (p) free(p);
}

//--------------------------------------------------------------------
// cigar_t
//--------------------------------------------------------------------
#define NUM_INITIAL_CIGAR_OPS 200

typedef struct cigar {
  int num_ops;
  int num_allocated_ops;
  uint32_t ops_array[NUM_INITIAL_CIGAR_OPS];
  uint32_t *ops_pointer;
  uint32_t *ops;
} cigar_t;

//--------------------------------------------------------------------

static inline void cigar_init(cigar_t *p) {
  p->num_ops = 0;
  p->num_allocated_ops = NUM_INITIAL_CIGAR_OPS;
  p->ops = p->ops_array;
  p->ops_pointer = NULL;
}

//--------------------------------------------------------------------

static inline void cigar_clean(cigar_t *p) {
  if (p->ops_pointer) {
    free(p->ops_pointer);
  }
}

//--------------------------------------------------------------------

static inline cigar_t *cigar_new(int value, int name) {
  cigar_t *p = (cigar_t *) malloc(sizeof(cigar_t));
  cigar_init(p);
  p->ops[0] = ((value << 8) | (name & 255));
  p->num_ops = 1;
  //  printf("+++++ cigar_new: p = %x\n", p);
  return p;
}

//--------------------------------------------------------------------

static inline cigar_t *cigar_new_empty() {
  cigar_t *p = (cigar_t *) malloc(sizeof(cigar_t));
  cigar_init(p);
  p->num_ops = 0;
  //  printf("+++++ cigar_new_empty: p = %x\n", p);
  return p;
}

//--------------------------------------------------------------------

static inline void cigar_free(cigar_t *p) {
  //  printf("---------- cigar_free: p = %x (%s)\n", p, cigar_to_string(p));
  //  printf("---------- cigar_free: p = %x\n", p);
  if (p) {
    cigar_clean(p);
    free(p);
  }
}

//--------------------------------------------------------------------

static inline void cigar_set_zero(cigar_t *p) {
  if (p) {
    p->num_ops = 0;
  }
}

//--------------------------------------------------------------------

static inline void cigar_get_op(int index, int *value, int *name, cigar_t *p) {
  //  printf("cigar_get_op: (index, num_ops) = (%i, %i) %s\n", 
  //	 index, p->num_ops, cigar_to_string(p));
  if (index >= p->num_ops) {
    printf("cigar_get_op: (index, num_ops) = (%i, %i)\n", index, p->num_ops);
    assert(index < p->num_ops);
  }
  *name = (p->ops[index] & 255);
  *value = (p->ops[index] >> 8);
}

//--------------------------------------------------------------------


static inline char *cigar_to_string(cigar_t *p) {
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

static inline char *cigar_to_M_string(int *num_mismatches, int *num_cigar_ops, cigar_t *p) {
  char *str = (char *) malloc(p->num_ops * 10);
  int name, value, num_ops = 0, num_m = 0, mis = 0;
  str[0] = 0;
  for (int i = 0; i < p->num_ops; i++) {
    cigar_get_op(i, &value, &name, p);
    if (name == '=') {
      num_m += value;
    } else if (name == 'X') {
      num_m += value;
      mis += value;
    } else {
      if (num_m > 0) {
	num_ops++;
	sprintf(str, "%s%iM", str, num_m);
	num_m = 0;
      }
      num_ops++;
      sprintf(str, "%s%i%c", str, value, name);
    }
  }
  if (num_m > 0) {
    num_ops++;
    sprintf(str, "%s%iM", str, num_m);
  }

  *num_mismatches = mis;
  *num_cigar_ops = num_ops;
  return str;
}

//--------------------------------------------------------------------

static inline void cigar_set_op(int index, int value, int name, cigar_t *p) {
  p->ops[index] = ((value << 8) | (name & 255));
}

//--------------------------------------------------------------------
/*
static inline void _cigar_append_op(int value, int name, cigar_t *p) {
  if ((p->num_ops + 1) == p->num_allocated_ops) {
    printf("_cigar_append_op: (num_ops, num_allocated_ops) = (%i, %i)", 
	   p->num_ops, p->num_allocated_ops);
    exit(-1);
  }
  cigar_set_op(p->num_ops, value, name, p);
  p->num_ops++;
}
*/
//--------------------------------------------------------------------

static inline void cigar_append_op(int value, int name, cigar_t *p) {
  if (p->num_ops == 0) {
    cigar_set_op(0, value, name, p);
    p->num_ops++;
  } else {
    int last_op_name, last_op_value;
    cigar_get_op(p->num_ops - 1, &last_op_value, &last_op_name, p);
    if (last_op_name == name) {
      cigar_set_op(p->num_ops - 1, last_op_value + value, name, p);      
    } else {
      
      if (p->num_allocated_ops <= p->num_ops) {
	p->num_allocated_ops += (2 * NUM_INITIAL_CIGAR_OPS);
	uint32_t *aux = malloc(p->num_allocated_ops * sizeof(uint32_t));
	memcpy(aux, p->ops, p->num_ops * sizeof(uint32_t));
	if (p->ops_pointer) {
	  free(p->ops_pointer);
	}
	p->ops_pointer = aux;
	p->ops = p->ops_pointer;
	//	printf("cigar_append_op: (num_ops, num_allocated_ops) = (%i, %i)\n", 
	//	       p->num_ops, p->num_allocated_ops);
	//	exit(-1);
      }

      cigar_set_op(p->num_ops, value, name, p);
      p->num_ops++;      
    }
  }
}

//--------------------------------------------------------------------

static inline void cigar_concat(cigar_t *src, cigar_t *dst) {
  //check if there is space in dst to copy src
  if (dst->num_allocated_ops < (dst->num_ops + src->num_ops)) {
    dst->num_allocated_ops += (src->num_ops + (2 * NUM_INITIAL_CIGAR_OPS));
    uint32_t *aux = malloc(dst->num_allocated_ops * sizeof(uint32_t));
    memcpy(aux, dst->ops, dst->num_ops * sizeof(uint32_t));
    if (dst->ops_pointer) {
      free(dst->ops_pointer);
    }
    dst->ops_pointer = aux;
    dst->ops = dst->ops_pointer;
   //    printf("cigar_concat: dst->num_allocated_ops = %i < (dst->num_ops = %i + src->num_ops = %i)\n", 
    //	   dst->num_allocated_ops, dst->num_ops, src->num_ops);
    //    exit(-1);
  }
  
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

static inline void cigar_copy(cigar_t *dst, cigar_t *src) {
  if (src->num_ops > 0) {
    //check if there is space in dst to copy src
    if (dst->num_allocated_ops < src->num_ops) {
      dst->num_allocated_ops = (src->num_ops + (2 * NUM_INITIAL_CIGAR_OPS));
      uint32_t *aux = malloc(dst->num_allocated_ops * sizeof(uint32_t));
      memcpy(aux, dst->ops, dst->num_ops * sizeof(uint32_t));
      if (dst->ops_pointer) {
	free(dst->ops_pointer);
      }
      dst->ops_pointer = aux;
      dst->ops = dst->ops_pointer;
      //      printf("cigar_copy: dst->num_allocated_ops = %i < src->num_ops = %i)\n", 
      //	   dst->num_allocated_ops, src->num_ops);
      //      exit(-1);
    }

    dst->num_ops = src->num_ops;
    memcpy(dst->ops, src->ops, src->num_ops * sizeof(uint32_t));
  } else {
    cigar_set_zero(dst);
  }
}

//--------------------------------------------------------------------

static inline void cigar_revcopy(cigar_t *dst, cigar_t *src) {
  if (src->num_ops > 0) {
    //check if there is space in dst to copy src
    if (dst->num_allocated_ops < src->num_ops) {
      dst->num_allocated_ops = (src->num_ops + (2 * NUM_INITIAL_CIGAR_OPS));
      uint32_t *aux = malloc(dst->num_allocated_ops * sizeof(uint32_t));
      memcpy(aux, dst->ops, dst->num_ops * sizeof(uint32_t));
      if (dst->ops_pointer) {
	free(dst->ops_pointer);
      }
      dst->ops_pointer = aux;
      dst->ops = dst->ops_pointer;
      //      printf("cigar_revcopy: dst->num_allocated_ops = %i < src->num_ops = %i)\n", 
      //	   dst->num_allocated_ops, src->num_ops);
      //      exit(-1);
    }

    dst->num_ops = src->num_ops;
    for (int i = 0, j = src->num_ops - 1; i < src->num_ops; i++, j--) {
      dst->ops[i] = src->ops[j];
    }
  }
}

//--------------------------------------------------------------------

static inline void cigar_rev(cigar_t *p) {
  if (p->num_ops > 0) {
    cigar_t aux;
    cigar_init(&aux);
    cigar_copy(&aux, p);
    for (int i = 0, j = p->num_ops - 1; i < p->num_ops; i++, j--) {
      p->ops[i] = aux.ops[j];
    }
    cigar_clean(&aux);
  }
}

//--------------------------------------------------------------------

static inline int cigar_get_length(cigar_t *p) {
  int op_value, op_name, num_ops, len = 0;
  if ((num_ops = p->num_ops) > 0) {
    for (int i = 0; i < num_ops; i++) {
      cigar_get_op(i, &op_value, &op_name, p);
      if (op_name == 'X' || op_name == '=' || op_name == 'I' || op_name == 'S' || op_name== 'M') {
	len += op_value;
      }
    }
  }
  return len;
}

//--------------------------------------------------------------------

static inline void cigar_ltrim_ops(int len, cigar_t *p) {
  int num_ops;
  if ((num_ops = p->num_ops) > 0) {
    if (len >= num_ops) {
      cigar_set_zero(p);
    } else {
      num_ops -= len;
      for (int i = 0, j = len; i < num_ops; i++, j++) {
	p->ops[i] = p->ops[j];
      }
      p->num_ops = num_ops;
      //	memcpy(p->ops, &p->ops[len], p->num_ops * sizeof(uint32_t));
    }
  }
}

//--------------------------------------------------------------------

static inline void cigar_rtrim_ops(int len, cigar_t *p) {
  int num_ops;
  if ((num_ops = p->num_ops) > 0) {
    if (len >= num_ops) {
      cigar_set_zero(p);
    } else {
      p->num_ops -= len;
    }
  }
}

//--------------------------------------------------------------------

static inline float cigar_compute_score(float match_score, float mismatch_penalty,
				 float open_gap_penalty, float extend_gap_penalty,
				 cigar_t *p) {
  int num_matches = 0, num_mismatches = 0;
  int num_open_gaps = 0, num_extend_gaps = 0;
  int op_value, op_name, num_ops = p->num_ops;

  for (int i = 0; i < num_ops; i++) {
    cigar_get_op(i, &op_value, &op_name, p);
    if (op_name == '=') {
      num_matches += op_value;
    } else if (op_name == 'X') {
      num_mismatches += op_value;
    } else if (op_name == 'I' || op_name == 'D') {
      num_open_gaps++;
      num_extend_gaps += (op_value - 1);
    }
  }

  return (num_matches * match_score) + (num_mismatches * mismatch_penalty) +
    (num_open_gaps * open_gap_penalty) + (num_extend_gaps * extend_gap_penalty);
}

				 

//--------------------------------------------------------------------
// seed_t
//--------------------------------------------------------------------

typedef struct seed {
  size_t read_start;
  size_t read_end;
  size_t genome_start;
  size_t genome_end;

  size_t suf_read_start;
  size_t suf_read_end;
  size_t suf_genome_start;
  size_t suf_genome_end;

  int strand;
  int chromosome_id;
  int num_mismatches;
  int num_open_gaps;
  int num_extend_gaps;

  cigar_t cigar;
} seed_t;

//--------------------------------------------------------------------

static inline seed_t *seed_new(size_t read_start, size_t read_end,
			size_t genome_start, size_t genome_end) {

  seed_t *p = (seed_t *) malloc(sizeof(seed_t));

  p->read_start = read_start;
  p->read_end = read_end;
  p->genome_start = genome_start;
  p->genome_end = genome_end;

  p->suf_read_start = read_start;
  p->suf_read_end = read_end;
  p->suf_genome_start = genome_start;
  p->suf_genome_end = genome_end;

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

//--------------------------------------------------------------------

static inline void seed_ltrim_read(int len, seed_t *p) {
  // look at left-side cigar
  cigar_t *cigar = &p->cigar;
  int op_value, op_name;
  int remain = len, q_trim = 0, r_trim = 0, trim_ops = 0, num_ops = cigar->num_ops;

  for (int i = 0; i < num_ops; i++) {
    cigar_get_op(i, &op_value, &op_name, cigar);

    if (op_name == '=' || op_name == 'X') {

      if (op_value >= remain) {
	q_trim += remain;
	r_trim += remain;
	cigar_set_op(i, op_value - remain, op_name, cigar);
	break;
      } else {
	q_trim += op_value;
	r_trim += op_value;
	remain -= op_value;
	trim_ops++;
      }

    } if (op_name == 'I') {

      if (op_value >= remain) {
	q_trim += remain;
	cigar_set_op(i, op_value - remain, op_name, cigar);
	break;
      } else {
	q_trim += op_value;
	remain -= op_value;
	trim_ops++;
      }

    } else if (op_name == 'D') {

      r_trim += op_value;
      trim_ops++;

    }
  }
  p->read_start += q_trim;
  p->genome_start += r_trim;
  if (trim_ops) {
    cigar_ltrim_ops(trim_ops, cigar);
  }
}

//--------------------------------------------------------------------

static inline void seed_rtrim_read(int len, seed_t *p) {
  // look at right-side cigar
  cigar_t *cigar = &p->cigar;
  int op_value, op_name;
  int remain = len, q_trim = 0, r_trim = 0, trim_ops = 0, num_ops = cigar->num_ops;

  for (int i = num_ops - 1; i >= 0; i--) {
    cigar_get_op(i, &op_value, &op_name, cigar);

    if (op_name == '=' || op_name == 'X') {

      if (op_value >= remain) {
	q_trim += remain;
	r_trim += remain;
	cigar_set_op(i, op_value - remain, op_name, cigar);
	break;
      } else {
	q_trim += op_value;
	r_trim += op_value;
	remain -= op_value;
	trim_ops++;
      }

    } if (op_name == 'I') {

      if (op_value >= remain) {
	q_trim += remain;
	cigar_set_op(i, op_value - remain, op_name, cigar);
	break;
      } else {
	q_trim += op_value;
	remain -= op_value;
	trim_ops++;
      }

    } else if (op_name == 'D') {

      r_trim += op_value;
      trim_ops++;

    }
  }
  p->read_end -= q_trim;
  p->genome_end -= r_trim;
  if (trim_ops) {
    cigar_rtrim_ops(trim_ops, cigar);
  }
}

//--------------------------------------------------------------------

static inline void seed_ltrim_genome(int len, seed_t *p) {
  // look at left-side cigar
  cigar_t *cigar = &p->cigar;
  int op_value, op_name;
  int remain = len, q_trim = 0, r_trim = 0, trim_ops = 0, num_ops = cigar->num_ops;

  for (int i = 0; i < num_ops; i++) {
    cigar_get_op(i, &op_value, &op_name, cigar);

    if (op_name == '=' || op_name == 'X') {

      if (op_value >= remain) {
	q_trim += remain;
	r_trim += remain;
	cigar_set_op(i, op_value - remain, op_name, cigar);
	break;
      } else {
	q_trim += op_value;
	r_trim += op_value;
	remain -= op_value;
	trim_ops++;
      }

    } if (op_name == 'D') {

      if (op_value >= remain) {
	q_trim += remain;
	cigar_set_op(i, op_value - remain, op_name, cigar);
	break;
      } else {
	q_trim += op_value;
	remain -= op_value;
	trim_ops++;
      }

    } else if (op_name == 'I') {

      r_trim += op_value;
      trim_ops++;

    }
  }
  p->read_start += q_trim;
  p->genome_start += r_trim;
  if (trim_ops) {
    cigar_ltrim_ops(trim_ops, cigar);
  }
}

//--------------------------------------------------------------------

static inline void seed_rtrim_genome(int len, seed_t *p) {
  // look at right-side cigar
  cigar_t *cigar = &p->cigar;
  int op_value, op_name;
  int remain = len, q_trim = 0, r_trim = 0, trim_ops = 0, num_ops = cigar->num_ops;

  for (int i = num_ops - 1; i >= 0; i--) {
    cigar_get_op(i, &op_value, &op_name, cigar);

    if (op_name == '=' || op_name == 'X') {

      if (op_value >= remain) {
	q_trim += remain;
	r_trim += remain;
	cigar_set_op(i, op_value - remain, op_name, cigar);
	break;
      } else {
	q_trim += op_value;
	r_trim += op_value;
	remain -= op_value;
	trim_ops++;
      }

    } if (op_name == 'D') {

      if (op_value >= remain) {
	q_trim += remain;
	cigar_set_op(i, op_value - remain, op_name, cigar);
	break;
      } else {
	q_trim += op_value;
	remain -= op_value;
	trim_ops++;
      }

    } else if (op_name == 'I') {

      r_trim += op_value;
      trim_ops++;

    }
  }
  p->read_end -= q_trim;
  p->genome_end -= r_trim;
  if (trim_ops) {
    cigar_rtrim_ops(trim_ops, cigar);
  }
}

//--------------------------------------------------------------------
// cigarset_t
//--------------------------------------------------------------------
typedef struct cigarset_info {
  int active;
  int overlap;
  cigar_t *cigar;
  seed_t *seed;
} cigarset_info_t;

static inline cigarset_info_t *cigarset_info_set(int active, int overlap,
					  cigar_t *cigar, seed_t *seed,
					  cigarset_info_t *p) {
  p->active = active;
  p->overlap = overlap;
  p->cigar = cigar;
  p->seed = seed;
}

//--------------------------------------------------------------------

typedef struct cigarset {
  int size;
  cigarset_info_t *info;
} cigarset_t;

//--------------------------------------------------------------------

static inline cigarset_t *cigarset_new(int size) {
  cigarset_t *p = (cigarset_t *) malloc(sizeof(cigarset_t));
  p->size = size;
  p->info = (cigarset_info_t *) malloc(size * sizeof(cigarset_info_t));
  return p;
}

//--------------------------------------------------------------------

static inline void cigarset_free(cigarset_t *p) {
  if (p) {
    if (p->info) free(p->info);
    free(p);
  }
}

//--------------------------------------------------------------------
// seed_cal_t
//--------------------------------------------------------------------

typedef struct seed_cal {
  size_t chromosome_id;
  short int strand;
  size_t start;
  size_t end;

  int AS; // SAM format flag
  int read_area;
  int invalid;

  int cigar_len;

  int num_mismatches;
  int num_open_gaps;
  int num_extend_gaps;

  float score;
  cigar_t cigar;

  fastq_read_t *read;
  linked_list_t *seed_list;
  cigarset_t *cigarset;
} seed_cal_t;

//--------------------------------------------------------------------

static inline seed_cal_t *seed_cal_new(const size_t chromosome_id,
				const short int strand,
				const size_t start,
				const size_t end,
				linked_list_t *seed_list) {

  seed_cal_t *p = (seed_cal_t *) malloc(sizeof(seed_cal_t));

  p->strand = strand;
  p->chromosome_id = chromosome_id;
  p->start = start;
  p->end = end;

  p->AS = 0;
  p->read_area = 0;
  p->invalid = 0;

  p->cigar_len = 0;

  p->num_mismatches = 0;
  p->num_open_gaps = 0;
  p->num_extend_gaps = 0;

  p->score = 0.0f;
  cigar_init(&p->cigar);

  p->read = NULL;
  p->seed_list = seed_list;
  p->cigarset = NULL;

  return p;
}

//--------------------------------------------------------------------

void seed_cal_free(seed_cal_t *p);

//--------------------------------------------------------------------

static inline void seed_cal_update_info(seed_cal_t *cal) {
  int read_area = 0, num_mismatches = 0, num_open_gaps = 0, num_extend_gaps = 0;
  seed_t *seed;

  for (linked_list_item_t *item = cal->seed_list->first; 
       item != NULL; item = item->next) {
    seed = item->item;
    num_mismatches += seed->num_mismatches;
    num_open_gaps += seed->num_open_gaps;
    num_extend_gaps += seed->num_extend_gaps;
    read_area += (seed->read_end - seed->read_start - num_mismatches - num_open_gaps - num_extend_gaps);
    
  }
  cal->num_mismatches = num_mismatches;
  cal->num_open_gaps = num_open_gaps;
  cal->num_extend_gaps = num_extend_gaps;
  cal->num_mismatches = num_mismatches;
  cal->read_area = read_area;
}

//--------------------------------------------------------------------

static inline void seed_cal_set_cigar_by_seed(seed_t *seed, seed_cal_t *cal) {
  cal->num_mismatches = seed->num_mismatches;
  cal->num_open_gaps = seed->num_open_gaps;
  cal->num_extend_gaps = seed->num_extend_gaps;
  cigar_concat(&seed->cigar, &cal->cigar);
}

//--------------------------------------------------------------------
void print_seed(char *msg, seed_t *s);

static inline void seed_cal_print(seed_cal_t *cal) {
  printf(" CAL (%c)[%lu:%lu-%lu] (%s, x:%i, og:%i, eg:%i) score = %0.2f: (read id %s)\n", 
	 (cal->strand == 0 ? '+' : '-'), 
	 cal->chromosome_id, cal->start, cal->end, cigar_to_string(&cal->cigar), cal->num_mismatches,
	 cal->num_open_gaps, cal->num_extend_gaps, cal->score,
	 cal->read->id);
  printf("\tSEEDS LIST:\n");
  if (cal->seed_list == NULL || cal->seed_list->size == 0) {
    printf("\t\tno seeds\n");
  } else {
    for (linked_list_item_t *item = cal->seed_list->first; 
	 item != NULL; item = item->next) {
      seed_t *seed = item->item;
      print_seed("\t\t", seed);
      //      printf(" [%lu|%lu - %lu|%lu] (%s, x:%i, og:%i, eg:%i)", seed->genome_start, seed->read_start, 
      //	     seed->read_end, seed->genome_end, cigar_to_string(&seed->cigar), seed->num_mismatches,
      //	     seed->num_open_gaps, seed->num_extend_gaps);
    }
    printf("\n");
  }
}

//--------------------------------------------------------------------
// utils
//--------------------------------------------------------------------

float get_max_score(array_list_t *cal_list);
int get_min_num_mismatches(array_list_t *cal_list);
int get_max_read_area(array_list_t *cal_list);

//--------------------------------------------------------------------

void filter_cals_by_min_read_area(int read_area, array_list_t **list);
void filter_cals_by_max_read_area(int read_area, array_list_t **list);
void filter_cals_by_max_score(float score, array_list_t **list);
void filter_cals_by_max_num_mismatches(int num_mismatches, array_list_t **list);
void filter_cals_by_pair_mode(int pair_mode, int pair_min_distance, int pair_max_distance, 
			      int num_lists, array_list_t **cal_lists);

//--------------------------------------------------------------------

void create_bam_alignments(array_list_t *cal_list, fastq_read_t *read, 
		       array_list_t *mapping_list);

void create_alignments(array_list_t *cal_list, fastq_read_t *read, 
		       int bam_format, array_list_t *mapping_list);

//--------------------------------------------------------------------

void display_suffix_mappings(int strand, size_t r_start, size_t suffix_len, 
			     size_t low, size_t high, sa_index3_t *sa_index);
void display_sequence(uint j, sa_index3_t *index, uint len);
char *get_subsequence(char *seq, size_t start, size_t len);
void display_cmp_sequences(fastq_read_t *read, sa_index3_t *sa_index);

//--------------------------------------------------------------------

void complete_pairs(sa_mapping_batch_t *batch);

//--------------------------------------------------------------------
//--------------------------------------------------------------------

#endif // _SA_DNA_COMMONS_H
