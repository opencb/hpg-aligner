#include <stdio.h>
#include <stdlib.h>

#include "commons/workflow_scheduler.h"
#include "bioformats/bam/bam_file.h"

#include "batch_writer.h"

#include "dna/commons.h"
#include "dna/doscadfun.h"
#include "sa/sa_index3.h"


// SA index3 new + Compressed Row Storage
// SA mapper based on workflow_scheduler
// commons-lib/commons/workflow_scheduler

#define MISMATCH_PERC 0.10f

#define MAX_NUM_MISMATCHES   4
#define MAX_NUM_SUFFIXES   2000

//--------------------------------------------------------------------

char *get_subsequence(char *seq, size_t start, size_t len) {
  char *subseq = (char *) malloc((len + 1) * sizeof(char));
  memcpy(subseq, seq + start, len);
  subseq[len] = 0;
  return subseq;
}

//--------------------------------------------------------------------

void display_cmp_sequences(fastq_read_t *read, char *revcomp_seq, sa_index3_t *sa_index) {
  size_t pos, chrom, strand;
  char *ref, *seq, *chrom_str, *aux, *p1, *p2;
  
  aux = strdup(read->id);
  p1 = strstr(aux, "_");
  *p1 = 0;
  chrom_str = strdup(aux);
  for (chrom = 0; chrom < sa_index->genome->num_chroms; chrom++) {
    if (strcmp(chrom_str, sa_index->genome->chrom_names[chrom]) == 0) {
      break;
    }
  }
  p2 = strstr(p1 + 1, "_");
  *p2 = 0;
  pos = atol(p1 + 1);
  
  p1 = strstr(p2 + 1, "_");
  p2 = strstr(p1 + 1, "_");
  *p2 = 0;
  strand = atoi(p1 + 1);
  
  free(aux);
  free(chrom_str);
  
  printf("\n\n======> %s\n", read->id);
  for (int i = 0; i < read->length; i++) {
    if (i % 10 == 0) {
      printf("%i", i / 10);
    } else {
      printf(" ");
    }
  }
  printf("\n");
  for (int i = 0; i < read->length; i++) {
    printf("%i", i % 10);
  }
  printf("\n");
  if (strand) {
    seq = revcomp_seq;
  } else {
   seq = read->sequence;
  }
  printf("%s\n", seq);
  ref = &sa_index->genome->S[pos + sa_index->genome->chrom_offsets[chrom] - 1];
  for (int i = 0; i < read->length; i++) {
    if (seq[i] == ref[i]) {
      printf("|");
    } else {
      printf("x");
    }
  }
  printf("\n");
  for (int i = 0; i < read->length; i++) {
    printf("%c", ref[i]);
  }
  printf("\n");
  //  printf("chrom = %i\n", chrom);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

typedef struct cigar {
  uint32_t ops[100];
  int num_ops;
} cigar_t;

inline cigar_t *cigar_new(int value, int name) {
  cigar_t *p = (cigar_t *) malloc(sizeof(cigar_t));
  p->ops[0] = ((value << 8) | (name & 255));
  p->num_ops = 1;
  //  printf("+++++ cigar_new: p = %x\n", p);
  return p;
}

inline cigar_t *cigar_new_empty() {
  cigar_t *p = (cigar_t *) malloc(sizeof(cigar_t));
  p->num_ops = 0;
  //  printf("+++++ cigar_new_empty: p = %x\n", p);
  return p;
}

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

inline void cigar_free(cigar_t *p) {
  //  printf("---------- cigar_free: p = %x (%s)\n", p, cigar_to_string(p));
  //  printf("---------- cigar_free: p = %x\n", p);
  if (p) free(p);
}

inline void cigar_init(cigar_t *p) {
  p->num_ops = 0;
}

inline void cigar_get_op(int index, int *value, int *name, cigar_t *p) {
  assert(index < p->num_ops);
  *name = (p->ops[index] & 255);
  *value = (p->ops[index] >> 8);
}

inline void cigar_set_op(int index, int value, int name, cigar_t *p) {
  p->ops[index] = ((value << 8) | (name & 255));
}

inline void _cigar_append_op(int value, int name, cigar_t *p) {
  cigar_set_op(p->num_ops, value, name, p);
  p->num_ops++;
}

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

//====================================================================
// STRUCTURES
//====================================================================

#ifdef _TIMING

#define FUNC_SEARCH_SUFFIX           0
#define FUNC_SEARCH_PREFIX           1
#define FUNC_SEARCH_SA               2
#define FUNC_CALS_FROM_EXACT_READ    3
#define FUNC_CALS_FROM_SUFFIXES      4
#define FUNC_INIT_CALS_FROM_SUFFIXES 5
#define FUNC_SET_POSITIONS           6
#define FUNC_SET_REF_SEQUENCE        7
#define FUNC_SKIP_SUFFIXES           8
#define FUNC_MINI_SW_RIGHT_SIDE      9
#define FUNC_REV_QUERY_REF          10
#define FUNC_MINI_SW_LEFT_SIDE      11
#define FUNC_SEED_REGION_NEW        12
#define FUNC_SEED_LIST_INSERT       13
#define FUNC_CAL_NEW                14
#define FUNC_CAL_MNG_INSERT         15
#define FUNC_CAL_MNG_TO_LIST        16
#define FUNC_FILTER_BY_READ_AREA    17
#define FUNC_PRE_SW                 18
#define FUNC_SW                     19
#define FUNC_POST_SW                20
#define FUNC_OTHER                  21
#define FUNC_CREATE_ALIGNMENTS      22


#define NUM_TIMING (FUNC_CREATE_ALIGNMENTS + 1)
double func_times[NUM_TIMING];

const char *func_names[NUM_TIMING] = {"search_suffix", 
				      "search_prefix",
				      "search_sa",
				      "generate_cals_from_exact_read",
				      "generate_cals_from_suffixes",
				      "init_cals_from_suffixes",
				      "set_positions",
				      "set_reference_sequence",
				      "skip_suffixes",
				      "mini_sw_right_side",
				      "rev_query_ref",
				      "mini_sw_left_side",
				      "seed_region_new",
				      "seed_list_insert",
				      "cal_new",
				      "cal_mng_insert",
				      "cal_mng_to_array_list",
				      "filter_by_read_area",
				      "sw_pre_processing",
				      "sw_execution",
				      "sw_post_processing",
				      "fill_middle_gaps",
				      "merge_gaps",
				      "fill_end_gaps",
				      "other functions",
				      "create_alignments"};

#endif

typedef struct sa_mapping_batch {
  int sam_format;
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
  p->sam_format = 0;
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
// sq_wf_batch
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
// sa_wf_input
//--------------------------------------------------------------------

typedef struct sa_wf_input {
  int sam_format;
  fastq_batch_reader_input_t *fq_reader_input;
  sa_wf_batch_t *wf_batch;
} sa_wf_input_t;

//--------------------------------------------------------------------

inline sa_wf_input_t *sa_wf_input_new(int sam_format,
				      fastq_batch_reader_input_t *fq_reader_input,
				      sa_wf_batch_t *wf_batch) {
  sa_wf_input_t *p = (sa_wf_input_t *) malloc(sizeof(sa_wf_input_t));
  p->sam_format = sam_format;
  p->fq_reader_input = fq_reader_input;
  p->wf_batch = wf_batch;
  return p;
}

//--------------------------------------------------------------------

inline void sa_wf_input_free(sa_wf_input_t *p) {
  if (p) free(p);
}

//====================================================================
// FUNCTIONS
//====================================================================

void *sa_fq_reader(void *data);
int sa_mapper(void *data);
int sa_bam_writer(void *data);

//====================================================================
// PRODUCER
//====================================================================

//--------------------------------------------------------------------
// sa fq reader
//--------------------------------------------------------------------

size_t fastq_fread_se_ex(array_list_t *reads, size_t num_reads, fastq_file_t *fq_file) {
  size_t count = 0;
  char *p;
  char header1[MAX_READ_ID_LENGTH];
  char sequence[MAX_READ_SEQUENCE_LENGTH];
  char header2[MAX_READ_ID_LENGTH];
  char qualities[MAX_READ_SEQUENCE_LENGTH];
  int header_length, sequence_length, quality_length;
  fastq_read_t *read;
  
  while (count < num_reads && fgets(header1, MAX_READ_ID_LENGTH, fq_file->fd) != NULL) {
    fgets(sequence, MAX_READ_SEQUENCE_LENGTH, fq_file->fd);
    fgets(header2, MAX_READ_ID_LENGTH, fq_file->fd);
    fgets(qualities, MAX_READ_SEQUENCE_LENGTH, fq_file->fd);
    
    header_length = strlen(header1);
    sequence_length = strlen(sequence);
    quality_length = strlen(qualities);
    
    // '\n' char is removed, but '\0' is left
    chomp_at(header1, header_length - 1);
    if ((p = strstr(header1, " ")) != NULL) {
      *p = 0;
    }
    chomp_at(sequence, sequence_length - 1);
    chomp_at(qualities, quality_length - 1);

    read = fastq_read_new(&header1[1], sequence, qualities);
    array_list_insert(read, reads);
    
    count++;
  }
  
  return count;
}

//--------------------------------------------------------------------

void *sa_fq_reader(void *input) {
  sa_wf_input_t *wf_input = (sa_wf_input_t *) input;

  sa_wf_batch_t *new_wf_batch = NULL;
  sa_wf_batch_t *curr_wf_batch = wf_input->wf_batch;
  
  fastq_batch_reader_input_t *fq_reader_input = wf_input->fq_reader_input;
  array_list_t *reads = array_list_new(fq_reader_input->batch_size, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  
  fastq_fread_se_ex(reads, fq_reader_input->batch_size, fq_reader_input->fq_file1);
  
  size_t num_reads = array_list_size(reads);
  
  if (num_reads == 0) {
    array_list_free(reads, (void *)fastq_read_free);
  } else {
    sa_mapping_batch_t *sa_mapping_batch = sa_mapping_batch_new(reads);
    sa_mapping_batch->sam_format = wf_input->sam_format;


    new_wf_batch = sa_wf_batch_new(curr_wf_batch->sa_index,
				   curr_wf_batch->writer_input, 
				   sa_mapping_batch);
  }
  return new_wf_batch;
}

//====================================================================
// CONSUMER
//====================================================================

//--------------------------------------------------------------------
// sa bam writer
//--------------------------------------------------------------------

int sa_bam_writer(void *data) {
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;
  
  sa_mapping_batch_t *mapping_batch = (sa_mapping_batch_t *) wf_batch->mapping_batch;
  if (mapping_batch == NULL) {
    printf("bam_writer1: error, NULL mapping batch\n");
    exit(-1);
  }
  int sam = mapping_batch->sam_format;

  if (sam) {
    bam1_t *bam1;
    array_list_t *mapping_list;
    FILE *bam_file = (FILE *) wf_batch->writer_input->bam_file;
    
    size_t num_reads, num_mappings;
    num_reads = mapping_batch->num_reads;
    for (size_t i = 0; i < num_reads; i++) {
      mapping_list = mapping_batch->mapping_lists[i];
      num_mappings = array_list_size(mapping_list);
      for (size_t j = 0; j < num_mappings; j++) {
	bam1 = (bam1_t *) array_list_get(j, mapping_list);
	fprintf(bam_file, "%s\t%i\t%lu\t%s\t%s\t%s\n", 
		bam1_qname(bam1), 
		bam1->core.tid,
		bam1->core.pos,
		convert_to_cigar_string(bam1_cigar(bam1), bam1->core.n_cigar),
		convert_to_sequence_string(bam1_seq(bam1), bam1->core.l_qseq),
		convert_to_sequence_string(bam1_seq(bam1), bam1->core.l_qseq)
		);
	bam_destroy1(bam1);	 
      }
      array_list_free(mapping_list, (void *) NULL);
    }
  } else {
    bam1_t *bam1;
    array_list_t *mapping_list;
    bam_file_t *bam_file = wf_batch->writer_input->bam_file;
    
    size_t num_reads, num_mappings;
    num_reads = mapping_batch->num_reads;
    for (size_t i = 0; i < num_reads; i++) {
      mapping_list = mapping_batch->mapping_lists[i];
      num_mappings = array_list_size(mapping_list);
      for (size_t j = 0; j < num_mappings; j++) {
	bam1 = (bam1_t *) array_list_get(j, mapping_list);
	//	if (j == 0) {
	  bam_fwrite(bam1, bam_file);
	  //	}
	bam_destroy1(bam1);	 
      }
      array_list_free(mapping_list, (void *) NULL);
    }
  }

  // free memory
  sa_mapping_batch_free(mapping_batch);

  if (wf_batch) sa_wf_batch_free(wf_batch);

  return 0;
}

//--------------------------------------------------------------------
#ifdef _VERBOSE
int num_dup_reads = 0;
int num_total_dup_reads = 0;
#endif

int sa_sam_writer(void *data) {
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;
  
  sa_mapping_batch_t *mapping_batch = (sa_mapping_batch_t *) wf_batch->mapping_batch;
  if (mapping_batch == NULL) {
    printf("bam_writer1: error, NULL mapping batch\n");
    exit(-1);
  }

  #ifdef _TIMING
  for (int i = 0; i < NUM_TIMING; i++) {
    func_times[i] += mapping_batch->func_times[i];
  }
  #endif

  int flag, pnext = 0, tlen = 0;
  char *rnext = "*", *optional_flags = "NM:i:3";

  fastq_read_t *read;
  array_list_t *read_list = mapping_batch->fq_reads;

  alignment_t *alig;
  array_list_t *mapping_list;
  FILE *out_file = (FILE *) wf_batch->writer_input->bam_file;

  sa_genome3_t *genome = wf_batch->sa_index->genome;

  size_t num_reads, num_mappings;
  num_reads = mapping_batch->num_reads;
  for (size_t i = 0; i < num_reads; i++) {
    read = (fastq_read_t *) array_list_get(i, read_list);
    mapping_list = mapping_batch->mapping_lists[i];
    num_mappings = array_list_size(mapping_list);
    #ifdef _VERBOSE
    if (num_mappings > 1) {
      num_dup_reads++;
      num_total_dup_reads += num_mappings;
    }
    #endif

    if (num_mappings > 0) {
      for (size_t j = 0; j < num_mappings; j++) {
	alig = (alignment_t *) array_list_get(j, mapping_list);

	flag = (alig->seq_strand ? 16 : 0);

	fprintf(out_file, "%s\t%i\t%s\t%lu\t%i\t%s\t%s\t%lu\t%i\t%s\t%s\t%s\n", 
		alig->query_name,
		flag,
		genome->chrom_names[alig->chromosome],
		alig->position + 1,
		alig->map_quality,
		alig->cigar,
		rnext,
		pnext,
		tlen,
		alig->sequence,
		alig->quality,
		optional_flags
		);
	alignment_free(alig);	 
      }
    } else {
      fprintf(out_file, "%s\t4\t*\t0\t0\t*\t*\t0\t0\t%s\t%s\n", 
	      read->id,
	      read->sequence,
	      read->quality
	      );
    }
    array_list_free(mapping_list, (void *) NULL);
  }

  // free memory
  sa_mapping_batch_free(mapping_batch);

  if (wf_batch) sa_wf_batch_free(wf_batch);

  return 0;
}

//====================================================================
// SA MAPPER
//====================================================================

const uint max_uint = 4294967295;

//--------------------------------------------------------------------
// structures for sa mapper
//--------------------------------------------------------------------
//const int MAX_NUM_MISMATCHES = 4;

typedef struct mapped_reg {
  int num_matches;
  int num_mismatches;
  int strand;
  int chromosome;
  int r_start;
  int r_end;
  size_t g_start;
  size_t g_end;
  int mismatches_pos[MAX_NUM_MISMATCHES];
} mapped_reg_t;

mapped_reg_t *mapped_reg_new(int r_start, int r_end, int strand, int chromosome, 
			     size_t g_start, size_t g_end, int num_mismatches) {
  mapped_reg_t *p = (mapped_reg_t *) malloc(sizeof(mapped_reg_t));
  p->r_start = r_start;
  p->r_end = r_end;
  p->strand = strand;
  p->chromosome = chromosome;
  p->g_start = g_start;
  p->g_end = g_end;
  p->num_mismatches = num_mismatches;
  p->num_matches = r_end - r_start - num_mismatches + 1;
  return  p;
}

void mapped_reg_free(mapped_reg_t *p) {
  if (p) free(p);
}

//--------------------------------------------------------------------
// utils
//--------------------------------------------------------------------

void display_sequence(uint j, sa_index3_t *index, uint len) {


  char *p = &index->genome->S[index->SA[j]];
  char chrom = index->CHROM[j];
  for (int i = 0; i < len; i++) {
    printf("%c", *p);
    p++;
  }
  printf("\t%u\t%s:%u\n", index->SA[j], 
	 index->genome->chrom_names[chrom], index->SA[j] - index->genome->chrom_offsets[chrom]);
}

//--------------------------------------------------------------------
// binary searches
//--------------------------------------------------------------------

uint binary_search_first(char* key, uint len, uint ulow, uint uhigh, 
			 uint skip, uint *matched, sa_index3_t *index) {
  // Finds  first string >= key assisted by binary search on suffix array,
  // but without using lcp.
  uint count = 0;
  uint i, j;
  
  uint max_matched = 0;
  uint best_low = max_uint;

  uint n = index->num_suffixes;
  char *S = index->genome->S;
  uint *SA = index->SA;

  //  uint mid;
  unsigned long mid, low = (unsigned long) ulow, high = (unsigned long) uhigh;

  while (low <= high) {
    mid = (low + high) >> 1;
    //    mid = ((unsigned long) low + (unsigned long) high) / 2;
    //    mid = (low >> 1 ) + (high >> 1);
    //mid = (low + high) >> 1; // mid = (low + high) / 2;
    //printf("(low, high) = (%u, %u) -> mid = %u\n", low, high, mid);
    j = SA[mid] + skip; // Position in s
    i = skip;       // Position in key
    printf("\tfirst: (low, high) = (%u, %u) -> mid = %u\n", low, high, mid);
    // Like strcmp
    while (S[j] == key[i] && key[i]) {
      j++;
      i++;
    } 
    if (i >= max_matched && mid < best_low) {
      max_matched = i;
      best_low = mid;
    }
    printf("\t\t");
    display_sequence(mid, index, len);
    if (key[i] == 0 || S[j] > key[i]) {
      high = mid - 1;
    } else {
      low = mid + 1;
    }
  }
  
  //  printf("\t>>> binary_search_first: low = %u (matching %u of %u)\n\n", low, i, len);
  //  if (i < max_matched) {
    printf("\t\t>>> retrive max. matched. best low = %u (matching %u of %u)\n\n", best_low, max_matched, len);
    *matched = max_matched;
    return (uint) best_low;
    //  } else {
    //    *matched = i;
    //    return (uint) low;
    //  }
}

//--------------------------------------------------------------------

uint binary_search_last(char* key, uint len, uint ulow, uint uhigh, 
			uint skip, uint matched, sa_index3_t *index) {
  // Finds last string <= key assisted by binary search on suffix array,
  // but without using lcp.
  uint i, j;
  
  uint n = index->num_suffixes;
  char *S = index->genome->S;
  uint *SA = index->SA;

  //  uint mid;
  unsigned long mid, low = (unsigned long) ulow, high = (unsigned long) uhigh;
  
  while (low <= high) {
    mid = (low + high) >> 1;
    //    mid = ((unsigned long) low + (unsigned long) high) / 2;
    //    mid = (low >> 1 ) + (high >> 1);
    //mid = (low + high) / 2;
    j = SA[mid] + skip; // Position in s
    i = skip;       // Position in key
    printf("\tlast: (low, high) = (%u, %u) -> mid = %u\n", low, high, mid);

    // Like strcmp
    while (S[j] == key[i] && key[i]) {
      j++; 
      i++; 
    } 
    printf("\t\t");
    display_sequence(mid, index, len);
    if (key[i] == 0 || S[j] < key[i]) {
      low = mid + 1;
    } else {
      high = mid - 1;
    }
    //    if (i == matched) { break; }
  }
  
  printf("\t>>> binary_search_last: high = %u (matching %u of %u)\n\n", high, i, len);
  return high;
}

//--------------------------------------------------------------------

size_t search_prefix(char *sequence, size_t *low, size_t *high, sa_index3_t *sa_index, int display) {
  size_t num_mappings = 0;

  char *seq;
  size_t value, row, col;
  uint ia, ia1, ia2, found_ia;
  uint ja, found_ja;
  uint a, a1, a2;

  //  printf("prefix: %s\n", sequence);
  //  display_prefix(sequence, sa_index->k_value);
  value = compute_prefix_value(sequence, sa_index->k_value);
  row = value >> 8;
  col = 255LLU & value;
  //  printf(" -> prefix value = %lu -> (row, col) = (%lu, %lu)\n", value, row, col); 
  
  ia1 = sa_index->IA[row];
  //  printf("\tIA[%lu] = %lu\n", row, ia1);
  if (ia1 == max_uint) {
    //    printf("\t\t\t----> ROW NOT FOUND !!!\n");
    return num_mappings;
  }
  
  size_t row2 = row + 1;
  while ((ia2 = sa_index->IA[row2]) == max_uint) {
    row2++;
    if (row2 >= sa_index->IA_items) {
      //      printf("\trow2 reaches IA_items limit (IA items = %lu)\n", sa_index->IA_items);
      row2 = sa_index->IA_items;
      ia2 = sa_index->A_items;
      break;
    }
  }
  //  printf("\tfrom IA[%lu] = %lu to IA[%lu] = %lu -> num. columns = %lu\n", row, ia1, row2, ia2, ia2 - ia1);
  
  found_ja = 0;

  int idx, size = ia2 - ia1;
  for (ia = ia1; ia < ia2; ia++) {
    a1 = sa_index->A[ia];
    ja = sa_index->JA[ia];
    //      printf("\t\tA[%lu] = %lu\t JA[%lu] = %lu\n", 
    //      	     ia, a1, ia, ja);
    if (sa_index->JA[ia] == col) {
      found_ja = 1;
      if (ia + 1 >= sa_index->A_items) {
	a2 = sa_index->num_suffixes;
      } else {
	a2 = sa_index->A[ia + 1];
      }
      break;
    }
  }

  if (found_ja) {
    num_mappings = a2 - a1;
    *low = a1;
    *high = a2;
  }

  return num_mappings;
}

//--------------------------------------------------------------------

size_t search_suffix(char *seq, uint len, sa_index3_t *sa_index, 
		     size_t *low, size_t *high, size_t *suffix_len
                     #ifdef _TIMING
		     , sa_mapping_batch_t *mapping_batch
                     #endif
		     ) {

  #ifdef _TIMING
  struct timeval stop, start;
  #endif

  int display = 1;
  char *ref, *query;
  size_t num_suffixes = 0;
  uint matched, max_matched = 0;

  *suffix_len = 0;

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif
  size_t num_prefixes = search_prefix(seq, low, high, sa_index, display);
  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_SEARCH_PREFIX] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif

  #ifdef _VERBOSE	  
  printf("\t\tnum. prefixes = %lu\n", num_prefixes);
  #endif


  if (num_prefixes && num_prefixes < MAX_NUM_SUFFIXES) {   

    //    *suffix_len = sa_index->k_value;
    //    return num_prefixes;
    
    #ifdef _TIMING
    gettimeofday(&start, NULL);
    #endif


    #ifdef _VERBOSE1	  
    {
      char *ss = get_subsequence(seq + sa_index->k_value, 0, 40);
      printf("\tquery:\t%s\n", ss);
      free(ss);
      for (size_t i = *low; i < *high; i++) {
	printf("\t%lu\t", i);
	ref = &sa_index->genome->S[sa_index->SA[i]] + sa_index->k_value;
	char *ss = get_subsequence(ref, 0, 40);
	printf("%s\n", ss);
	free(ss);
      }
    }
    #endif


    size_t first = *low, last = *low;

    if (num_prefixes == 1) {
      query = seq + sa_index->k_value;
      ref = &sa_index->genome->S[sa_index->SA[*low]] + sa_index->k_value;
      matched = 0;
      while (query[matched] == ref[matched]) {
	matched++;
      }
      *high = *low;
      *suffix_len = matched + sa_index->k_value;
      return num_prefixes;
    }

    /*
    size_t mid, i, l, r;
    // first binary search -> first
    #ifdef _VERBOSE	  
    printf("low = %lu, high = %lu\n", *low, *high);
    #endif

    l = *low;
    r = *high;
    while (l <= r) {
      mid = (l + r) / 2;
      query = seq + sa_index->k_value;
      ref = &sa_index->genome->S[sa_index->SA[mid]] + sa_index->k_value;
      i = 0;
      #ifdef _VERBOSE	  
      printf("\t(l, r) = (%lu, %lu) -> mid = %lu\n", l, r, mid);
      #endif
      while (query[i] && query[i] == ref[i]) {
	i++;
      }
      if (i > max_matched) max_matched = i;
      if (query[i] == 0 || ref[i] > query[i]) {
	r = mid - 1;
        #ifdef _VERBOSE	  
	printf("\t\tref > query -> r = mid - 1 = %lu (matched = %i)\n", r, i);
	char *ss = get_subsequence(query, 0, i + 1);
	printf("\t\tquery: %s\n", ss);
	free(ss);
	ss = get_subsequence(ref, 0, i + 1);
	printf("\t\tref. : %s\n", ss);
	free(ss);
	#endif
      } else {
	l = mid + 1;
        #ifdef _VERBOSE	  
	printf("\t\tref <= query -> l = mid + 1 = %lu (matched = %i)\n", l, i);
	char *ss = get_subsequence(query, 0, i + 1);
	printf("\t\tquery: %s\n", ss);
	free(ss);
	ss = get_subsequence(ref, 0, i + 1);
	printf("\t\tref. : %s\n", ss);
	free(ss);
	#endif
      }
    }
    if (max_matched) {
      first = l;
      #ifdef _VERBOSE
      printf("\t\t\t---> first = %lu (max. matched = %i)\n", first, max_matched);
      #endif
      
      r = *high;
      while (l <= r) {
	mid = (l + r) / 2;
	query = seq + sa_index->k_value;
	ref = &sa_index->genome->S[sa_index->SA[mid]] + sa_index->k_value;

//	  printf("(l, r) = (%lu, %lu) -> mid = %lu\n", l, r, mid);
//	  char *ss = get_subsequence(query, 0, max_matched);
//	  printf("\tquery: %s\n", ss);
//	  free(ss);
//	  ss = get_subsequence(ref, 0, max_matched);
//	  printf("\tref. : %s\n", ss);
//	  free(ss);

	if (strncmp(ref, query, max_matched) <= 0) {
	  l = mid + 1;
	  //printf("\t\t ref <= query -> l = mid + 1 = %lu\n", l);
	} else {
	  r = mid - 1;
	  //printf("\t\t ref > query -> r = mid - 1 = %lu\n", r);
	}
      }
      
      last = r;
      #ifdef _VERBOSE
      printf("\t\t\t---> last = %lu\n", last);
      #endif
    } else {
      first = *low;
      last = first + num_prefixes - 1;
    }

    first = *low;
    last = *low;
    max_matched = 0;
*/
    for (size_t i = *low; i < *high; i++) {
      query = seq + sa_index->k_value;
      ref = &sa_index->genome->S[sa_index->SA[i]] + sa_index->k_value;
      matched = 0;
      while (query[matched] == ref[matched]) {
	matched++;
      }
      if (matched > max_matched) {
	first = i;
	last = i;
	max_matched = matched;
	//	break;
      } else if (matched == max_matched) {
	last = i;
      } else {
	break;
      }
    }
    /*
    printf("----> sequential method: (low, high) = (%lu, %lu) : (first, last) = (%lu, %lu) -> num. suffixes = %lu (max. matched = %i)\n", 
	   *low, *high, first, last, last - first + 1, max_matched);
    exit(-1);
    */

    if (first <= last) {
      *low = first;
      *high = last;
      *suffix_len = max_matched + sa_index->k_value;
      num_suffixes = last - first + 1;
    }

    #ifdef _TIMING
    gettimeofday(&stop, NULL);
    mapping_batch->func_times[FUNC_SEARCH_SA] += 
      ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
    #endif
  }
  
  return num_suffixes;
}

//--------------------------------------------------------------------
/*
alignment_t *create_alignment(fastq_read_t *read, char *revcomp, mapped_reg_t *region, int num_items) {

  int distance = region->num_mismatches;
  int AS = (region->num_matches * 5) - (distance * 4);

  size_t map_len, num_cigar_ops = 1;

  // sequence and quality
  map_len = region->r_end - region->r_start + 1;
  char *sequence = (char *) malloc(map_len + 1);
  memcpy(sequence, &read->sequence[region->r_start], map_len);
  sequence[map_len] = 0;
  
  char *quality = (char *) malloc(map_len + 1);
  memcpy(quality, &read->quality[region->r_start], map_len);
  quality[map_len] = 0;  

  // cigar
  char *cigar = (char *) calloc(100, sizeof(char));
  if (region->r_start > 0) {
    num_cigar_ops++;
    sprintf(cigar, "%s%iH", cigar, region->r_start);
  } 
  sprintf(cigar, "%s%iM", cigar, map_len);
  if (region->r_end < read->length - 1) {
    num_cigar_ops++;
    sprintf(cigar, "%s%iH", cigar, read->length - region->r_end - 1);
  } 
  
  //  printf("seq  : %s\n", sequence);
  //  printf("qual : %s\n", quality);
  //  printf("cigar: %s\n", cigar);
  //  printf("(strand, chrom, pos) = (%lu, %lu, %lu)\n", region->strand, region->chromosome, region->g_start);

  // set optional fields
  int optional_fields_length = 100;
  char *optional_fields = (char *) calloc(optional_fields_length, sizeof(char));
  
  char *p = optional_fields;
  
  sprintf(p, "ASi");
  p += 3;
  memcpy(p, &AS, sizeof(int));
  p += sizeof(int);
  
  sprintf(p, "NHi");
  p += 3;
  memcpy(p, &num_items, sizeof(int));
  p += sizeof(int);
  
  sprintf(p, "NMi");
  p += 3;
  memcpy(p, &distance, sizeof(int));
  p += sizeof(int);
  //      *p = '\0';
  optional_fields_length = p - optional_fields;
  
  // create the alignment
  alignment_t *alignment = alignment_new();
  alignment_init_single_end(strdup(read->id), sequence, quality, 
			    region->strand, 
			    region->chromosome, 
			    region->g_start - 1,
			    cigar, num_cigar_ops, (AS * 254) / (read->length * 5), 1, (num_items > 1),
			    optional_fields_length, optional_fields, 0, alignment);  
  return alignment;
}

//--------------------------------------------------------------------

size_t extend_seeds(char *r_seq, uint read_length, int strand, uint seed_pos, 
		    size_t low, size_t high, sa_index3_t *sa_index, array_list_t *region_list) {
  
  int chrom, min_num_matches, num_matches;
  int num_left_mismatches, left_mismatches_pos[MAX_NUM_MISMATCHES];
  int num_right_mismatches, right_mismatches_pos[MAX_NUM_MISMATCHES];
  
  size_t r_left, r_right, g_left, g_right; // read and genomic positions

  int found;
  mapped_reg_t *reg;
  size_t num_regions;

  char *S = sa_index->genome->S;
  min_num_matches = read_length - (read_length * 0.25);
    
  // for each prefix
  for (size_t a = low; a < high; a++) {
    
    chrom = sa_index->CHROM[a];
    
    // extend to the left side
    num_left_mismatches = 0;
    r_left = seed_pos;
    g_left = sa_index->SA[a];
    if (seed_pos > 0) {
      for (int k = seed_pos - 1; k >= 0; k--) {
	r_left--;
	g_left--;
	//	   printf("%i (%c, %c)\n", k, r_seq[r_left], S[g_left]);
	if (S[g_left] != r_seq[r_left]) {
	  if (num_left_mismatches + 1 == MAX_NUM_MISMATCHES) {
	    r_left++;
	    g_left++;
	    break;
	  } else {
	    left_mismatches_pos[num_left_mismatches++] = r_left;
	    //	       printf("\t\t\tmismatch at %lu\n", r_left);
	  }
	}
      }
      //      printf("\t\t\t\textended to left: from %lu to %lu (total mismatches %lu)\n", 
      //	     seed_pos - 1, r_left, num_left_mismatches);
    }

    // extend to the right side
    num_right_mismatches = 0;
    r_right = seed_pos + sa_index->k_value - 1;
    g_right = sa_index->SA[a] + sa_index->k_value - 1;
    for (int k = r_right + 1; k < read_length; k++) {
      r_right++;
      g_right++;
      //	 printf("%i (%c, %c)\n", k, r_seq[r_right], S[g_right]);
      if (S[g_right] != r_seq[r_right]) {
	if (num_right_mismatches + 1 >= MAX_NUM_MISMATCHES) {
	  r_right--;
	  g_right--;
	  break;
	} else {
	  right_mismatches_pos[num_right_mismatches++] = r_right;
	  //	     printf("\t\t\tmismatch at %lu\n", r_right);
	}
      }
    }
    //    printf("\t\t\t\textended to right: from %lu to %lu (total mismatches %lu)\n", 
    //	   seed_pos + sa_index->k_value, r_right, num_right_mismatches);
    
    num_matches = r_right - r_left - num_left_mismatches - num_right_mismatches + 1;
    //    printf("\t\t\t\t\t%i-%i (matches: %lu (min. %lu), mismatches: %lu, at",
    //	   r_left, r_right, num_matches, min_num_matches, num_left_mismatches + num_right_mismatches);
    //    printf(")\n");
    
    // insert in the list the best mapped region (max. num. matches)
    if ((num_regions = array_list_size(region_list)) > 0) {
      reg = (mapped_reg_t *) array_list_get(0, region_list);
      if (num_matches > reg->num_matches) {
	array_list_clear(region_list, (void *) mapped_reg_free);
	mapped_reg_t *mapped_reg = mapped_reg_new(r_left, r_right, strand, chrom,
						  g_left - sa_index->genome->chrom_offsets[chrom] + 1,
						  g_right - sa_index->genome->chrom_offsets[chrom] + 1,
						  num_left_mismatches + num_right_mismatches);
	array_list_insert(mapped_reg, region_list);
      } else if (num_matches == reg->num_matches) {
	g_left = g_left - sa_index->genome->chrom_offsets[chrom] + 1;
	g_right = g_right - sa_index->genome->chrom_offsets[chrom] + 1;
	for (size_t r = 0; r < num_regions; r++) {
	  reg = (mapped_reg_t *) array_list_get(r, region_list);
	  if (reg->g_start == g_left && reg->chromosome == chrom) {
	    found = 1;
	    break;
	  }
	}
	if (!found) {
	  mapped_reg_t *mapped_reg = mapped_reg_new(r_left, r_right, strand, chrom,
						    g_left, g_right,
						    num_left_mismatches + num_right_mismatches);
	  array_list_insert(mapped_reg, region_list);
	}
      }
    } else {
      mapped_reg_t *mapped_reg = mapped_reg_new(r_left, r_right, strand, chrom,
						g_left - sa_index->genome->chrom_offsets[chrom] + 1,
						g_right - sa_index->genome->chrom_offsets[chrom] + 1,
						num_left_mismatches + num_right_mismatches);
      array_list_insert(mapped_reg, region_list);
    }
  }

  return array_list_size(region_list);
}
*/
//--------------------------------------------------------------------

void print_seed_region(char *msg, seed_region_t *s) {
  printf("%s%c:%i[%lu|%lu - %lu|%lu] (cigar: %s, num. mismatches = %i)\n",  msg, (s->strand == 0 ? '+' : '-'),
	 s->chromosome_id, s->genome_start, s->read_start, s->read_end, s->genome_end,
	 cigar_to_string(s->info), s->num_mismatches);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

typedef struct cal_mng {
  int min_read_area;
  int max_read_area;
  int read_length;
  int num_chroms;
  size_t low_prefix[MAX_NUM_SUFFIXES];
  size_t high_prefix[MAX_NUM_SUFFIXES];
  size_t low_suffix[MAX_NUM_SUFFIXES];
  size_t high_suffix[MAX_NUM_SUFFIXES];
  linked_list_t **cals_lists;
} cal_mng_t;

//--------------------------------------------------------------------

cal_mng_t * cal_mng_new(sa_genome3_t *genome) {

  int num_chroms = genome->num_chroms;

  linked_list_t **cals_lists = (linked_list_t **) malloc (sizeof(linked_list_t *) * num_chroms);
  for (unsigned int i = 0; i < num_chroms; i++) {
    cals_lists[i] = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
  }

  cal_mng_t *p = (cal_mng_t *) calloc(1, sizeof(cal_mng_t));
  p->read_length = 10;
  p->min_read_area = 100;
  p->max_read_area = 0;
  p->num_chroms = num_chroms;
  p->cals_lists = cals_lists;
  return p;
}

//--------------------------------------------------------------------

void cal_mng_free(cal_mng_t *p) {
  if (p) {
    if (p->cals_lists) {
      for (unsigned int i = 0; i < p->num_chroms; i++) {
	if (p->cals_lists[i]) {
	  linked_list_free(p->cals_lists[i], cal_free);
	}
      }
      free(p->cals_lists);
    }
    free(p);
  }
}

//--------------------------------------------------------------------

void cal_mng_clear(cal_mng_t *p) {
  if (p) {
    if (p->cals_lists) {
      for (unsigned int i = 0; i < p->num_chroms; i++) {
	if (p->cals_lists[i]) {
	  linked_list_clear(p->cals_lists[i], cal_free);
	}
      }
    }
  }
}

//--------------------------------------------------------------------

void cal_mng_update(seed_region_t *seed, cal_mng_t *p) {
  int is_used = 0;
  if (p->cals_lists) {
    cal_t *cal;
    seed_region_t *s_last;
    linked_list_t *sr_list;
    linked_list_t *cal_list = p->cals_lists[seed->chromosome_id];
    if (cal_list) {
      #ifdef _VERBOSE
      printf("\t\t\tinsert this seed to the CAL manager:\n");
      print_seed_region("\t\t\t", seed);
      #endif

      //      if (cal->read_area < p->min_read_area) {
	//	printf("****************** set read min. area from %i to %i\n", p->min_read_area, cal->read_area);
      //	p->min_read_area = cal->read_area;
      //      }
      if (linked_list_size(cal_list) <= 0) {
	// create CAL and insert it into the CAL manager
	is_used = 1;
	sr_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	linked_list_insert(seed, sr_list);
	cal = cal_new(seed->chromosome_id, seed->strand, seed->genome_start, seed->genome_end, 1, sr_list, NULL);
	cal->read_area = seed->read_end - seed->read_start + 1;
	cal->num_mismatches = seed->num_mismatches;
	linked_list_insert(cal, cal_list);
      } else {
	// insert (by order)
	linked_list_iterator_t* itr = linked_list_iterator_new(cal_list);
	cal_t *item = (cal_t *) linked_list_iterator_curr(itr);
	while (item != NULL) {
	  #ifdef _VERBOSE1
	  printf("---> merging with this CAL?\n");
	  cal_print(item);
	  #endif
	  //	  assert(cal->end > item->start);
	  s_last = linked_list_get_last(item->sr_list);
	  if (seed->read_start - s_last->read_end < 100 && 
	      seed->genome_start - s_last->genome_end < 100) {
	    is_used = 1;
	    linked_list_insert_last(seed, item->sr_list);
	    item->end = seed->genome_end;
	    item->read_area += (seed->read_end - seed->read_start + 1);
	    item->num_mismatches += seed->num_mismatches;

	    //	    cal_free(cal);
            #ifdef _VERBOSE
	    printf("---> yes, merging CAL, result:\n");
	    cal_print(item);
	    #endif
	    break;
	  } else {
	    if (seed->genome_end < item->start) {
	      // create CAL and insert it into the CAL manager
	      is_used = 1;
	      sr_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	      linked_list_insert(seed, sr_list);
	      cal = cal_new(seed->chromosome_id, seed->strand, seed->genome_start, seed->genome_end, 
			    1, sr_list, NULL);
	      cal->read_area = seed->read_end - seed->read_start + 1;
	      cal->num_mismatches = seed->num_mismatches;
	      linked_list_iterator_insert(cal, itr);
	      linked_list_iterator_prev(itr);
	      break;
	    }
	  }
	  //continue loop...
	  linked_list_iterator_next(itr);
	  item = linked_list_iterator_curr(itr);
	}
	if (item == NULL) {
	  // create CAL and insert it into the CAL manager
	  is_used = 1;
	  sr_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	  linked_list_insert(seed, sr_list);
	  cal = cal_new(seed->chromosome_id, seed->strand, seed->genome_start, seed->genome_end, 1, sr_list, NULL);
	  cal->read_area = seed->read_end - seed->read_start + 1;
	  cal->num_mismatches = seed->num_mismatches;
	  linked_list_insert_last(cal, cal_list);
	}
	linked_list_iterator_free(itr);
      }
    }
  }
  if (!is_used) {
    //    seed_region_free(seed);
  }
}

//--------------------------------------------------------------------

int cal_mng_find(int chrom, size_t start, size_t end, cal_mng_t *p) {
  #ifdef _VERBOSE1
  printf("\t\t***** searching CAL: chrom %i: %lu-%lu\n", chrom, start, end);
  #endif
  int found_cal = 0;
  if (p->cals_lists) {
    linked_list_t *cal_list = p->cals_lists[chrom];
    if (cal_list) {
      cal_t *cal;
      for (linked_list_item_t *item = cal_list->first; 
	   item != NULL; 
	   item = item->next) {
	cal = item->item;
	#ifdef _VERBOSE1
	printf("\t\t\t***** searching CAL: comparing %lu-%lu to cal item %lu-%lu\n", 
	       start, end, cal->start, cal->end);
	#endif
	if (cal->start <= start && cal->end >= end) {
	  found_cal = 1;
	  break;
	} 
	if (cal->start > end) {
	  break;
	}
      }
    }
  }
  #ifdef _VERBOSE1
  printf("\t\t\t\t***** searching CAL: found_cal = %i\n", found_cal);
  #endif
  return found_cal;
}

//--------------------------------------------------------------------

void cal_mng_to_array_list(int read_area, array_list_t *out_list, cal_mng_t *p) {
  cal_t *cal;
  linked_list_iterator_t itr;

  if (p->cals_lists) {
    linked_list_t *cal_list;
    for (unsigned int i = 0; i < p->num_chroms; i++) {
      cal_list = p->cals_lists[i];
      while (cal = (cal_t *) linked_list_remove_last(cal_list)) {
	if ((cal->end - cal->start) >= read_area) {
	  array_list_insert(cal, out_list);
	} else {
	  cal_free(cal);
	}
	/*
	if (p->min_read_area <= read_area && cal->read_area <= read_area) {
	  array_list_insert(cal, out_list);
	} else {
	  cal_free(cal);
	}
	*/
      }
    }
  }
}

//--------------------------------------------------------------------

void cal_mng_select_best(int read_area, array_list_t *valid_list, array_list_t *invalid_list, 
			 cal_mng_t *p) {
  cal_t *cal;
  linked_list_iterator_t itr;

  if (p->cals_lists) {
    linked_list_t *cal_list;
    for (unsigned int i = 0; i < p->num_chroms; i++) {
      cal_list = p->cals_lists[i];
      while (cal = (cal_t *) linked_list_remove_last(cal_list)) {
	if (p->min_read_area <= read_area && cal->read_area <= read_area) {
	  array_list_insert(cal, valid_list);
	} else {
	  array_list_insert(cal, invalid_list);
	}
      }
    }
  }
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

typedef struct seed_mng {
  int num_strands;
  int num_chroms;
  int min_cal_size;
  linked_list_t ***cals_lists;
} seed_mng_t;

//--------------------------------------------------------------------

seed_mng_t * seed_mng_new(sa_genome3_t *genome, int min_cal_size) {

  int num_strands = 2;
  int num_chroms = genome->num_chroms;

  linked_list_t ***cals_lists = (linked_list_t ***) malloc(sizeof(linked_list_t **) * num_strands);

  for (unsigned int i = 0; i < num_strands; i++) {
    cals_lists[i] = (linked_list_t **) malloc (sizeof(linked_list_t *) * num_chroms);
    for (unsigned int j = 0; j < num_chroms; j++) {
      cals_lists[i][j] = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    }
  }

  seed_mng_t *p = (seed_mng_t *) calloc(1, sizeof(seed_mng_t));
  p->num_strands = num_strands;
  p->num_chroms = num_chroms;
  p->min_cal_size = min_cal_size;
  p->cals_lists = cals_lists;
  return p;
}

//--------------------------------------------------------------------

void seed_mng_free(seed_mng_t *p) {
  if (p) {
    if (p->cals_lists) {
      for (unsigned int i = 0; i < p->num_strands; i++) {
	if (p->cals_lists[i]) {
	  for (unsigned int j = 0; j < p->num_chroms; j++) {
	    if (p->cals_lists[i][j]) {
	      linked_list_free(p->cals_lists[i][j], short_cal_free);
	    }
	  }
	  free(p->cals_lists[i]);
	}
      }
      free(p->cals_lists);
    }
    free(p);
  }
}

//--------------------------------------------------------------------

void seed_mng_clear(seed_mng_t *p) {
  if (p) {
    if (p->cals_lists) {
      for (unsigned int i = 0; i < p->num_strands; i++) {
	if (p->cals_lists[i]) {
	  for (unsigned int j = 0; j < p->num_chroms; j++) {
	    if (p->cals_lists[i][j]) {
	      linked_list_clear(p->cals_lists[i][j], short_cal_free);
	    }
	  }
	}
      }
    }
  }
}

//--------------------------------------------------------------------

void seed_mng_update(int num_suffixes, size_t low, size_t high,
		     int suffix_len, size_t read_pos, size_t read_length,
		     int strand, int id, sa_index3_t *sa_index, 
		     seed_mng_t *seed_mng) {

  size_t r_start, r_end, g_start, g_end;
  int chrom;
  region_t *region;

  // for each suffix
  //printf("(low, high) = (%lu, %lu)\n", low, high);
  for (size_t suff = low; suff <= high; suff++) {
    r_start = read_pos;
    r_end = r_start + suffix_len - 1;
    chrom = sa_index->CHROM[suff];
    g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];
    g_end = g_start + suffix_len - 1;
    
    region = region_bwt_new(chrom, strand, g_start, g_end, r_start, r_end, suffix_len, id);

    //printf("Insert (%c)[%i:%lu|%i-%i|%lu]\n", strand == 0 ? '+' : ' -', chrom, 
    //	   g_start, r_start, r_end, g_end);

    my_cp_list_append_linked_list(seed_mng->cals_lists[strand][chrom], region, read_length, 1000);
    region_bwt_free(region);
  }
}

//--------------------------------------------------------------------

void seed_mng_generate_cals(array_list_t *cal_list, sa_index3_t *sa_index, seed_mng_t *seed_mng) {
  seed_region_t *seed_region, *s;
  linked_list_iterator_t itr;
  
  cal_t *cal;
  short_cal_t *short_cal;
  linked_list_item_t *list_item_cal;
  
  uint num_strands = seed_mng->num_strands;
  uint num_chroms = seed_mng->num_chroms;
  uint min_cal_size = seed_mng->min_cal_size;
  
  // create CAL list from CAL manager structure   
  for (uint j = 0; j < num_chroms; j++) {
    for (uint i = 0; i < num_strands; i++) {
      linked_list_iterator_init(seed_mng->cals_lists[i][j], &itr);
      list_item_cal = linked_list_iterator_list_item_curr(&itr);
      while (list_item_cal != NULL) {
	short_cal = (short_cal_t *) list_item_cal->item;
	if (short_cal->end - short_cal->start + 1 >= min_cal_size) {
	  linked_list_t *list_aux = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	  while (s = (seed_region_t *)linked_list_remove_last(short_cal->sr_list)) {
	    //TODO: Change all parameters to seed_region_t
	    append_seed_region_linked_list(list_aux,
					   s->read_start, s->read_end, 
					   s->genome_start, s->genome_end, 
					   s->id);	    
	    seed_region_free(s);	    
	  }
	  cal = cal_new(j, i, short_cal->start, 
			short_cal->end, short_cal->num_seeds, 
			list_aux, short_cal->sr_duplicate_list);
	  array_list_insert(cal, cal_list);
	  // for debugging purporses
	  //cal_print(cal);
	  
	  //	  short_cal->sr_duplicate_list = NULL;
	  short_cal->sr_duplicate_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
	}
	
	linked_list_iterator_next(&itr);
	list_item_cal = linked_list_iterator_list_item_curr(&itr);
      }
    }
  }
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

void filter_valid_cals(array_list_t **list) {
  array_list_t *cal_list = *list;
  uint num_cals = array_list_size(cal_list);

  uint min_seeds, max_seeds;
  cal_t *cal;

  int founds[num_cals], found = 0;
  for (size_t j = 0; j < num_cals; j++) {
    founds[j] = 0;
    cal = array_list_get(j, cal_list);
    //    printf("\tcal %i of %i: sr_list size = %i (cal->num_seeds = %i) %i:%lu-%lu\n", 
    //	   j, num_cals, cal->sr_list->size, cal->num_seeds,
    //	   cal->chromosome_id, cal->start, cal->end);
    if (cal->sr_list->size > 0) {
      int start = 0;
      size_t genome_start = 0;
      int first = 1;
      for (linked_list_item_t *list_item = cal->sr_list->first; list_item != NULL; list_item = list_item->next) {
	seed_region_t *s = list_item->item;
	
	//	printf("\t\t:: star %lu > %lu s->read_start\n", start, s->read_start);
	if (start > s->read_start || s->read_start >= s->read_end) {
	  //	  printf("\t\t\t:: remove\n");
	  found++;
	  founds[j] = 1;
	}
	
	first = 0;
	start = s->read_end + 1;
	genome_start = s->genome_end + 1;
      }
    } else {
      // CAL without seeds is marked to be removed
      found++;
      founds[j] = 1;
    }
  }

  if (found) {
    min_seeds = 100000;
    max_seeds = 0;
    array_list_t *new_cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    for (size_t j = 0; j < num_cals; j++) {
      if (!founds[j]) {
	cal = array_list_get(j, list);
	cal->num_seeds = cal->sr_list->size;
	if (cal->num_seeds > max_seeds) max_seeds = cal->num_seeds;
	if (cal->num_seeds < min_seeds) min_seeds = cal->num_seeds;
	array_list_insert(cal, new_cal_list);
	array_list_set(j, NULL, cal_list);
      }
    }
    array_list_free(cal_list, (void *) cal_free);
    *list = new_cal_list;
  }
}

//--------------------------------------------------------------------

void filter_best_cals(int min_perc, fastq_read_t *read, char *revcomp_seq,
		      sa_index3_t *sa_index, array_list_t **list) {
  array_list_t *cal_list = *list;
  uint num_mismatches, num_regions, num_cals = array_list_size(cal_list);

  char *ref = sa_index->genome->S;
  char *query;

  cal_t *cal;
  seed_region_t *s, *prev_s, *curr_s;

  int invalid[num_cals], num_invalid = 0, matches;
  int read_pos, chrom;

  int max_num_mismatches = read->length - read->length * min_perc / 100;
  int min_num_mismatches = 100000;

  size_t g_pos;

  for (size_t j = 0; j < num_cals; j++) {
    matches = 0;
    invalid[j] = 0;
    cal = array_list_get(j, cal_list);
    cal->read_area = 0;

    #ifdef _VERBOSE
    //    printf("filter_best, cal:\n");
    cal_print(cal);
    #endif

    query = (cal->strand ? revcomp_seq : read->sequence);
    chrom = cal->chromosome_id;

    num_regions = cal->sr_list->size;
    if (num_regions == 0) {
      // must we check the sr_duplicated_list ??
      num_invalid++;
      invalid[j] = 1;
    } else {
      num_mismatches = 0;

      // intermediated gaps
      linked_list_item_t *prev_item = cal->sr_list->first; 
      prev_s = prev_item->item;
      for (linked_list_item_t *curr_item = prev_item->next; 
	   curr_item != NULL; 
	   curr_item = curr_item->next) {
	curr_s = curr_item->item;
	g_pos = sa_index->genome->chrom_offsets[chrom] + prev_s->genome_end + 1;
	for (read_pos = prev_s->read_end + 1; read_pos < curr_s->read_start; read_pos++) {
	  //printf("read pos. %lu %c - %c %lu genome pos.\n", read_pos, query[read_pos], ref[g_pos], g_pos);
	  if (query[read_pos] != ref[g_pos]) {
	    num_mismatches++;
	    if (num_mismatches > max_num_mismatches) {
	      // CAL with bad 'coverage' is marked to be removed
	      num_invalid++;
	      invalid[j] = 1;
	      break;
	    }
	  }
	  g_pos++;
	}
	prev_item = curr_item;
	prev_s = curr_s;
      }

      if (!invalid[j]) {
	
	// first gap
	s = linked_list_get_first(cal->sr_list);
	if (s->read_start > 0) {
	  g_pos = sa_index->genome->chrom_offsets[chrom] + s->genome_start - 1;
	  for (read_pos = s->read_start - 1; read_pos >= 0; read_pos--) {
	    //	  printf("read pos. %lu %c - %c %lu genome pos.\n", read_pos, query[read_pos], ref[g_pos], g_pos);
	    if (query[read_pos] != ref[g_pos]) {
	      num_mismatches++;
	      if (num_mismatches > max_num_mismatches) {
		// CAL with bad 'coverage' is marked to be removed
		num_invalid++;
		invalid[j] = 1;
		break;
	      }
	    }
	    g_pos--;
	  }
	}
	// last gap
	s = linked_list_get_last(cal->sr_list);
	if (!invalid[j] && s->read_end < read->length - 1) {
	  g_pos = sa_index->genome->chrom_offsets[chrom] + s->genome_end + 1;
	  for (read_pos = s->read_end + 1; read_pos < read->length; read_pos++) {
	    //	  printf("read pos. %lu %c - %c %lu genome pos.\n", read_pos, query[read_pos], ref[g_pos], g_pos);
	    if (query[read_pos] != ref[g_pos]) {
	      num_mismatches++;
	      if (num_mismatches > max_num_mismatches) {
		// CAL with bad 'coverage' is marked to be removed
		num_invalid++;
		invalid[j] = 1;
		break;
	      }
	    }
	    g_pos++;
	  }
	}
      }
      cal->read_area = num_mismatches;
      if (num_mismatches < min_num_mismatches) {
	min_num_mismatches = num_mismatches;
      }
      #ifdef _VERBOSE	  
      printf("\t\t ---> num. mismatches = %i (max. %i) -> is valid ? %i\n", 
	     num_mismatches, max_num_mismatches, !invalid[j]);
      #endif	  
    }
    
    /*
    //    printf("\tcal %i of %i: sr_list size = %i (cal->num_seeds = %i) %i:%lu-%lu\n", 
    //	   j, num_cals, cal->sr_list->size, cal->num_seeds,
    //	   cal->chromosome_id, cal->start, cal->end);
    if (cal->sr_list->size > 1) {
      continue;
    } else if (cal->sr_list->size == 1) {
      int start = 0;
      size_t genome_start = 0;
      int first = 1;
      for (linked_list_item_t *list_item = cal->sr_list->first; list_item != NULL; list_item = list_item->next) {
	seed_region_t *s = list_item->item;
	matches += (s->read_end - s->read_start + 1);
      }
      if (matches * 100 / read_len < min_perc) {
	// CAL with bad 'coverage' is marked to be removed
	num_invalid++;
	invalid[j] = 1;
      }
    } else {
      // CAL without seeds is marked to be removed
      num_invalid++;
      invalid[j] = 1;
    }
    */
  }

  //  if (num_invalid) {
    //  if (num_invalid != num_cals && num_invalid > 20) {
    array_list_t *new_cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    for (size_t j = 0; j < num_cals; j++) {
      if (!invalid[j]) {
	//	printf("---> read_area = %i, min_num_mismatches = %i\n", cal->read_area, min_num_mismatches);
	cal = array_list_get(j, cal_list);
	if (cal->read_area == min_num_mismatches) {
	  array_list_insert(cal, new_cal_list);
	  array_list_set(j, NULL, cal_list);
	}
      }
    }
    array_list_free(cal_list, (void *) cal_free);
    *list = new_cal_list;
    //  }
}

//--------------------------------------------------------------------

void display_seed(seed_region_t *s, int chrom, sa_index3_t *sa_index) {
  size_t start = s->genome_start;// + 1;
  size_t end = s->genome_end;// + 1;
  size_t len = end - start + 1;
  char *seq = sa_genome_get_sequence(chrom, start, end, sa_index->genome);
  printf("\tseed: [%i|%i - %i|%i] %s (len = %i)\n", 
	 s->genome_start, s->read_start, s->read_end, s->genome_end, seq, strlen(seq));
  free(seq);
}

//--------------------------------------------------------------------

size_t get_distance(char *seq1, char *seq2, int len) {
  size_t distance = 0;
  for (size_t i = 0; i < len; i++) {
    if (seq1[i] != seq2[i]) distance++;
  }
  #ifdef _VERBOSE
  printf("\t\tdistance = %lu\n", distance);
  printf("\t\t\tquery: %s\n", seq1);
  printf("\t\t\tref. : %s\n", seq2);
  #endif
  return distance;
}

//--------------------------------------------------------------------

#define SINGLE_FLANK 0 //2
#define DOUBLE_FLANK 0 //4

void fill_cal_gaps(fastq_read_t *read, char *revcomp_seq,
		   sa_index3_t *sa_index, array_list_t *cal_list) {
  uint num_mismatches, num_regions, num_cals = array_list_size(cal_list);
  
  char *seq, *ref, *query;

  cal_t *cal;
  seed_region_t *s, *new_s, *prev_s, *curr_s;
  linked_list_iterator_t* itr;

  cigar_op_t *cigar_op;
  cigar_code_t *cigar_code;

  size_t start, end, read_len = read->length;
  size_t gap_read_start, gap_read_end, gap_read_len;
  size_t gap_genome_start, gap_genome_end, gap_genome_len;
  size_t min_gap = sa_index->k_value;

  int distance, min_distance, left_flank, right_flank;
  sw_prepare_t *sw_prepare;
  array_list_t *sw_prepare_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  size_t g_pos;

  for (size_t j = 0; j < num_cals; j++) {
    //    matches = 0;
    //    invalid[j] = 0;
    cal = array_list_get(j, cal_list);
    cal->read_area = 0;

    #ifdef _VERBOSE
    printf("fill_cal_gaps, cal %i of %i:\n", j, num_cals);
    cal_print(cal);
    #endif

    prev_s = NULL;
    itr = linked_list_iterator_new(cal->sr_list);
    s = (seed_region_t *) linked_list_iterator_curr(itr);
    while (s != NULL) {

      #ifdef _VERBOSE
      display_seed(s, cal->chromosome_id, sa_index);
      #endif

      // set the cigar for the current region
      gap_read_len = s->read_end - s->read_start + 1;
      cigar_code = cigar_code_new();
      cigar_code_append_op(cigar_op_new(gap_read_len, 'M'), cigar_code);
      s->info = (void *) cigar_code;
      
      cigar_code = NULL;
      sw_prepare = NULL;
      
      if ((prev_s == NULL && s->read_start != 0) || (prev_s != NULL)) {
	num_mismatches = 0;
	if (prev_s == NULL) {
	  // gap at the first position
	  gap_read_start = 0;
	  gap_read_end = s->read_start - 1;
	  gap_read_len = gap_read_end - gap_read_start + 1;
	  
	  gap_genome_start = s->genome_start - s->read_start;
	  gap_genome_end = s->genome_start - 1;
	  gap_genome_len = gap_genome_end - gap_genome_start + 1;
	  
	  cal->start = gap_genome_start;
	  
	  assert(gap_read_len != 0);
	  assert(gap_genome_len != 0);
	  
	  if (gap_read_len > min_gap) {
	    // the gap is too big, may be there's another CAL to cover it
            #ifdef _VERBOSE
            printf("\tthe left-most gap is too big (%lu, max. size = %lu)\n", gap_read_len, min_gap);
            #endif
	    cigar_code = cigar_code_new();
	    cigar_code_append_op(cigar_op_new(gap_read_len, 'H'), cigar_code);	      
	  } else {
	    // the gap is small, we can process it now !
            #ifdef _VERBOSE
            printf("\tthe left-most gap is small, it will be processed now\n");
            #endif
	    left_flank = 0;
	    right_flank = DOUBLE_FLANK;
	  }
	} else {
	  // gap in a middle position
	  //assert(prev_s->read_end < s->read_start);
	  if (prev_s->read_end >= s->read_start) {
	    printf("\t\t\tassert failed: read %s: prev_s->read_end (%lu) >= (%lu) s->read_start\n", 
		   read->id, prev_s->read_end, s->read_start);
	  }
	 
	  if (prev_s->read_end == s->read_start) {
	    gap_read_len = 0;
	  } else {
	    gap_read_start = prev_s->read_end + 1;
	    gap_read_end = s->read_start - 1;
	    gap_read_len = gap_read_end - gap_read_start + 1;
	  }	    

	  gap_genome_start = prev_s->genome_end + 1;
	  gap_genome_end = s->genome_start - 1;
	  gap_genome_len = gap_genome_end - gap_genome_start + 1;
	  
	  #ifdef _VERBOSE
	  printf("gap (read, genome) = (%i, %i)\n", gap_read_len, gap_genome_len);
	  #endif
	  
	  if (gap_genome_len == 0) { printf("#@#: %s\n", read->id); }
	  assert(gap_genome_len != 0);
	  
	  if (gap_read_len == 0) {
	    // there's a deletion just between two consecutives seeds
	    cigar_code = (cigar_code_t *)prev_s->info;
	    
	    cigar_code_append_op(cigar_op_new(gap_genome_len, 'D'), cigar_code);
	    cigar_code->distance += gap_genome_len;
	    
	    cigar_code_append_op(cigar_op_new(s->read_end - s->read_start + 1, 'M'), cigar_code);
	    cigar_code->distance += ((cigar_code_t *)s->info)->distance;
	    
	    prev_s->read_end = s->read_end;
	    prev_s->genome_end = s->genome_end;
	    
            #ifdef _VERBOSE
	    printf("prev cigar = %s\n", new_cigar_code_string((cigar_code_t *)prev_s->info));
            #endif
	    
	    // continue loop...
	    linked_list_iterator_remove(itr);
	    s = linked_list_iterator_curr(itr);
	    continue;
	  }
	
	  left_flank = SINGLE_FLANK;
	  right_flank = SINGLE_FLANK;
	} // end of middle gap

	if (!cigar_code) {
	  // we have to try to fill this gap and get a cigar
	  seq = (cal->strand ? revcomp_seq : read->sequence);
	  ref = sa_genome_get_sequence(cal->chromosome_id, gap_genome_start, gap_genome_end, sa_index->genome);
	  if (gap_read_len == gap_genome_len) {
	    // first, compare char by char
	    distance = get_distance(&seq[gap_read_start], ref, gap_read_len);
	    min_distance = ceil(0.3 * gap_read_len); // 30%
	    if (distance <= min_distance) {
	      cigar_code = cigar_code_new();
	      cigar_code_append_op(cigar_op_new(gap_read_len, 'M'), cigar_code);
	      cigar_code_inc_distance(distance, cigar_code);
	      free(ref);
	    }
	  }

	  if (!cigar_code) {
	    // otherwise, prepare sequqences to run SW
	    seq = get_subsequence(seq, gap_read_start, gap_read_len);
	    sw_prepare = sw_prepare_new(seq, ref, 0, 0, (prev_s == NULL ? FIRST_SW : MIDDLE_SW));
	    sw_prepare->seed_region = s;
	    sw_prepare->cal = cal;
	    sw_prepare->read = read;
	    array_list_insert(sw_prepare, sw_prepare_list);
	    
            #ifdef _VERBOSE
	    display_seed(s, cal->chromosome_id, sa_index);
	    #endif
	  }
	}
      
	// insert gap in the list
	new_s = seed_region_new(gap_read_start, gap_read_end, gap_genome_start, gap_genome_end, 0);
	new_s->info = (void *) cigar_code;
	linked_list_iterator_insert(new_s, itr);
	
	if (sw_prepare) {
	  sw_prepare->seed_region = new_s;
	  sw_prepare->cal = cal;
	  sw_prepare->read = read;
	}
      }
      
      // continue loop...
      prev_s = s;
      linked_list_iterator_next(itr);
      s = linked_list_iterator_curr(itr);
    } // end of while 
      
    // check for a gap at the last position
    sw_prepare = NULL;
    if (prev_s != NULL && prev_s->read_end < read_len - 1) { 
      cigar_code = NULL;
      
      // gap at the last position
      gap_read_start = prev_s->read_end + 1;
      gap_read_end = read_len - 1;
      gap_read_len = gap_read_end - gap_read_start + 1;
      assert(gap_read_len != 0);
      
      gap_genome_len = gap_read_len;
      gap_genome_start = prev_s->genome_end + 1;
      gap_genome_end = gap_genome_start + gap_genome_len - 1;
      assert(gap_genome_len != 0);
      
      cal->end = gap_genome_end;
      
      if (gap_read_len > min_gap) {
	// the gap is too big, may be there's another CAL to cover it
        #ifdef _VERBOSE
	printf("\tthe right-most gap is too big (%lu, max. size = %lu)\n", gap_read_len, min_gap);
        #endif
	cigar_code = cigar_code_new();
	cigar_code_append_op(cigar_op_new(gap_read_len, 'H'), cigar_code);	      
      } else {
	// the gap is small, we can process it now !
        #ifdef _VERBOSE
	printf("\tthe right-most gap is small, it will be processed now\n");
        #endif

	// we have to try to fill this gap and get a cigar
	// first, compare char by char
	seq = (cal->strand ? revcomp_seq : read->sequence);
	ref = sa_genome_get_sequence(cal->chromosome_id, gap_genome_start, gap_genome_end, sa_index->genome);
	
	distance = get_distance(&seq[gap_read_start], ref, gap_read_len);
	min_distance = ceil(0.3 * gap_read_len); // 30%
	if (distance <= min_distance) {
	  cigar_code = cigar_code_new();
	  cigar_code_append_op(cigar_op_new(gap_read_len, 'M'), cigar_code);
	  cigar_code_inc_distance(distance, cigar_code);
	  free(ref);
	} else {
	  // otherwise, prepare sequences to run SW
	  seq = get_subsequence(seq, gap_read_start, gap_read_len);
	  sw_prepare = sw_prepare_new(seq, ref, 0, 0, LAST_SW);
	  sw_prepare->seed_region = prev_s;
	  sw_prepare->cal = cal;
	  sw_prepare->read = read;
	  array_list_insert(sw_prepare, sw_prepare_list);
	  
          #ifdef _VERBOSE
	  display_seed(prev_s, cal->chromosome_id, sa_index);
	  #endif
	}
      }
      
      // insert gap in the list
      new_s = seed_region_new(gap_read_start, gap_read_end, gap_genome_start, gap_genome_end, 0);
      new_s->info = (void *) cigar_code;
      linked_list_insert_last(new_s, cal->sr_list);
      
      if (sw_prepare) {
	sw_prepare->seed_region = new_s;
	sw_prepare->cal = cal;
	sw_prepare->read = read;
      }
    }
    linked_list_iterator_free(itr);      
  }

  // run Smith-Waterman
  size_t sw_count = array_list_size(sw_prepare_list); 
  #ifdef _VERBOSE
  printf("num. sw to do = %lu\n", sw_count);
  #endif
 if (sw_count <= 0) {
    // no Smith-Waterman to run
   array_list_free(sw_prepare_list, (void *) NULL);
    return;
  }

  sw_optarg_t sw_optarg;
  sw_optarg_init(10, 0.5, 5, -4, &sw_optarg);

  char *q[sw_count], *r[sw_count];
  for (int i = 0; i < sw_count; i++) {
    sw_prepare = array_list_get(i, sw_prepare_list);
    q[i] = sw_prepare->query;
    r[i] = sw_prepare->ref;
    #ifdef _VERBOSE
    printf("\t\t%i: query: %s\n", i, q[i]);
    printf("\t\t%i: ref. : %s\n", i, r[i]);
    printf("\n");
    #endif
  }
  
  sw_multi_output_t *sw_output = sw_multi_output_new(sw_count);
  smith_waterman_mqmr(q, r, sw_count, &sw_optarg, 1, sw_output);
  #ifdef _VERBOSE
  sw_multi_output_save(sw_count, sw_output, stdout);
  #endif

  // process Smith-Waterman output
  for (int i = 0; i < sw_count; i++) {
    sw_prepare = array_list_get(i, sw_prepare_list);
    s = sw_prepare->seed_region;

    int read_gap_len = s->read_end - s->read_start + 1;
    int genome_gap_len = s->genome_end - s->genome_start + 1;

    int read_gap_len_ex = read_gap_len + sw_prepare->left_flank + sw_prepare->right_flank;
    int genome_gap_len_ex = genome_gap_len + sw_prepare->left_flank + sw_prepare->right_flank;

    cigar_code_t *cigar_c = generate_cigar_code(sw_output->query_map_p[i], sw_output->ref_map_p[i],
						strlen(sw_output->query_map_p[i]), sw_output->query_start_p[i],
						sw_output->ref_start_p[i], read_gap_len, genome_gap_len,
						&distance, sw_prepare->ref_type);

    #ifdef _VERBOSE
    printf("\tpost-processing gap (read %lu-%lu, genome %lu-%lu) = (%i, %i): read %s\n", 
	   s->read_start, s->read_end, s->genome_start, s->genome_end,
	   read_gap_len, genome_gap_len, sw_prepare->read->id);
    display_seed(s, sw_prepare->cal->chromosome_id, sa_index);
    printf("\t\tflanks (left, right) = (%i, %i)\n", sw_prepare->left_flank, sw_prepare->right_flank);
    printf("\t\tquery : %s\n", sw_prepare->query);
    printf("\t\tref   : %s\n", sw_prepare->ref);
    printf("\t\tmquery: %s (start %i)\n", sw_output->query_map_p[i], sw_output->query_start_p[i]);
    printf("\t\tmref  : %s (start %i)\n", sw_output->ref_map_p[i], sw_output->ref_start_p[i]);
    printf("\t\t\tscore : %0.2f, cigar: %s (distance = %i)\n", 
	   sw_output->score_p[i], new_cigar_code_string(cigar_c), distance);
    #endif

    cigar_op = cigar_code_get_op(0, cigar_c);
    if (cigar_op) {
      if (cigar_op->name == 'H') {
	if (sw_output->ref_start_p[i] == 0) { 
	  cigar_op->name = 'I';
	} else {
	  cigar_op->name = 'M';
	}
      } else if (cigar_op->name == '=') cigar_op->name = 'M';
    }
    
    cigar_op = cigar_code_get_last_op(cigar_c);
    if (cigar_op && cigar_op->name == 'H') cigar_op->name = 'I';

    #ifdef _VERBOSE
    printf("\t\t\tgap_read_len = %i, cigar_code_length (%s) = %i\n", 
	   read_gap_len, new_cigar_code_string(cigar_c), cigar_code_nt_length(cigar_c));
    #endif
    if (read_gap_len != cigar_code_nt_length(cigar_c)) {
      printf("\t\t\tassert failed: read %s: gap_read_len = %i, cigar_code_length (%s) = %i\n", 
	     read->id, read_gap_len, new_cigar_code_string(cigar_c), cigar_code_nt_length(cigar_c));
      exit(-1);      
    }
    //    assert(read_gap_len == cigar_code_nt_length(cigar_c));

    // and now set the cigar for this gap
    s->info = (void *) cigar_c;

    // free
    sw_prepare_free(sw_prepare);
  }

  // free memory
  sw_multi_output_free(sw_output);
  array_list_free(sw_prepare_list, (void *) NULL);
}

//--------------------------------------------------------------------

void merge_cal_seeds(array_list_t *cal_list) {
  seed_region_t *s, *s_first;
  cigar_code_t *cigar_code, *cigar_code_prev;
  cigar_op_t *cigar_op, *cigar_op_prev;
  int num_ops;
  int op;

  linked_list_item_t *list_item;
  linked_list_iterator_t itr;

  cal_t *cal;
  uint num_cals = array_list_size(cal_list);

  for (size_t i = 0; i < num_cals; i++) {
    cal = array_list_get(i, cal_list);
    linked_list_iterator_init(cal->sr_list, &itr);
    
    s_first = linked_list_iterator_curr(&itr);      

    if (s_first) {
      cigar_code_prev = (cigar_code_t *)s_first->info;
      s = linked_list_iterator_next(&itr);
      while (s) {
	cigar_code = (cigar_code_t *)s->info;
	if (cigar_code) {
	  num_ops = array_list_size(cigar_code->ops);
	  for (op = 0, cigar_op = array_list_get(op, cigar_code->ops); 
	       op < num_ops;
	       op++, cigar_op = array_list_get(op, cigar_code->ops)) {
	    array_list_set(op, NULL, cigar_code->ops); 
	    cigar_code_append_op(cigar_op, cigar_code_prev);	    
	  }
	  cigar_code_prev->distance += cigar_code->distance;
	  cigar_code_free(cigar_code);
	} 
	
	s_first->read_end = s->read_end;
	s_first->genome_end = s->genome_end;

	seed_region_free(s);

	linked_list_iterator_remove(&itr);	
	s = linked_list_iterator_curr(&itr);
      }
      cal->info = (void *)cigar_code_prev;
    } else {
      cal->info = NULL;
    }
  }
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

void concat_cigar_op(int number, char name, cigar_code_t *dst) {
  #ifdef _VERBOSE
  printf("\tconcating %i%c to %s...\n", number, name, new_cigar_code_string(dst));
  #endif
  if (cigar_code_get_num_ops(dst) > 0) {
    cigar_op_t *op = cigar_code_get_last_op(dst);
    if (name == 'M' && op->name == 'M') {
      op->number += number;
    } else {
      cigar_code_append_new_op(number, name, dst);
    }
  } else {
    cigar_code_append_new_op(number, name, dst);
  }
  #ifdef _VERBOSE
  printf("\t\tconcating done: %s\n", new_cigar_code_string(dst));
  #endif
}

//--------------------------------------------------------------------

void concat_cigar_ops(cigar_code_t *src, int start, int end, cigar_code_t *dst) {
  cigar_op_t *op;
  #ifdef _VERBOSE
  printf("concating %s (%i-%i) to %s...\n", new_cigar_code_string(src), start, end, new_cigar_code_string(dst));
  #endif
  
  for (int i = start; i <= end; i++) {
    op = cigar_code_get_op(i, src);
    concat_cigar_op(op->number, op->name, dst);
  }

  dst->distance += src->distance;
}

//--------------------------------------------------------------------

void fill_cal_end_gaps(fastq_read_t *read, char *revcomp_seq,
		       sa_index3_t *sa_index, array_list_t *cal_list) {
  
  cal_t *cal;

  seed_region_t *s;

  cigar_op_t *cigar_op, *first_cigar_op, *last_cigar_op;
  cigar_code_t *cigar_code, *cal_cigar_code;

  size_t gap_read_start, gap_read_end, gap_read_len;
  size_t gap_genome_start, gap_genome_end, gap_genome_len;

  int min_distance, distance, min_H, mode;
  sw_prepare_t *sw_prepare;
  array_list_t *sw_prepare_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  
  char *seq, *ref, *query;

  min_H = 5;

  int num_cigar_ops;
  size_t read_len = read->length;
  size_t num_cals = array_list_size(cal_list);

  // processing each CAL from this read, checking the initial and
  // ending gaps to initialize query and ref. sequences to Smith-Waterman
  for(size_t i = 0; i < num_cals; i++) {

    // get cal
    cal = array_list_get(i, cal_list);
    if (cal->sr_list->size == 0) continue;

    //    cal->info = (void *) cigar_code_new();

    sw_prepare = NULL;
    s = (seed_region_t *) linked_list_get_first(cal->sr_list);
    cigar_code = (cigar_code_t *) s->info;

    num_cigar_ops = cigar_code_get_num_ops(cigar_code);

    #ifdef _VERBOSE
    printf("fill_cal_end_gaps, cal %i\n", i);
    printf("\t\tCAL %i of %i (strand %i), sr_list size = %i, cigar = %s (num. ops = %i, distance = %i)\n", 
	   i, num_cals, cal->strand, cal->sr_list->size, 
	   new_cigar_code_string(cigar_code), num_cigar_ops, cigar_code->distance);
    #endif
     
    if (num_cigar_ops == 1) {
      cal->info = (void *) cigar_code;
      s->info = (void *) NULL;
      continue;
    }

    cal_cigar_code = cigar_code_new();

    first_cigar_op = cigar_code_get_first_op(cigar_code);
    last_cigar_op = cigar_code_get_last_op(cigar_code);
    assert(first_cigar_op != last_cigar_op);

    // first gap, (looking for hard-clipping)
    if (first_cigar_op->name == 'H' && first_cigar_op->number > min_H) {
      // fill this gap: first, nt by nt, otherwise, prepare SW
      gap_read_start = 0;
      gap_read_end = first_cigar_op->number - 1;
      gap_read_len = gap_read_end - gap_read_start + 1;
      gap_genome_start = s->genome_start;
      gap_genome_end = gap_genome_start + first_cigar_op->number - 1;
      gap_genome_len = gap_genome_end - gap_genome_start + 1;

      // get query and ref sequences, revcomp if necessary
      seq = (cal->strand ? revcomp_seq : read->sequence);
      ref = sa_genome_get_sequence(cal->chromosome_id, gap_genome_start, gap_genome_end, sa_index->genome);

      distance = get_distance(&seq[gap_read_start], ref, gap_read_len);
      min_distance = ceil(0.3 * gap_read_len); // 30%

      #ifdef _VERBOSE
      printf("\t\tfirst gap: %i%c (distance = %i, max. distance = %i)\n", 
	     first_cigar_op->number, first_cigar_op->name, distance, min_distance);
      #endif

      if (distance <= min_distance) {
	// concat from this cigar op. (replacing H by M) to the last one
	cigar_code_append_new_op(gap_read_len, 'M', cal_cigar_code);
	cigar_code_inc_distance(distance, cal_cigar_code);
	concat_cigar_ops(cigar_code, 1, num_cigar_ops - 2, cal_cigar_code);
	free(ref);
      } else {
	seq = get_subsequence(seq, gap_read_start, gap_read_len);
	sw_prepare = sw_prepare_new(seq, ref, 0, 0, FIRST_SW);
	sw_prepare->seed_region = s;
	sw_prepare->cal = cal;
	sw_prepare->read = read;
	array_list_insert(sw_prepare, sw_prepare_list);

        #ifdef _VERBOSE
	display_seed(s, cal->chromosome_id, sa_index);
	#endif
      }
    } else {
      // concat cigar ops., from the first one to the last one
      concat_cigar_ops(cigar_code, 0, num_cigar_ops - 2, cal_cigar_code);
    }

    // last gap, (looking for hard-clipping)
    if (last_cigar_op->name == 'H' && last_cigar_op->number > min_H) {
      // fill this gap: first, nt by nt, otherwise, prepare SW
      gap_read_start = read_len - last_cigar_op->number;
      gap_read_end = read_len - 1;
      gap_read_len = gap_read_end - gap_read_start + 1;
      gap_genome_end = s->genome_end;
      gap_genome_start = gap_genome_end - last_cigar_op->number + 1;
      gap_genome_len = gap_genome_end - gap_genome_start + 1;

      // get query and ref sequences, revcomp if necessary
      seq = (cal->strand ? revcomp_seq : read->sequence);
      ref = sa_genome_get_sequence(cal->chromosome_id, gap_genome_start, gap_genome_end, sa_index->genome);

      distance = get_distance(&seq[gap_read_start], ref, gap_read_len);
      min_distance = ceil(0.3 * gap_read_len); // 30%

      #ifdef _VERBOSE
      printf("\t\tlast gap: %i%c (distance = %i, max. distance = %i)\n", 
	     last_cigar_op->number, last_cigar_op->name, distance, min_distance);
      #endif

      if (distance <= min_distance) {
	// concat from this cigar op. (replacing H by M) if the cal cigar is not null,
	// otherwise, save this last cigar op. in the cal cigar (same effect as concating)
	concat_cigar_op(gap_read_len, 'M', cal_cigar_code);
	cigar_code_inc_distance(distance, cal_cigar_code);
	free(ref);
      } else {
	seq = get_subsequence(seq, gap_read_start, gap_read_len);
	sw_prepare = sw_prepare_new(seq, ref, 0, 0, LAST_SW);
	sw_prepare->seed_region = s;
	sw_prepare->cal = cal;
	sw_prepare->read = read;
	array_list_insert(sw_prepare, sw_prepare_list);

        #ifdef _VERBOSE
	display_seed(s, cal->chromosome_id, sa_index);
	#endif
      }
    } else {
      // concat this last cigar op. if the cal cigar is not null,
      // otherwise, save this last cigar op. in the cal cigar (same effect as concating)
      concat_cigar_op(last_cigar_op->number, last_cigar_op->name, cal_cigar_code);
    }

    cal->info = cal_cigar_code;

    /*
    for (int k = 0; k < 2; k++) {
      mode = MIDDLE_SW;
      if (k == 0) {
	if ((cigar_op = cigar_code_get_op(0, cigar_code)) &&
	    cigar_op->name == 'H' && cigar_op->number > min_H) {

	  #ifdef _VERBOSE
	  printf("\t\tfirst sw: %i%c\n", cigar_op->number, cigar_op->name);
	  #endif
	  
	  mode = FIRST_SW;
	  gap_read_start = 0;
	  gap_read_end = cigar_op->number - 1;
	  gap_genome_start = s->genome_start;
	  gap_genome_end = gap_genome_start + cigar_op->number - 1;
	}
      } else {
	if ((cigar_op = cigar_code_get_last_op(cigar_code)) &&
	    cigar_op->name == 'H' && cigar_op->number > min_H) {
	  #ifdef _VERBOSE
	  printf("\t\tlast sw: %i%c\n", cigar_op->number, cigar_op->name);
	  #endif
	  
	  mode = LAST_SW;
	  gap_read_start = read_len - cigar_op->number;
	  gap_read_end = read_len - 1;
	  gap_genome_end = s->genome_end;
	  gap_genome_start = gap_genome_end - cigar_op->number + 1;
	}
      }
      
      if (mode == MIDDLE_SW) continue;

      gap_read_len = gap_read_end - gap_read_start + 1;
      gap_genome_len = gap_genome_end - gap_genome_start + 1;
      #ifdef _VERBOSE
      printf("\t\t(gap_read_len, gap_genome_len) = (%lu, %lu)\n", gap_read_len, gap_genome_len);
      printf("\t\t\t%lu|%lu-%lu|%lu\n", gap_genome_start, gap_read_start, gap_read_end, gap_genome_end);
      #endif

      // get query and ref sequences, revcomp if necessary
      seq = (cal->strand ? revcomp_seq : read->sequence);
      ref = sa_genome_get_sequence(cal->chromosome_id, gap_genome_start, gap_genome_end, sa_index->genome);

      distance = get_distance(&seq[gap_read_start], ref, gap_read_len);
      min_distance = ceil(0.3 * gap_read_len); // 30%
      if (distance <= min_distance) {
	cigar_op->name = 'M';
	cigar_op->number += gap_read_len;
	cigar_code->distance += distance;
	free(ref);
      } else {
	seq = get_subsequence(seq, gap_read_start, gap_read_len);
	sw_prepare = sw_prepare_new(seq, ref, 0, 0, mode);
	sw_prepare->seed_region = s;
	sw_prepare->cal = cal;
	sw_prepare->read = read;
	array_list_insert(sw_prepare, sw_prepare_list);

        #ifdef _VERBOSE
	display_seed(s, cal->chromosome_id, sa_index);
	#endif
      }
    }
    */
  }

  // run Smith-Waterman
  size_t sw_count = array_list_size(sw_prepare_list);
  #ifdef _VERBOSE
  printf("num. sw to do = %lu\n", sw_count);
  #endif
  if (sw_count <= 0) {
    // no Smith-Waterman to run
    array_list_free(sw_prepare_list, (void *) NULL);
    return;
  }

  sw_optarg_t sw_optarg;
  sw_optarg_init(10, 0.5, 5, -4, &sw_optarg);

  char *q[sw_count], *r[sw_count];
  for (int i = 0; i < sw_count; i++) {
    sw_prepare = array_list_get(i, sw_prepare_list);
    q[i] = sw_prepare->query;
    r[i] = sw_prepare->ref;
  }

  sw_multi_output_t *sw_output = sw_multi_output_new(sw_count);
  smith_waterman_mqmr(q, r, sw_count, &sw_optarg, 1, sw_output);
  #ifdef _VERBOSE
  sw_multi_output_save(sw_count, sw_output, stdout);
  #endif

  // process Smith-Waterman output
  for (int i = 0; i < sw_count; i++) {
    sw_prepare = array_list_get(i, sw_prepare_list);
    s = sw_prepare->seed_region;
    cal = sw_prepare->cal;

    int read_gap_len = strlen(sw_prepare->query); //s->read_end - s->read_start + 1;
    int genome_gap_len = strlen(sw_prepare->ref); //s->genome_end - s->genome_start + 1;

    //    int read_gap_len_ex = read_gap_len + sw_prepare->left_flank + sw_prepare->right_flank;
    //    int genome_gap_len_ex = genome_gap_len + sw_prepare->left_flank + sw_prepare->right_flank;

    cigar_code_t *cigar_c = generate_cigar_code(sw_output->query_map_p[i], sw_output->ref_map_p[i],
						strlen(sw_output->query_map_p[i]), sw_output->query_start_p[i],
						sw_output->ref_start_p[i], read_gap_len, genome_gap_len,
						&distance, sw_prepare->ref_type);
    #ifdef _VERBOSE
    printf("\tpost-processing gap (read %lu-%lu, genome %lu-%lu) = (%i, %i): read %s\n", 
	   s->read_start, s->read_end, s->genome_start, s->genome_end,
	   read_gap_len, genome_gap_len, sw_prepare->read->id);
    display_seed(s, sw_prepare->cal->chromosome_id, sa_index);
    printf("\t\tflanks (left, right) = (%i, %i)\n", sw_prepare->left_flank, sw_prepare->right_flank);
    printf("\t\tquery : %s\n", sw_prepare->query);
    printf("\t\tref   : %s\n", sw_prepare->ref);
    printf("\t\tmquery: %s (start %i)\n", sw_output->query_map_p[i], sw_output->query_start_p[i]);
    printf("\t\tmref  : %s (start %i)\n", sw_output->ref_map_p[i], sw_output->ref_start_p[i]);
    printf("\t\t\tscore : %0.2f, cigar: %s (distance = %i)\n", 
	   sw_output->score_p[i], new_cigar_code_string(cigar_c), distance);
    #endif

    cal_cigar_code = (cigar_code_t *) cal->info;
    num_cigar_ops = cigar_code_get_num_ops(cal_cigar_code);
    if (sw_prepare->ref_type == FIRST_SW) {
      // if num. cigar ops > 0, cal_cigar_code contains the cigar for the last gap
      // the concat middle cigar ops and then, the cigar from the last gap
      cigar_code = s->info;
      concat_cigar_ops(cigar_code, 1, cigar_code_get_num_ops(cigar_code) - 2, cigar_c);
      if (num_cigar_ops > 0) {
	// this means that in cal_cigar_code we saved the cigar for the last gap,
	// we must concat it
	concat_cigar_ops(cal_cigar_code, 0, num_cigar_ops - 1, cigar_c);
      }
      cigar_code_free(cal_cigar_code);
      cal->info = (cigar_code_t *) cigar_c;
      #ifdef _VERBOSE
      printf("\t\t\t\tLAST_SW: final cigar: %s (distance = %i)\n", 
	     new_cigar_code_string(cigar_c), cigar_c->distance);
      #endif
    } else if (sw_prepare->ref_type == LAST_SW) {
      num_cigar_ops = cigar_code_get_num_ops(cigar_c);
      if (cigar_code_get_last_op(cigar_c)->name == 'D') {
	num_cigar_ops--;
      }
      concat_cigar_ops(cigar_c, 0, num_cigar_ops - 1, cal_cigar_code);
      cigar_code_free(cigar_c);      
      #ifdef _VERBOSE
      printf("\t\t\t\tLAST_SW: final cigar: %s (distance = %i)\n", 
	     new_cigar_code_string(cal_cigar_code), cal_cigar_code->distance);
      #endif
    }

    /*
    cigar_op = cigar_code_get_op(0, cigar_c);
    if (cigar_op) {
      if (cigar_op->name == 'H') {
	if (sw_output->ref_start_p[i] == 0) { 
	  cigar_op->name = 'I';
	} else {
	  cigar_op->name = 'M';
	}
      } else if (cigar_op->name == '=') cigar_op->name = 'M';
    }
    
    cigar_op = cigar_code_get_last_op(cigar_c);
    if (cigar_op && cigar_op->name == 'H') cigar_op->name = 'I';

    #ifdef _VERBOSE
    printf("\t\t\tgap_read_len = %i, cigar_code_length (%s) = %i\n", 
	   read_gap_len, new_cigar_code_string(cigar_c), cigar_code_nt_length(cigar_c));
    #endif
    if (read_gap_len != cigar_code_nt_length(cigar_c)) {
      printf("\t\t\tassert failed: read %s: gap_read_len = %i, cigar_code_length (%s) = %i\n", 
	     read->id, read_gap_len, new_cigar_code_string(cigar_c), cigar_code_nt_length(cigar_c));
      exit(-1);      
    }
    //    assert(read_gap_len == cigar_code_nt_length(cigar_c));

    // and now set the cigar for this gap
    if (s->info) cigar_code_free(s->info);
    s->info = (void *) cigar_c;
    */


    // free
    sw_prepare_free(sw_prepare);
  }  

  // free memory
  sw_multi_output_free(sw_output);
  array_list_free(sw_prepare_list, (void *) NULL);
}

//--------------------------------------------------------------------

void filter_invalid_cals(fastq_read_t *read, array_list_t **list) {
  // filter-incoherent CALs
  cal_t *cal;
  array_list_t *cal_list = *list;
  
  size_t num_cals = array_list_size(cal_list);
  int founds[num_cals], found = 0;

  for (size_t j = 0; j < num_cals; j++) {
    founds[j] = 0;
    cal = array_list_get(j, cal_list);
    if (cal->sr_list->size > 0) {
      int start = 0;
      size_t genome_start = 0;
      int first = 1;
      for (linked_list_item_t *list_item = cal->sr_list->first; 
	   list_item != NULL; 
	   list_item = list_item->next) {
	seed_region_t *s = list_item->item;
	
	if (start > s->read_start || s->read_start >= s->read_end) {
	  found++;
	  founds[j] = 1;
	}

	if (!first && 
	    ((s->genome_start < genome_start) || 
	     (s->genome_start - genome_start) > 2*read->length)) {
	  found++;
	  founds[j] = 1;
	}
	
	first = 0;
	start = s->read_end + 1;
	genome_start = s->genome_end + 1;
      }
    } else {
      found++;
      founds[j] = 1;
    }
  }

  if (found) {
    array_list_t *new_cal_list = array_list_new(num_cals - found + 1, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    for (size_t j = 0; j < num_cals; j++) {
      if (!founds[j]) {
	cal = array_list_get(j, cal_list);
	array_list_insert(cal, new_cal_list);
	array_list_set(j, NULL, cal_list);
      }
    }
    array_list_free(cal_list, (void *) cal_free);
    *list = new_cal_list;
  }
}

//--------------------------------------------------------------------

void filter_cals_by_min_read_area(int read_area, array_list_t **list) {
  cal_t *cal;
  array_list_t *cal_list = *list;
  size_t num_cals = array_list_size(cal_list);
  array_list_t *new_cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  for (size_t j = 0; j < num_cals; j++) {
    cal = array_list_get(j, cal_list);
    if (cal->read_area <= read_area) {
      array_list_insert(cal, new_cal_list);
      array_list_set(j, NULL, cal_list);
    }
  }
  array_list_free(cal_list, (void *) cal_free);
  *list = new_cal_list;
}

//--------------------------------------------------------------------

int get_min_num_mismatches(array_list_t *cal_list) {
  cal_t *cal;
  size_t num_cals = array_list_size(cal_list);
  int min_num_mismatches = 100000;
  for (size_t j = 0; j < num_cals; j++) {
    cal = array_list_get(j, cal_list);
    if (cal->num_mismatches < min_num_mismatches) {
      min_num_mismatches = cal->num_mismatches;
    }
  }
  return min_num_mismatches;
}

int get_max_read_area(array_list_t *cal_list) {
  cal_t *cal;
  size_t num_cals = array_list_size(cal_list);
  int max_read_area = 0;
  for (size_t j = 0; j < num_cals; j++) {
    cal = array_list_get(j, cal_list);
    if (cal->read_area > max_read_area) {
      max_read_area = cal->read_area;
    }
  }

  return max_read_area;
}

void filter_cals_by_max_read_area(int read_area, array_list_t **list) {
  cal_t *cal;
  array_list_t *cal_list = *list;
  size_t num_cals = array_list_size(cal_list);
  array_list_t *new_cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  for (size_t j = 0; j < num_cals; j++) {
    cal = array_list_get(j, cal_list);
    if (cal->read_area >= read_area) {
      array_list_insert(cal, new_cal_list);
      array_list_set(j, NULL, cal_list);
    }
  }
  array_list_free(cal_list, (void *) cal_free);
  *list = new_cal_list;
}

//--------------------------------------------------------------------

void filter_cals_by_max_num_mismatches(int num_mismatches, array_list_t **list) {
  cal_t *cal;
  array_list_t *cal_list = *list;
  size_t num_cals = array_list_size(cal_list);
  array_list_t *new_cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  for (size_t j = 0; j < num_cals; j++) {
    cal = array_list_get(j, cal_list);
    if (cal->num_mismatches <= num_mismatches) {
      array_list_insert(cal, new_cal_list);
      array_list_set(j, NULL, cal_list);
    }
  }
  array_list_free(cal_list, (void *) cal_free);
  *list = new_cal_list;
}

//--------------------------------------------------------------------

void filter_cals_by_score(int min_perc, fastq_read_t *read, char *revcomp_seq,
			 sa_index3_t *sa_index, array_list_t **list) {

  cal_t *cal;
  array_list_t *cal_list = *list;
  
  seed_region_t *s;
  cigar_code_t *cigar_code;

  size_t read_len = read->length;
  size_t num_cals = array_list_size(cal_list);

  int min_distance = read_len - read_len * min_perc / 100;

  int invalid[num_cals], num_invalid = 0;

  // processing each CAL from this read, checking the initial and
  // ending gaps to initialize query and ref. sequences to Smith-Waterman
  for(size_t i = 0; i < num_cals; i++) {
    // get cal
    cal = array_list_get(i, cal_list);
    if (cal->sr_list->size == 0) {
      num_invalid++;
      invalid[i] = 1;
    } else {    
      cigar_code = (cigar_code_t *) cal->info;
      if (cigar_code->distance > min_distance) {
	num_invalid++;
	invalid[i] = 1;
      }
      cigar_code_free(cigar_code);
    }
  }

  if (num_invalid) {
    array_list_t *new_cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    for (size_t j = 0; j < num_cals; j++) {
      if (!invalid[j]) {
	cal = array_list_get(j, cal_list);
	array_list_insert(cal, new_cal_list);
	array_list_set(j, NULL, cal_list);
      }
    }
    array_list_free(cal_list, (void *) cal_free);
    *list = new_cal_list;
  }
}

//--------------------------------------------------------------------

void create_alignments(array_list_t *cal_list, fastq_read_t *read, 
		       array_list_t *mapping_list) {

  // CAL
  cal_t *cal;
  uint num_cals = array_list_size(cal_list);

  // alignments
  char cigar[1000];
  int num_cigar_ops;
  alignment_t *alignment;

  if (num_cals <= 0) {
    // no CALs -> no alignment
    return;
  }

  int AS;
  size_t i, pos;
  linked_list_item_t *list_item; 
  seed_region_t *s_first, *s_last;

  for (i = 0; i < num_cals; i++) {
    cigar[0] = 0;
    cal = array_list_get(i, cal_list);

    #ifdef _VERBOSE	  
    printf("--> CAL #%i (cigar %s):\n", i, cigar_to_string(cal->info));
    cal_print(cal);
    #endif

    // computing aligments
    AS = read->length * 5;
    s_first = linked_list_get_first(cal->sr_list);
    s_last = linked_list_get_last(cal->sr_list);

    // sanity checking
    //    assert(s_first != NULL && s_last != NULL);

    if (s_first != NULL && s_last != NULL) {
      //      printf("-----> cigar_code = %x\n", s_first->info);
      //      printf("-----> cigar_code = %x\n", s_last->info);
      /*
      if (num_cals > 1) {
	printf("\t\t%s\t%i of %i\t[%i|%i - %i|%i]\n", 
	       read->id, i, num_cals, 
	       s_first->genome_start, s_first->read_start, s_first->read_end, s_first->genome_end);
	/*
	cigar_code_t *cigar = (cigar_code_t *) cal->info;
	printf("\t\t%s\t%i of %i\t[%i|%i - %i|%i] (cigar %s, distance = %i)\n", 
	       read->id, i, num_cals, 
	       s_first->genome_start, s_first->read_start, s_first->read_end, s_first->genome_end,
	       (cigar != NULL ? new_cigar_code_string(cigar) : "none"), 
	       (cigar != NULL ? cigar->distance : 0));
	*/
      //      }

      num_cigar_ops = 0;
      pos = cal->start;
      if (s_first->read_start > 0) {
	pos -= s_first->read_start;
	//sprintf(cigar, "%s%iS", cigar, s_first->read_start);
	//      num_cigar_ops++;
      }
      sprintf(cigar, "%s%iM", cigar, read->length);
      num_cigar_ops++;
      /*
      s_first->read_start = 0; // be carefull !!! to remove...now !!!
      sprintf(cigar, "%s%iM", cigar, s_last->read_end - s_first->read_start + 1);
      num_cigar_ops++;
      if (s_last->read_end < read->length - 1) {
	sprintf(cigar, "%s%iS", cigar, read->length - s_last->read_end - 1);
	num_cigar_ops++;
      }
      */
      alignment = alignment_new();	       
      alignment_init_single_end(strdup(read->id), strdup(read->sequence), strdup(read->quality), 
				cal->strand, cal->chromosome_id, 
				pos,
				strdup(cigar), num_cigar_ops, (AS * 254) / (read->length * 5), 1, (num_cals > 1),
				0, 0, alignment);  
      
      array_list_insert(alignment, mapping_list);
      //      array_list_insert(convert_to_bam(alignment, 33), mapping_list);
      
      // free memory
      //      alignment_free(alignment);
    } else if (cal->sr_duplicate_list && cal->sr_duplicate_list->size < 20) {
      /*
      // sr_list is empty, check sr_duplicate_list
      for (linked_list_item_t *list_item = cal->sr_duplicate_list->first; 
	   list_item != NULL; 
	   list_item = list_item->next) {
	seed_region_t *s = list_item->item;
	if (s->read_end - s->read_start + 1 == read->length) {
	  sprintf(cigar, "%iM", read->length);
	  num_cigar_ops = 1;
	  pos = s->genome_start;
	  
	  alignment = alignment_new();	       
	  alignment_init_single_end(strdup(read->id), strdup(read->sequence), strdup(read->quality), 
				    cal->strand, cal->chromosome_id, 
				    pos,
				    strdup(cigar), num_cigar_ops, (AS * 254) / (read->length * 5), 1, (num_cals > 1),
				    0, 0, 0, alignment);  
	  
	  array_list_insert(alignment, mapping_list);
	  //	  array_list_insert(convert_to_bam(alignment, 33), mapping_list);
	  printf("\t\t%s\t[%i|%i - %i|%i]\n", read->id, s->genome_start, s->read_start, s->read_end, s->genome_end);
	  
	  // free memory
	  //	  alignment_free(alignment);
	}
      */
    }
    
    if (cal->info) cigar_free(cal->info);
    cal_free(cal);
  }
}

//--------------------------------------------------------------------

void display_suffix_mappings(int strand, size_t r_start, size_t suffix_len, 
			     size_t low, size_t high, sa_index3_t *sa_index) {
  int chrom;
  size_t r_end, g_start, g_end;
  for (size_t suff = low; suff <= high; suff++) {
    r_end = r_start + suffix_len - 1;
    chrom = sa_index->CHROM[suff];
    g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];
    g_end = g_start + suffix_len - 1;
    printf("\t\t[%lu|%lu-%lu|%lu] %c chrom %s\n",
	   g_start, r_start, r_end, g_end, (strand == 0 ? '+' : '-'), 
	   sa_index->genome->chrom_names[chrom]);
  }
}

//--------------------------------------------------------------------

void generate_cals_from_exact_read(int strand, fastq_read_t *read, char *revcomp,
				   size_t low, size_t high, sa_index3_t *sa_index, 
				   array_list_t *cal_list) {
  
  size_t g_start, g_end;
  int chrom;

  cal_t *cal;
  cigar_t *cigar;
  seed_region_t *seed_region;
  linked_list_t *sr_list;
  linked_list_t *sr_duplicate_list;
  
  for (size_t suff = low; suff <= high; suff++) {
    chrom = sa_index->CHROM[suff];
    g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];
    g_end = g_start + read->length - 1;
    sr_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    sr_duplicate_list = linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    seed_region = seed_region_new(0, read->length - 1, g_start, g_end, 0);
    cigar = cigar_new(read->length, 'M');
    seed_region->info = (void *) cigar;
    linked_list_insert(seed_region, sr_list);
    cal = cal_new(chrom, strand, g_start, g_end, 1, sr_list, sr_duplicate_list);
    cal->read_area = read->length;
    cal->num_mismatches = 0;
    array_list_insert(cal, cal_list);
  }
}

//--------------------------------------------------------------------

int generate_cals_from_suffixes(int strand, fastq_read_t *read, char *revcomp,
				int read_pos, int suffix_len, size_t low, size_t high, 
				sa_index3_t *sa_index, cal_mng_t *cal_mng
                                #ifdef _TIMING
				, sa_mapping_batch_t *mapping_batch
                                #endif
				) {

  #ifdef _TIMING
  struct timeval stop, start;
  #endif

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif


  size_t r_start_suf, r_end_suf, g_start_suf, g_end_suf;
  size_t r_start, r_end, r_len, g_start, g_end, g_len;
  int found_cal, chrom, diff, max_map_len = 0, num_mismatches = 0;

  float score;
  alig_out_t alig_out;

  cigar_t *cig, *cigar;
  seed_region_t *seed;

  char *g_seq, *r_seq;
  r_seq = (strand ? revcomp : read->sequence);
  
  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_INIT_CALS_FROM_SUFFIXES] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif

  for (size_t suff = low; suff <= high; suff++) {
    #ifdef _TIMING
    gettimeofday(&start, NULL);
    #endif
    chrom = sa_index->CHROM[suff];

    num_mismatches = 0;
    seed = NULL;

    // extend suffix to right side
    r_start = read_pos;
    r_end = r_start + suffix_len - 1;
    r_len = read->length - r_end - 1;
    
    g_start = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom];
    g_end = g_start + suffix_len - 1;
    g_len = r_len + 5;

    #ifdef _TIMING
    gettimeofday(&stop, NULL);
    mapping_batch->func_times[FUNC_SET_POSITIONS] += 
      ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
    #endif

    // skip suffixes
    #ifdef _TIMING
    gettimeofday(&start, NULL);
    #endif
    found_cal = cal_mng_find(chrom, g_start, g_end, cal_mng);
    #ifdef _TIMING
    gettimeofday(&stop, NULL);
    mapping_batch->func_times[FUNC_SKIP_SUFFIXES] += 
      ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
    #endif
    // found cal...next suffix
    if (found_cal) continue;

    g_start_suf = g_start;
    g_end_suf = g_end;
    r_start_suf = r_start;
    r_end_suf = r_end;

    #ifdef _VERBOSE
    printf("\t\tsuffix at [%lu|%lu-%lu|%lu] %c chrom %s\n",
	   g_start, r_start, r_end, g_end, (strand == 0 ? '+' : '-'), 
	   sa_index->genome->chrom_names[chrom]);
    #endif
    if (r_end >= read->length - 10) {
      // create seed
      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      seed = seed_region_new(r_start, r_end, g_start, g_end, 0);

      // fill the mini-gap
      diff = read->length - seed->read_end - 1;
      g_seq = &sa_index->genome->S[seed->genome_end + sa_index->genome->chrom_offsets[chrom] + 1];
      for (size_t k1 = seed->read_end + 1, k2 = 0; k1 < read->length; k1++, k2++) {
	if (r_seq[k1] != g_seq[k2]) {
	  num_mismatches++;
	}
      }
      seed->read_end += diff;
      seed->genome_end += diff;

      cigar = cigar_new(seed->read_end - seed->read_start + 1, 'M');
      seed->info = (void *) cigar;
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_SEED_REGION_NEW] += 
	  ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif
      #ifdef _VERBOSE
      printf("\t\t\t creating seed (r_end >= read->length - 3)\n");
      print_seed_region("\t\t\t", seed);
      #endif  
    } else {
      // extend to the right side
      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      g_seq = &sa_index->genome->S[g_end + sa_index->genome->chrom_offsets[chrom] + 1];
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_SET_REF_SEQUENCE] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif

      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      score = doscadfun(&r_seq[r_end + 1], r_len, g_seq, g_len, MISMATCH_PERC,
			&alig_out);
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_MINI_SW_RIGHT_SIDE] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif
      if (score > 0 || (alig_out.map_len1 + suffix_len) > 20) {
	num_mismatches = alig_out.mismatch + alig_out.gap_open + alig_out.gap_extend;
	
	// create seed
        #ifdef _TIMING
	gettimeofday(&start, NULL);
        #endif
	seed = seed_region_new(r_start, r_end + alig_out.map_len1, 
			       g_start, g_end + alig_out.map_len2, 0);

	// fill the mini-gap
	diff = read->length - seed->read_end - 1;
	if (diff < 10) {
	  g_seq = &sa_index->genome->S[seed->genome_end + sa_index->genome->chrom_offsets[chrom] + 1];
	  for (size_t k1 = seed->read_end + 1, k2 = 0; k1 < read->length; k1++, k2++) {
	    if (r_seq[k1] != g_seq[k2]) {
	      num_mismatches++;
	    }
	  }
	  seed->read_end += diff;
	  seed->genome_end += diff;
	}

	// it must be improved !!!
	cigar = cigar_new(seed->read_end - seed->read_start + 1, 'M');
	seed->info = (void *) cigar;
        #ifdef _TIMING
	gettimeofday(&stop, NULL);
	mapping_batch->func_times[FUNC_SEED_REGION_NEW] += 
	  ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
        #endif
        #ifdef _VERBOSE
        printf("\t\t\tcreating seed when extending to RIGHT side: score > 0 || (alig_out.map_len1 + suffix_len) > 20\n");
        print_seed_region("\t\t\t", seed);
        #endif
      }
    }

    if (read_pos > 0) {
      // extend suffix to left side
      r_start = 0;
      r_end = read_pos - 1;
      r_len = read_pos;

      g_len = r_len + 5;
      g_end = sa_index->SA[suff] - sa_index->genome->chrom_offsets[chrom] - 1;
      g_start = g_end - g_len;
      //      printf("(g_start, g_end) = (%lu, %lu)\n", g_start, g_end);
      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      g_seq = &sa_index->genome->S[g_start + sa_index->genome->chrom_offsets[chrom] + 1];
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_SET_REF_SEQUENCE] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif

      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      score = doscadfun_inv(r_seq, r_len, g_seq, g_len, MISMATCH_PERC,
			    &alig_out);
			    
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_MINI_SW_LEFT_SIDE] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif
      if (seed) {
	if (score > 0) {
	  // update seed
	  num_mismatches += alig_out.mismatch + alig_out.gap_open + alig_out.gap_extend;
	  seed->read_start -= alig_out.map_len1;
	  seed->genome_start -= alig_out.map_len2;

	  cig = cigar_new(alig_out.map_len1, 'M');
	  cigar_concat(seed->info, cig);
	  cigar_free(seed->info);
	  seed->info = cig;

	  // fill the mini-gap
	  if (seed->read_start > 0 && seed->read_start < 10) {
	    g_seq = &sa_index->genome->S[seed->genome_start + sa_index->genome->chrom_offsets[chrom] - seed->read_start];
	    for (size_t k1 = 0, k2 = 0; k1 < seed->read_start; k1++, k2++) {
	      if (r_seq[k1] != g_seq[k2]) {
		num_mismatches++;
	      }
	    }
	    cig = cigar_new(seed->read_start, 'M');
	    cigar_concat(seed->info, cig);
	    cigar_free(seed->info);
	    seed->info = cig;

	    seed->genome_start -= seed->read_start;
	    seed->read_start = 0;
	  }
	}
      } else {
	if (score > 0 || (alig_out.map_len1 + suffix_len) > 20) {
	  num_mismatches += alig_out.mismatch + alig_out.gap_open + alig_out.gap_extend;
	  // create seed
          #ifdef _TIMING
	  gettimeofday(&start, NULL);
          #endif
	  r_end -= (alig_out.match + alig_out.mismatch - 1);
	  g_end -= (alig_out.match + alig_out.mismatch - 1);
	  seed = seed_region_new(r_end, r_end_suf, g_end, g_end_suf, 0);
	  // it must be improved !!!
	  cig = cigar_new(seed->read_end - seed->read_start + 1, 'M');
	  seed->info = cig;

	  // fill the mini-gap
	  if (seed->read_start < 10) {
	    g_seq = &sa_index->genome->S[seed->genome_start + sa_index->genome->chrom_offsets[chrom] - seed->read_start];
	    for (size_t k1 = 0, k2 = 0; k1 < seed->read_start; k1++, k2++) {
	      if (r_seq[k1] != g_seq[k2]) {
		num_mismatches++;
	      }
	    }
	    cig = cigar_new(seed->read_start, 'M');
	    cigar_concat(seed->info, cig);
	    cigar_free(seed->info);
	    seed->info = cig;

	    seed->genome_start -= seed->read_start;
	    seed->read_start = 0;
	  }

	  cigar = cigar_new(seed->read_end - seed->read_start + 1, 'M');
	  seed->info = (void *) cigar;
          #ifdef _TIMING
	  gettimeofday(&stop, NULL);
	  mapping_batch->func_times[FUNC_SEED_REGION_NEW] += 
	    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
          #endif
          #ifdef _VERBOSE
	  printf("\t\t\t creating seed when extending to LEFT side: score > 0 || (alig_out.map_len1 + suffix_len) > 50\n");
	  print_seed_region("\t\t\t", seed);
          #endif
	}
	//      } else {
	//	num_mismatches += r_len;
      }
    }

    // update CAL manager with this seed
    if (seed) {
      seed->strand = strand;
      seed->chromosome_id = chrom;
      seed->num_mismatches = num_mismatches;
      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      cal_mng_update(seed, cal_mng);
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_CAL_MNG_INSERT] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif
    }
  }

  //  return max_map_len + suffix_len + 1;
  return (sa_index->k_value / 2);
  //  return (sa_index->k_value * 4);
}

//--------------------------------------------------------------------

array_list_t *step_one(fastq_read_t *read, char *revcomp_seq,
		       sa_mapping_batch_t *mapping_batch, 
		       sa_index3_t *sa_index, cal_mng_t *cal_mng) {

  #ifdef _TIMING
  struct timeval stop, start;
  #endif

  int suffix_len, num_suffixes;
  char *r_seq = read->sequence;

  size_t low, high;

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif

  array_list_t *cal_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
  cal_mng->min_read_area = read->length;
  cal_mng->read_length = read->length;

  int max_read_area;// = read->length * MISMATCH_PERC;
  int read_pos, read_end_pos, read_inc = sa_index->k_value / 2;
 
  #ifdef _VERBOSE	  
  printf("\n\n====>>>> STEP ONE <<<<====\n");
  display_cmp_sequences(read, revcomp_seq, sa_index);
  #endif

  // fill in the CAL manager structure
  read_end_pos = read->length - sa_index->k_value;
  //    memset(saved_pos, 0, sizeof(saved_pos));

  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_OTHER] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif

  // first step, searching mappings in both strands
  // distnce between seeds >= prefix value (sa_index->k_value)
  for (int strand = 0; strand < 2; strand++) {
    #ifdef _VERBOSE	  
    printf("=======> STRAND %c\n", (strand == 0 ? '+' : '-'));
    #endif

    for (read_pos = 0; read_pos < read_end_pos; )  {	
      // save this position and search suffixes from this read position
      //	saved_pos[strand][read_pos] = 1;
      #ifdef _VERBOSE	  
      printf("\tread pos. = %lu\n", read_pos);
      #endif

      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      num_suffixes = search_suffix(&r_seq[read_pos], sa_index->k_value, sa_index, 
				   &low, &high, &suffix_len
                                   #ifdef _TIMING
				   , mapping_batch
                                   #endif
				   );
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_SEARCH_SUFFIX] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif
      
      #ifdef _VERBOSE	  
      printf("\t\tnum. suffixes = %lu (suffix length = %lu)\n", num_suffixes, suffix_len);
      #endif
      if (num_suffixes < MAX_NUM_SUFFIXES && suffix_len) {
        #ifdef _VERBOSE	  
	//display_suffix_mappings(strand, read_pos, suffix_len, low, high, sa_index);
        #endif 
	
	// exact search
	if (suffix_len == read->length) {
          #ifdef _TIMING
	  gettimeofday(&start, NULL);
          #endif
	  generate_cals_from_exact_read(strand, read, revcomp_seq,
					low, high, sa_index, cal_list);
          #ifdef _TIMING
	  gettimeofday(&stop, NULL);
	  mapping_batch->func_times[FUNC_CALS_FROM_EXACT_READ] += 
	    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
          #endif
	  break;
	} else {
          #ifdef _TIMING
	  gettimeofday(&start, NULL);
          #endif
	  read_pos += generate_cals_from_suffixes(strand, read, revcomp_seq,
						  read_pos, suffix_len, low, high, sa_index, cal_mng
                                                  #ifdef _TIMING
						  , mapping_batch
                                                  #endif
						  );
          #ifdef _TIMING
	  gettimeofday(&stop, NULL);
	  mapping_batch->func_times[FUNC_CALS_FROM_SUFFIXES] += 
	    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
          #endif
	}
      } else {
	read_pos += read_inc;
      }
    } // end of for read_pos
    
      // update cal list from cal manager
    #ifdef _TIMING
    gettimeofday(&start, NULL);
    #endif
    cal_mng_to_array_list((read->length / 3), cal_list, cal_mng);
    #ifdef _TIMING
    gettimeofday(&stop, NULL);
    mapping_batch->func_times[FUNC_CAL_MNG_TO_LIST] += 
      ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
    #endif
    
    // next, - strand
    r_seq = revcomp_seq;
  } // end of for strand
  
  //  printf("**************** filter min_read_area: = %i, num_cals = %i\n", 
  //	 cal_mng->min_read_area, array_list_size(cal_list));
  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif
  max_read_area = get_max_read_area(cal_list);
  if (max_read_area > read->length) max_read_area = read->length;
  filter_cals_by_max_read_area(max_read_area, &cal_list);
  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_FILTER_BY_READ_AREA] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif

  #ifdef _VERBOSE
  printf("\t***> max_read_area = %i\n", max_read_area);
  #endif
  
  return cal_list;
}

//--------------------------------------------------------------------

void step_two(fastq_read_t *read, char *revcomp_seq,
	      sa_mapping_batch_t *mapping_batch, 
	      sa_index3_t *sa_index, cal_mng_t *cal_mng,
	      array_list_t *valid_cal_list, array_list_t *invalid_cal_list) {
  #ifdef _TIMING
  struct timeval stop, start;
  #endif

  int suffix_len, num_suffixes;
  char *r_seq = read->sequence;

  size_t low, high;

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif

  cal_mng->min_read_area = 100;

  int max_read_area = read->length * MISMATCH_PERC;
  int read_pos, read_end_pos, read_inc = sa_index->k_value / 2;
 
  #ifdef _VERBOSE	  
  printf("\n\n====>>>> STEP TWO <<<<====\n");
  display_cmp_sequences(read, revcomp_seq, sa_index);
  #endif

  // fill in the CAL manager structure
  read_end_pos = read->length - sa_index->k_value;
  //    memset(saved_pos, 0, sizeof(saved_pos));

  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_OTHER] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif

  // first step, searching mappings in both strands
  // distnce between seeds >= prefix value (sa_index->k_value)
  for (int strand = 0; strand < 2; strand++) {
    #ifdef _VERBOSE	  
    printf("=======> STRAND %c\n", (strand == 0 ? '+' : '-'));
    #endif

    for (read_pos = 0; read_pos < read_end_pos; )  {	
      // save this position and search suffixes from this read position
      //	saved_pos[strand][read_pos] = 1;
      #ifdef _VERBOSE	  
      printf("\tread pos. = %lu\n", read_pos);
      #endif

      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      num_suffixes = search_prefix(&r_seq[read_pos], &low, &high, sa_index, 0);
      suffix_len = sa_index->k_value;
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_SEARCH_PREFIX] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif
      
      #ifdef _VERBOSE	  
      printf("\t\tnum. suffixes = %lu (suffix length = %lu)\n", num_suffixes, suffix_len);
      #endif
      if (num_suffixes && num_suffixes < MAX_NUM_SUFFIXES && suffix_len) {
        #ifdef _VERBOSE	  
	//display_suffix_mappings(strand, read_pos, suffix_len, low, high, sa_index);
        #endif 
	
	assert(suffix_len != read->length);

        #ifdef _TIMING
	gettimeofday(&start, NULL);
        #endif
	read_pos += generate_cals_from_suffixes(strand, read, revcomp_seq,
						read_pos, suffix_len, low, high, sa_index, cal_mng
                                                #ifdef _TIMING
						, mapping_batch
                                                #endif
						);
        #ifdef _TIMING
	gettimeofday(&stop, NULL);
	mapping_batch->func_times[FUNC_CALS_FROM_SUFFIXES] += 
	  ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
        #endif
      } else {
	read_pos += read_inc;
      }
    } // end of for read_pos
    
      // update cal list from cal manager
    #ifdef _TIMING
    gettimeofday(&start, NULL);
    #endif
    cal_mng_select_best(max_read_area, valid_cal_list, invalid_cal_list, cal_mng);
    #ifdef _TIMING
    gettimeofday(&stop, NULL);
    mapping_batch->func_times[FUNC_CAL_MNG_TO_LIST] += 
      ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
    #endif
    
    // next, - strand
    r_seq = revcomp_seq;
  } // end of for strand
  /*  
  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif
  filter_cals_by_max_read_area(cal_mng->min_read_area, &cal_list);
  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_FILTER_BY_READ_AREA] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif
  */
}

//--------------------------------------------------------------------

typedef struct cigarset {
  int num_cigars;
  int *active;
  cigar_t **cigars;
} cigarset_t;

cigarset_t *cigarset_new(int num_cigars) {
  cigarset_t *p = (cigarset_t *) malloc(sizeof(cigarset_t));
  p->num_cigars = num_cigars;
  p->active = (int *) calloc(num_cigars, sizeof(int));
  p->cigars = (cigar_t **) malloc(num_cigars * sizeof(cigar_t*));
  return p;
}

void cigarset_free(cigarset_t *p) {
  if (p) {
    if (p->active) free(p->active);
    if (p->cigars) free(p->cigars);
    free(p);
  }
}

//--------------------------------------------------------------------

void step_three(fastq_read_t *read, char *revcomp_seq, 
		sa_mapping_batch_t *mapping_batch, sa_index3_t *sa_index, 
		array_list_t *cal_list) {
  size_t seed_count, num_seeds, num_cals = array_list_size(cal_list);

  #ifdef _VERBOSE	  
  printf("\n\n====>>>> STEP THREE <<<<====\n");
  #endif

  #ifdef _TIMING
  struct timeval stop, start;
  #endif

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif

  char *seq, *ref;
  seed_region_t *prev_seed, *seed;
  int gap_read_len, gap_genome_len;
  size_t gap_read_start, gap_read_end;
  size_t gap_genome_start, gap_genome_end;

  linked_list_item_t *item;
  sw_prepare_t *sw_prepare;
  array_list_t *sw_prepare_list = array_list_new(1000, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  cal_t *cal;
  cigar_t *cigar;
  cigarset_t *cigarset;

  for (int i = 0; i < num_cals; i++) {
    cal = array_list_get(i, cal_list);
    #ifdef _VERBOSE
    cal_print(cal);
    #endif

    // if not seeds, then next cal
    num_seeds = cal->sr_list->size;
    if (num_seeds <= 0) continue;

    // cal cigar
    cigarset = cigarset_new(num_seeds * 2 + 1);
    cal->info = (void *) cigarset;

    // first seed
    seed = linked_list_get_first(cal->sr_list);

    // to remove
    if (((int) seed->read_start) < 0) seed->read_start = 0;

    if (seed->read_start > 0) {
      gap_genome_start = seed->genome_start - seed->read_start - 1;// - 10;
      gap_genome_end = seed->genome_start - 1;
      ref = sa_genome_get_sequence(cal->chromosome_id, gap_genome_start, gap_genome_end, sa_index->genome);
      
      seq = get_subsequence((cal->strand ? revcomp_seq : read->sequence), 
			    0, seed->read_start + 1);
      
      sw_prepare = sw_prepare_new(seq, ref, 0, 0, FIRST_SW);
      sw_prepare->seed_region = 0;
      sw_prepare->cal = cal;
      sw_prepare->read = read;
      array_list_insert(sw_prepare, sw_prepare_list);

      cigarset->active[0] = 1;
    }
    // seeds at the middle positions
    prev_seed = seed;
    cigarset->active[1] = 1;
    cigarset->cigars[1] = seed->info;
    seed->info = NULL;
    seed_count = 0;
    for (item = cal->sr_list->first->next; item != NULL; item = item->next) {
      seed_count++;
      seed = item->item;

      gap_read_start = prev_seed->read_end + 1;
      gap_read_end = seed->read_start - 1;
      gap_read_len = gap_read_end - gap_read_start + 1;

      gap_genome_start = prev_seed->genome_end + 1;
      gap_genome_end = seed->genome_start - 1;
      gap_genome_len = gap_genome_end - gap_genome_start + 1;

      if (gap_read_len <= 0) {
	// deletion
	cigarset->active[seed_count * 2] = 1;
	cigarset->cigars[seed_count * 2] = cigar_new(gap_genome_len, 'D');	
	prev_seed = seed;
	continue;
      }

      if (gap_genome_len <= 0) {
	// insertion
	cigarset->active[seed_count * 2] = 1;
	cigarset->cigars[seed_count * 2] = cigar_new(gap_genome_len, 'I');	
	prev_seed = seed;
	continue;
      }

      #ifdef _VERBOSE1
      print_seed_region("", prev_seed);
      print_seed_region("", seed);
      printf("read id = %s\n", read->id);
      printf("genome gap (start, end) = (%lu, %lu), len = %i\n", gap_genome_start, gap_genome_end, gap_genome_len);
      printf("read gap (start, end) = (%lu, %lu), len = %i\n", gap_read_start, gap_read_end, gap_read_len);
      exit(-1);
      #endif

      assert(gap_read_len > 0);
      seq = get_subsequence((cal->strand ? revcomp_seq : read->sequence), 
			    gap_read_start, gap_read_len);

      assert(gap_genome_len > 0);
      ref = sa_genome_get_sequence(cal->chromosome_id, gap_genome_start, gap_genome_end, sa_index->genome);
      
      
      sw_prepare = sw_prepare_new(seq, ref, 0, 0, MIDDLE_SW);
      sw_prepare->seed_region = seed_count * 2;
      sw_prepare->cal = cal;
      sw_prepare->read = read;
      array_list_insert(sw_prepare, sw_prepare_list);
      cigarset->active[seed_count * 2] = 1;

      prev_seed = seed;
      cigarset->active[seed_count * 2 + 1] = 1;
      cigarset->cigars[seed_count * 2 + 1] = seed->info;
      seed->info = NULL;
    }
    // last seed
    seed = linked_list_get_last(cal->sr_list);
    if (seed->read_end < read->length - 1) {
      gap_genome_start = seed->genome_end + 1;
      gap_genome_end = gap_genome_start + (read->length - seed->read_end);// + 10;
      ref = sa_genome_get_sequence(cal->chromosome_id, gap_genome_start, gap_genome_end, sa_index->genome);

      seq = get_subsequence((cal->strand ? revcomp_seq : read->sequence), 
			    seed->read_end + 1, read->length - seed->read_end - 1);

      sw_prepare = sw_prepare_new(seq, ref, 0, 0, LAST_SW);
      sw_prepare->seed_region = num_seeds * 2;
      sw_prepare->cal = cal;
      sw_prepare->read = read;
      array_list_insert(sw_prepare, sw_prepare_list);

      cigarset->active[num_seeds * 2] = 1;
    }
  }

  // apply smith-waterman
  sw_optarg_t sw_optarg;
  sw_optarg_init(10, 0.5, 5, -4, &sw_optarg);
  
  size_t sw_count = array_list_size(sw_prepare_list); 
  char *q[sw_count], *r[sw_count];
  for (int i = 0; i < sw_count; i++) {
    sw_prepare = array_list_get(i, sw_prepare_list);
    q[i] = sw_prepare->query;
    r[i] = sw_prepare->ref;
    #ifdef _VERBOSE
    printf("\t\t%i: query: %s\n", i, q[i]);
    printf("\t\t%i: ref. : %s\n", i, r[i]);
    printf("\n");
    #endif
  }
  
  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_PRE_SW] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif


  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif

  sw_multi_output_t *sw_output = sw_multi_output_new(sw_count);
  smith_waterman_mqmr(q, r, sw_count, &sw_optarg, 1, sw_output);
  #ifdef _VERBOSE
  sw_multi_output_save(sw_count, sw_output, stdout);
  #endif

  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_SW] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif

  // process Smith-Waterman output
  int op_name, op_value, diff, len;
  for (int i = 0; i < sw_count; i++) {
    sw_prepare = array_list_get(i, sw_prepare_list);
    cal = sw_prepare->cal;
    seed = sw_prepare->seed_region;
    
    cigarset = cal->info;
    cigar = cigar_new_empty();
    op_value = 0;
    op_name = 'M';
    if (sw_output->query_start_p[i] > 0 && sw_output->ref_start_p[i] > 0) {
      diff = sw_output->query_start_p[i] - sw_output->ref_start_p[i];
      if (diff < 0) {
	//	cigar_append_op(abs(diff), 'D', cigar);
	//	op_value = sw_output->query_start_p[i];
	//	op_name = 'M';
      } else if (diff > 0) {
	cigar_append_op(abs(diff), 'I', cigar);
	op_value = sw_output->query_start_p[i];
	op_name = 'M';
      } else {
	op_value = sw_output->query_start_p[i];
	op_name = 'M';
      }
      if (op_value) {
	seq = &q[i][abs(diff)];
	ref = &r[i][abs(diff)];
	for(int j = 0; j < op_value; j++) {
	  if (seq[j] != ref[j]) {
	    cal->num_mismatches++;	
	  }
	}
      }
    } else if (sw_output->query_start_p[i] > 0) {
      op_value = sw_output->query_start_p[i];
      op_name = 'I';
    } else {
      //      op_value = sw_output->ref_start_p[i];
      //      op_name = 'D';
    }

    len = strlen(sw_output->query_map_p[i]);
    for(int j = 0; j < len; j++) {
      if (sw_output->query_map_p[i][j] == '-') {
	// deletion (in the query)
	if (op_name != 'D' && op_value > 0) {
	  cigar_append_op(op_value, op_name, cigar);
	  op_value = 0;
	  op_name = 'D';
	}
	op_value++;
	cal->num_mismatches++;	
      } else if (sw_output->ref_map_p[i][j] == '-') {
	// insertion (in the query)
	if (op_name != 'I' && op_value > 0) {
	  cigar_append_op(op_value, op_name, cigar);
	  op_value = 0;
	  op_name = 'I';
	}
	op_value++;
	cal->num_mismatches++;	
      } else {
	if (op_name != 'M' && op_value > 0) {
	  cigar_append_op(op_value, op_name, cigar);
          #ifdef _VERBOSE
	  printf("****************** (value, name) = (%i, %c) : op. change: current cigar: %s\n", 
		 op_value, op_name, cigar_to_string(cigar));
          #endif
	  op_value = 0;
	  op_name = 'M';
	}
	op_value++;
	if (sw_output->query_map_p[i][j] != sw_output->ref_map_p[i][j]) {
	  cal->num_mismatches++;	
	}
      }
    }
    if (op_value > 0) {
      cigar_append_op(op_value, op_name, cigar);
    }
    cigarset->active[(int)sw_prepare->seed_region] = 1;
    cigarset->cigars[(int)sw_prepare->seed_region] = cigar;
    #ifdef _VERBOSE
    printf("************** sw_count: %i, cigar for gap %i: %s\n", 
	   i, sw_prepare->seed_region, cigar_to_string(cigar));
    #endif

    // free
    sw_prepare_free(sw_prepare);
  }

  for (int i = 0; i < num_cals; i++) {
    cal = array_list_get(i, cal_list);
    cigar = cigar_new_empty();
    cigarset = cal->info;
    for (int j = 0; j < cigarset->num_cigars; j++) {
      if (cigarset->active[j]) {
	#ifdef _VERBOSE
	printf("************** gap %i -> concat cigar %s into %s\n",
	       j, cigar_to_string(cigarset->cigars[j]), cigar_to_string(cigar));
	#endif
	cigar_concat(cigarset->cigars[j], cigar);
	cigar_free(cigarset->cigars[j]);
      }
    }
    cigarset_free(cigarset);
    cal->info = cigar;
  }
  // free memory
  sw_multi_output_free(sw_output);
  array_list_free(sw_prepare_list, (void *) NULL);

  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_POST_SW] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif
}

//--------------------------------------------------------------------

void filter_min_read_area_cals(array_list_t **list) {
  array_list_t *cal_list = *list;

  size_t num_cals = array_list_size(cal_list);
  if (num_cals <= 0) return;

  cal_t *cal;
  int min_read_area = 10000, invalid[num_cals], num_invalid = 0;
  
  printf(">>>> filter_min_read_area_cals:\n");
  for(size_t i = 0; i < num_cals; i++) {
    cal = array_list_get(i, cal_list);

    cal_print(cal);

    if (cal->read_area < min_read_area) {
      min_read_area = cal->read_area;
    }
  }
  printf(">>>> filter_min_read_area_cals: min_read_area = %i\n", min_read_area);

  for(size_t i = 0; i < num_cals; i++) {
    invalid[i] = 0;
    cal = array_list_get(i, cal_list);
    if (cal->read_area > min_read_area) {
      num_invalid++;
      invalid[i] = 1;
    }
  }

  if (num_invalid) {
    array_list_t *new_cal_list = array_list_new(num_cals - num_invalid, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
    for (size_t j = 0; j < num_cals; j++) {
      cal = array_list_get(j, cal_list);
      if (!invalid[j]) {
	array_list_insert(cal, new_cal_list);
	array_list_set(j, NULL, cal_list);
      }
    }
    array_list_free(cal_list, (void *) cal_free);
    *list = new_cal_list;
  }
}

//--------------------------------------------------------------------
// sa mapper
//--------------------------------------------------------------------

int sa_mapper(void *data) {

  #ifdef _TIMING
  struct timeval stop, start;
  #endif

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif
  
  sa_wf_batch_t *wf_batch = (sa_wf_batch_t *) data;
  
  sa_mapping_batch_t *mapping_batch = wf_batch->mapping_batch;
  sa_index3_t *sa_index = (sa_index3_t *) wf_batch->sa_index;
 
  size_t num_reads = mapping_batch->num_reads;
  int max_read_area, min_num_mismatches;
  
  // CAL management
  size_t num_cals;
  cal_t *cal;
  cal_mng_t *cal_mng;
  array_list_t *cal_list;

  fastq_read_t *read;

  //  int saved_pos[2][1024];
  // TODO !!! 20 = min. cal size
  uint min_cal_size = 20;
  cal_mng = cal_mng_new(sa_index->genome);
  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_OTHER] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif

  // for each read
  for (int i = 0; i < num_reads; i++) {
    read = array_list_get(i, mapping_batch->fq_reads);
    //printf("read id = %s\n", read->id);
    //max_num_mismatches = read->length * MISMATCH_PERC;

    // step one:: extend using mini-sw from suffix
    cal_list = step_one(read, mapping_batch->revcomp_seqs[i], mapping_batch, sa_index, cal_mng);

    if (array_list_size(cal_list) > 0) {
      min_num_mismatches = get_min_num_mismatches(cal_list);
      #ifdef _VERBOSE
      printf("\t*** before SW> min_num_mismatches = %i\n", min_num_mismatches);
      #endif
      
      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      filter_cals_by_max_num_mismatches(min_num_mismatches, &cal_list);
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_FILTER_BY_READ_AREA] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif

      // step three (Smith-Waterman to fill in the gaps)
      step_three(read, mapping_batch->revcomp_seqs[i], mapping_batch, sa_index, cal_list);

      min_num_mismatches = get_min_num_mismatches(cal_list);
      #ifdef _VERBOSE
      printf("\t*** after SW> min_num_mismatches = %i\n", min_num_mismatches);
      #endif
      
      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      filter_cals_by_max_num_mismatches(min_num_mismatches, &cal_list);
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_FILTER_BY_READ_AREA] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif
    }

    //    printf("******************* cal list size = %i\n", array_list_size(cal_list));
    //    exit(-1);
    /*
    if ((num_cals = array_list_size(cal_list)) <= 0) {
      // step two: extend using mini-sw from prefix
      array_list_t *invalid_cal_list = array_list_new(MAX_CALS, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);
      step_two(read, mapping_batch->revcomp_seqs[i], mapping_batch, sa_index, cal_mng, 
	       cal_list, invalid_cal_list);

      //      printf("********************** after step_two num. valid = %i, invalid = %i\n", 
      //	     array_list_size(cal_list), array_list_size(invalid_cal_list));

      if (array_list_size(cal_list) > 0) {
	array_list_free(invalid_cal_list, (void *) cal_free);
      } else {
	// step three: apply smith-waterman
	array_list_free(cal_list, (void *) NULL);	
	//step_three(read, mapping_batch->revcomp_seqs[i], mapping_batch, sa_index, invalid_cal_list);
	cal_list = invalid_cal_list;

        #ifdef _TIMING
	gettimeofday(&start, NULL);
        #endif
	filter_cals_by_max_read_area(max_read_area, &cal_list);
        #ifdef _TIMING
	gettimeofday(&stop, NULL);
	mapping_batch->func_times[FUNC_FILTER_BY_READ_AREA] += 
	  ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
        #endif
      }
    */

    /*
    if (array_list_size(cal_list) > 0) {
      #ifdef _TIMING
      gettimeofday(&start, NULL);
      #endif
      filter_min_read_area_cals(&cal_list);
      #ifdef _TIMING
      gettimeofday(&stop, NULL);
      mapping_batch->func_times[FUNC_FILTER_BY_READ_AREA] += 
	((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
      #endif
    }
    */
    // fill in gaps on CALs and create alignments structures
    #ifdef _TIMING
    gettimeofday(&start, NULL);
    #endif
    create_alignments(cal_list, read, mapping_batch->mapping_lists[i]);
    #ifdef _TIMING
    gettimeofday(&stop, NULL);
    mapping_batch->func_times[FUNC_CREATE_ALIGNMENTS] += 
      ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
    #endif
    
    // free cal list and clear seed manager for next read
    array_list_free(cal_list, (void *) NULL);
  } // end of for reads

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif
  cal_mng_free(cal_mng);
  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  mapping_batch->func_times[FUNC_OTHER] += 
    ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif
  
  return -1;
}

//--------------------------------------------------------------------

bam_header_t *create_bam_header(sa_genome3_t *genome) {

  bam_header_t *bam_header = (bam_header_t *) calloc(1, sizeof(bam_header_t));

  int num_targets = genome->num_chroms;

  bam_header->n_targets = num_targets;
  bam_header->target_name = (char **) calloc(num_targets, sizeof(char *));
  bam_header->target_len = (uint32_t*) calloc(num_targets, sizeof(uint32_t));
  for (int i = 0; i < num_targets; i++) {
    bam_header->target_name[i] = strdup(genome->chrom_names[i]);
    bam_header->target_len[i] = genome->chrom_lengths[i] + 1;
  }
  bam_header->text = strdup("@PG\tID:HPG-Aligner\tVN:1.0\n");
  bam_header->l_text = strlen(bam_header->text);

  return bam_header;
}

//--------------------------------------------------------------------

void *write_sam_header(sa_genome3_t *genome, FILE *f) {
  fprintf(f, "@PG\tID:HPG-Aligner\tVN:2.0\n");
  for (int i = 0; i < genome->num_chroms; i++) {
    fprintf(f, "@SQ\tSN:%s\tLN:%lu\n", genome->chrom_names[i], genome->chrom_lengths[i] + 1);
  }
}

//--------------------------------------------------------------------
// main 
//--------------------------------------------------------------------

void dna_aligner(options_t *options) {
  if (argc != 6) {
    printf("Usage: %s <sa-dirname> <fastq-filename> <output-filename-bam-or-sam> <batch-size> <num-threads>\n", argv[0]);
    exit(-1);
  }

  #ifdef _TIMING
  for (int i = 0; i < NUM_TIMING; i++) {
    func_times[i] = 0;
  }
  #endif

  char *sa_dirname = argv[1];
  char *fastq_filename = argv[2];
  char *bam_filename = argv[3];
  int batch_size = atoi(argv[4]);
  int num_threads = atoi(argv[5]);
  int sam = 0;

  genome_t *genome = NULL;

  // load SA index
  struct timeval stop, start;
  printf("\n");
  printf("Loading SA tables...\n");
  gettimeofday(&start, NULL);
  sa_index3_t *sa_index = sa_index3_new(sa_dirname);
  gettimeofday(&stop, NULL);
  sa_index3_display(sa_index);
  printf("End of loading SA tables in %0.2f s. Done!!\n", 
	 (stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  

  // preparing input FastQ file
  fastq_batch_reader_input_t reader_input;
  fastq_batch_reader_input_init(fastq_filename, NULL, 
				0, batch_size, 
				NULL, &reader_input);
  
  reader_input.fq_file1 = fastq_fopen(fastq_filename);
  reader_input.fq_file2 = NULL;
  
  // preparing output BAM file
  batch_writer_input_t writer_input;
  batch_writer_input_init(bam_filename, NULL, NULL, NULL, NULL, &writer_input);
  
  if (strstr(bam_filename, ".sam")) {
    sam = 1;
    writer_input.bam_file = (bam_file_t *) fopen(bam_filename, "w");    
    write_sam_header(sa_index->genome, writer_input.bam_file);
  } else {
    sam = 0;
    bam_header_t *bam_header = create_bam_header(sa_index->genome);
    writer_input.bam_file = bam_fopen_mode(bam_filename, bam_header, "w");
    bam_fwrite_header(bam_header, writer_input.bam_file);
  }
  
  //--------------------------------------------------------------------------------------
  // workflow management
  //
  sa_wf_batch_t *wf_batch = sa_wf_batch_new((void *)sa_index, &writer_input, NULL);
  sa_wf_input_t *wf_input = sa_wf_input_new(sam, &reader_input, wf_batch);

  // create and initialize workflow
  workflow_t *wf = workflow_new();
  
  workflow_stage_function_t stage_functions[] = {sa_mapper};
  char *stage_labels[] = {"SA mapper"};
  workflow_set_stages(1, &stage_functions, stage_labels, wf);
  
  // optional producer and consumer functions
  workflow_set_producer(sa_fq_reader, "FastQ reader", wf);
  workflow_set_consumer(sa_sam_writer, "SAM writer", wf);

  workflow_run_with(num_threads, wf_input, wf);
  
  printf("----------------------------------------------\n");
  workflow_display_timing(wf);
  printf("----------------------------------------------\n");

  #ifdef _TIMING
  char func_name[1024];
  double total_func_times = 0;
  for (int i = 0; i < NUM_TIMING; i++) {
    if (i != FUNC_SEARCH_PREFIX && i != FUNC_SEARCH_SA 
	&& i < FUNC_INIT_CALS_FROM_SUFFIXES || i > FUNC_CAL_MNG_INSERT) {
      total_func_times += func_times[i];
    }
  }
  printf("Timing in seconds:\n");
  for (int i = 0; i < NUM_TIMING; i++) {
    if (i == FUNC_SEARCH_PREFIX || i == FUNC_SEARCH_SA ||
	(i >= FUNC_INIT_CALS_FROM_SUFFIXES && i <= FUNC_CAL_MNG_INSERT)) {
      printf("\t");
    }
    printf("\t%0.2f %%\t%0.4f\tof %0.4f\t%s\n", 
	   100.0 * func_times[i] / total_func_times, func_times[i], total_func_times, func_names[i]);
  }
  #endif

  //  printf("Total num. mappings: %u\n", total_num_mappings);

  // free memory
  workflow_free(wf);
  sa_wf_input_free(wf_input);
  sa_wf_batch_free(wf_batch);
  //
  // end of workflow management
  //--------------------------------------------------------------------------------------

  // free memory
  if (sa_index) sa_index3_free(sa_index);
  
  //closing files
  fastq_fclose(reader_input.fq_file1);
  if (sam) {
    fclose((FILE *) writer_input.bam_file);
  } else {
    bam_fclose(writer_input.bam_file);
  }

  #ifdef _VERBOSE
  printf("*********> num_dup_reads = %i, num_total_dup_reads = %i\n", 
	 num_dup_reads, num_total_dup_reads);
  #endif

}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
