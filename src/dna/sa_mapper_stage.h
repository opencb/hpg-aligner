#ifndef _SA_MAPPER_STAGE_H
#define _SA_MAPPER_STAGE_H

#include "sa/sa_search.h"
#include "dna/sa_dna_commons.h"
#include "dna/doscadfun.h"

//--------------------------------------------------------------------
// cal_mng_t struct
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

cal_mng_t * cal_mng_new(sa_genome3_t *genome);
void cal_mng_free(cal_mng_t *p);
void cal_mng_simple_free(cal_mng_t *p);
void cal_mng_simple_clear(cal_mng_t *p);
void cal_mng_clear(cal_mng_t *p);

void cal_mng_update(seed_t *seed, fastq_read_t *read, cal_mng_t *p);
int cal_mng_find(int strand, int chrom, size_t start, size_t end, cal_mng_t *p);

void cal_mng_to_array_list(int read_area, array_list_t *out_list, cal_mng_t *p);
void cal_mng_select_best(int read_area, array_list_t *valid_list, 
			 array_list_t *invalid_list, cal_mng_t *p);

//--------------------------------------------------------------------
// sa mapper
//--------------------------------------------------------------------

int sa_single_mapper(void *data);
int sa_pair_mapper(void *data);

array_list_t *step_one(fastq_read_t *read, char *revcomp_seq,
		       sa_mapping_batch_t *mapping_batch, 
		       sa_index3_t *sa_index, cal_mng_t *cal_mng);

//--------------------------------------------------------------------
//--------------------------------------------------------------------
#endif // _SA_MAPPER_STAGE_H
