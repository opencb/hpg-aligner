#ifndef _SA_MAPPER_STAGE_H
#define _SA_MAPPER_STAGE_H

#include "adapter.h"

#include "sa/sa_search.h"
#include "dna/sa/sa_dna_commons.h"
#include "dna/doscadfun.h"

#include "dna/suffix_mng.h"
#include "dna/cal_mng.h"

//--------------------------------------------------------------------
// sa mapper
//--------------------------------------------------------------------

int sa_single_mapper(void *data);
int sa_pair_mapper(void *data);

array_list_t *step_one(fastq_read_t *read, char *revcomp_seq,
		       sa_mapping_batch_t *mapping_batch, 
		       sa_index3_t *sa_index, cal_mng_t *cal_mng);

//--------------------------------------------------------------------

void fastq_read_revcomp(fastq_read_t *read);

//--------------------------------------------------------------------
//--------------------------------------------------------------------
#endif // _SA_MAPPER_STAGE_H
