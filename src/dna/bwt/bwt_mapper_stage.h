#ifndef _BWT_MAPPER_STAGE_H
#define _BWT_MAPPER_STAGE_H

#include "adapter.h"

#include "sa/sa_search.h"
#include "dna/sa_dna_commons.h"
#include "dna/doscadfun.h"

#include "dna/suffix_mng.h"

//--------------------------------------------------------------------
// bwt mapper
//--------------------------------------------------------------------

int bwt_single_mapper(void *data);
int bwt_pair_mapper(void *data);

//--------------------------------------------------------------------

void fastq_read_revcomp(fastq_read_t *read);

//--------------------------------------------------------------------
//--------------------------------------------------------------------
#endif // _BWT_MAPPER_STAGE_H
