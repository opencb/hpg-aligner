#ifndef _SA_IO_STAGES_H
#define _SA_IO_STAGES_H

#include "sa/sa_index3.h"

#include "batch_writer.h"

#include "dna/sa_dna_commons.h"

//--------------------------------------------------------------------

void *sa_fq_reader(void *input);
void *sa_bam_reader_single(void *input);
void *sa_bam_reader_pairend(void *input);
void *sa_bam_reader_unmapped(void *input);

//--------------------------------------------------------------------

int sa_sam_writer(void *data);
void write_sam_header(sa_genome3_t *genome, FILE *f);

//--------------------------------------------------------------------

bam_header_t *create_bam_header(sa_genome3_t *genome);
int sa_bam_writer(void *data);

//--------------------------------------------------------------------
//--------------------------------------------------------------------

#endif // _SA_IO_STAGES_H
