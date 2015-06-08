#ifndef _INDEX_H
#define _INDEX_H

#include "options.h"

#include "bwt_server.h"
#include "sa/sa_index3.h"

//#include "dna/sa/sa_io_stages.h"

//--------------------------------------------------------------------

#define SA_INDEX  0
#define BWT_INDEX 1

//--------------------------------------------------------------------
// index_t struct
//--------------------------------------------------------------------

typedef struct index {
  int mode;
  sa_index3_t *sa_index;
  bwt_index_t *bwt_index;
  genome_t *genome;
} index_t;

index_t *index_new(char *dirname);
void index_free(index_t *index);

bam_header_t *index_create_bam_header(options_t *options, index_t *index);
void index_write_sam_header(options_t *options, FILE *file, index_t *index);

//--------------------------------------------------------------------
//--------------------------------------------------------------------
#endif // _INDEX_H
