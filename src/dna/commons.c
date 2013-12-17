#include "commons.h"

//--------------------------------------------------------------------
// structs
//--------------------------------------------------------------------

wf_batch_t *wf_batch_new(void *index, void *optarg,
			 batch_writer_input_t *writer_input,
			 mapping_batch_t *mapping_batch) {
  wf_batch_t *p = (wf_batch_t *) malloc(sizeof(wf_batch_t));
  p->index = index;
  p->optarg = optarg;
  p->writer_input = writer_input;
  p->mapping_batch = mapping_batch;
  return p;
}

void wf_batch_free(wf_batch_t *p) {
  if (p) free(p);
}

//--------------------------------------------------------------------

wf_in_t *wf_in_new(fastq_batch_reader_input_t *fq_reader_input,
		   wf_batch_t *wf_batch) {
  wf_in_t *p = (wf_in_t *) malloc(sizeof(wf_in_t));
  p->fq_reader_input = fq_reader_input;
  p->wf_batch = wf_batch;
  return p;
}

void wf_in_free(wf_in_t *p) {
  if (p) free(p);
}


//--------------------------------------------------------------------
//--------------------------------------------------------------------
