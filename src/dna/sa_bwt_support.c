#include "dna/sa_bwt_support.h"

//--------------------------------------------------------------------
// generic functions to support SA or BWT search
//--------------------------------------------------------------------

void search_support_free(search_support_t *p) {
  if (p) {
    if (p->dirname) free(p->dirname);
    free(p);
  }
}

//--------------------------------------------------------------------
// SA search support
//--------------------------------------------------------------------

int sa_load_index(void *data) {
  search_support_t *p = (search_support_t *) data;
  if (p) {
    p->index = (void *) sa_index3_new(p->dirname);
    p->genome = ((sa_index3_t *) p->index)->genome;
  }
}

//--------------------------------------------------------------------

int sa_free_index(void *data) {
  search_support_t *p = (search_support_t *) data;
  
  if (p && p->index) {
      sa_index3_free((sa_index3_t *) p->index);
  }
}

//--------------------------------------------------------------------

int sa_get_seed_size(search_support_t *p) {
  return ((sa_index3_t *) p->index)->k_value;
}

int search_support_get_seed_size(search_support_t *p) {
  return p->get_seed_size(p)
}

//--------------------------------------------------------------------

search_support_t *sa_search_support_new(char *dirname) {
  search_support_t *p = (search_support_t *) calloc(1, sizeof(search_support_t));

  p->dirname = strdup(dirname);
  p->index = NULL;
  p->genome = NULL;

  p->load_index = sa_load_index;
  p->free_index = sa_free_index;
  p->get_seed_size = sa_get_seed_size;

  return p;
}

//--------------------------------------------------------------------
// BWT search support
//--------------------------------------------------------------------

int bwt_load_index(void *data) {
  search_support_t *p = (search_support_t *) data;
  if (p) {
    p->index = (void *) bwt_index_new(p->dirname, false);
    p->genome = (void *) genome_new("dna_compression.bin", p->dirname, BWT_MODE); 
  }
}

//--------------------------------------------------------------------

int bwt_free_index(void *data) {
  search_support_t *p = (search_support_t *) data;
  
  if (p) {
    if (p->index) bwt_index_free((bwt_index_t *) p->index);
    if (p->genome) genome_free((genome_t *) p->genome);
  }
}

//--------------------------------------------------------------------

search_support_t *bwt_search_support_new(char *dirname) {
  search_support_t *p = (search_support_t *) calloc(1, sizeof(search_support_t));

  p->dirname = strdup(dirname);
  p->index = NULL;
  p->genome = NULL;

  p->load_index = bwt_load_index;
  p->free_index = bwt_free_index;

  return p;
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
