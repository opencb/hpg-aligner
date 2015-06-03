#ifndef SA_BWT_SUPPORT_H
#define SA_BWT_SUPPORT_H

#include <stdio.h>
#include <stdlib.h>

#include "bwt_server.h"

#include "sa/sa_index3.h"
#include "sa/sa_search.h"

//--------------------------------------------------------------------
// generic functions to support SA or BWT search
//--------------------------------------------------------------------

typedef int (*generic_function_t) (void *);

typedef struct search_support {
  char *dirname;
  void *index;  
  void *genome;
  generic_function_t *load_index;
  generic_function_t *free_index;
} search_support_t;

void search_support_free(search_support_t *p);

//--------------------------------------------------------------------
// SA search support
//--------------------------------------------------------------------

int sa_load_index(void *data);
int sa_free_index(void *data);

search_support_t *sa_search_support_new(char *dirname);

//--------------------------------------------------------------------
// BWT search support
//--------------------------------------------------------------------

int bwt_load_index(void *data);
int bwt_free_index(void *data);

search_support_t *bwt_search_support_new(char *dirname);

//--------------------------------------------------------------------
//--------------------------------------------------------------------
#endif // SA_BWT_SUPPORT_H
