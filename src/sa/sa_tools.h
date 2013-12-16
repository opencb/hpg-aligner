#ifndef SA_TOOLS_H
#define SA_TOOLS_H

#include <stdio.h>

#include "sa_tools.h"

//--------------------------------------------------------------------------------------

typedef unsigned int uint;

extern size_t PREFIX_TABLE_K_VALUE;
extern size_t PREFIX_TABLE_NT_VALUE[256];

//--------------------------------------------------------------------------------------

char *read_s(char *filename, uint *len);
void compute_sa(uint *sa, uint num_sufixes);
void compute_lcp(char *s, uint *sa, uint *lcp, uint num_sufixes);
void compute_child(uint *lcp, int *child, uint num_sufixes);

//--------------------------------------------------------------------------------------

size_t compute_prefix_value(char *prefix, int len);

//--------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------

#endif // SA_TOOLS_H
