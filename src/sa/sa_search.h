#ifndef _SA_SEARCH_H
#define _SA_SEARCH_H

#include "sa/sa_tools.h"
#include "sa/sa_index3.h"

//--------------------------------------------------------------------

size_t search_prefix(char *sequence, size_t *low, size_t *high, 
		     sa_index3_t *sa_index, int display);

size_t search_suffix(char *seq, uint len, int max_num_suffixes,
		     sa_index3_t *sa_index, 
		     size_t *low, size_t *high, size_t *suffix_len
                     #ifdef _TIMING
		     , double *prefix_time, double *suffix_time
                     #endif
		     );

//--------------------------------------------------------------------
//--------------------------------------------------------------------
#endif // _SA_SEARCH_H

