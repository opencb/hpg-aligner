#include "sa/sa_search.h"

//--------------------------------------------------------------------

size_t search_prefix(char *sequence, size_t *low, size_t *high, 
		     sa_index3_t *sa_index, int display) {
  size_t num_mappings = 0;

  char *seq;
  size_t value, row, col;
  uint ia, ia1, ia2, found_ia;
  uint ja, found_ja;
  uint a, a1, a2;

  //  printf("prefix: %s\n", sequence);
  //  display_prefix(sequence, sa_index->k_value);
  value = compute_prefix_value(sequence, sa_index->k_value);
  row = value >> 8;
  col = 255LLU & value;
  //  printf(" -> prefix value = %lu -> (row, col) = (%lu, %lu)\n", value, row, col); 
  
  ia1 = sa_index->IA[row];
  //  printf("\tIA[%lu] = %lu\n", row, ia1);
  if (ia1 == max_uint) {
    //    printf("\t\t\t----> ROW NOT FOUND !!!\n");
    return num_mappings;
  }
  
  size_t row2 = row + 1;
  //if (row2 >= sa_index->IA_items) {
  //row2 = sa_index->IA_items;
  //ia2 = sa_index->A_items;
  //} else {
    while ((ia2 = sa_index->IA[row2]) == max_uint) {
      row2++;
      if (row2 >= sa_index->IA_items) {
	//      printf("\trow2 reaches IA_items limit (IA items = %lu)\n", sa_index->IA_items);
	row2 = sa_index->IA_items;
	ia2 = sa_index->A_items;
	break;
      }
    }
    //}
  //  printf("\tfrom IA[%lu] = %lu to IA[%lu] = %lu -> num. columns = %lu\n", row, ia1, row2, ia2, ia2 - ia1);
  
  found_ja = 0;

  int idx, size = ia2 - ia1;
  for (ia = ia1; ia < ia2; ia++) {
    a1 = sa_index->A[ia];
    ja = sa_index->JA[ia];
    //      printf("\t\tA[%lu] = %lu\t JA[%lu] = %lu\n", 
    //      	     ia, a1, ia, ja);
    if (sa_index->JA[ia] == col) {
      found_ja = 1;
      if (ia + 1 >= sa_index->A_items) {
	a2 = sa_index->num_suffixes;
      } else {
	a2 = sa_index->A[ia + 1];
      }
      break;
    }
  }

  if (found_ja) {
    num_mappings = a2 - a1;
    *low = a1;
    *high = a2;
  }

  return num_mappings;
}

//--------------------------------------------------------------------

size_t search_suffix(char *seq, uint len, int max_num_suffixes,
		     sa_index3_t *sa_index, 
		     size_t *low, size_t *high, size_t *suffix_len
                     #ifdef _TIMING
		     , double *prefix_time, double *suffix_time
                     #endif
		     ) {
  #ifdef _TIMING
  struct timeval stop, start;
  #endif

  int display = 1;
  char *ref, *query;
  size_t num_suffixes = 0;
  uint matched, max_matched = 0;

  *suffix_len = 0;

  #ifdef _TIMING
  gettimeofday(&start, NULL);
  #endif
  size_t num_prefixes = search_prefix(seq, low, high, sa_index, display);
  #ifdef _TIMING
  gettimeofday(&stop, NULL);
  *prefix_time = ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
  #endif

  #ifdef _VERBOSE	  
  printf("\t\tnum. prefixes = %lu\n", num_prefixes);
  #endif

  if (num_prefixes && num_prefixes < max_num_suffixes) {   

    //    *suffix_len = sa_index->k_value;
    //    return num_prefixes;
    
    #ifdef _TIMING
    gettimeofday(&start, NULL);
    #endif


    #ifdef _VERBOSE1	  
    {
      char *ss = get_subsequence(seq + sa_index->k_value, 0, 40);
      printf("\tquery:\t%s\n", ss);
      free(ss);
      for (size_t i = *low; i < *high; i++) {
	printf("\t%lu\t", i);
	ref = &sa_index->genome->S[sa_index->SA[i]] + sa_index->k_value;
	char *ss = get_subsequence(ref, 0, 40);
	printf("%s\n", ss);
	free(ss);
      }
    }
    #endif


    size_t first = *low, last = *low;

    if (num_prefixes == 1) {
      query = seq + sa_index->k_value;
      ref = &sa_index->genome->S[sa_index->SA[*low]] + sa_index->k_value;
      matched = 0;
      while (query[matched] == ref[matched]) {
	matched++;
      }
      *high = *low;
      *suffix_len = matched + sa_index->k_value;
      return num_prefixes;
    }

    
    for (size_t i = *low; i < *high; i++) {
      query = seq + sa_index->k_value;
      ref = &sa_index->genome->S[sa_index->SA[i]] + sa_index->k_value;
      matched = 0;
      while (query[matched] == ref[matched]) {
	matched++;
      }
      if (matched > max_matched) {
	first = i;
	last = i;
	max_matched = matched;
	//	break;
      } else if (matched == max_matched) {
	last = i;
      } else {
	break;
      }
    }
    

    if (first <= last) {
      *low = first;
      *high = last;
      *suffix_len = max_matched + sa_index->k_value;
      num_suffixes = last - first + 1;
    }

    #ifdef _TIMING
    gettimeofday(&stop, NULL);
    *suffix_time = ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
    #endif
  }
  
  return num_suffixes;
}


//--------------------------------------------------------------------
//--------------------------------------------------------------------
