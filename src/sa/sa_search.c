#include "sa/sa_search.h"

//--------------------------------------------------------------------

size_t search_prefix(char *sequence, size_t *low, size_t *high, 
		     sa_index3_t *sa_index, int display) {
  size_t num_mappings = 0;


  size_t value, row, col;
  uint ia, ia1, ia2;
  uint  found_ja;
  uint a1, a2;

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
  if (row2 >= sa_index->IA_items) {
    row2 = sa_index->IA_items;
    ia2 = sa_index->A_items;
  } else {
    while ((ia2 = sa_index->IA[row2]) == max_uint) {
      row2++;
      if (row2 >= sa_index->IA_items) {
	//      printf("\trow2 reaches IA_items limit (IA items = %lu)\n", sa_index->IA_items);
	row2 = sa_index->IA_items;
	ia2 = sa_index->A_items;
	break;
      }
    }
  }
  //}
  //  printf("\tfrom IA[%lu] = %lu to IA[%lu] = %lu -> num. columns = %lu\n", row, ia1, row2, ia2, ia2 - ia1);
  
  found_ja = 0;


  for (ia = ia1; ia < ia2; ia++) {
    a1 = sa_index->A[ia];
    //ja = sa_index->JA[ia];
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

  *suffix_len = 0;
  num_suffixes = num_prefixes;

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
      num_suffixes = num_prefixes;
    } else {
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
    }

    #ifdef _TIMING
    gettimeofday(&stop, NULL);
    *suffix_time = ((stop.tv_sec - start.tv_sec) + (stop.tv_usec - start.tv_usec) / 1000000.0f);  
    #endif
  }

  
  return num_suffixes;

}

//--------------------------------------------------------------------
      //char *ref_x = &sa_index->genome->S[sa_index->SA[*low]] + sa_index->k_value;

      ////////////////////////////////////////
      //-----------------------------------------------------------------------------------
      /*
      size_t s, e;
      
      s = start_g;
      e = end_g;
      
      size_t group_start = s/4;
      size_t group_end   = e/4;
      
      size_t nucleotide_start = s%4;
      size_t nucleotide_end   = e%4;
      
      unsigned int seq_pos = 0;
      char *str_tmp;
      
      printf("%lu - %lu, %lu - %lu, %lu, %lu\n", group_start, group_end, nucleotide_start, nucleotide_end, strlen(seq), end_g - start_g);
      
      //printf("Get sequence 1\n");
      str_tmp = genome->code_table[genome->X[group_start]];
      for(unsigned int i = nucleotide_start; i < 4; i++){
	sequence[seq_pos++] = str_tmp[i];
      }
      group_start++;
      
      //printf("Get sequence 2\n");
      while (group_start < group_end) {
	str_tmp = genome->code_table[genome->X[group_start]];
	sequence[seq_pos++] = str_tmp[0];
	sequence[seq_pos++] = str_tmp[1];
	sequence[seq_pos++] = str_tmp[2];
	sequence[seq_pos++] = str_tmp[3];
	group_start++;
      }
      
      //printf("Get sequence 3\n");
      if (group_start <= group_end) {  
	str_tmp = genome->code_table[genome->X[group_start]];
	for(unsigned int i = 0; i <= nucleotide_end; i++){
	  sequence[seq_pos++] = str_tmp[i];
	}
      }
      
      printf("%i\n", seq_pos);
      sequence[seq_pos] = '\0';

      //-----------------------------------------------------------------------------------

      char *ref_n = &sa_index->genome->S[sa_index->SA[*low]] + sa_index->k_value;
      char r[1024];
      strncpy(r, ref_n, strlen(seq) + 1);
      r[strlen(seq) + 1] = '\0';

      int chrom;
      for (int c = 0; c < genome->num_chromosomes; c++) {
	if (start_g >= genome->chr_offset[c]) {
	  chrom = c;
	}
      }
      
      size_t g_start = start_g - sa_index->genome->chrom_offsets[chrom];
      
      if (strcmp(sequence, r) != 0) {
	printf("[%i]%lu - %lu : %s\n", chrom, g_start, g_start + strlen(seq), sequence);
	printf("[%i]%lu - %lu : %s\n", chrom, g_start, g_start + strlen(seq), r);
	exit(-1);
      }
      */
      ///////////////////////////////////////////

size_t search_suffix_rna(char *seq, uint len, int max_num_suffixes,
			 sa_index3_t *sa_index,  size_t *low, size_t *high, 
			 size_t *suffix_len, genome_t *genome) {
  
  int display = 1;
  char *query;
  size_t num_suffixes = 0;
  uint matched, max_matched = 0;
  size_t num_prefixes = search_prefix(seq, low, high, sa_index, display);
  
  *suffix_len = 0;
  num_suffixes = num_prefixes;
  
  if (num_prefixes && num_prefixes < max_num_suffixes) {   
    
    size_t first = *low, last = *low;
    
    if (num_prefixes == 1) {
      query = seq + sa_index->k_value;
      
      //char ref[strlen(seq)];      
      size_t start_g = sa_index->SA[*low] + sa_index->k_value;
      size_t end_g = start_g + strlen(seq);
      //genome_read_sequence_sa(ref, &start_g, &end_g, genome);

      size_t group_start = start_g/4;
      size_t group_end   = end_g/4;

      size_t nucleotide_start = start_g%4;
      size_t nucleotide_end   = end_g%4;

      unsigned int seq_pos = 0;
      char *str_tmp;

      matched = 0;
      int mismatch = 0;
      str_tmp = genome->code_table[genome->X[group_start]];
      for(unsigned int i = nucleotide_start; i < 4; i++){
	if (str_tmp[i] == query[matched]) { matched++; }
	else { mismatch = 1; break; }
      }
      group_start++;

      if (!mismatch) {	
	while (group_start < group_end) {
	  str_tmp = genome->code_table[genome->X[group_start]];
	  if (str_tmp[0] == query[matched]) { matched++; }
	  else { mismatch = 1; break; }
	  if (str_tmp[1] == query[matched]) { matched++; }
	  else { mismatch = 1; break; }
	  if (str_tmp[2] == query[matched]) { matched++; }
	  else { mismatch = 1; break; }
	  if (str_tmp[3] == query[matched]) { matched++; }
	  else { mismatch = 1; break; }
	  
	  group_start++;
	}
      }

      if (!mismatch) {
	if (group_start <= group_end) {  
	  str_tmp = genome->code_table[genome->X[group_start]];
	  for(unsigned int i = 0; i <= nucleotide_end; i++){
	    if (str_tmp[i] == query[matched]) { matched++; }
	    else { mismatch = 1; break; }
	  }
	}
      }

      //ref = &sa_index->genome->S[sa_index->SA[*low]] + sa_index->k_value;      
      //matched = 0;
      //while (query[matched] == ref[matched]) {
      //matched++;
      //}

      *high = *low;
      *suffix_len = matched + sa_index->k_value;
      num_suffixes = num_prefixes;

    } else {
      for (size_t i = *low; i < *high; i++) {
	query = seq + sa_index->k_value;

	//char ref[strlen(seq)];
	size_t start_g = sa_index->SA[i] + sa_index->k_value;
	size_t end_g = start_g + strlen(seq);
	
	size_t group_start = start_g/4;
	size_t group_end   = end_g/4;
	
	size_t nucleotide_start = start_g%4;
	size_t nucleotide_end   = end_g%4;

	unsigned int seq_pos = 0;
	//genome_read_sequence_sa(ref, &start_g, &end_g, genome);

	matched = 0;
	int mismatch = 0;
	char *str_tmp = genome->code_table[genome->X[group_start]];
	for(unsigned int i = nucleotide_start; i < 4; i++){
	  if (str_tmp[i] == query[matched]) { matched++; }
	  else { mismatch = 1; break; }
	}
	group_start++;
	
	if (!mismatch) {	
	  while (group_start < group_end) {
	    str_tmp = genome->code_table[genome->X[group_start]];
	    if (str_tmp[0] == query[matched]) { matched++; }
	    else { mismatch = 1; break; }
	    if (str_tmp[1] == query[matched]) { matched++; }
	    else { mismatch = 1; break; }
	    if (str_tmp[2] == query[matched]) { matched++; }
	    else { mismatch = 1; break; }
	    if (str_tmp[3] == query[matched]) { matched++; }
	    else { mismatch = 1; break; }
	    
	    group_start++;
	  }
	}
	
	if (!mismatch) {
	  if (group_start <= group_end) {  
	    str_tmp = genome->code_table[genome->X[group_start]];
	    for(unsigned int i = 0; i <= nucleotide_end; i++){
	      if (str_tmp[i] == query[matched]) { matched++; }
	      else { mismatch = 1; break; }
	    }
	  }
	}

	//matched = 0;
	//while (query[matched] == ref[matched]) {
	//matched++;
	//}

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
    }
  }
  

  
  return num_suffixes;
  
}


//--------------------------------------------------------------------
//--------------------------------------------------------------------

