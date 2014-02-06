#include "doscadfun.h"

//--------------------------------------------------------------------
//--------------------------------------------------------------------

inline void set_map_trace(int *map_counter, char *st1_map, char st1_value, 
			  char *nt_map, char nt_value, char *st2_map, char st2_value) {
  int count = *map_counter;
  st1_map[count] = st1_value;
  nt_map[count] = nt_value;
  st2_map[count] = st2_value;
  *map_counter = (++count);
}

//--------------------------------------------------------------------

float doscadfun(char* st1, int ll1, char* st2, int ll2, 
		float mismatch_perc, alig_out_t *out) {
  int  i, j;
  char ch;
  float  faux;
  int aborto=0;
  int num_errors=0;
  float score;

  int map_len1, map_len2, match, mism, gap1, gapmas, st_m_len;

  // variables to recover mapped data when errors > MAX_NUM_ERRORS
  int first_i, first_j, first_match, first_mism, first_gap1, first_gapmas, first_st_map_len;
  cigar_t first_cigar;

  // init output
  alig_out_init(out);

#ifdef _VERBOSE
  int max_len = (ll1 > ll2 ? ll1 : ll2), map_counter = 0;
  char st1_map[max_len * 2], st2_map[max_len * 2], nt_map[max_len * 2];
#endif
  
  // init recovery variables
  first_i=0; first_j=0;
  first_match=0; first_mism=0;
  first_gap1=0; first_gapmas=0;
  first_st_map_len=0;
  cigar_init(&first_cigar);
  
  // init variables
  match=0; mism=0; gap1=0; gapmas=0; st_m_len=0;
  	
  if (mismatch_perc > 0.0f) {
    int max_num_mismatches = round(ll1 * mismatch_perc);
    if (ll1 < 10) max_num_mismatches *= 2;
    for (i = 0; i < ll1; i++) {
      if (st1[i] == st2[i]) {
	//(*match)++;
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[i]);
#endif
	match++;
      } else {
	//(*mism)++;
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[i]);
#endif
	mism++;
	if (mism > max_num_mismatches) break;
      }
    }
    
    if (mism <= max_num_mismatches) {
      //*map_len1 = ll1;
      //*map_len2 = ll1;
      map_len1 = ll1;
      map_len2 = ll1;
      st_m_len = ll1;
      //
      //definicion de la funcion::::::::::::::::::::::
      //out_alig_set(int map_l1, int map_l2, int mm, int mism,
      //                     int gap1, int gapmas, st_m_len, out_alig_t *out) 		  
      alig_out_set(map_len1, map_len2, match, mism, gap1, gapmas, st_m_len, out);
      cigar_append_op(map_len1, 'M', &out->cigar);
      
      //score = (float) (((*match) * 5.0f)-((*mism) * 4.0f));
      score = (float) ((match * 5.0f)-(mism * 4.0f));
      
#ifdef _VERBOSE
      st1_map[map_counter] = 0;
      nt_map[map_counter] = 0;
      st2_map[map_counter] = 0;
      
      printf("\n\t%s\n\t%s\n\t%s\n", st1_map, nt_map, st2_map);
      
      //  printf("\t\t\tscore = %0.2f (mapped %i of %i): matches = %d, mismatches = %d, open-gaps = %d, extended-gaps = %d\n", 
      //	 score, *map_len1, ll1, *match, *mism, *gap1, *gapmas);
      printf("\tscore = %0.2f (mapped %i of %i): matches = %d, mismatches = %d, open-gaps = %d, extended-gaps = %d, cigar = %s\n", 
	     score, map_len1, ll1, match, mism, gap1, gapmas, cigar_to_string(&out->cigar));	 
#endif
      
      
      
      return score;
    }
    //---> initialize again.
    //*match=0; *mism=0; score = 0.0f;
    match=0; mism=0; score = 0.0f;
    
#ifdef _VERBOSE
    map_counter=0;//===>> fundamental para resetear los strings
#endif
  }
  /*  
      #ifdef _VERBOSE
      printf("\t\t\tS1=%s : %4d\n", st1, ll1);
      char *ref = get_subsequence(st2, 0, ll2 - 1);
      printf("\t\t\tS2=%s : %4d\n", ref, ll2);
      free(ref);
      #endif
  */
  num_errors=0;
  for (i=0, j=0; (i<(ll1-2)) && (j<ll2); i++, j++){
    if (st1[i]==st2[j]){
#ifdef _VERBOSE
      set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[j]);
#endif
      //(*match)++;
      match++;
      cigar_append_op(1, 'M', &out->cigar);
      num_errors = 0; //--> reset a la cuenta de num_errors!!!!
    } 
    else {//---> que cerrare antes del for
      if (num_errors == MAX_NUM_ERRORS){
		
#ifdef _VERBOSE
	printf("\t\t\t---- abort at (%i, %i) first error at (%i, %i, %s) -- %d consecutive errors --\n", 
	       i, j, first_i, first_j, cigar_to_string(&first_cigar), num_errors);
#endif
	alig_out_set(first_i, first_j, first_match, first_mism, first_gap1, first_gapmas, first_st_map_len, out);
	cigar_copy(&out->cigar, &first_cigar);
	   
	return (-1.0f);
      } // to avoid concatenating more than MAX_NUM_ERRORS consecutive. (27-XI-2013)
      else {
	num_errors++;
	if (num_errors == 1) {
	  first_i = i;
	  first_j = j;
	  
	  first_match  = match;
	  first_mism   = mism;
	  first_gap1   = gap1;
	  first_gapmas = gapmas;

	  //------> NEW: add to 11-XII-2013
#ifdef _VERBOSE
	  first_st_map_len= map_counter;
#endif
	  cigar_copy(&first_cigar, &out->cigar);
	}
      }
      //----------------------->> sigo: hay pocos num_errors consecutivos.  
      if ((st1[i+1]==st2[j+1]) && (st1[i+2]==st2[j+2])) { //--> mismatch de 1nt 
#ifdef _VERBOSE  
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]); 
	set_map_trace(&map_counter, st1_map, st1[i+1], nt_map, '|', st2_map, st2[j+1]);  
	set_map_trace(&map_counter, st1_map, st1[i+2], nt_map, '|', st2_map, st2[j+2]);
#endif
	//(*mism)++; i++; j++;
	//(*match)++;    
	//(*match)++;    i++; j++;
	mism++; match+=2;
	cigar_append_op(3, 'M', &out->cigar);	
	i+=2; j+=2;
      }
      //--> busco "x | x | | "
      else if((st1[i]!=st2[j]) && (st1[i+1]==st2[j+1]) && (st1[i+2]!=st2[j+2]) && 
	      (st1[i+3]==st2[j+3]) && (st1[i+4]==st2[j+4])) {
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i+1], nt_map, '|', st2_map, st2[j+1]);
	set_map_trace(&map_counter, st1_map, st1[i+2], nt_map, 'x', st2_map, st2[j+2]);
	set_map_trace(&map_counter, st1_map, st1[i+3], nt_map, '|', st2_map, st2[j+3]);
	set_map_trace(&map_counter, st1_map, st1[i+4], nt_map, '|', st2_map, st2[j+4]);
#endif
	//(*mism)++; 
	//(*match)++; i++; j++;    
	//(*mism)++;  i++; j++;
	//(*match)++; i++; j++;
	//(*match)++; i++; j++;
	mism+=2; match+=3; 
	cigar_append_op(5, 'M', &out->cigar);	
	i+=4; j+=4;
      }
      //--> busco gaps de 1nt
      else if((st1[i]==st2[j+1]) && (st1[i+1]==st2[j+2])) {
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[j+1]);
	set_map_trace(&map_counter, st1_map, st1[i+1], nt_map, '|', st2_map, st2[j+2]);
#endif
	//(*match)++; 
	//(*match)++; i++; j++;
	//(*gap1)++; j++;
	match+=2; gap1++; 
	cigar_append_op(1, 'D', &out->cigar);	
	cigar_append_op(2, 'M', &out->cigar);	
	i++; j+=2;
      }
      else if((st1[i+1]==st2[j]) && (st1[i+2]==st2[j+1])) { 
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, ' ', st2_map, '-');
	set_map_trace(&map_counter, st1_map, st1[i+1], nt_map, '|', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i+2], nt_map, '|', st2_map, st2[j+1]);
#endif
	//(*match)++; 
	//(*match)++; i++; j++;
	//(*gap1)++;   i++;
	match+=2; gap1++; 
	cigar_append_op(1, 'I', &out->cigar);	
	cigar_append_op(2, 'M', &out->cigar);	
	i+=2; j++;
      }
      else if ((st1[i+2]==st2[j+2])  && (st1[i+3]==st2[j+3])) { //--> mismatch de 2nt 
	
	if (i>ll1-4){//--> dos mismatch NO caben
#ifdef _VERBOSE
	  printf("\t\t\t---- break at (%i, %i) --------\n", i, j);
#endif
	  break;
	}//---> para no salirme de la read en las comparaciones.
	
#ifdef _VERBOSE   
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i+1], nt_map, 'x', st2_map, st2[j+1]);
	set_map_trace(&map_counter, st1_map, st1[i+2], nt_map, '|', st2_map, st2[j+2]);
	set_map_trace(&map_counter, st1_map, st1[i+3], nt_map, '|', st2_map, st2[j+3]);
#endif
	//(*mism)++; i++; j++;
	//(*mism)++; i++; j++;
	//(*match)++;    
	//(*match)++;    i++; j++;
	mism+=2; match+=2; 
	cigar_append_op(4, 'M', &out->cigar);	
	i+=3; j+=3;  
      }
      //--> busco gaps de 2nt
      else if((st1[i]==st2[j+2]) && (st1[i+1]==st2[j+3])) {
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j+1]);
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[j+2]);
	set_map_trace(&map_counter, st1_map, st1[i+1], nt_map, '|', st2_map, st2[j+3]);
#endif
	//(*gap1)++; j++; 
	//(*gapmas)++; j++;
	//(*match)++;
	//(*match)++;    i++; j++; 
	match+=2; gap1++; gapmas++; 
	cigar_append_op(2, 'D', &out->cigar);	
	cigar_append_op(2, 'M', &out->cigar);	
	i++; j+=3;
      }
      else if((st1[i+2]==st2[j]) && (st1[i+3]==st2[j+1])) { 
	
	if (i>ll1-4){//--> dos gaps arriba NO caben
#ifdef _VERBOSE
	  printf("\t\t\t---- break at (%i, %i) --------\n", i, j);
#endif
          break;
	}//--->
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, ' ', st2_map, '-');
	set_map_trace(&map_counter, st1_map, st1[i+1], nt_map, ' ', st2_map, '-');
	set_map_trace(&map_counter, st1_map, st1[i+2], nt_map, '|', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i+3], nt_map, '|', st2_map, st2[j+1]);
#endif
	//(*gap1)++; i++;
	//(*gapmas)++; i++;
	//(*match)++;
	//(*match)++;    i++; j++;
	match+=2; gap1++; gapmas++; 
	cigar_append_op(2, 'I', &out->cigar);	
	cigar_append_op(2, 'M', &out->cigar);	
	i+=3; j++;
      }
      //--> NUEVO: busco mut+gaps de 1nt
      else if((st1[i+1]==st2[j+2]) && (st1[i+2]==st2[j+3])) {
#ifdef _VERBOSE
	//------>> prioridad "x-"
	//set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]);
	//set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j+1]);
	
	//------>> prioridad "-x"
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j+1]);
	
	set_map_trace(&map_counter, st1_map, st1[i+1], nt_map, '|', st2_map, st2[j+2]);
	set_map_trace(&map_counter, st1_map, st1[i+2], nt_map, '|', st2_map, st2[j+3]);
#endif
	//(*mism)++; i++; j++; 
	//(*gap1)++; j++;
	//(*match)++;
	//(*match)++;    i++; j++; 
	mism++; match+=2; gap1++; 
	cigar_append_op(1, 'D', &out->cigar);	
	cigar_append_op(3, 'M', &out->cigar);	
	i+=2; j+=3;
      }
      else if((st1[i+2]==st2[j+1]) && (st1[i+3]==st2[j+2])) { 
	
	if (i>ll1-4){//--> mismatch+gap arriba NO caben
#ifdef _VERBOSE
	  printf("\t\t\t---- break at (%i, %i) --------\n", i, j);
#endif
	  break;
	}//--->
	
#ifdef _VERBOSE
	//------>> prioridad "x-"
	//set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]);
	//set_map_trace(&map_counter, st1_map, st1[i+1], nt_map, ' ', st2_map, '-');
	//------>> prioridad "-x"
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, ' ', st2_map, '-');
	set_map_trace(&map_counter, st1_map, st1[i+1], nt_map, 'x', st2_map, st2[j]);
	  
	set_map_trace(&map_counter, st1_map, st1[i+2], nt_map, '|', st2_map, st2[j+1]);
	set_map_trace(&map_counter, st1_map, st1[i+3], nt_map, '|', st2_map, st2[j+2]);
#endif
	//(*mism)++; i++; j++;
	//(*gap1)++; i++;
	//(*match)++;
	//(*match)++;    i++; j++;
	mism++; match+=2; gap1++;
	cigar_append_op(1, 'I', &out->cigar);	
	cigar_append_op(3, 'M', &out->cigar);	
	i+=3; j+=2;	  
      }                                  
      //--> busco gaps de 3nt
      else if ((st1[i]==st2[j+3]) && (st1[i+1]==st2[j+4])) {
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j+1]);
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j+2]);
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[j+3]);
	set_map_trace(&map_counter, st1_map, st1[i+1], nt_map, '|', st2_map, st2[j+4]);
#endif
	//(*gap1)++; (*gapmas)++; (*gapmas)++; j++; j++; j++;
	//(*match)++;                              
	//(*match)++;    i++; j++; 
	match+=2; gap1++; gapmas+=2; 
	cigar_append_op(3, 'D', &out->cigar);	
	cigar_append_op(2, 'M', &out->cigar);	
	i++; j+=4;
      }
      else if((st1[i+3]==st2[j]) && (st1[i+4]==st2[j+1])) { 
	if (i>ll1-5){//--> tres GAPS arriba NO caben
#ifdef _VERBOSE
	  printf("\t\t\t---- 3*break at (%i, %i) --------\n", i, j);
#endif
          break;
	}
		  
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, ' ', st2_map, '-');
	set_map_trace(&map_counter, st1_map, st1[i+1], nt_map, ' ', st2_map, '-');
	set_map_trace(&map_counter, st1_map, st1[i+2], nt_map, ' ', st2_map, '-');
	set_map_trace(&map_counter, st1_map, st1[i+3], nt_map, '|', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i+4], nt_map, '|', st2_map, st2[j+1]);
#endif
	//(*gap1)++; (*gapmas)++; (*gapmas)++; i++; i++; i++;
	//(*match)++; 
	//(*match)++;    i++; j++;
	match+=2; gap1++; gapmas+=2; 
	cigar_append_op(3, 'I', &out->cigar);	
	cigar_append_op(2, 'M', &out->cigar);	
	i+=4; j++;	  
      }    
      else if ((st1[i+3]==st2[j+3]) && (st1[i+4]==st2[j+4])) { //--> mismatch de 3nt 
#ifdef _VERBOSE
	  
	if (i>ll1-5){//--> tres GAPS arriba NO caben
#ifdef _VERBOSE
	  printf("\t\t\t---- 3*break at (%i, %i) --------\n", i, j);
#endif
          break;
	}
	  
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]);  
	set_map_trace(&map_counter, st1_map, st1[i+1], nt_map, 'x', st2_map, st2[j+1]);
	set_map_trace(&map_counter, st1_map, st1[i+2], nt_map, 'x', st2_map, st2[j+2]);
	set_map_trace(&map_counter, st1_map, st1[i+3], nt_map, '|', st2_map, st2[j+3]);
	set_map_trace(&map_counter, st1_map, st1[i+4], nt_map, '|', st2_map, st2[j+4]);
#endif
	//(*mism)++; i++; j++;
	//(*mism)++; i++; j++;
	//(*mism)++; i++; j++;
	//(*match)++; 
	//(*match)++; i++; j++;
	mism+=3; match+=2; 
	cigar_append_op(5, 'M', &out->cigar);	
	i+=4; j+=4;
	  
      }
      //--> busco "x | x | x | | "
     
      else if ((st1[i]!=st2[j]) && (st1[i+1]==st2[j+1]) && (st1[i+2]!=st2[j+2]) && (st1[i+3]==st2[j+3]) 
	       && (st1[i+4]!=st2[j+4]) && (st1[i+5]==st2[j+5]) && (st1[i+6]==st2[j+6])) {
		
	if (i>ll1-7){//--> tres mismatch intercalados NO caben
#ifdef _VERBOSE
	  printf("\t\t\t---- break at (%i, %i) --------\n", i, j);
#endif
          break;
	}
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i+1], nt_map, '|', st2_map, st2[j+1]);
	set_map_trace(&map_counter, st1_map, st1[i+2], nt_map, 'x', st2_map, st2[j+2]);
	set_map_trace(&map_counter, st1_map, st1[i+3], nt_map, '|', st2_map, st2[j+3]);
	set_map_trace(&map_counter, st1_map, st1[i+4], nt_map, 'x', st2_map, st2[j+4]);
	set_map_trace(&map_counter, st1_map, st1[i+5], nt_map, '|', st2_map, st2[j+5]);
	set_map_trace(&map_counter, st1_map, st1[i+6], nt_map, '|', st2_map, st2[j+6]);
#endif
	//(*mism)++; 
	//(*match)++; i++; j++;    
	//(*mism)++;  i++; j++;
	//(*match)++; i++; j++;
	//(*mism)++;  i++; j++;
	//(*match)++; i++; j++;
	//(*match)++; i++; j++;
	mism+=3; match+=4; 
	cigar_append_op(7, 'M', &out->cigar);	
	i+=6; j+=6;
		
      }

      else {
	//--> nada de nada: corto la comparacion!!!!  
	if (i + DEL_FINAL >= ll1) {//--> estoy a DEL_FINAL nt del final ... NO ABORTO!!
	  break;
	}
#ifdef _VERBOSE
	printf("\t\t\t------abort at (%i, %i)-------\n", i, j);
	/*
	  st1_map[map_counter] = 0;
	  nt_map[map_counter] = 0;
	  st2_map[map_counter] = 0;
	  printf("\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n", st1_map, nt_map, st2_map);
	  score= (float) (((*match)*5.0)-((*mism)*4.0));
	  score= score - (*gap1)*10.0;
	  score= score - (*gapmas)*0.5;
	  printf("\t\t\tscore = %0.2f, matches = %d, mismatches = %d, open-gaps = %d, extended-gaps = %d\n", 
	  score, *match, *mism, *gap1, *gapmas);
	*/
#endif
	//      exit(-1);
	// aborto= 1;
	//*map_len1 = i;
	//*map_len2 = j;
	  
	map_len1 = i;
	map_len2 = j;
	  	  
	if (num_errors>0)
	  st_m_len = first_st_map_len;
#ifdef _VERBOSE
	else st_m_len= map_counter;		
#endif
      
	alig_out_set(map_len1, map_len2, match, mism, gap1, gapmas, st_m_len, out);
      
	return (-1.0f);
	// i=ll1; j=ll2; //--> salgo llevando los índices al final.
      }  
    }//--> del else del primer IF  
    //getchar();
  }//----------> fin del bucle_for
  
  //-------> Comprobamos los dos últimos nt.

  //==================== DOS ULTIMOS NUCLEOTIDOS  ========================
  if(i==(ll1-2)){ 
    if (st1[i]==st2[j]){ //--> coincide 1
#ifdef _VERBOSE
      set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[j]);
#endif
      //(*match)++;
      match++;
      cigar_append_op(1, 'M', &out->cigar);	
      i++; j++;
      if (st1[i]==st2[j]){ //--> coincide 2
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[j]);
#endif
	//(*match)++;
	match++;
	cigar_append_op(1, 'M', &out->cigar);	
	i++; j++;
      }
      else{ //--> NO coincide 2: meto Mismatch
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]);
#endif
	//(*mism)++;
	mism++;
	cigar_append_op(1, 'M', &out->cigar);	
	i++; j++;
      }
    } 
    else {
      if ((st1[i+1]==st2[j+1]))  { //--> mach en 2nt y mismacht en 1nt
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i+1], nt_map, '|', st2_map, st2[j+1]);
#endif
	//(*mism)++;
	//(*match)++;
	//i++; j++; i++; j++;
	mism++; match++;
	cigar_append_op(2, 'M', &out->cigar);	
	i+=2; j+=2;			   
			   
      }
      else if ((st1[i]==st2[j+1]) && (st1[i+1]==st2[j+2])) { //--> GAP en 1nt y 2*match
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[j+1]);
	set_map_trace(&map_counter, st1_map, st1[i+1], nt_map, '|', st2_map, st2[j+2]);
#endif
	//(*gap1)++; j++;
	//(*match)++; (*match)++;
	//i++; j++; i++; j++;
	gap1++; match+=2;
	cigar_append_op(1, 'D', &out->cigar);	
	cigar_append_op(2, 'M', &out->cigar);	
	i+=2; j+=3;	
      }
      else  if ((st1[i]==st2[j+2]) && (st1[i+1]==st2[j+3])) { //--> 2*GAP y 2*match
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j+1]);
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[j+2]);
	set_map_trace(&map_counter, st1_map, st1[i+1], nt_map, '|', st2_map, st2[j+3]);
#endif
	//(*gap1)++; (*gapmas)++; j++; j++;
	//(*match)++; (*match)++;
	//i++; j++; i++; j++;
	gap1++; gapmas++; match+=2;
	cigar_append_op(2, 'D', &out->cigar);	
	cigar_append_op(2, 'M', &out->cigar);	
	i+=2; j+=4;	
      }
      else if ((st1[i]==st2[j+3]) && (st1[i+1]==st2[j+4])) { //--> 3*GAP y 2*match
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j+1]);
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j+2]);
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[j+3]);
	set_map_trace(&map_counter, st1_map, st1[i+1], nt_map, '|', st2_map, st2[j+4]);
#endif
	//(*gap1)++; (*gapmas)++; (*gapmas)++; j++; j++; j++;
	//(*match)++; (*match)++;
	//i++; j++; i++; j++;
	gap1++; gapmas++; match+=2;
	cigar_append_op(3, 'D', &out->cigar);	
	cigar_append_op(2, 'M', &out->cigar);	
	i+=2; j+=5;	
      }
      else {
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i+1], nt_map, 'x', st2_map, st2[j+1]);
#endif
	//(*mism)++; (*mism)++; 
	//i++; j++; i++; j++;
	mism+=2;
	cigar_append_op(2, 'M', &out->cigar);	
	i+=2; j+=2;	
      }
    }      
  }//--> del if() de los dos ultimos

  if(i==(ll1-1)){ //--> esoy en el últimpo nt!!
    if (st1[i]==st2[j]){ //--> coincide 1
#ifdef _VERBOSE
      set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[j]);
#endif
      //(*match)++;
      match++;
      cigar_append_op(1, 'M', &out->cigar);	
      i++; j++;
    }
    else{ //--> NO coincide 2: meto Mismatch
#ifdef _VERBOSE
      set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]);
#endif
      //(*mism)++;
      mism++;
      cigar_append_op(1, 'M', &out->cigar);	
      i++; j++;
    }
  }   
			
  //*map_len1 = i;
  //*map_len2 = j;

  map_len1 = i;
  map_len2 = j;

  //score= (float) (((*match)*5.0)-((*mism)*4.0));
  //score= score - (*gap1)*10.0;
  //score= score - (*gapmas)*0.5;
  
  score= (float) ((match*5.0)-(mism*4.0));
  score= score - (gap1*10.0);
  score= score - (gapmas*0.5);
  
#ifdef _VERBOSE
  st1_map[map_counter] = 0;
  nt_map[map_counter] = 0;
  st2_map[map_counter] = 0;
  
  printf("\n\t%s\n\t%s\n\t%s\n", st1_map, nt_map, st2_map);
  printf("\t\t\tscore = %0.2f (query: mapped %i of %i, ref: mapped %i of %i): matches = %d, mismatches = %d, open-gaps = %d, extended-gaps = %d, cigar = %s\n", 
	 score, map_len1, ll1, map_len2, ll2, match, mism, gap1, gapmas, cigar_to_string(&out->cigar));
#endif

  map_len1 = i;
  map_len2 = j;
#ifdef _VERBOSE
  st_m_len = map_counter;	
#endif

  //getchar();

  alig_out_set(map_len1, map_len2, match, mism, gap1, gapmas, st_m_len, out);

  
  // getchar();
  return (score);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------

float doscadfun_inv(char* st1, int ll1, char* st2, int ll2, 
		    float mismatch_perc, alig_out_t *out) {
  int  i, j;
  char ch;
  float  faux;
  int aborto=0;
  int num_errors=0;
  float score;
  
  int map_len1, map_len2, match, mism, gap1, gapmas, st_m_len;
	
  int first_i, first_j, first_match, first_mism, first_gap1, first_gapmas, first_st_map_len;
  cigar_t first_cigar;

  // init output
  alig_out_init(out);
  
#ifdef _VERBOSE
  printf("\t\t\tdoscadfun_inv: (len st1, len st2) = (%i, %i)\n", ll1, ll2); 
  int max_len = (ll1 > ll2 ? ll1 : ll2), map_counter = 0;
  char st1_map[max_len * 2], st2_map[max_len * 2], nt_map[max_len * 2];
#endif
  
  // init recovery variables
  first_i=0; first_j=0;
  first_match=0; first_mism=0;
  first_gap1=0; first_gapmas=0;
  first_st_map_len=0;
  cigar_init(&first_cigar);
  
  // init variables
  match=0; mism=0; gap1=0; gapmas=0;;

  if (mismatch_perc > 0.0f) {
    int max_num_mismatches = round(ll1 * mismatch_perc);
    if (ll1 < 10) max_num_mismatches *= 2;
    char *st22 = &st2[ll2 - ll1];
    for (i = ll1 - 1; i >= 0; i--) {
      if (st1[i] == st22[i]) {
	//(*match)++;
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st22[i]);
#endif
	match++;
      } else {
	//(*mism)++;
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st22[i]);
#endif
	mism++;
	if (mism > max_num_mismatches) break;
      }
    }

    if (mism <= max_num_mismatches) {
      map_len1 = ll1;
      map_len2 = ll1;
      st_m_len = ll1;
      
      alig_out_set(map_len1, map_len2, match, mism, gap1, gapmas, st_m_len, out);
      cigar_append_op(map_len1, 'M', &out->cigar);
      
      //score = (float) (((*match) * 5.0f)-((*mism) * 4.0f));
      score = (float) ((match * 5.0f)-(mism * 4.0f));
#ifdef _VERBOSE
      st1_map[map_counter] = 0;
      nt_map[map_counter] = 0;
      st2_map[map_counter] = 0;
      printf("\t\t\tmap_len1 = %i\n", map_len1);
      printf("\t\t\tmap_len2 = %i\n", map_len2);
      printf("\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n", st1_map, nt_map, st2_map);
      printf("\t\t\tscore = %0.2f (mapped %i of %i): matches = %d, mismatches = %d, open-gaps = %d, extended-gaps = %d -> cigar = %s\n", 
	     score, map_len1, ll1, match, mism, gap1, gapmas, cigar_to_string(&out->cigar));
#endif
      
      return score;
    }
    //*match=0; *mism=0; score = 0.0f;
    match=0; mism=0; score = 0.0f;
#ifdef _VERBOSE
    map_counter=0; //===>> fundamental para resetear los strings
#endif
  }
  
  /*  
      #ifdef _VERBOSE
      printf("\t\t\tS1=%s : %4d\n", st1, ll1);
      char *ref = get_subsequence(st2, 0, ll2 - 1);
      printf("\t\t\tS2=%s : %4d\n", ref, ll2);
      free(ref);
      #endif
  */
  
  num_errors = 0;
  for (i=(ll1-1), j=(ll2-1); (i>=2) && (j>=0); i--, j--){
    if (st1[i]==st2[j]){
#ifdef _VERBOSE
      set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[j]);
#endif
      //(*match)++;
      match++;
      cigar_append_op(1, 'M', &out->cigar);
      num_errors = 0; //--> reset a la cuenta de num_errors!!!!
    } 
    else {//---> que cerrare antes del for
      if (num_errors==MAX_NUM_ERRORS){
	
#ifdef _VERBOSE
	printf("\t\t\t---- abort at (%i, %i) first error at (%i, %i) --> %d consecutive errors --\n", 
	       i, j, first_i, first_j, num_errors);
#endif
	
	alig_out_set((ll1 - 1) - first_i, (ll2 - 1) - first_j, first_match, first_mism, 
		     first_gap1, first_gapmas, first_st_map_len, out);
	cigar_revcopy(&out->cigar, &first_cigar);
	
	return (-1.0f);
      }//---> para no concatenar mas de MAX_NUM_ERRORS consecutivos. (27-XI-2013)
      else {
	num_errors++;
	if (num_errors == 1) {
	  first_i = i;
	  first_j = j;
	  
	  first_match = match;
	  first_mism  = mism;
	  first_gap1  = gap1;
	  first_gapmas= gapmas;
	  //------> NEW: add to 11-XII-2013
#ifdef _VERBOSE
	  first_st_map_len= map_counter;
#endif
	  cigar_copy(&first_cigar, &out->cigar);
	}
      }
      //----------------------->> sigo: hay pocos num_errors consecutivos.  
      if ((st1[i-1]==st2[j-1]) && (st1[i-2]==st2[j-2])) { //--> mismatch de 1nt 
#ifdef _VERBOSE   
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]); 
	set_map_trace(&map_counter, st1_map, st1[i-1], nt_map, '|', st2_map, st2[j-1]); 
	set_map_trace(&map_counter, st1_map, st1[i-2], nt_map, '|', st2_map, st2[j-2]);
#endif
	//(*mism)++; i--; j--;
	//(*match)++;    
	//(*match)++; i--; j--;
	mism++; match+=2;
	cigar_append_op(3, 'M', &out->cigar);
	i-=2; j-=2;
      }
      
      //--> busco "x | x | |"
      else if((st1[i]!=st2[j]) && (st1[i-1]==st2[j-1]) && (st1[i-2]!=st2[j-2]) && 
	      (st1[i-3]==st2[j-3]) && (st1[i-4]==st2[j-4])) {
	
	if (i<4){
#ifdef _VERBOSE
	  printf("\t\t\t---- break at (%i, %i) --------\n", (ll1-1)-i, (ll2-1)-j);
#endif
	  break;
	}//---> para no salirme de la read en las comparaciones.				
	
	
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i-1], nt_map, '|', st2_map, st2[j-1]);
	set_map_trace(&map_counter, st1_map, st1[i-2], nt_map, 'x', st2_map, st2[j-2]);
	set_map_trace(&map_counter, st1_map, st1[i-3], nt_map, '|', st2_map, st2[j-3]);
	set_map_trace(&map_counter, st1_map, st1[i-4], nt_map, '|', st2_map, st2[j-4]);
#endif
	//(*mism)++; 
	//(*match)++; i--; j--;    
	//(*mism)++;  i--; j--;
	//(*match)++; i--; j--;
	//(*match)++; i--; j--;
	mism+=2; match+=3; 
	cigar_append_op(5, 'M', &out->cigar);	
	i-=4; j-=4;
      }
      
      //--> busco gaps de 1nt
      else if((st1[i]==st2[j-1]) && (st1[i-1]==st2[j-2])) {
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[j-1]);
	set_map_trace(&map_counter, st1_map, st1[i-1], nt_map, '|', st2_map, st2[j-2]);
#endif
	//(*match)++; 
	//(*match)++; i--; j--;
	//(*gap1)++; j--;
	match+=2; gap1++; 
	cigar_append_op(1, 'D', &out->cigar);	
	cigar_append_op(2, 'M', &out->cigar);	
	i--; j-=2;
      }
      else if((st1[i-1]==st2[j]) && (st1[i-2]==st2[j-1])) { 
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, ' ', st2_map, '-');
	set_map_trace(&map_counter, st1_map, st1[i-1], nt_map, '|', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i-2], nt_map, '|', st2_map, st2[j-1]);
#endif
	//(*match)++; 
	//(*match)++; i--; j--;
	//(*gap1)++;   i--; 
	match+=2; gap1++; 
	cigar_append_op(1, 'I', &out->cigar);	
	cigar_append_op(2, 'M', &out->cigar);	
	i-=2; j--;
      }
      
      else if ((st1[i-2]==st2[j-2])  && (st1[i-3]==st2[j-3])) { //--> mismatch de 2nt 
	
	//if (i>ll1-4){//--> dos mismatch NO caben
	if (i<3){//--> dos mismatch NO caben
#ifdef _VERBOSE
	  printf("\t\t\t---- break at (%i, %i) --------\n", (ll1-1)-i, (ll2-1)-j);
#endif
          break;
	}//---> para no salirme de la read en las comparaciones.	
	
#ifdef _VERBOSE  
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i-1], nt_map, 'x', st2_map, st2[j-1]);
	set_map_trace(&map_counter, st1_map, st1[i-2], nt_map, '|', st2_map, st2[j-2]);
	set_map_trace(&map_counter, st1_map, st1[i-3], nt_map, '|', st2_map, st2[j-3]);
#endif
	//(*mism)++; i--; j--;
	//(*mism)++; i--; j--;
	//(*match)++;    
	//(*match)++;    i--; j--;
	mism+=2; match+=2; 
	cigar_append_op(4, 'M', &out->cigar);	
	i-=3; j-=3;
      }
      
      //--> busco gaps de 2nt
      else if((st1[i]==st2[j-2]) && (st1[i-1]==st2[j-3])) {
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j-1]);
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[j-2]);
	set_map_trace(&map_counter, st1_map, st1[i-1], nt_map, '|', st2_map, st2[j-3]);
#endif
	//(*gap1)++; j--; 
	//(*gapmas)++; j--;
	//(*match)++;
	//(*match)++;    i--; j--; 
	match+=2; gap1++; gapmas++; 
	cigar_append_op(2, 'D', &out->cigar);	
	cigar_append_op(2, 'M', &out->cigar);	
	i--; j-=3;
      }
      else if((st1[i-2]==st2[j]) && (st1[i-3]==st2[j-1])) { 
	
	// if (i>ll1-4){//--> dos gaps arriba NO caben
	if (i<3){//--> dos gaps arriba NO caben
#ifdef _VERBOSE
	  printf("\t\t\t---- break at (%i, %i) --------\n", (ll1-1)-i, (ll2-1)-j);
#endif
          break;
	}//--->
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, ' ', st2_map, '-');
	set_map_trace(&map_counter, st1_map, st1[i-1], nt_map, ' ', st2_map, '-');
	set_map_trace(&map_counter, st1_map, st1[i-2], nt_map, '|', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i-3], nt_map, '|', st2_map, st2[j-1]);
#endif
	//(*gap1)++; i--;
	//(*gapmas)++; i--;
	//(*match)++;
	//(*match)++; i--; j--;
	match+=2; gap1++; gapmas++; 
	cigar_append_op(2, 'I', &out->cigar);	
	cigar_append_op(2, 'M', &out->cigar);	
	i-=3; j--;
      }
      
      //--> NUEVO: busco mut+gaps de 1nt
      else if((st1[i-1]==st2[j-2]) && (st1[i-2]==st2[j-3])) {
#ifdef _VERBOSE
	//------>> prioridad "-x" en la directa. En la inversa es "x-"
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j-1]);
	set_map_trace(&map_counter, st1_map, st1[i-1], nt_map, '|', st2_map, st2[j-2]);
	set_map_trace(&map_counter, st1_map, st1[i-2], nt_map, '|', st2_map, st2[j-3]);
#endif
	//(*mism)++; i--; j--; 
	//(*gap1)++; j--;
	//(*match)++;
	//(*match)++;    i--; j--; 
	mism++; match+=2; gap1++; 
	cigar_append_op(1, 'M', &out->cigar);	
	cigar_append_op(1, 'D', &out->cigar);	
	cigar_append_op(2, 'M', &out->cigar);	
	i-=2; j-=3;
      }
      else if((st1[i-2]==st2[j-1]) && (st1[i-3]==st2[j-2])) { 
	
	//if (i>ll1-4){//--> mismatch+gap arriba NO caben
	if (i<3){//--> mismatch+gap arriba NO caben
#ifdef _VERBOSE
	  printf("\t\t\t---- break at (%i, %i) --------\n", (ll1-1)-i, (ll2-1)-j);
#endif
          break;
	}//--->
#ifdef _VERBOSE
	//------>> prioridad "-x" en la directa. En la inversa es "x-"
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i-1], nt_map, ' ', st2_map, '-');
	set_map_trace(&map_counter, st1_map, st1[i-2], nt_map, '|', st2_map, st2[j-1]);
	set_map_trace(&map_counter, st1_map, st1[i-3], nt_map, '|', st2_map, st2[j-2]);
#endif
	//(*mism)++; i--; j--;
	//(*gap1)++; i--;
	//(*match)++;
	//(*match)++;    i--; j--;
	mism++; match+=2; gap1++;
	cigar_append_op(1, 'M', &out->cigar);	
	cigar_append_op(1, 'I', &out->cigar);	
	cigar_append_op(2, 'M', &out->cigar);	
	i-=3; j-=2;
      }                                  
      //--> busco gaps de 3nt
      else if ((st1[i]==st2[j-3]) && (st1[i-1]==st2[j-4])) {
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j-1]);
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j-2]);
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[j-3]);
	set_map_trace(&map_counter, st1_map, st1[i-1], nt_map, '|', st2_map, st2[j-4]);
#endif
	//(*gap1)++; (*gapmas)++; (*gapmas)++; j--; j--; j--;
	//(*match)++;                              
	//(*match)++;    i--; j--; 
	match+=2; gap1++; gapmas+=2; 
	cigar_append_op(3, 'D', &out->cigar);	
	cigar_append_op(2, 'M', &out->cigar);	
	i--; j-=4;
      }
      else if((st1[i-3]==st2[j]) && (st1[i-4]==st2[j-1])) { 
	
	//if (i>ll1-5){//--> tres GAPS arriba NO caben
	if (i<4){//--> tres GAPS arriba NO caben
#ifdef _VERBOSE
	  printf("\t\t\t---- break at (%i, %i) --------\n", (ll1-1)-i, (ll2-1)-j);
#endif
          break;
	}
	
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, ' ', st2_map, '-');
	set_map_trace(&map_counter, st1_map, st1[i-1], nt_map, ' ', st2_map, '-');
	set_map_trace(&map_counter, st1_map, st1[i-2], nt_map, ' ', st2_map, '-');
	set_map_trace(&map_counter, st1_map, st1[i-3], nt_map, '|', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i-4], nt_map, '|', st2_map, st2[j-1]);
#endif
	//(*gap1)++; (*gapmas)++; (*gapmas)++; i--; i--; i--;
	//(*match)++; 
	//(*match)++;    i--; j--;
	match+=2; gap1++; gapmas+=2; 
	cigar_append_op(3, 'I', &out->cigar);	
	cigar_append_op(2, 'M', &out->cigar);	
	i-=4; j--;
      }    
      else if ((st1[i-3]==st2[j-3]) && (st1[i-4]==st2[j-4])) { //--> mismatch de 3nt 
	
	//if (i>ll1-5){//--> tres mismatch NO caben
	if (i<4){//--> tres GAPS arriba NO caben
#ifdef _VERBOSE
	  printf("\t\t\t---- break at (%i, %i) --------\n", (ll1-1)-i, (ll2-1)-j);
#endif
          break;
	}
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]);  
	set_map_trace(&map_counter, st1_map, st1[i-1], nt_map, 'x', st2_map, st2[j-1]);
	set_map_trace(&map_counter, st1_map, st1[i-2], nt_map, 'x', st2_map, st2[j-2]);
	set_map_trace(&map_counter, st1_map, st1[i-3], nt_map, '|', st2_map, st2[j-3]);
	set_map_trace(&map_counter, st1_map, st1[i-4], nt_map, '|', st2_map, st2[j-4]);
#endif
	//(*mism)++; i--; j--;
	//(*mism)++; i--; j--;
	//(*mism)++; i--; j--;
	//(*match)++; 
	//(*match)++; i--; j--;
	mism+=3; match+=2; 
	cigar_append_op(5, 'M', &out->cigar);	
	i-=4; j-=4;
      }
      //--> busco "x | x | x | |"
      
      else if ((st1[i]!=st2[j]) && (st1[i-1]==st2[j-1]) && (st1[i-2]!=st2[j-2]) && (st1[i-3]==st2[j-3]) 
	       && (st1[i-4]!=st2[j-4]) && (st1[i-5]==st2[j-5]) && (st1[i-6]==st2[j-6])) {
	//if (i>ll1-7){//--> tres mismatch intercalados NO caben
	if (i<6){//--> tres mismatch intercalados NO caben
#ifdef _VERBOSE
	  printf("\t\t\t---- break at (%i, %i) --------\n", (ll1-1)-i, (ll2-1)-j);
#endif
	  break;
	}
	
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i-1], nt_map, '|', st2_map, st2[j-1]);
	set_map_trace(&map_counter, st1_map, st1[i-2], nt_map, 'x', st2_map, st2[j-2]);
	set_map_trace(&map_counter, st1_map, st1[i-3], nt_map, '|', st2_map, st2[j-3]);
	set_map_trace(&map_counter, st1_map, st1[i-4], nt_map, 'x', st2_map, st2[j-4]);
	set_map_trace(&map_counter, st1_map, st1[i-5], nt_map, '|', st2_map, st2[j-5]);
	set_map_trace(&map_counter, st1_map, st1[i-6], nt_map, '|', st2_map, st2[j-6]);
#endif
	//(*mism)++; 
	//(*match)++; i--; j--;    
	//(*mism)++;  i--; j--;
	//(*match)++; i--; j--;
	//(*mism)++;  i--; j--;
	//(*match)++; i--; j--;
	//(*match)++; i--; j--;
	mism+=3; match+=4; 
	cigar_append_op(7, 'M', &out->cigar);	
	i-=6; j-=6;
      }
      
      else {
	//--> nada de nada: corto la comparacion!!!!  
	if (i - DEL_FINAL <= 0) {//--> estoy a DEL_FINAL nt del final ... NO ABORTO!!
	  //	  printf("\t\t\t------break (close to the end), at (%i, %i)-------\n", (ll1-1)-i, (ll2-1)-j);
	  break;
	}
#ifdef _VERBOSE
	printf("\t\t\t------abort at (%i, %i)-------\n", (ll1-1)-i, (ll2-1)-j);
	
	st1_map[map_counter] = 0;
	nt_map[map_counter] = 0;
	st2_map[map_counter] = 0;
	printf("\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n", st1_map, nt_map, st2_map);
	score= (float) ((match*5.0)-(mism*4.0));
	score= score - (gap1*10.0);
	score= score - (gapmas*0.5);
	printf("\t\t\tscore = %0.2f, matches = %d, mismatches = %d, open-gaps = %d, extended-gaps = %d\n", 
	       score, match, mism, gap1, gapmas);
	
	//	getchar();
	
#endif
	map_len1 = (ll1)-i;
	map_len2 = (ll2)-j;
	
	if (num_errors>0)
	  st_m_len= first_st_map_len;
#ifdef _VERBOSE
	else st_m_len= map_counter;
#endif
    
	alig_out_set(map_len1, map_len2, match, mism, gap1, gapmas, st_m_len, out);
	cigar_rev(&out->cigar);
	
	return (-1.0f);
	// i=ll1; j=ll2; //--> salgo llevando los índices al final.
      }  
    }//--> del else del primer IF  
    //getchar();
  }//----------> fin del bucle_for
  
  //-------> Comprobamos los dos últimos nt.

  //==========  DOS ULTIMOS NUCLEOTIDOS ========================

  if(i==1){ 
    if (st1[i]==st2[j]){ //--> coincide 1
#ifdef _VERBOSE
      set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[j]);
#endif
      //(*match)++;
      match++;
      cigar_append_op(1, 'M', &out->cigar);	
      i--; j--;
      if (st1[i]==st2[j]){ //--> coincide 2
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[j]);
#endif
	//(*match)++;
	match++;
	cigar_append_op(1, 'M', &out->cigar);	
	i--; j--;
      }
      else{ //--> NO coincide 2: meto Mismatch
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]);
#endif
	//(*mism)++;
	mism++;
	cigar_append_op(1, 'M', &out->cigar);	
	i--; j--;
      }
    } 
    else {
      if ((st1[i-1]==st2[j-1]))  { //--> mach en 2nt y mismacht en 1nt
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i-1], nt_map, '|', st2_map, st2[j-1]);
#endif
	//(*mism)++;
	//(*match)++;
	mism++; match++;
	cigar_append_op(2, 'M', &out->cigar);	
	i-=2; j-=2;	
      }
      else if ((st1[i]==st2[j-1]) && (st1[i-1]==st2[j-2])) { //--> GAP en 1nt y 2*match
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[j-1]);
	set_map_trace(&map_counter, st1_map, st1[i-1], nt_map, '|', st2_map, st2[j-2]);
#endif
	//(*gap1)++; j--;
	//(*match)++; (*match)++;
	//i--; j--; i--; j--;
	gap1++; match+=2;
	cigar_append_op(1, 'D', &out->cigar);	
	cigar_append_op(2, 'M', &out->cigar);	
	i-=2; j-=2;
      }
      else  if ((st1[i]==st2[j-2]) && (st1[i-1]==st2[j-3])) { //--> 2*GAP y 2*match
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j-1]);
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[j-2]);
	set_map_trace(&map_counter, st1_map, st1[i-1], nt_map, '|', st2_map, st2[j-3]);
#endif
	//(*gap1)++; (*gapmas)++; j--; j--;
	//(*match)++; (*match)++;
	//i--; j--; i--; j--;
	gap1++; gapmas++; match+=2;
	cigar_append_op(2, 'D', &out->cigar);	
	cigar_append_op(2, 'M', &out->cigar);	
	i-=2; j-=4;
      }
      else if ((st1[i]==st2[j-3]) && (st1[i-1]==st2[j-4])) { //--> 3*GAP y 2*match
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j-1]);
	set_map_trace(&map_counter, st1_map, '-', nt_map, ' ', st2_map, st2[j-2]);
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[j-3]);
	set_map_trace(&map_counter, st1_map, st1[i-1], nt_map, '|', st2_map, st2[j-4]);
#endif
	//(*gap1)++; (*gapmas)++; (*gapmas)++; j--; j--; j--;
	//(*match)++; (*match)++;
	//i--; j--; i--; j--;
	gap1++; gapmas+=2; match+=2;
	cigar_append_op(3, 'D', &out->cigar);	
	cigar_append_op(2, 'M', &out->cigar);	
	i-=2; j-=5;
      }
      else {
#ifdef _VERBOSE
	set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]);
	set_map_trace(&map_counter, st1_map, st1[i-1], nt_map, 'x', st2_map, st2[j-1]);
#endif
	//(*mism)++; (*mism)++; 
	//i--; j--; i--; j--;
	mism+=2;
	cigar_append_op(2, 'M', &out->cigar);	
	i-=2; j-=2;
      }
    }      
  }//--> del if() de los dos ultimos
  
  if(i==0){ //--> esoy en el últimpo nt!!
    if (st1[i]==st2[j]){ //--> coincide 1
#ifdef _VERBOSE
      set_map_trace(&map_counter, st1_map, st1[i], nt_map, '|', st2_map, st2[j]);
#endif
      //(*match)++;
      match++;
      cigar_append_op(1, 'M', &out->cigar);	
      i--; j--;
    }
    else{ //--> NO coincide 2: meto Mismatch
#ifdef _VERBOSE
      set_map_trace(&map_counter, st1_map, st1[i], nt_map, 'x', st2_map, st2[j]);
#endif
      //(*mism)++;
      mism++;
      cigar_append_op(1, 'M', &out->cigar);	
      i--; j--;
    }
  }   
  
  //*map_len1 = ll1 - i - 1;
  //*map_len2 = ll2 - j - 1;
  
  map_len1 = ll1 - i - 1;
  map_len2 = ll2 - j - 1;
  
  //score= (float) (((*match)*5.0)-((*mism)*4.0));
  //score= score - (*gap1)*10.0;
  //score= score - (*gapmas)*0.5;
  
  score= (float) ((match*5.0)-(mism*4.0));
  score= score - (gap1*10.0);
  score= score - (gapmas*0.5);
  
  
#ifdef _VERBOSE
  st1_map[map_counter] = 0;
  nt_map[map_counter] = 0;
  st2_map[map_counter] = 0;
  printf("\t\t\tmap_len1 = %i (i = %i, ll1 = %i)\n", map_len1, i, ll1);
  printf("\t\t\tmap_len2 = %i (j = %i, ll2 = %i)\n", map_len2, j, ll2);
  printf("\t\t\t%s\n\t\t\t%s\n\t\t\t%s\n", st1_map, nt_map, st2_map);
  printf("\t\t\tscore = %0.2f (query: mapped %i of %i, ref: mapped %i of %i): matches = %d, mismatches = %d, open-gaps = %d, extended-gaps = %d, cigar = %s\n", 
	 score, map_len1, ll1, map_len2, ll2, match, mism, gap1, gapmas, cigar_to_string(&out->cigar));
  
  //  if (score == 314 && strncmp(st2_map, "AGAGC", 5) == 0) {
  //    exit(-1);
  //  }
#endif
  
  //*****************************************
  /* codigo alternativo para ver el mapeo al reves de lo que se hace
     #ifdef _VERBOSE
     for (i=(map_counter-1); i>=0; i--)//--> bucle inverso
     printf("%c", st1_map[i]);
     printf("\n");	  
     for (i=(map_counter-1); i>=0; i--)//--> bucle inverso
     printf("%c", nt_map[i]);
     printf("\n");		  
     for (i=(map_counter-1); i>=0; i--)//--> bucle inverso
     printf("%c", st2_map[i]);	  
     printf("\n");
     #endif
     // */
  
  alig_out_set(map_len1, map_len2, match, mism, gap1, gapmas, st_m_len, out);
  cigar_rev(&out->cigar);
  
  // getchar();
  return (score);
}

//--------------------------------------------------------------------
//--------------------------------------------------------------------
