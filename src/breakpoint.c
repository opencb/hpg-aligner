#include "breakpoint.h"
#include "rna/rna_splice.h"

//#include "extrae_user_events.h" 

#define CHUNK_SHIFT   10

//5 //10
//--------------------------------------------------------------------------------------

array_list_t *breakpoint_list = NULL;

//--------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------

cigar_op_t *cigar_op_new(int number, char name) {
  cigar_op_t *p = (cigar_op_t *) malloc(sizeof(cigar_op_t));
  p->number = number;
  p->name = name;
  return p;
}

void cigar_op_free(cigar_op_t *p) {
  if (p) free(p);
}

//--------------------------------------------------------------------------------------
//
//--------------------------------------------------------------------------------------

cigar_code_t *cigar_code_new() {
  cigar_code_t *p = (cigar_code_t *) calloc(1, sizeof(cigar_code_t));

  p->distance = 0;
  p->cigar_str = NULL;
  p->ops = array_list_new(20, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  return p;
}

cigar_code_t *cigar_code_dup(cigar_code_t *c) {
  cigar_code_t *p = cigar_code_new();
  size_t num_ops = array_list_size(c->ops);
  
  for (int i = 0; i < num_ops; i++) {
    cigar_op_t *op = array_list_get(i, c->ops);
    array_list_insert(cigar_op_new(op->number, op->name), p->ops);
  }

  return p;

}

cigar_code_t *cigar_code_new_by_string(char *cigar_str) {
  cigar_code_t *p = cigar_code_new();
  
  int cigar_len = strlen(cigar_str);

  //p->cigar_str = strdup(cigar_str);
  //p->num_allocated_ops = (cigar_len > 10 ? cigar_len / 2 : cigar_len);
  //p->ops = (cigar_op_t*) calloc(p->num_allocated_ops, sizeof(cigar_op_t));


  int c = 0;
  char op;
  char op_value[1024];

  for (int j = 0; j < cigar_len; j++) {
    op = cigar_str[j];
    if (op < 58) {
      op_value[c++] = op;
    } else {
      //op_value[c++] = op;
      op_value[c] = '\0';
      cigar_code_append_new_op(atoi(op_value), op, p);
      c = 0;     
    }
  }
  
  return p;
}

//--------------------------------------------------------------------------------------

void cigar_code_free(cigar_code_t* p) {
  if (p) {
    //    if (p->ops) array_list_free(p->ops, (void *) cigar_op_free);
    if (p->ops) array_list_free(p->ops, (void *) NULL);
    if (p->cigar_str) free(p->cigar_str);
    free(p);
  }
}

//--------------------------------------------------------------------------------------

void cigar_code_merge(cigar_code_t *p, cigar_code_t *merge_p) {
  cigar_op_t *op;
  for (int i = 0; i < array_list_size(merge_p->ops); i++) {
    op = array_list_get(i, merge_p->ops);
    array_list_insert(op, p->ops);
  }
  
}

//--------------------------------------------------------------------------------------

void cigar_code_update(cigar_code_t *p) {
  cigar_op_t *prev_op, *curr_op;
  size_t j, num_ops = array_list_size(p->ops);
  prev_op = array_list_get(0, p->ops);
  j = 0;
  for (size_t i = 1; i < num_ops; i++) {
    curr_op = array_list_get(i, p->ops);
    if (prev_op->name == curr_op->name) {
      prev_op->number += curr_op->number;
    } else {
      j++;
      prev_op = array_list_get(j, p->ops);
      prev_op->name = curr_op->name;
      prev_op->number = curr_op->number;
    }
  }
  for (size_t i = num_ops - 1; i > j; i--) {
    array_list_remove_at(i, p->ops);
  }
}

//--------------------------------------------------------------------------------------
int cigar_code_get_num_ops(cigar_code_t *p) {
  int num = 0;
  if (p && p->ops) { 
    return array_list_size(p->ops);
  }
  return num;
}

//--------------------------------------------------------------------------------------

cigar_op_t *cigar_code_get_first_op(cigar_code_t *p) {
  if (!p) { return NULL; }  
  return array_list_get(0, p->ops);
}

cigar_op_t *cigar_code_get_op(int index, cigar_code_t *p) {
  int num_ops = cigar_code_get_num_ops(p);
  if (num_ops > 0 && index < num_ops) {
    return array_list_get(index, p->ops);
  }
  return NULL;
}

cigar_op_t *cigar_code_get_last_op(cigar_code_t *p) {
  int num_ops = cigar_code_get_num_ops(p);

  return array_list_get(num_ops - 1, p->ops);
    //}
  //return NULL;
}

//--------------------------------------------------------------------------------------

void cigar_code_append_op(cigar_op_t *op, cigar_code_t *p) {
  if (p && p->ops && op) {
    cigar_op_t *last = cigar_code_get_last_op(p);
    if (last && last->name == op->name) {
      last->number += op->number;
      cigar_op_free(op);
    } else {
      array_list_insert(op, p->ops);
    }
  }
}

//--------------------------------------------------------------------------------------

void cigar_code_insert_first_op(cigar_op_t *op, cigar_code_t *p) {
  if (p && p->ops && op) {
    cigar_op_t *first = cigar_code_get_first_op(p);

    if (first != NULL && first->name == op->name) {
      first->number += op->number;
    } else {
      if (array_list_size(p->ops) == 0) {
	array_list_insert(op, p->ops);
      } else {
	array_list_insert_at(0, op, p->ops);
      }
    }
  }

}

//--------------------------------------------------------------------------------------

void cigar_code_append_new_op(int value, char name, cigar_code_t *p) {
  cigar_code_append_op(cigar_op_new(value, name), p);
}

//--------------------------------------------------------------------------------------

void cigar_code_inc_distance(int distance, cigar_code_t *p) {
  if (p && p->ops) {
    p->distance += distance;
  }
}

//--------------------------------------------------------------------------------------

char *new_cigar_code_string(cigar_code_t *p) {
  
  if (!p) { return NULL; }

  if (p->cigar_str) {
    free(p->cigar_str);
  }

  int num_ops = cigar_code_get_num_ops(p);
  //char *str = (char *)calloc(num_ops*10, sizeof(char));

  char str[2048] = "\0";

  if (num_ops == 0) { return 0; }
  
  cigar_op_t *op;
  
  for (int i = 0; i < num_ops; i++) {
    op = array_list_get(i, p->ops);
    sprintf(str, "%s%i%c", str, op->number, op->name);
  }

  //p->cigar_str = strdup(str);

  //printf("CIGAR: %s\n", p->cigar_str);

  return strdup(str);

}

//--------------------------------------------------------------------------------------

int cigar_match_coverage(cigar_code_t *p) {
  int coverage = 0;

  if (p) {    
    size_t num_ops = array_list_size(p->ops);
    for (size_t i = 0; i < num_ops; i++) {
      cigar_op_t *cigar_op = array_list_get(i, p->ops);
      if (cigar_op->name == 'M' || cigar_op->name == '=') {
	coverage += cigar_op->number;
      }
    }
  }

  return coverage;
}


int cigar_read_coverage(cigar_code_t *p) {
  int coverage = 0;

  if (p) {    
    size_t num_ops = array_list_size(p->ops);
    for (size_t i = 0; i < num_ops; i++) {
      cigar_op_t *cigar_op = array_list_get(i, p->ops);
      if (cigar_op->name == 'M' || cigar_op->name == 'I') {
	coverage += cigar_op->number;
      }
    }
  }

  return coverage;
}

int cigar_genome_coverage(cigar_code_t *p) {
  int coverage = 0;

  if (p) {    
    size_t num_ops = array_list_size(p->ops);
    for (size_t i = 0; i < num_ops; i++) {
      cigar_op_t *cigar_op = array_list_get(i, p->ops);
      if (cigar_op->name == 'M' || cigar_op->name == 'D') {
	coverage += cigar_op->number;
      }
    }
  }

  return coverage;
}

int cigar_code_nt_length(cigar_code_t *p) {
  if (!p) {
    return 0;
  }

  int len = 0;
  int num_ops = array_list_size(p->ops);

  cigar_op_t *op;
  for (int i = 0; i < num_ops; i++) {
    op = array_list_get(i, p->ops);
    if (op->name == 'M' || op->name == 'I' || op->name == '=') {
      len += op->number;
    }
  }

  return len;
}

//--------------------------------------------------------------------------------------

int cigar_code_validate(int read_length, cigar_code_t *p) {
  if (!p || !p->ops) { return 0; }

  int cigar_len = 0;
  for (int i = 0; i < array_list_size(p->ops); i++) {
    cigar_op_t *op = array_list_get(i, p->ops);
    if (op->number <= 0) { /*printf("ERROR CIGAR %s\n", new_cigar_code_string(p));*/ return 0; }
    if (op->name == 'M' || op->name == 'I') {
      cigar_len += op->number;
    }
  }
  //printf("cigar len = %i\n", cigar_len);
  return read_length == cigar_len;

}

int cigar_code_validate_(fastq_read_t *fq_read, cigar_code_t *p) {
  if (!p || !p->ops) { return 0; }
  int read_length = fq_read->length;

  int cigar_len = 0;
  for (int i = 0; i < array_list_size(p->ops); i++) {
    cigar_op_t *op = array_list_get(i, p->ops);
    if (op->number <= 0) { 
      /*printf("ERROR CIGAR %s: %s\n", fq_read->id, new_cigar_code_string(p));*/ 
      return 0;
    }
    if (op->name == 'M' || 
	op->name == 'I' || 
	op->name == 'S' ||
	op->name == 'H') {
      cigar_len += op->number;
    }
  }
  //printf("cigar len = %i\n", cigar_len);
  return read_length == cigar_len;

}

//--------------------------------------------------------------------------------------

void cigar_code_print(cigar_code_t *cigar_code) {
  for (int i = 0; i < cigar_code_get_num_ops(cigar_code); i++) {
    cigar_op_t * op = array_list_get(i, cigar_code->ops);
    printf("%i%c", op->number, op->name);
  }
  printf("\n");
}

//--------------------------------------------------------------------------------------

int cigar_code_score(cigar_code_t *cigar_code, int read_length) {
  int num_M = 0, num_D = 0, num_I = 0, t_D = 0, t_I = 0;

  for (int c = 0; c < array_list_size(cigar_code->ops); c++) {
    cigar_op_t *op = array_list_get(c, cigar_code->ops);
    if (op->name == 'M') {
      num_M += op->number;
    } else if (op->name == 'D') {
      num_D += op->number;
      t_D++;
    } else if (op->name == 'I') {
      num_I += op->number;
      t_I++;
    }
  }      
  	
  int tot_DI = num_D + num_I;
  int mismatches = cigar_code->distance - (tot_DI);

  if (mismatches < 0) {
    mismatches = 0;
  }

  int matches = num_M - mismatches;
  int score_m = (matches * 1) - ((mismatches + num_D + num_I)* 0.4);//mismatches + num_D + num_I == cigar_code->distance
  int read_span = (int)score_m * 100 / read_length;

  //printf("score = %i, distance = %i\n", read_span, cigar_code->distance);

  return read_span;
}

//--------------------------------------------------------------------------------------

void cigar_code_delete_nt(int nt, int direction, cigar_code_t *cigar_code) {
  int refresh = nt;
  int pos;
  cigar_op_t *op;
  int num_ops = cigar_code_get_num_ops(cigar_code);
  if (direction == 1) {
    pos = num_ops - 1;
    while (refresh > 0) {
      if (pos < 0) { break; }
      op = array_list_get(pos, cigar_code->ops);
      if (op->name != 'M' && op->name != 'I') {
	pos--;
	continue;
      } 
      if (op->number > refresh) {
	op->number -= refresh;
	assert(op->number > 0);
	break;
      } else {
	op = array_list_remove_at(pos--, cigar_code->ops);
	refresh -= op->number;
	cigar_op_free(op);
      }
    }
  } else {
    pos = 0;
    while (refresh > 0) {
      if (pos >= num_ops) { break; }
      op = array_list_get(pos, cigar_code->ops);
      if (op->name != 'M' && op->name != 'I') {
	pos++;
	continue;
      }
      if (op->number > refresh) {
	op->number -= refresh;
	assert(op->number > 0);
	break;
      }  else {
	pos++;
	refresh -= op->number;
      }
    }

    for (int i = pos - 1; i >= 0; i--) {
      cigar_op_t *op = array_list_remove_at(i, cigar_code->ops);
      cigar_op_free(op);
    }
  }

}

//--------------------------------------------------------------------------------------
/*
cigar_code_t *cigar_code_merge_sp(cigar_code_t *cc_left,
				  cigar_code_t *cc_middle, 
				  cigar_code_t *cc_right,
				  int l_flank, int r_flank) {
  
  cigar_code_t *cigar_code = cigar_code_new();
  cigar_op_t *op;
  int num_ops = cigar_code_get_num_ops(cc_left);
  assert(num_ops >= 2);

  printf("l_flank = %i, r_flank = %i\n", l_flank, r_flank);

  cigar_code_print(cc_middle);
  cigar_code_print(cc_right);

  //Delete and free last operations 'H'
  op = array_list_remove_at(cigar_code_get_num_ops(cc_left) - 1, cc_left->ops);
  assert(op->name == 'H');
  cigar_op_free(op);

  //Refresh left ops <------|
  //Remove flank, only 'I' & 'M' OPs
 
  printf("CIGAR LEFT REFRESH: \n");
  cigar_code_print(cc_left);

  //Merge left cigars operations
  num_ops = cigar_code_get_num_ops(cc_left);
  for (int i = 0; i < num_ops; i++) {
    op = array_list_get(i, cc_left->ops);
    cigar_code_append_op(op, cigar_code);
  }

  //Now, merge middle cigar operations
  for (int i = 0; i < cigar_code_get_num_ops(cc_middle); i++) {
    op = array_list_get(i, cc_middle->ops);
    cigar_code_append_op(op, cigar_code);    
  }

  printf("CIGAR MIDDLE REFRESH: \n");
  cigar_code_print(cigar_code);


  //Delete and free first operations 'H'
  op = array_list_remove_at(0, cc_right->ops);
  assert(op->name == 'H');
  cigar_op_free(op);

  //Refresh right ops |------>
  //Remove flank, only 'I' & 'M' OPs
  num_ops = cigar_code_get_num_ops(cc_right);
  refresh = r_flank;
  pos = 0;
  printf("CIGAR RIGHT REFRESH (pos value %i): \n", pos);
  cigar_code_print(cc_right);

  //Delete and free operations in right cigar
  for (int i = 0; i < pos; i++) {
    op = array_list_get(i, cc_right->ops);
    printf("\tFree %i%c\n", op->number, op->name);
    cigar_op_free(op);
  }

  //Merge right cigars operations
  num_ops = cigar_code_get_num_ops(cc_right);
  for (int i = pos; i < num_ops; i++) {
    op = array_list_get(i, cc_right->ops);
    cigar_code_append_op(op, cigar_code);
  }
  
  cigar_code_free(cc_left);
  cigar_code_free(cc_middle);
  cigar_code_free(cc_right);

  printf("FINAL CIGAR...\n");
  cigar_code_print(cigar_code);

  return cigar_code;

}
**/
//--------------------------------------------------------------------------------------

float cigar_code_get_score(int read_len, cigar_code_t *p) {
  //float ret = 0.0f;
  //int cigar_len = cigar_code_nt_length(p);
  //int distance = (abs(read_len - cigar_len) * 2) + p->distance;
  //ret = read_len - distance;
  //  LOG_DEBUG_F("score = %0.2f (distance = %i, read len = %i, cigar len = %i)\n", 
  //	      ret, p->distance, read_len, cigar_len);
  //if (ret < 0.0f) {
  //ret = 0;
  //    LOG_FATAL_F("score is negative %0.2f (distance = %i)\n", ret, p->distance);
  //}  
  
  int match_counts = 0;
  
  if (p) {
    match_counts = cigar_match_coverage(p);    
    match_counts -= p->distance; 
  }

  return 1.0f * match_counts / read_len;
}

//--------------------------------------------------------------------------------------

void init_cigar_string(cigar_code_t *p) {
  return;

  /*if (p->cigar_str) {
    free(p->cigar_str);
  }

  int num_ops = cigar_code_get_num_ops(p);
  char *str = malloc(num_ops * 5 * sizeof(char));
  *str = 0;

  cigar_op_t *op;
  for (int i = 0; i < num_ops; i++) {
    op = array_list_get(i, p->ops);
    sprintf(str, "%s%i%c", str, op->number, op->name);
  }

  p->cigar_str = str;*/
}

//--------------------------------------------------------------------------------------

cigar_code_t *generate_cigar_code(char *query_map, char *ref_map, unsigned int map_len,
				  unsigned int query_start, unsigned int ref_start,
				  unsigned int query_len, unsigned int ref_len,
				  int *distance, int ref_type) {
  
  cigar_code_t *p = cigar_code_new();



  unsigned char status;
  unsigned char transition;
  short int cigar_soft;
  short int value = 0;
  unsigned int number_op = 0;
  char operation;
  unsigned int perfect = 0;  
  unsigned int deletions_tot = 0;
  unsigned int insertions_tot = 0;
  unsigned int map_ref_len;
  unsigned int map_seq_len;
  unsigned int last_h, last_h_aux;
  int dist = 0;

  //printf("### OrigSeqLen(%d)::startSeq::%d : %s\n", query_len, query_start, query_map);
  //printf("### OrigRefLen(%d)::startRef::%d : %s LenMap(%d)\n", ref_len, ref_start, ref_map, map_len);
  
  // hard clipping start

  if (query_start > 0) {
    if (ref_type == FIRST_SW) {
      //Normal Case
      if (query_start <= 5) {
	cigar_code_append_op(cigar_op_new(query_start, 'M'), p);
      } else {
	cigar_code_append_op(cigar_op_new(query_start, 'H'), p);
	dist += query_start;
      }
    } else {
      //Middle or last ref
      if (ref_start == 0) {
	cigar_code_append_op(cigar_op_new(query_start, 'I'), p);
	dist += query_start;
      } else {
	if (ref_start == query_start) {
	  cigar_code_append_op(cigar_op_new(query_start, 'M'), p);
	} else {
	  if (ref_start > query_start) {
	    cigar_code_append_op(cigar_op_new(ref_start - query_start, 'D'), p);
	    cigar_code_append_op(cigar_op_new(query_start, 'M'), p);
	    dist += (ref_start - query_start);
	  } else {
	    cigar_code_append_op(cigar_op_new(query_start - ref_start, 'I'), p);
	    cigar_code_append_op(cigar_op_new(ref_start, 'M'), p);
	    dist += (query_start - ref_start);
	  } 
	}
      }
    }
  } else if (ref_start > 0) {
    if (ref_type != FIRST_SW) {
      cigar_code_append_op(cigar_op_new(ref_start, 'D'), p);
      dist += ref_start;
    } 
  }
  


  // first Status
  if (query_map[0] != '-' && ref_map[0] != '-') {
    status = CIGAR_MATCH_MISMATCH;
    // soft clipping
    cigar_soft = 0;
    while ((ref_map[cigar_soft] != '-') && (query_map[cigar_soft] != '-') && 
	   (ref_map[cigar_soft] != query_map[cigar_soft])) {
      cigar_soft++;
      value++;
    }
    if (value > 0) {
      cigar_code_append_op(cigar_op_new(value, 'S'), p);
    } 
  } else if (query_map[0] == '-') {
    if (ref_map[0] == '-') {
      status = CIGAR_PADDING;
    } else {
      status = CIGAR_DELETION;
    }
  } else if(ref_map[0] == '-') {
    status = CIGAR_INSERTION;
  }
  
  for (int i = value; i < map_len; i++) {
    // transition
    if (query_map[i] != '-' && ref_map[i] != '-') {
      transition = CIGAR_MATCH_MISMATCH;
      if (query_map[i] == ref_map[i]) {
        perfect++;
      } else {
	dist++;
      }
    } else if(query_map[i] == '-') {
      if (ref_map[i] == '-') {
        transition = CIGAR_PADDING;
      } else {
        transition = CIGAR_DELETION;
        deletions_tot++;
	dist++;
      }
    } else if (ref_map[i] == '-') {
      transition = CIGAR_INSERTION;
      insertions_tot++;
      dist++;
    }
    
    if (transition != status) {
      // insert operation in cigar string
      operation = select_op(status);
      cigar_code_append_op(cigar_op_new(number_op, operation), p);
      number_op = 1;
      status = transition;
    } else {
      number_op++;
    }
  }
  
  /*if ((map_len == perfect) && (perfect == query_len)) {
    status = CIGAR_PERFECT_MATCH;
    }*/
  
  operation = select_op(status);

  // hard and Soft clipped end
  if (status == CIGAR_MATCH_MISMATCH) {
    cigar_soft = map_len - 1;
    value = 0;
    while (cigar_soft >= 0 && 
	   (ref_map[cigar_soft] != '-') && (query_map[cigar_soft] != '-') && 
	   (query_map[cigar_soft] != ref_map[cigar_soft])){
      cigar_soft--;
      value++;
      //printf("(Soft %c!=%c)", output_p->mapped_ref_p[i][cigar_soft], output_p->mapped_seq_p[i][cigar_soft]);
    }
    
    cigar_code_append_op(cigar_op_new(number_op - value, operation), p);
    
    if (value > 0) {
      number_op -= value;
      cigar_code_append_op(cigar_op_new(value, 'S'), p);
    }
  } else {
    cigar_code_append_op(cigar_op_new(number_op, operation), p);
  }

  //printf("%d+%d < %d\n", length - deletions_tot, start_seq, seq_orig_len);
  //last_h = ((map_len - deletions_tot) + query_start);
  //if (last_h < query_len) {
  //cigar_code_append_op(cigar_op_new(query_len - last_h, 'H'), p);
  //}
  //printf("IN-->SW CIGAR %s\n", new_cigar_code_string(p));
  //printf("deletions_tot = %i, insertions_tot = %i\n", deletions_tot, insertions_tot);

  map_seq_len  = ((map_len - deletions_tot) + query_start);
  map_ref_len  = ((map_len - insertions_tot) + ref_start);

  //printf("query_start = %i, ref_start = %i, map_seq_len = %i, map_ref_len = %i, query_len = %i, ref_len = %i, map_len = %i\n", 
  //	 query_start, ref_start, map_seq_len, map_ref_len, query_len, ref_len, map_len);


  if (map_seq_len < query_len) {
    last_h = query_len - map_seq_len;
    //printf("last_h = %i\n", last_h);
    if (ref_type == LAST_SW) {
      //Normal Case
      //cigar_code_append_op(cigar_op_new(query_start, 'H'), p);
      if (query_start <= 5) {
	cigar_code_append_op(cigar_op_new(last_h, 'M'), p);
      } else {
	cigar_code_append_op(cigar_op_new(last_h, 'H'), p);
	dist += last_h;
      }

    } else {
      //Middle or first ref
      if (map_ref_len == ref_len) {
	cigar_code_append_op(cigar_op_new(last_h, 'I'), p);
	dist += last_h;
      } else {
	last_h_aux = ref_len - map_ref_len;
	//printf("last_h_aux = %i\n", last_h_aux);
	if (last_h_aux == last_h) {
	  cigar_code_append_op(cigar_op_new(last_h, 'M'), p);
	} else {	  
	  if (last_h_aux > last_h) {
	    cigar_code_append_op(cigar_op_new(last_h_aux - last_h, 'D'), p);
	    cigar_code_append_op(cigar_op_new(last_h, 'M'), p);
	    dist += (last_h_aux - last_h);
	  } else {
	    cigar_code_append_op(cigar_op_new(last_h - last_h_aux, 'I'), p);
	    cigar_code_append_op(cigar_op_new(last_h_aux, 'M'), p);
	    dist += (last_h - last_h_aux);
	  } 
	}
      }
    }
  } else if (map_ref_len < ref_len) {
    if (ref_type != LAST_SW) {
      dist += (ref_len - map_ref_len);
      cigar_code_append_op(cigar_op_new(ref_len - map_ref_len, 'D'), p);
    }
  }
  
  //printf("%d-%d\n", length, *number_op_tot);
  *distance = dist;

  init_cigar_string(p);
  p->distance = dist;

  return p;
}


//--------------------------------------------------------------------------------------
//        M E T A E X O N   S T R U C T U R E S   I M P L E M E N T A T I O N
//--------------------------------------------------------------------------------------

metaexon_t *metaexon_new(size_t start, size_t end) {
  metaexon_t *metaexon = (metaexon_t *)malloc(sizeof(metaexon_t));

  metaexon->start  = start;
  metaexon->end    = end;
  metaexon->left_closed = 0;
  metaexon->right_closed = 0;  

  metaexon->left_breaks = array_list_new(10, 1.25f, 
					 COLLECTION_MODE_ASYNCHRONIZED);
  metaexon->right_breaks = array_list_new(10, 1.25f, 
					 COLLECTION_MODE_ASYNCHRONIZED);

  return metaexon;

}

//-----------------------------------------------------------------------------

void metaexon_merge_breaks(void *source, void *target) {  
  //return;
  metaexon_t *meta_source = source;
  metaexon_t *meta_target = target;

  if (meta_source->left_closed) {
    meta_target->left_closed = 1;
    for (int i = 0; i < array_list_size(meta_source->left_breaks); i++) {
      void *info = array_list_get(i, meta_source->left_breaks);
      int insert = 1;
      for (int i = 0; i < array_list_size(meta_target->left_breaks); i++) {
	void *aux_info = array_list_get(i, meta_target->left_breaks);
	if (aux_info == info) {
	  insert = 0;
	  break;
	}
      }
      if (insert) {
	array_list_insert(info, meta_target->left_breaks);
      }
    }
  }
  
  if (meta_source->right_closed) {
    meta_target->right_closed = 1;
    for (int i = 0; i < array_list_size(meta_source->right_breaks); i++) {
      void *info = array_list_get(i, meta_source->right_breaks);
      int insert = 1;
      for (int i = 0; i < array_list_size(meta_target->right_breaks); i++) {
	void *aux_info = array_list_get(i, meta_target->right_breaks);
	if (aux_info == info) {
	  insert = 0;
	  break;
	}
      }
      if (insert) {
	array_list_insert(info, meta_target->right_breaks);
      }
    }
  }
  
}

int metaexon_insert_break(void *info, int type, metaexon_t *metaexon, int db_type) {
  //printf("INSERT POSITION %lu \n", ((avl_node_t *)info)->position);
  if (type == METAEXON_LEFT_END) {
    //printf("CLOSE LEFT META-EXON\n");
    metaexon->left_closed = 1;
    for (int i = 0; i < array_list_size(metaexon->left_breaks); i++) {
      void *aux_info = array_list_get(i, metaexon->left_breaks);
      if (aux_info == info) {
	return 0;
      } 
    }
    array_list_insert(info, metaexon->left_breaks);
    /*if (db_type) {
      for (int i = 0; i < array_list_size(metaexon->left_breaks); i++) {
	void *aux_info = array_list_get(i, metaexon->left_breaks);	
	//printf("START POSITION METAEXON L: %lu \n ",	  
	//     ((avl_node_t *)aux_info)->position);
      }
      //printf("NEW POSITION L: %lu \n",
      //     ((avl_node_t *)info)->position);
      return -1;
    }*/
    return 1;
  } else if (type == METAEXON_RIGHT_END) {
    //printf("CLOSE RIGHT META-EXON\n");
    metaexon->right_closed = 1;
    for (int i = 0; i < array_list_size(metaexon->right_breaks); i++) {
      void *aux_info = array_list_get(i, metaexon->right_breaks);
      if (aux_info == info) {
	return 0;
      }
    }
    array_list_insert(info, metaexon->right_breaks);
    /*if (db_type) {
      for (int i = 0; i < array_list_size(metaexon->right_breaks); i++) {
	void *aux_info = array_list_get(i, metaexon->right_breaks);	
	//printf("START POSITION METAEXON R: %lu \n ",
	//     ((avl_node_t *)aux_info)->position);
      }
      //printf("NEW POSITION R: %lu \n",
      //     ((avl_node_t *)info)->position);
      return -1;    
    }*/
    return 1;
  }  
  return 0;
}

//-----------------------------------------------------------------------------

void metaexon_free(metaexon_t *metaexon) {
  array_list_free(metaexon->left_breaks, NULL);
  array_list_free(metaexon->right_breaks, NULL);

  free(metaexon);
}

void metaexon_free_callback(void *data) {
  metaexon_t *metaexon = (metaexon_t *)data;
  
  array_list_free(metaexon->left_breaks, NULL);
  array_list_free(metaexon->right_breaks, NULL);
  
  free(metaexon);
  
}


void metaexon_set_callback(void *data, size_t start, size_t end) {
  metaexon_t *metaexon = (metaexon_t *)data;
  
  metaexon->start = start;
  metaexon->end   = end;  
}



//-----------------------------------------------------------------------------

linked_list_item_t *__metaexon_insert(linked_list_t* list_p,
				      linked_list_item_t *list_item,
				      size_t start, 
				      size_t end, 
				      size_t max_distance,
				      array_list_t *delete_items) {
  unsigned char actualization = 0;
  metaexon_t *item, *item_aux, *new_item_p;

  
  
  //size_t start = new_start;
  //size_t end   = new_end;
  
  linked_list_iterator_t* itr = linked_list_iterator_new_by_item(list_p, list_item);
  linked_list_item_t *target_item, *item_free;

  if (linked_list_size(list_p) <= 0) {
    item = metaexon_new(start, end);
    linked_list_insert(item, list_p);
    target_item = list_p->first;
    //final_list_item = list_p->first;
  } else {
    item = (metaexon_t *)linked_list_iterator_curr(itr);
    while (item != NULL) {
      if (start < item->start) {
	if (end + max_distance < item->start) {
	  /*********************************************
	   *    Case 1: New item insert before item.   *
           *                                           *
           *        new item     item                  *
           *       |-------| |--------|                *
           ********************************************/
	  new_item_p = metaexon_new(start, end);
	  linked_list_iterator_insert(new_item_p, itr);
	  linked_list_iterator_prev(itr);
	  //target_item = linked_list_iterator_list_item_curr(itr);
	  //metaexon = new_item_p;
	} else {
	  /********************************************
           *  Case 2: Actualization item start        *
           *           new item                       *
           *          |-------|   item                *
           *                   |--------|             *                            
           ********************************************/
	  item->start = start;
	  //metaexon = item;
	  if (end > item->end) {
	    /**************************************************
             *  Case 3: Actualization item start and item end *
             *          new item                              *
             *         |------------|                         *
             *              item                              *    
             *           |--------|                           *                                    
             **************************************************/
	    item->end = end;
	    actualization = 1;
	    //metaexon = item;
	  }
	  // metaexon = item;
	}
	target_item = linked_list_iterator_list_item_curr(itr);
	break;
      } else {
	if (end <= item->end) {
	  /**************************************************                                       
           *  Case 4: The new item don't insert in the list *                             
           *              item                              * 
           *         |-------------|                        * 
           *             new item                           * 
           *            |--------|                          * 
           **************************************************/
	  //metaexon = item;
	  target_item = linked_list_iterator_list_item_curr(itr);
	  break;
	} else if (item->end + max_distance >= start) {
	  /********************************************                                              
           *  Case 5: Actualization item end          *
           *            item                          *                                              
           *          |-------| new item              *                                            
           *                 |--------|               *                                              
           ********************************************/
	  item->end = end;
	  actualization = 1;
	  target_item = linked_list_iterator_list_item_curr(itr);
	  //metaexon = item;
	  break;
	}
      } // end else
      //continue loop...
      item = linked_list_iterator_next(itr);
      //item = linked_list_iterator_curr(itr);      

    } // end while

    if (item == NULL) {
      /******************************************************* 
       * Case 5: Insert new item at the end of the list      * 
       *                 item    new item                    * 
       *              |-------| |--------|                   *    
       *******************************************************/
      new_item_p = metaexon_new(start, end);
      linked_list_insert_last(new_item_p, list_p);
      target_item = list_p->last;
      //metaexon = new_item_p;
    }

    //printf("Insert OK! and now actualization\n");
    if (actualization == 1) {
      linked_list_iterator_next(itr);
      item_aux = linked_list_iterator_curr(itr);
      
      while (item_aux != NULL) {
	if (item->end + max_distance < item_aux->start) {
	  break;
	} else {
	  if (item->end < item_aux->end) {
	    item->end = item_aux->end;
	  }
	  item_free = linked_list_iterator_remove_2(itr);
	  if (item_free) {
	    array_list_insert(item_free, delete_items);
	  }
	  //if (item_free) { metaexon_free(item_free); }
	  item_aux = linked_list_iterator_curr(itr);
	}                                                                                             
      }
    }

  }//end first else

  linked_list_iterator_free(itr);

  return target_item;
}

//-----------------------------------------------------------------------------

metaexons_t *metaexons_new(unsigned int num_chromosomes, size_t *chr_size) {
  metaexons_t *metaexons = (metaexons_t *)malloc(sizeof(metaexons_t));
  //unsigned int num_chromosomes = genome->num_chromosomes;
  size_t num_chunks;
  size_t tot_chunks = 0;

  metaexons->num_chromosomes =  num_chromosomes;
  //metaexons->chunk_size      =  10;
  metaexons->num_chunks      =  (size_t *)calloc(num_chromosomes, sizeof(size_t));
  metaexons->mutex           =  (pthread_mutex_t *)calloc(num_chromosomes, sizeof(pthread_mutex_t));
  metaexons->metaexons_list  =  (linked_list_t **)calloc(num_chromosomes, sizeof(linked_list_t *));
  metaexons->bypass_pointer  =  (metaexon_pair_t **)calloc(num_chromosomes, sizeof(metaexon_pair_t *));

  //metaexons->metaexons_x  =  (linked_list_t ***)calloc(num_chromosomes, sizeof(linked_list_t **));

  //(linked_list_t *)calloc(num_chromosomes, sizeof(linked_list_t));
  for (unsigned int i = 0; i < num_chromosomes; i++) {
    num_chunks = chr_size[i] >> CHUNK_SHIFT;
    if (chr_size[i] % (int)pow(2, CHUNK_SHIFT) > 0) {
      num_chunks++;
    }
    metaexons->metaexons_list[i] =  linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    metaexons->bypass_pointer[i] = (metaexon_pair_t *)calloc(num_chunks, sizeof(metaexon_pair_t));

    //metaexons->metaexons_x[i]  =  (linked_list_t **)calloc(num_chunks, sizeof(linked_list_t *));
    //for (unsigned int j = 0; j < num_chunks; j++) {
    //metaexons->bypass_pointer[i][j] = (linked_list_item_t **)calloc(2, sizeof(linked_list_item_t *));
    //metaexons->metaexons_x[i][j]  =  linked_list_new(COLLECTION_MODE_ASYNCHRONIZED);
    //}

    metaexons->num_chunks[i]     = num_chunks;
    tot_chunks                  += num_chunks;

    pthread_mutex_init(&metaexons->mutex[i], NULL);

  }
  
  //printf("Tot chunks = %i \n", tot_chunks);

  return metaexons;

}

/*

metaexons_t *metaexons_new(unsigned int num_chromosomes, size_t *chr_size) {
  metaexons_t *metaexons = (metaexons_t *)malloc(sizeof(metaexons_t));
  //unsigned int num_chromosomes = genome->num_chromosomes;
  size_t num_chunks;
  size_t tot_chunks = 0;

  metaexons->num_chromosomes =  num_chromosomes;
  metaexons->metaexons_list  =  (skip_list_t **)calloc(num_chromosomes, sizeof(skip_list_t *));
  metaexons->mutex           =  (pthread_mutex_t *)calloc(num_chromosomes, sizeof(pthread_mutex_t));

  for (unsigned int i = 0; i < num_chromosomes; i++) {

    metaexons->metaexons_list[i] =  skip_list_new(10,
						  COLLECTION_MODE_ASYNCHRONIZED,
						  metaexon_free_callback,
						  metaexon_set_callback,
						  metaexon_merge_breaks,
						  NULL,
						  NULL);
    pthread_mutex_init(&metaexons->mutex[i], NULL);
  }
  
  
  return metaexons;

}
*/


void metaexons_free(metaexons_t *metaexons) {
  for (int i = 0; i < metaexons->num_chromosomes; i++) {
    linked_list_free(metaexons->metaexons_list[i], (void *)metaexon_free);
    free(metaexons->bypass_pointer[i]);
  }

  free(metaexons->num_chunks);
  free(metaexons->mutex);
  free(metaexons->metaexons_list);
  free(metaexons->bypass_pointer);
  
  free(metaexons);

}


/*
void metaexons_free(metaexons_t *metaexons) {
  if (metaexons) {
    for (int i = 0; i < metaexons->num_chromosomes; i++) {
      skip_list_free(metaexons->metaexons_list[i]);
    }
    
    free(metaexons->metaexons_list);
    free(metaexons->mutex);
    free(metaexons);
  }
}
*/



int metaexon_search(unsigned int strand, 
		    unsigned int chromosome,
		    size_t start, size_t end,
		    metaexon_t **metaexon_found,	    
		    metaexons_t *metaexons) {
  /*
  extern size_t search_calls;
  extern size_t insert_calls;
  extern pthread_mutex_t mutex_calls;
  
  struct timeval t_end, t_start;
  double time = 0.0;
  extern double time_search;
  extern double time_insert;

  start_timer(t_start);
  */
  int chunk_start = start >> CHUNK_SHIFT; //metaexons->chunk_size;
  int chunk_end   = end   >> CHUNK_SHIFT; //metaexons->chunk_size;
  linked_list_iterator_t itr;

  *metaexon_found = NULL;
  
  //printf("Search [%i:%lu-%lu], [%i-%i]\n", chromosome, start, end,
  //	 chunk_start, chunk_end);

  //pthread_mutex_lock(&metaexons->mutex[chromosome]);

  //Search the first reference
  int chk = chunk_start;
  linked_list_item_t *list_item;
  metaexon_t *metaexon;

  while (chk <= chunk_end && 
	 !metaexons->bypass_pointer[chromosome][chk].first) { chk++; }

  if (chk <= chunk_end && metaexons->bypass_pointer[chromosome][chk].first) {
    list_item = metaexons->bypass_pointer[chromosome][chk].first;

    linked_list_iterator_init_by_item(metaexons->metaexons_list[chromosome], list_item, &itr);   

    metaexon = linked_list_iterator_curr(&itr); 
   
    while (metaexon) {
      if (start <= metaexon->end && end >= metaexon->start) {
	//Found! We have good luck!
	*metaexon_found = metaexon;	
      } else if (start <= metaexon->end) {
	break;
      }
      metaexon = linked_list_iterator_next(&itr);
    }
  }
  
  //stop_timer(t_start, t_end, time);

  //pthread_mutex_unlock(&metaexons->mutex[chromosome]);  

  //linked_list_iterator_free(itr);
  /*
  pthread_mutex_lock(&(mutex_calls));
  search_calls++;
  time_search += time;
  pthread_mutex_unlock(&(mutex_calls));
  */
  return *metaexon_found == NULL ? 0 : 1;
  
}


/*
int metaexon_search(unsigned int strand, 
		    unsigned int chromosome,
		    size_t start, size_t end,
		    metaexon_t **metaexon_found,	    
		    metaexons_t *metaexons) {

  extern size_t search_calls;
  extern size_t insert_calls;
  extern pthread_mutex_t mutex_calls;
  
  struct timeval t_end, t_start;
  double time = 0.0;
  extern double time_search;
  extern double time_insert;

  //Extrae_event(6000019, 5);
  
  skip_list_item_t *skip_item;
  start_timer(t_start); 
  skip_item = skip_list_get_first_overlapped_item(metaexons->metaexons_list[chromosome], 
						  start, end);
  stop_timer(t_start, t_end, time);

  *metaexon_found = (skip_item) ? (skip_item)->data : NULL;
  
  //Extrae_event(6000019, 0);   

  pthread_mutex_lock(&(mutex_calls));
  search_calls++;
  time_search += time;
  pthread_mutex_unlock(&(mutex_calls));
  
  return *metaexon_found == NULL ? 0 : 1;
  
}
*/

/*
int metaexon_insert(unsigned int strand, unsigned int chromosome,
		     size_t start, size_t end, int min_intron_size, 
		     unsigned char type, void *info_break, 
		     metaexons_t *metaexons) {
  //Calculating start chunk and end chunk
  //assert(chromosome <= 25);

  if (end < start) { 
    printf("META-ERR: %lu - %lu\n", start, end); 
    exit(-1);
  }

  //assert(start < end);

  int ck_start = start >> CHUNK_SHIFT;
  int ck_end   = end  >> CHUNK_SHIFT;
  int new_ck_start, new_ck_end;

  int ck;
  metaexon_t *metaexon, *metaexon_ref;
  int db_type = strand;
  linked_list_item_t *list_item = NULL, *list_item_ref = NULL;
  int loop_start, loop_end;
  array_list_t *delete_items = array_list_new(50, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  extern double time_timer0;
  extern double time_timer1;
  extern double time_timer2;
  extern double time_timer3;

  struct timeval time_free_s, time_free_e;
  //-------->Timer0
  //start_timer(time_free_s);
  pthread_mutex_lock(&metaexons->mutex[chromosome]);
  //stop_timer(time_free_s, time_free_e, time_timer0);
  //-------->End Timer0

  //-------->Timer1
  //start_timer(time_free_s);
  for (ck = ck_start; ck <= ck_end; ck++) {
    if (metaexons->bypass_pointer[chromosome][ck].first) {
      //You have found a reference. Insert metaexon
      //printf("chr:%i, chunk :%i \n", chromosome, ck);
      list_item     = metaexons->bypass_pointer[chromosome][ck].first;
      list_item_ref = __metaexon_insert(metaexons->metaexons_list[chromosome],
					list_item, start, end,
					min_intron_size,
					delete_items);
      metaexon_ref = (metaexon_t *)list_item_ref->item;
      new_ck_start = metaexon_ref->start >> CHUNK_SHIFT;
      new_ck_end   = metaexon_ref->end >> CHUNK_SHIFT;

      loop_start = new_ck_start;
      loop_end   = new_ck_end;
      if (metaexons->bypass_pointer[chromosome][new_ck_start].first) {
	loop_start++;
	//printf("\tStar chunk actualization %i\n", new_ck_start);
	
	list_item = metaexons->bypass_pointer[chromosome][new_ck_start].first;
	assert(list_item != NULL);
	metaexon  = (metaexon_t *)list_item->item;
	if (metaexon_ref->start <= metaexon->start) {
	  metaexons->bypass_pointer[chromosome][new_ck_start].first = list_item_ref;
	}
	
	list_item = metaexons->bypass_pointer[chromosome][new_ck_start].last;
	assert(list_item != NULL);
	metaexon  = (metaexon_t *)list_item->item;
	if (metaexon_ref->end >= metaexon->end) {
	  metaexons->bypass_pointer[chromosome][new_ck_start].last = list_item_ref;
	}	
      }
      
      if (new_ck_start != new_ck_end && metaexons->bypass_pointer[chromosome][new_ck_end].first) {
	loop_end--;
	//printf("\tStar end actualization %i\n", new_ck_end);
	list_item = metaexons->bypass_pointer[chromosome][new_ck_end].first;
	assert(list_item);
	metaexon  = (metaexon_t *)list_item->item;

	if (metaexon_ref->start <= metaexon->start) {
	  metaexons->bypass_pointer[chromosome][new_ck_end].first = list_item_ref;
	}
	
	list_item = metaexons->bypass_pointer[chromosome][new_ck_end].last;
	assert(list_item);
	metaexon  = (metaexon_t *)list_item->item;
	//printf(" %lu >= %lu ??\n", metaexon_ref->start, metaexon->start);
	if (metaexon_ref->end >= metaexon->end) {
	  metaexons->bypass_pointer[chromosome][new_ck_end].last = list_item_ref;
	}	
      }
      
      for (int ck2 = loop_start; ck2 <= loop_end; ck2++) {
	metaexons->bypass_pointer[chromosome][ck2].first = list_item_ref;
	metaexons->bypass_pointer[chromosome][ck2].last = list_item_ref;
      }      
      break;
    }
  } //end for

  //stop_timer(time_free_s, time_free_e, time_timer1);
  //-------->End Timer1

  //-------->Timer2
  //start_timer(time_free_s);
  if (!list_item) {
    //printf("List item not insert :( ...\n");
    int ck_s = ck_start - 1;
    int ck_e = ck_end + 1;
    for (int i = 0; i < metaexons->num_chunks[chromosome]; i++) {
      if (ck_s - i >= 0 && 
	  metaexons->bypass_pointer[chromosome][ck_s - i].last ) {
	//printf("chunk 2: %i \n", ck_s - i);
	list_item     = metaexons->bypass_pointer[chromosome][ck_s - i].last; 
	list_item_ref = __metaexon_insert(metaexons->metaexons_list[chromosome],
					  list_item, start, end,
					  min_intron_size,
					  delete_items);      	
	break;
      }
      if (ck_e + i < metaexons->num_chunks[chromosome] &&
	  metaexons->bypass_pointer[chromosome][ck_e + i].first) {
	//printf("chunk 2: %i \n", ck_e + i);
	list_item     = metaexons->bypass_pointer[chromosome][ck_e + i].first;
	list_item_ref = __metaexon_insert(metaexons->metaexons_list[chromosome],
					  list_item, start, end,
					  min_intron_size,
					  delete_items);
	break;
      }
    }

    if (!list_item) {
      //printf("chunk 3: NO\n");
      list_item_ref = __metaexon_insert(metaexons->metaexons_list[chromosome],
					list_item, start, end,
					min_intron_size,
					delete_items);
    }

    metaexon_ref = (metaexon_t *)list_item_ref->item;
    new_ck_start = metaexon_ref->start >> CHUNK_SHIFT;
    new_ck_end   = metaexon_ref->end >> CHUNK_SHIFT;

    for (int ck = new_ck_start; ck <= new_ck_end; ck++) {
      metaexons->bypass_pointer[chromosome][ck].first = list_item_ref;
      metaexons->bypass_pointer[chromosome][ck].last  = list_item_ref;
    }

  }
  //stop_timer(time_free_s, time_free_e, time_timer2);
  //-------->End Timer2

  if (info_break) {
    metaexon_insert_break(info_break, type, metaexon_ref, db_type);    
  }

  pthread_mutex_unlock(&metaexons->mutex[chromosome]);

  //-------->Timer3
  //start_timer(time_free_s);
  for (int i = array_list_size(delete_items) - 1; i >= 0; i--) {
    list_item = array_list_get(i, delete_items);
    linked_list_item_free(list_item, (void (*)(void *))metaexon_free);
  }

  array_list_free(delete_items, (void *)NULL);
  //stop_timer(time_free_s, time_free_e, time_timer3);
  //-------->End Timer3

  return db_type;

}
*/

/*
int metaexon_insert(unsigned int strand, unsigned int chromosome,
		    size_t start, size_t end, int min_intron_size, 
		    unsigned char type, void *info_break, 
		    metaexons_t *metaexons) {

  extern size_t search_calls;
  extern size_t insert_calls;
  extern pthread_mutex_t mutex_calls;
  
  struct timeval t_stop, t_start;
  double time = 0.0;
  extern double time_search;
  extern double time_insert;

  metaexon_t *metaexon;
  
  metaexon = metaexon_new(start, end);

  //Extrae_event(6000019, 4);   
  pthread_mutex_lock(&(metaexons->mutex[chromosome]));
  //Extrae_event(6000019, 0);   
  
  //Extrae_event(6000019, 6);   
  start_timer(t_start); 
  skip_list_insert_item(metaexons->metaexons_list[chromosome], metaexon, start, end, min_intron_size);
  stop_timer(t_start, t_stop, time);
  //Extrae_event(6000019, 0);

  if (info_break) {
    metaexon_insert_break(info_break, type, metaexon, strand);
  }

  pthread_mutex_unlock(&(metaexons->mutex[chromosome]));


  pthread_mutex_lock(&(mutex_calls));
  insert_calls++;
  time_insert += time;
  pthread_mutex_unlock(&(mutex_calls));

  return 1;

}
*/


int metaexon_insert(unsigned int strand, unsigned int chromosome,
		     size_t start, size_t end, int min_intron_size, 
		     unsigned char type, void *info_break, 
		     metaexons_t *metaexons) {

    //assert(chromosome <= 25);
  /*
  extern size_t search_calls;
  extern size_t insert_calls;
  extern pthread_mutex_t mutex_calls;
  
  struct timeval t_stop, t_start;
  double time = 0.0;
  extern double time_search;
  extern double time_insert;
  */
  if (end < start) {
    printf("META-ERR: %lu - %lu\n", start, end);
    exit(-1);
  }

    //assert(start < end);

    int ck_start = start >> CHUNK_SHIFT;
    int ck_end = end >> CHUNK_SHIFT;
    int ck;
    metaexon_t *metaexon, *metaexon_ref;
    int db_type = strand;
    linked_list_item_t *list_item = NULL, *list_item_ref = NULL;
    int loop_start, loop_end;
    array_list_t *delete_items = array_list_new(50, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

    pthread_mutex_lock(&metaexons->mutex[chromosome]);

    //start_timer(t_start);

    for (ck = ck_start; ck <= ck_end; ck++) {
      if (metaexons->bypass_pointer[chromosome][ck].first) {
	//You have found a reference. Insert metaexon
	//printf("chr:%i, chunk :%i \n", chromosome, ck);
	list_item = metaexons->bypass_pointer[chromosome][ck].first;
	list_item_ref = __metaexon_insert(metaexons->metaexons_list[chromosome],
					  list_item, start, end,
					  min_intron_size,
					  delete_items);
	metaexon_ref = (metaexon_t *)list_item_ref->item;
	ck_start = metaexon_ref->start >> CHUNK_SHIFT;
	ck_end = metaexon_ref->end >> CHUNK_SHIFT;

	loop_start = ck_start;
	loop_end = ck_end;
	if (metaexons->bypass_pointer[chromosome][ck_start].first) {
	  loop_start++;
	  //printf("\tStar chunk actualization %i\n", ck_start);

	  list_item = metaexons->bypass_pointer[chromosome][ck_start].first;
	  assert(list_item != NULL);
	  metaexon = (metaexon_t *)list_item->item;
	  if (metaexon_ref->start <= metaexon->start) {
	    metaexons->bypass_pointer[chromosome][ck_start].first = list_item_ref;
	  }

	  list_item = metaexons->bypass_pointer[chromosome][ck_start].last;
	  assert(list_item != NULL);
	  metaexon = (metaexon_t *)list_item->item;
	  if (metaexon_ref->end >= metaexon->end) {
	    metaexons->bypass_pointer[chromosome][ck_start].last = list_item_ref;
	  }
	}
      
	if (ck_start != ck_end && metaexons->bypass_pointer[chromosome][ck_end].first) {
	  loop_end--;
	  //printf("\tStar end actualization %i\n", ck_end);
	  list_item = metaexons->bypass_pointer[chromosome][ck_end].first;
	  assert(list_item);
	  metaexon = (metaexon_t *)list_item->item;

	  if (metaexon_ref->start <= metaexon->start) {
	    metaexons->bypass_pointer[chromosome][ck_end].first = list_item_ref;
	  }

	  list_item = metaexons->bypass_pointer[chromosome][ck_end].last;
	  assert(list_item);
	  metaexon = (metaexon_t *)list_item->item;
	  //printf(" %lu >= %lu ??\n", metaexon_ref->start, metaexon->start);
	  if (metaexon_ref->end >= metaexon->end) {
	    metaexons->bypass_pointer[chromosome][ck_end].last = list_item_ref;
	  }
	}
      
	for (int ck2 = loop_start; ck2 <= loop_end; ck2++) {
	  metaexons->bypass_pointer[chromosome][ck2].first = list_item_ref;
	  metaexons->bypass_pointer[chromosome][ck2].last = list_item_ref;
	}
      
	break;
      }
    }
  
    if (!list_item) {
      //printf("List item not insert :( ...\n");
      int ck_s = ck_start - 1;
      int ck_e = ck_end + 1;
      for (int i = 0; i < metaexons->num_chunks[chromosome]; i++) {
	if (ck_s - i >= 0 &&
	    metaexons->bypass_pointer[chromosome][ck_s - i].last ) {
	  //printf("chunk 2: %i \n", ck_s - i);
	  list_item = metaexons->bypass_pointer[chromosome][ck_s - i].last;
	  list_item_ref = __metaexon_insert(metaexons->metaexons_list[chromosome],
					    list_item, start, end,
					    min_intron_size,
					    delete_items);
	  break;
	}
	if (ck_e + i < metaexons->num_chunks[chromosome] &&
	    metaexons->bypass_pointer[chromosome][ck_e + i].first) {
	  //printf("chunk 2: %i \n", ck_e + i);
	  list_item = metaexons->bypass_pointer[chromosome][ck_e + i].first;
	  list_item_ref = __metaexon_insert(metaexons->metaexons_list[chromosome],
					    list_item, start, end,
					    min_intron_size,
					    delete_items);
	  break;
	}
      }

      if (!list_item) {
	//printf("chunk 3: NO\n");
	list_item_ref = __metaexon_insert(metaexons->metaexons_list[chromosome],
					  list_item, start, end,
					  min_intron_size,
					  delete_items);
      }

      metaexon_ref = (metaexon_t *)list_item_ref->item;
      //ck_start = metaexon_ref->start >> CHUNK_SHIFT;
      //ck_end = metaexon_ref->end >> CHUNK_SHIFT;

      for (int ck = ck_start; ck <= ck_end; ck++) {
	metaexons->bypass_pointer[chromosome][ck].first = list_item_ref;
	metaexons->bypass_pointer[chromosome][ck].last = list_item_ref;
      }

    }

    if (info_break) {
      metaexon_insert_break(info_break, type, metaexon_ref, db_type);
    }

    pthread_mutex_unlock(&metaexons->mutex[chromosome]);

    for (int i = array_list_size(delete_items) - 1; i >= 0; i--) {
      list_item = array_list_get(i, delete_items);
      //metaexon = (metaexon_t *)list_item->item;
      //metaexon_merge_breaks(metaexon, metaexon_ref);
      
      linked_list_item_free(list_item, (void *)metaexon_free);
    }

    array_list_free(delete_items, (void *)NULL);

    /*
    stop_timer(t_start, t_stop, time);

    pthread_mutex_lock(&(mutex_calls));
    insert_calls++;
    time_insert += time;
    pthread_mutex_unlock(&(mutex_calls));
    */

    return db_type;

}
