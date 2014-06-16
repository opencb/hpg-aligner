#include "rna_splice.h"

extern size_t junction_id;
size_t cannonical_sp = 0;
size_t semi_cannonical_sp = 0;
size_t total_splice = 0;


//--------------------------------------------------------------------------------

intron_t *intron_new(unsigned char strand, unsigned int chromosome, size_t start, size_t end) {
  intron_t *intron = (intron_t *)malloc(sizeof(intron_t));

  intron->strand = strand;
  intron->chromosome = chromosome;
  intron->start = start;
  intron->end = end;

  return intron;

}

//--------------------------------------------------------------------------------

void intron_free(intron_t *intron) {
  free(intron);
}

//--------------------------------------------------------------------------------

void avl_end_process(cp_avlnode *node) {
  if (node->left) { avl_end_process(node->left); }
  avl_node_t *avl_node = (avl_node_t *)node->value;
  end_data_t *data = avl_node->data;
  for (size_t i = 0; i < array_list_size(data->list_starts); i++) {
    size_t end_sp = (size_t)array_list_get(i, data->list_starts);
    printf("%i-%i\n", end_sp, avl_node->position);
  }
  if (node->right) { avl_end_process(node->right); }
  return;
}


//--------------------------------------------------------------------------------

void avl_process(cp_avlnode *node) {
  if (node->left) { avl_process(node->left); }
  avl_node_t *avl_node = (avl_node_t *)node->value;
  start_data_t *data = avl_node->data;
  for (int i = 0; i < array_list_size(data->list_ends); i++) {
    splice_end_t *end_sp = array_list_get(i, data->list_ends);
    printf("%i-%i\n", avl_node->position, end_sp->end);
  }
  if (node->right) { avl_process(node->right); }
  return;
}


static inline void insert_ends(unsigned char strand, start_data_t *start_data, size_t start, array_list_t *intron_list) {
  //Node is in the range
  int num_ends = array_list_size(start_data->list_ends);

  for(size_t i = 0; i < num_ends; i++) {
    splice_end_t *end_sp = array_list_get(i, start_data->list_ends);
    intron_t *new_intron = intron_new(strand, 0, start - 1, end_sp->end + 1);
    array_list_insert(new_intron, intron_list);
  }

}

static inline void insert_starts(unsigned char strand, end_data_t *end_data, size_t end, array_list_t *intron_list) {
  //Node is in the range
  int num_starts = array_list_size(end_data->list_starts);

  for(size_t i = 0; i < num_starts; i++) {
    size_t start = (size_t)array_list_get(i, end_data->list_starts);
    intron_t *new_intron = intron_new(strand, 0, start - 1, end + 1);
    array_list_insert(new_intron, intron_list);
  }

}


array_list_t *search_candidate_sp_avl(unsigned char strand, size_t lim_start, 
				      size_t lim_end, array_list_t *intron_list, 
				      cp_avlnode *node, unsigned char type_search) {
  avl_node_t *avl_node;
  size_t position;

  if (node) {
    avl_node = node->value;
    position = avl_node->position;
    if (position >= lim_start && position <= lim_end) {
      if (type_search == SEARCH_STARTS) {
	insert_starts(strand, avl_node->data, position, intron_list);
      } else {
	insert_ends(strand, avl_node->data, position, intron_list);
      }
      if (node->left)  { search_candidate_sp_avl(strand, lim_start, lim_end, intron_list, node->left, type_search);  }
      if (node->right) { search_candidate_sp_avl(strand, lim_start, lim_end, intron_list, node->right, type_search); }    
    } else {
      if (position <= lim_start) {
	if (node->right) { search_candidate_sp_avl(strand, lim_start, lim_end, intron_list, node->right, type_search); }
      } else {
	if (node->left) { search_candidate_sp_avl(strand, lim_start, lim_end, intron_list, node->left, type_search); }
      }
    }
  }

  return intron_list;

}

//--------------------------------------------------------------------------------

void load_intron_file(genome_t *genome, char* intron_filename, avls_list_t *avls) {
  FILE *fd_intron = fopen(intron_filename, "r");
  size_t max_size = 2048;
  size_t line_size;
  char line[max_size];
  char value[max_size];
  int pos, value_pos;
  size_t chr, start, end;
  unsigned char strand;
  int found;
  node_element_splice_t *node;

  if (fd_intron == NULL) {
    return;
  }
  
  while (!feof(fd_intron)) {
    fgets(line, max_size, fd_intron);
    line_size = strlen(line);
    pos = 0;
    value_pos = 0;
    
    while (line[pos] != '\t' && pos < line_size) {
      value[value_pos++] = line[pos++];
    }
    value[value_pos] = '\0';

    found = 0;
    for (chr = 0; chr < genome->num_chromosomes; chr++) {
      if (strcmp(value, genome->chr_name[chr]) == 0) { found = 1; break; }
    }

    if (!found) { continue; }
    
    pos++;
    value_pos = 0;
    while (line[pos] != '\t' && pos < line_size) {
      value[value_pos++] = line[pos++];
    }
    value[value_pos] = '\0';
    start = atoi(value);

    pos++;
    value_pos = 0;
    while (line[pos] != '\t' && pos < line_size) {
      value[value_pos++] = line[pos++];
    }
    value[value_pos] = '\0';
    end = atoi(value);

    pos++;
    value_pos = 0;
    while (line[pos] != '\t' && pos < line_size) {
      value[value_pos++] = line[pos++];
    }
    value[value_pos] = '\0';
    if (strcmp(value, "-1") == 0 || strcmp(value, "-")) {
      strand = 1;
    } else { 
      strand = 0; 
    }

    allocate_start_node(chr, strand, start, 
			end, start, end, FROM_FILE, 
			-1, NULL, NULL, NULL, avls);
  }
}

//--------------------------------------------------------------------------------

void splice_end_type_new(char type_sp, char *splice_nt, splice_end_t *splice_end) {
  switch(type_sp) {
  case UNKNOWN_SPLICE:
    splice_end->splice_nt = strdup(splice_nt);
    break;
  case GT_AG_SPLICE:
    splice_end->splice_nt = strdup("GT-AG");
    break;
  case CT_AC_SPLICE:
    splice_end->splice_nt = strdup("CT-AC");
    break;
  case AT_AC_SPLICE:
    splice_end->splice_nt = strdup("AT-AC");
    break;
  case GT_AT_SPLICE:
    splice_end->splice_nt = strdup("GT-AT");
    break;
  case GC_AG_SPLICE:
    splice_end->splice_nt = strdup("GC-AG");
    break;
  case CT_GC_SPLICE:
    splice_end->splice_nt = strdup("CT-GC");
    break;
  default:
    splice_end->splice_nt = strdup("NONE");
    break;
  }
}

splice_end_t *splice_end_new(size_t end, size_t end_extend, unsigned char type_orig, char type_sp, char *splice_nt) {
  splice_end_t *splice_end = (splice_end_t *)malloc(sizeof(splice_end_t));
  if (type_orig == FROM_READ) {
    splice_end->reads_number = 1;
  } else {
    splice_end->reads_number = 0;
    //printf("From file\n");
  }
  splice_end->end = end;
  splice_end->end_extend = end_extend;
  splice_end->origin = type_orig;
  //splice_end->type_sp = type_sp;
  splice_end_type_new(type_sp, splice_nt, splice_end);

  return splice_end;
}

void splice_end_free(splice_end_t *splice_end) { 
  if (splice_end->splice_nt) { free(splice_end->splice_nt); }
  free(splice_end);
}

//--------------------------------------------------------------------------------

end_data_t *end_data_new() {
  end_data_t *end_data = (end_data_t *)malloc(sizeof(end_data_t));
  end_data->list_starts = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  return end_data;
}

void end_data_free(end_data_t *data) {
  array_list_free(data->list_starts, NULL);
  free(data);
}

//--------------------------------------------------------------------------------

start_data_t *start_data_new() {
  start_data_t *start_data = (start_data_t *)malloc(sizeof(start_data_t));
  start_data->list_ends = array_list_new(10, 1.25f, COLLECTION_MODE_ASYNCHRONIZED);

  return start_data;
}

void start_data_free(start_data_t *data) {
  array_list_free(data->list_ends, (void *)splice_end_free);
  free(data);
}

//--------------------------------------------------------------------------------

avl_node_t *avl_node_end_new(size_t position) {
  end_data_t *end_data = end_data_new();
  avl_node_t *node = (avl_node_t *)malloc(sizeof(avl_node_t));

  node->position = position;
  node->data = end_data;

  return node;
}

void avl_node_end_free(avl_node_t *avl_node) {
  end_data_free((end_data_t *)avl_node->data);
  free(avl_node);
}


avl_node_t *avl_node_new(size_t position) {
  start_data_t *start_data = start_data_new();
  avl_node_t *node = (avl_node_t *)malloc(sizeof(avl_node_t));

  node->position = position;
  node->data = start_data;

  return node;
}

void avl_node_free(avl_node_t *avl_node) {
  start_data_free((start_data_t *)avl_node->data);
  free(avl_node);
}

int avl_node_compare(avl_node_t* a, size_t b) {
  if(a->position == b) { 
    return 0;
  }else if(a->position < b) {
   return -1;
  }else{
    return 1;
  }
}


//--------------------------------------------------------------------------------

void allocate_end_splice(size_t end, size_t end_extend, int type_orig, char type_sp, start_data_t *data, char *splice_nt) {
  int num_ends = array_list_size(data->list_ends);
  splice_end_t *splice_end;

  //printf("Insert SP:%i\n", type_sp);
  for (int i = 0; i < num_ends; i++) {    
    splice_end = array_list_get(i, data->list_ends);
    if (splice_end->end == end) { 
      if (splice_end->end_extend < end_extend) { 
	splice_end->end_extend = end_extend;
      }
      if (splice_end->splice_nt) {
	  free(splice_end->splice_nt);
      }
      splice_end_type_new(type_sp, splice_nt, splice_end);
      
      if (type_orig != FROM_FILE) { splice_end->reads_number++; }
      return;
    }
  }

  splice_end = splice_end_new(end, end_extend, type_orig, type_sp, splice_nt);
  array_list_insert(splice_end, data->list_ends);
  return;
}

void allocate_start_splice(size_t start, end_data_t *data) {
  int num_starts = array_list_size(data->list_starts);
  size_t start_list;

  for (int i = 0; i < num_starts; i++) {
    start_list = (size_t)array_list_get(i, data->list_starts);
    if (start == start_list) { return; }    
  }

  array_list_insert((void *)start, data->list_starts);

  return;
}


//--------------------------------------------------------------------------------

void allocate_start_node(unsigned int chromosome, unsigned char strand, 
			 size_t start_extend, size_t end_extend, 
			 size_t start, size_t end, unsigned char type_orig,
			 char type_sp, char *splice_nt, 
			 avl_node_t **ref_node_start, avl_node_t **ref_node_end, 
			 avls_list_t *avls_list) {

  pthread_mutex_lock(&(avls_list->avls[strand][chromosome].mutex));
  //printf("Insert %lu - %lu\n", start, end);

  if (start > end) { 
    //fprintf(stderr, "ERROR [%i:%i]!!! START END AVL %lu vs %lu", 
    //	    strand, chromosome, start, end); 
    //exit(-1); 
    pthread_mutex_unlock(&(avls_list->avls[strand][chromosome].mutex));

    return;
  }

  cp_avltree *avl = avls_list->avls[strand][chromosome].avl;
  avl_node_t *node_start, *node_end;
  start_data_t *start_data;

  node_start = (avl_node_t *)cp_avltree_get(avl, (void *)start);
  if(node_start == NULL) {
    //printf("\tNot Exist S\n");
    node_start = cp_avltree_insert(avl, (void *)start, (void *)start);
    start_data = (start_data_t *)node_start->data;
    start_data->start_extend = start_extend;
  } else {
    //printf("\tExist S\n");
    start_data = (start_data_t *)node_start->data;
    if (start_data->start_extend > start_extend) {
    start_data->start_extend = start_extend;
    }
  }

  allocate_end_splice(end, end_extend, type_orig, type_sp, start_data, splice_nt);

  //For Extra speed we insert all ends in the other avl 
  avl = avls_list->ends_avls[strand][chromosome].avl;
  node_end = (avl_node_t *)cp_avltree_get(avl, (void *)end);

  if(node_end == NULL) {
    //printf("\tNot Exist E\n");
    node_end = cp_avltree_insert(avl, (void *)end, (void *)end);
  }  //else {
  //printf("\tExist E\n");
  //}

  allocate_start_splice(start, (end_data_t *)node_end->data);

  if (ref_node_start != NULL && ref_node_end != NULL) {
    *ref_node_start = node_start;
    *ref_node_end = node_end;
  }

  pthread_mutex_unlock(&(avls_list->avls[strand][chromosome].mutex));

}

//--------------------------------------------------------------------------------

avls_list_t* avls_list_new(size_t num_chromosomes) {
  avls_list_t* avls_list = (avls_list_t *)malloc(sizeof(avls_list_t));
  size_t i, st;
  avls_list->avls = (avl_tree_t **)malloc(sizeof(avl_tree_t *)*NUM_STRANDS);

  for(st = 0; st < NUM_STRANDS; st++) {
    avls_list->avls[st] = (avl_tree_t *)malloc(sizeof(avl_tree_t)*num_chromosomes);
    for(i = 0; i < num_chromosomes; i++) {
      avls_list->avls[st][i].avl = cp_avltree_create_by_option(COLLECTION_MODE_NOSYNC |
								COLLECTION_MODE_COPY |
								COLLECTION_MODE_DEEP,
								(cp_compare_fn) avl_node_compare,
								(cp_copy_fn) avl_node_new,
								(cp_destructor_fn)avl_node_free,
								(cp_copy_fn) avl_node_new,
								(cp_destructor_fn)avl_node_free);
      pthread_mutex_init(&(avls_list->avls[st][i].mutex), NULL);
    }
  }

  //For Extra Speed when search splice junctions supracals
  avls_list->ends_avls = (avl_tree_t **)malloc(sizeof(avl_tree_t *)*NUM_STRANDS);

  for(st = 0; st < NUM_STRANDS; st++) {
    avls_list->ends_avls[st] = (avl_tree_t *)malloc(sizeof(avl_tree_t)*num_chromosomes);
    for(i = 0; i < num_chromosomes; i++) {
      avls_list->ends_avls[st][i].avl = cp_avltree_create_by_option(COLLECTION_MODE_NOSYNC |
                                                                       COLLECTION_MODE_COPY |
                                                                       COLLECTION_MODE_DEEP,
                                                                       (cp_compare_fn) avl_node_compare,
                                                                       (cp_copy_fn) avl_node_end_new,
                                                                       (cp_destructor_fn)avl_node_end_free,
                                                                       (cp_copy_fn) avl_node_end_new,
                                                                       (cp_destructor_fn)avl_node_end_free);
      pthread_mutex_init(&(avls_list->ends_avls[st][i].mutex), NULL);
    }
  }

  return avls_list;
}

//--------------------------------------------------------------------------------


allocate_buffers_t* process_avlnode_ends(avl_node_t *node_val, unsigned char st,
					 unsigned int chromosome, allocate_buffers_t *allocate_batches) {
  int i;
  char strand[2] = {'+', '-'};
  list_item_t* item_p = NULL;
  unsigned int bytes_exact, bytes_extend;
  start_data_t *start_data = node_val->data;
  size_t num_ends = array_list_size(start_data->list_ends);
  const int write_size = 5000000;

  allocate_batches->write_exact_sp;

  for (i = 0; i < num_ends; i++) {
    if ((allocate_batches->write_exact_sp->size + 100) > write_size) {
      fwrite((char *)allocate_batches->write_exact_sp->buffer_p, allocate_batches->write_exact_sp->size, 1, allocate_batches->fd_exact);
      //fwrite((char *)allocate_batches->write_extend_sp->buffer_p, allocate_batches->write_extend_sp->size, 1, allocate_batches->fd_extend);
      allocate_batches->write_exact_sp->size = 0;
      //allocate_batches->write_extend_sp->size = 0;
    } 
    splice_end_t *end_sp = array_list_get(i, start_data->list_ends);
    if (end_sp->reads_number) {
      //printf("%lu - %lu\n", node_val->position, end_sp->end);
      bytes_exact = pack_junction(chromosome, st, 
				  node_val->position, end_sp->end, 
				  junction_id, end_sp->reads_number, 
				  end_sp->splice_nt,
				  &(((char *)allocate_batches->write_exact_sp->buffer_p)[allocate_batches->write_exact_sp->size]));
      
      /*bytes_extend = pack_junction(chromosome, st, start_data->start_extend, 
				   end_sp->end_extend, junction_id, end_sp->reads_number, 
				   end_sp->splice_nt,
				   &(((char *)allocate_batches->write_extend_sp->buffer_p)[allocate_batches->write_extend_sp->size])); 
      */
      allocate_batches->write_exact_sp->size += bytes_exact;
      //allocate_batches->write_extend_sp->size += bytes_extend;
      
      total_splice += end_sp->reads_number;
      junction_id++;
      //if (end_sp->type_sp >= AT_AC_SPLICE) {
      //semi_cannonical_sp++;
      //} else {
      //cannonical_sp++;
      //}
    }
  }

  return allocate_batches;

}

allocate_buffers_t * process_avlnode(cp_avlnode *node, unsigned char st, 
				     unsigned int chromosome, allocate_buffers_t *allocate_batches) {
  if (node->left) { 
    allocate_batches = process_avlnode(node->left, st, chromosome, allocate_batches);
  }
  
  allocate_batches = process_avlnode_ends((avl_node_t *)node->value, st, 
					  chromosome, allocate_batches);
  
  if (node->right) {
    allocate_batches = process_avlnode(node->right, st, chromosome, allocate_batches);
  }
  
  return allocate_batches;

}

//--------------------------------------------------------------------------------


void write_chromosome_avls( char *extend_sp, char *exact_sp, 
			    size_t num_chromosomes, avls_list_t *avls_list) {
  int c, chr;
  unsigned char st;
  allocate_buffers_t *allocate_batches = (allocate_buffers_t *)malloc(sizeof(allocate_buffers_t));
  const int write_size = 5000000;

  FILE *fd_exact = fopen(exact_sp, "w");  
  if (!fd_exact) { 
    printf("Imposible to create FILE: %s\n", exact_sp);
    exit(-1);
  }

  FILE *fd_extend = NULL;//fopen(extend_sp, "w");
  //if (!fd_extend) { 
  //printf("Imposible to create FILE: %s\n", exact_sp);
  //exit(-1);
  //}

  allocate_batches->fd_exact = fd_exact;
  //allocate_batches->fd_extend = fd_extend;

  allocate_batches->write_exact_sp  = write_batch_new(write_size, SPLICE_EXACT_FLAG);
  //allocate_batches->write_extend_sp  = write_batch_new(write_size, SPLICE_EXTEND_FLAG);

  for(st = 0; st < NUM_STRANDS; st++) {      
    for(c = 0; c < num_chromosomes; c++) {
      if(avls_list->avls[st][c].avl->root != NULL) {
	chr = c + 1;	
	//printf("Chromosome %i(%i)\n", chr, st);
	allocate_batches = process_avlnode(avls_list->avls[st][c].avl->root, st, chr, allocate_batches);
	
	if(allocate_batches->write_exact_sp != NULL) {
	  if(allocate_batches->write_exact_sp->size > 0) {	  
	    fwrite((char *)allocate_batches->write_exact_sp->buffer_p, allocate_batches->write_exact_sp->size, 1, allocate_batches->fd_exact);
	    //fwrite((char *)allocate_batches->write_extend_sp->buffer_p, allocate_batches->write_extend_sp->size, 1, allocate_batches->fd_extend);
	    allocate_batches->write_exact_sp->size = 0;
	    //allocate_batches->write_extend_sp->size = 0;
	  }
	}
      } //end IF chromosome splice not NULL
      cp_avltree_destroy(avls_list->avls[st][c].avl);
      cp_avltree_destroy(avls_list->ends_avls[st][c].avl);
    }
    free(avls_list->avls[st]);
    free(avls_list->ends_avls[st]);
  }

  free(avls_list->avls);
  free(avls_list->ends_avls);

  free(avls_list);
  //fclose(fd_extend);
  fclose(fd_exact);
  
  write_batch_free(allocate_batches->write_exact_sp);
  //write_batch_free(allocate_batches->write_extend_sp);

  free(allocate_batches);

  //basic_statistics_sp_init(total_splice, cannonical_sp, semi_cannonical_sp, basic_st);
    
}

/*

int search_end_splice(node_element_splice_t *node, size_t end) {
  for (int i = 0; i < node->number_allocate_ends; i++) {
    //printf("%i == %i, %i == %i\n", node->allocate_ends[i]->end, end, node->allocate_ends[i]->strand, strand);
    if (node->allocate_ends[i]->end == end) { 
      return 1;
    }
  }
  return 0;
}


array_list_t *search_end_splice_avl(size_t read_length, unsigned char strand, 
				    size_t chromosome, size_t end, 
				    cp_avlnode *node, array_list_t *intron_list) {
  size_t start_limit = 0;
  size_t end_limit = 0;

  if (!node) {  return NULL; }

  if (node->left) {
    search_end_splice_avl(read_length, strand, chromosome, end, node->left, intron_list);
  }
  
  if (node) {
    //printf("Node study\n");
    node_element_splice_t *ends_array = node->value;
    start_limit = ends_array->splice_start;
    for(size_t i = 0; i < ends_array->number_allocate_ends; i++) {
      end_limit = ends_array->allocate_ends[i]->end;
      if (((end_limit - read_length) <= end ) && (end <= (end_limit + read_length))) {
	//Found Start. Now we create one CAL for each end in the Node	
	intron_t *new_intron = intron_new(start_limit - 1, end_limit + 1, strand);
	array_list_insert(new_intron, intron_list);
      }
    }
  }

  if (node->right) {
    search_end_splice_avl(read_length, strand, chromosome, end, node->right, intron_list);
  }

  return intron_list;
}


array_list_t *search_start_splice_avl(size_t read_length, unsigned char strand, 
				      size_t chromosome, size_t start, 
				      cp_avlnode *node, array_list_t *intron_list) {
  size_t start_limit = 0;
  node_element_splice_t *ends_array;

  if (!node) {  return NULL; }
  
  if (node->left) {
    search_start_splice_avl(read_length, strand, chromosome, start, node->left, intron_list);
  }

  if (node) {
    ends_array = node->value;
    start_limit = ends_array->splice_start;
  }
  //printf("Compare AVL start = %llu with start = %llu \n", start_limit, start);
  if (((start_limit - read_length) <= start ) && (start <= (start_limit + read_length))) {
    //Found Start. Now we create one CAL for each end in the Node

    for(size_t i = 0; i < ends_array->number_allocate_ends; i++) {
      intron_t *new_intron = intron_new(start_limit - 1, ends_array->allocate_ends[i]->end + 1, strand);
      array_list_insert(new_intron, intron_list);
    }
  }

  if (node->right) {
    search_start_splice_avl(read_length, strand, chromosome, start, node->right, intron_list);
  }

  return intron_list;
}




int node_compare(node_element_splice_t* a, size_t b) {
  if(a->splice_start == b){ 
    return 0;
  }else if(a->splice_start < b){
   return -1;
  }else{
    return 1;
  }
}

node_element_splice_t* node_copy(size_t b) {
 node_element_splice_t *node = (node_element_splice_t *)malloc(sizeof(node_element_splice_t));

 //assert(node != NULL);
 
 node->maximum_allocate_ends = 10;
 node->number_allocate_ends = 0;
 node->allocate_ends = (splice_end_t **)malloc(node->maximum_allocate_ends * sizeof(splice_end_t *));
 if (node->allocate_ends == NULL){ exit(-1); }
 
 node->splice_start = b;

 return node;
}

void node_free(node_element_splice_t* a) {  
  for (unsigned int i = 0; i < a->number_allocate_ends; i++) {
    free_splice_end(a->allocate_ends[i]);
  }

  free(a->allocate_ends);
  free(a);
}


node_element_splice_t* insert_end_splice(splice_end_t *splice_end_p, node_element_splice_t *element_p){
  element_p->allocate_ends[element_p->number_allocate_ends] = splice_end_p;
  element_p->number_allocate_ends++;

  if (element_p->number_allocate_ends >= element_p->maximum_allocate_ends) {
    element_p->maximum_allocate_ends = element_p->maximum_allocate_ends * 2; 
    element_p->allocate_ends = (splice_end_t **) realloc (element_p->allocate_ends, element_p->maximum_allocate_ends * sizeof (splice_end_t *));
    if (element_p->allocate_ends == NULL) { exit(-1); }
  }

  return element_p;
}


node_element_splice_t* search_and_insert_end_splice(unsigned int chromosome,
						    size_t end, size_t splice_start, 
						    size_t splice_end, int type_orig, 
						    int type_sp, node_element_splice_t *element_p) {
  unsigned int i;
  
  if (element_p->splice_start_extend > splice_start) {
    element_p->splice_start_extend = splice_start;
  }	  
  
  for (i = 0; i < element_p->number_allocate_ends; i++) {
    if (element_p->allocate_ends[i]->end == end) { 

      element_p->allocate_ends[i]->reads_number++;
      if (type_sp >= 0) {
	element_p->allocate_ends[i]->type_sp = type_sp;
      }
      if (element_p->allocate_ends[i]->splice_end_extend < splice_end) {
	element_p->allocate_ends[i]->splice_end_extend = splice_end;
      }  
	  
      return element_p;
    }
  }
  
  splice_end_t *splice_end_p = new_splice_end(end, type_orig, type_sp, splice_end);
  
  return insert_end_splice(splice_end_p, element_p);
}


allocate_splice_elements_t** new_allocate_splice_elements(size_t nchromosomes){
  allocate_splice_elements_t** chromosomes_avls = (allocate_splice_elements_t **)malloc(sizeof(allocate_splice_elements_t *)*NUM_STRANDS);
  size_t i, st;
  for(st = 0; st < NUM_STRANDS; st++) {
    chromosomes_avls[st] = (allocate_splice_elements_t *)malloc(sizeof(allocate_splice_elements_t)*nchromosomes);
    for(i = 0; i < nchromosomes; i++) {
      chromosomes_avls[st][i].avl_splice = cp_avltree_create_by_option(COLLECTION_MODE_NOSYNC | 
								       COLLECTION_MODE_COPY   |
								       COLLECTION_MODE_DEEP, 
								       (cp_compare_fn) node_compare, 
								       (cp_copy_fn) node_copy, 
								       (cp_destructor_fn)node_free, 
								       (cp_copy_fn) node_copy, 
								       (cp_destructor_fn)node_free);
      
      if(chromosomes_avls[st][i].avl_splice == NULL) { exit(-1); }
      pthread_mutex_init(&(chromosomes_avls[st][i].mutex), NULL);
    }
  }
  return chromosomes_avls;
}


allocate_splice_elements_t* init_allocate_splice_elements(allocate_splice_elements_t* chromosomes_avls_p, size_t nchromosomes){
   int i;
   for(i = 0; i < nchromosomes; i++){
     chromosomes_avls_p[i].avl_splice = cp_avltree_create_by_option(COLLECTION_MODE_NOSYNC | 
								    COLLECTION_MODE_COPY   |
								    COLLECTION_MODE_DEEP, 
								    (cp_compare_fn) node_compare, 
								    (cp_copy_fn) node_copy, 
								    (cp_destructor_fn)node_free, 
								    (cp_copy_fn) node_copy, 
								    (cp_destructor_fn)node_free);
     
     if(chromosomes_avls_p[i].avl_splice == NULL) {exit(-1);}
     pthread_mutex_init(&(chromosomes_avls_p[i].mutex), NULL);
   }
   
   return chromosomes_avls_p;
}

allocate_splice_elements_t* allocate_new_splice(unsigned int chromosome, unsigned char strand, 
						size_t end, size_t start, 
						size_t splice_start, size_t splice_end, int type_orig,
						int type_sp, allocate_splice_elements_t **chromosome_avls_p) {
  node_element_splice_t *node;

  node = (node_element_splice_t *)cp_avltree_get(chromosome_avls_p[strand][chromosome].avl_splice, (void *)start);

  if(node == NULL) {
    node = cp_avltree_insert(chromosome_avls_p[strand][chromosome].avl_splice, (void *)start, (void *)start);
    node->splice_start_extend = splice_start;
  }

  node = search_and_insert_end_splice(chromosome, end, splice_start, splice_end, type_orig, type_sp, node);

  return chromosome_avls_p;
}

splice_end_t* new_splice_end(size_t end, int type_orig, int type_sp, size_t splice_end) {

  splice_end_t* splice_end_p = (splice_end_t *)malloc(sizeof(splice_end_t));

  if(splice_end_p == NULL) {exit(-1);}
  
  splice_end_p->end = end;
  splice_end_p->splice_end_extend = splice_end;
  splice_end_p->type_sp = type_sp;

  if (type_orig == FROM_READ) {
    splice_end_p->reads_number = 1;
  } else {
    splice_end_p->reads_number = 0;
  }

  return splice_end_p;
}


void free_splice_end(splice_end_t *splice_end_p){
  free(splice_end_p);
}


allocate_buffers_t * process_avlnode_in_order(cp_avlnode *node, unsigned char st, unsigned int chromosome, 
					      list_t* write_list_p, unsigned int write_size,   allocate_buffers_t *allocate_batches){  
  if (node->left) { 
    allocate_batches = process_avlnode_in_order(node->left, st, chromosome, write_list_p, write_size, allocate_batches);
  }
  
  allocate_batches = process_avlnode_ends_in_order((node_element_splice_t *)node->value, st, chromosome, write_list_p, write_size,
						   allocate_batches);
  
  if (node->right) { 
    allocate_batches = process_avlnode_in_order(node->right, st, chromosome,  write_list_p, write_size, allocate_batches);
  }
  
  return allocate_batches;
  // return exact_splice_write_p;

}

allocate_buffers_t* process_avlnode_ends_in_order(node_element_splice_t *node, unsigned char st, unsigned int chromosome,
						  list_t* write_list_p, unsigned int write_size, 
						  allocate_buffers_t *allocate_batches) {
  int i;
  char strand[2] = {'+', '-'};
  list_item_t* item_p = NULL;
  unsigned int bytes_exact, bytes_extend;
  allocate_batches->write_exact_sp;
  //  write_batch_t* extend_splice_write_p = write_batch_new(write_size, SPLICE_EXTEND_FLAG);

  //printf("----------------->%i\n", chromosome);
  for(i = 0; i < node->number_allocate_ends; i++){
    if(( allocate_batches->write_exact_sp->size + 100) > write_size) {
      //item_p = list_item_new(0, WRITE_ITEM,  allocate_batches->write_exact_sp);
      //list_insert_item(item_p, write_list_p);
      //allocate_batches->write_exact_sp = write_batch_new(write_size, SPLICE_EXACT_FLAG);
      fwrite((char *)allocate_batches->write_exact_sp->buffer_p, allocate_batches->write_exact_sp->size, 1, allocate_batches->fd_exact);
      fwrite((char *)allocate_batches->write_extend_sp->buffer_p, allocate_batches->write_extend_sp->size, 1, allocate_batches->fd_extend);
      allocate_batches->write_exact_sp->size = 0;
      allocate_batches->write_extend_sp->size = 0;
    } 
    
    if (node->allocate_ends[i]->reads_number) {
      bytes_exact = pack_junction(chromosome, st, 
				  node->splice_start, node->allocate_ends[i]->end, 
				  junction_id, node->allocate_ends[i]->reads_number, 
				  node->allocate_ends[i]->type_sp,
				  &(((char *)allocate_batches->write_exact_sp->buffer_p)[allocate_batches->write_exact_sp->size]));
      
      bytes_extend = pack_junction(chromosome, st, node->splice_start_extend, 
				   node->allocate_ends[i]->splice_end_extend, junction_id, node->allocate_ends[i]->reads_number, 
				   node->allocate_ends[i]->type_sp,
				   &(((char *)allocate_batches->write_extend_sp->buffer_p)[allocate_batches->write_extend_sp->size])); 
      
      allocate_batches->write_exact_sp->size += bytes_exact;
      allocate_batches->write_extend_sp->size += bytes_extend;
      
      total_splice += node->allocate_ends[i]->reads_number;
      junction_id++;
      if (node->allocate_ends[i]->type_sp >= AT_AC_SPLICE) {
	semi_cannonical_sp++;
      } else {
	cannonical_sp++;
      }
    }
  }
  return allocate_batches;
  //return exact_splice_write_p;
}

void chromosome_avls_free(allocate_splice_elements_t **chromosome_avls, size_t nchromosomes) {
  for(int st = 0; st < NUM_STRANDS; st++) {
    for(int c = 0; c < nchromosomes; c++){
      cp_avltree_destroy(chromosome_avls[st][c].avl_splice);
    }
    free(chromosome_avls[st]);
  }
  free(chromosome_avls);
}

void write_chromosome_avls(allocate_splice_elements_t **chromosome_avls, 
			   list_t* write_list_p, char *extend_sp, char *exact_sp, 
			   unsigned int write_size, size_t nchromosomes) {
  int c, chr;
  unsigned char st;
  allocate_buffers_t *allocate_batches = (allocate_buffers_t *)malloc(sizeof(allocate_buffers_t));
  //write_batch_t *exact_splice_write_p;
  //write_batch_t *extend_splice_write_p;
  //printf("Open splice\n");
  FILE *fd_exact = fopen(exact_sp, "w");  
  if (!fd_exact) { 
    printf("Imposible to create FILE: %s\n", exact_sp);
    exit(-1);
  }

  FILE *fd_extend = fopen(extend_sp, "w");
  if (!fd_exact) { 
    printf("Imposible to create FILE: %s\n", exact_sp);
    exit(-1);
  }

  //  printf("Open splice\n");

  allocate_batches->fd_exact = fd_exact;
  allocate_batches->fd_extend = fd_extend;

  allocate_batches->write_exact_sp  = write_batch_new(write_size, SPLICE_EXACT_FLAG);
  allocate_batches->write_extend_sp  = write_batch_new(write_size, SPLICE_EXTEND_FLAG);
      
  for(st = 0; st < NUM_STRANDS; st++) {
    for(c = 0; c < nchromosomes; c++){
      //printf("Chromosome %i:\n", c);
      if(chromosome_avls[st][c].avl_splice->root != NULL) {
	chr = c + 1;
	//printf("\tYes\n");
	//allocate_batches->write_extend_sp  = write_batch_new(1000, SPLICE_EXTEND_FLAG);
	
	allocate_batches = process_avlnode_in_order(chromosome_avls[st][c].avl_splice->root, st, chr, write_list_p, write_size, allocate_batches);
      
      //exact_splice_write_p = allocate_batches->write_exact_sp;
      //extend_splice_write_p = allocate_batches->write_extend_sp;
      
	if(allocate_batches->write_extend_sp != NULL) {
	  if(allocate_batches->write_extend_sp->size > 0) {
	    //item_p = list_item_new(0, WRITE_ITEM, exact_splice_write_p);
	    //list_insert_item(item_p, write_list_p);
	    //} else {
	    //write_batch_free(exact_splice_write_p);
	    //}
	    fwrite((char *)allocate_batches->write_exact_sp->buffer_p, allocate_batches->write_exact_sp->size, 1, allocate_batches->fd_exact);
	    fwrite((char *)allocate_batches->write_extend_sp->buffer_p, allocate_batches->write_extend_sp->size, 1, allocate_batches->fd_extend);
	    allocate_batches->write_exact_sp->size = 0;
	    allocate_batches->write_extend_sp->size = 0;
	  }
	}
      } //end IF chromosome splice not NULL
      cp_avltree_destroy(chromosome_avls[st][c].avl_splice);
    }
    free(chromosome_avls[st]);
  }
  free(chromosome_avls);
  fclose(fd_extend);
  fclose(fd_exact);
  
  write_batch_free(allocate_batches->write_exact_sp);
  write_batch_free(allocate_batches->write_extend_sp);
  free(allocate_batches);
  basic_statistics_sp_init(total_splice, cannonical_sp, semi_cannonical_sp, basic_st);
  
  
}

*/
