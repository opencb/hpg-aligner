#include "array_list_bs.h"

/**
 * array_list_bs items functions
 */

//------------------------------------------------------------------------------------

array_list_bs_t* array_list_bs_new(size_t initial_capacity, float realloc_factor, int SYNC_MODE) {
  array_list_bs_t *array_list_bs_p = (array_list_bs_t*) malloc(sizeof(array_list_bs_t));
  array_list_bs_p->capacity = initial_capacity;
  array_list_bs_p->size = 0;
  array_list_bs_p->realloc_factor = realloc_factor;
  array_list_bs_p->mode = SYNC_MODE;
  array_list_bs_p->compare_fn = compare_items_bs;

  array_list_bs_p->items = (metil_data_t *) malloc(initial_capacity * sizeof(metil_data_t));
  
  pthread_mutex_init(&(array_list_bs_p->lock), NULL);
  
  return array_list_bs_p;
}

//------------------------------------------------------------------------------------

int array_list_bs_clear(array_list_bs_t *array_list_bs_p) {
  if(array_list_bs_p != NULL) {
    for(size_t i=0; i < array_list_bs_p->size; i++) {
      if (array_list_bs_p->items[i].query_name) free(array_list_bs_p->items[i].query_name);
    }
    // Set default parameters
    array_list_bs_p->size = 0;
    return 1;
  }
  return 0;
}

//------------------------------------------------------------------------------------

void array_list_bs_free(array_list_bs_t *array_list_bs_p) {
  if(array_list_bs_p != NULL) {
    array_list_bs_clear(array_list_bs_p);
    free(array_list_bs_p->items);
    free(array_list_bs_p);
  }
}

//------------------------------------------------------------------------------------

size_t array_list_bs_capacity(array_list_bs_t *array_list_bs_p) {
  if(array_list_bs_p != NULL) {
    return array_list_bs_p->capacity;
  }
  return ULONG_MAX;
}

//------------------------------------------------------------------------------------

size_t array_list_bs_size(array_list_bs_t *array_list_bs_p) {
  size_t length;
  
  if(array_list_bs_p != NULL) {
    if(array_list_bs_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
      pthread_mutex_lock(&array_list_bs_p->lock);
    }
    
    length = array_list_bs_p->size;
    
    if(array_list_bs_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
      pthread_mutex_unlock(&array_list_bs_p->lock);
    }
  }else{
    length = ULONG_MAX;
  }
  
  return length;
}

//------------------------------------------------------------------------------------

int array_list_bs_insert(array_list_bs_t *array_list_bs_p, char* query_name, char status,
	                 int chromosome, size_t start, char context, int strand, int zone) {
  if(array_list_bs_p != NULL) {
    if(array_list_bs_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
      pthread_mutex_lock(&array_list_bs_p->lock);
    }
    
    // Check if array_list_bs is full
    if(array_list_bs_p->size >= array_list_bs_p->capacity) {
      //printf("reallocate\n");
      array_list_bs_p = reallocate_bs(array_list_bs_p, 0);
    }
    
    array_list_bs_p->items[array_list_bs_p->size].query_name = strdup(query_name);
    array_list_bs_p->items[array_list_bs_p->size].status = status;
    array_list_bs_p->items[array_list_bs_p->size].chromosome = chromosome;
    array_list_bs_p->items[array_list_bs_p->size].start = start;
    array_list_bs_p->items[array_list_bs_p->size].context = context;
    array_list_bs_p->items[array_list_bs_p->size].strand = strand;
    array_list_bs_p->items[array_list_bs_p->size].zone = zone;
    
    array_list_bs_p->size++;
    if(array_list_bs_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
      pthread_mutex_unlock(&array_list_bs_p->lock);
    }
    return 1;
  }
  return 0;
}

//------------------------------------------------------------------------------------

metil_data_t * array_list_bs_get(size_t index, array_list_bs_t *array_list_bs_p) {
  metil_data_t *item_p;
  if(array_list_bs_p != NULL && index >= 0 && index <= array_list_bs_p->size) {
    if(array_list_bs_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
      pthread_mutex_lock(&array_list_bs_p->lock);
    }
    
    item_p = &(array_list_bs_p->items[index]);
    
    if(array_list_bs_p->mode == COLLECTION_MODE_SYNCHRONIZED) {
      pthread_mutex_unlock(&array_list_bs_p->lock);
    }
    
  }else{
    item_p = NULL;
  }
  
  return item_p;
}

//------------------------------------------------------------------------------------

array_list_bs_t *reallocate_bs(array_list_bs_t * array_list_bs_p, size_t inc_size) {
  // Capacity is increased in factor.
  size_t new_capacity;// = array_list_bs_t->capacity;
  if(!inc_size) {
    new_capacity = (int)ceil((float)array_list_bs_p->capacity * array_list_bs_p->realloc_factor);
  }else {
    new_capacity = array_list_bs_p->capacity + inc_size;
  }
  // Realloc items with the new capacity. Size remains equals.
  metil_data_t *items_aux = (metil_data_t*) realloc(array_list_bs_p->items, new_capacity * sizeof(metil_data_t));
  if(items_aux != NULL) {
    array_list_bs_p->items = items_aux;
    array_list_bs_p->capacity = new_capacity;
  }else {
    LOG_ERROR("Error in reallocate_bs");
  }
  return array_list_bs_p;
}

//------------------------------------------------------------------------------------

int compare_items_bs(const void *item1, const void *item2) {
  return item1 != item2;
}

//------------------------------------------------------------------------------------

int compare_bs(const void *a, const void *b) {
  return strcmp((char*)a, (char*)b);
}

//------------------------------------------------------------------------------------

