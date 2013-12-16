#ifndef ARRAY_LIST_BS_H
#define ARRAY_LIST_BS_H

#include <stdio.h>
#include <math.h>
#include <limits.h>

#include "commons/string_utils.h"
#include "commons/log.h"

#include "containers/containers.h"
#include "containers/cprops/hashtable.h"


//#define COLLECTION_MODE_SYNCHRONIZED        1
//#define COLLECTION_MODE_ASYNCHRONIZED       2

//=====================================================
// structures
//=====================================================

typedef struct metil_data {
  char*  query_name;
  char   status;
  int    chromosome;
  size_t start;
  char   context;
  int    strand;
  int    zone;
} metil_data_t;


typedef struct array_list_bs {
  size_t capacity;
  size_t size;
  float realloc_factor;
  int mode;
  int flag;

  int (*compare_fn)(const void *, const void *);

  pthread_mutex_t lock;
  pthread_cond_t condition;

  metil_data_t *items;
} array_list_bs_t;

//------------------------------------------------------------------------------------

/**
 * array_list_bs functions
 */

array_list_bs_t* array_list_bs_new(size_t initial_capacity, float realloc_factor, int SYNC_MODE);

void array_list_bs_free(array_list_bs_t* list_p);

//------------------------------------------------------------------------------------

size_t array_list_bs_capacity(array_list_bs_t *array_list_bs_p);

size_t array_list_bs_size(array_list_bs_t *array_list_bs_p);

size_t array_list_bs_index_of(void *item, array_list_bs_t *array_list_bs_p);

int array_list_bs_contains(void* item, array_list_bs_t *array_list_bs_p);

int array_list_bs_clear(array_list_bs_t* array_list_bs_p);

//------------------------------------------------------------------------------------

int array_list_bs_insert(array_list_bs_t *array_list_bs_p, char* query_name, char status,
	                 int chromosome, size_t start, char   context, int strand, int zone);

int array_list_bs_insert_at(size_t index, void* item_p, array_list_bs_t *array_list_bs_p);

int array_list_bs_insert_all(void** item_p, size_t num_items, array_list_bs_t *array_list_bs_p);

int array_list_bs_insert_all_at(size_t index, void** item_p, size_t num_items, array_list_bs_t* array_list_bs_p);

//------------------------------------------------------------------------------------

void* array_list_bs_remove(void *item, array_list_bs_t *array_list_bs_p);

void* array_list_bs_remove_at(size_t index, array_list_bs_t *array_list_bs_p);

void** array_list_bs_remove_range(size_t start, size_t end, array_list_bs_t* array_list_bs_p);

//------------------------------------------------------------------------------------

metil_data_t * array_list_bs_get(size_t index, array_list_bs_t *array_list_bs_p);

array_list_bs_t* array_list_bs_sublist(size_t start, size_t end, array_list_bs_t *array_list_bs_p, array_list_bs_t *sublist);

void* array_list_bs_set(size_t index, void* new_item, array_list_bs_t *array_list_bs_p);

//------------------------------------------------------------------------------------

void array_list_bs_print(array_list_bs_t *array_list_bs_p);

static array_list_bs_t *reallocate_bs(array_list_bs_t * array_list_bs_p, size_t inc_size);

static int compare_items_bs(const void *item1, const void *item2);

int array_list_bs_swap(const int pos1, const int pos2, array_list_bs_t *array_list_bs_p);

//------------------------------------------------------------------------------------

void array_list_bs_set_flag(int flag, array_list_bs_t *array_list_bs_p);
int array_list_bs_get_flag(array_list_bs_t *array_list_bs_p);

//------------------------------------------------------------------------------------

/**
*  @brief Returns an arraylist with the unique elements of the given array list
*  @param orig pointer to the origin array list
*  @param compare callback function for comparison 
*  @param[in,out] dest pointer to the destination array list
*  @return pointer to the destination array list
*  
*  Returns an arraylist with the unique elements of the given array list
*/
array_list_bs_t* array_list_bs_unique(array_list_bs_t *orig, int (*compare)(const void *a, const void *b), array_list_bs_t *dest);

//------------------------------------------------------------------------------------

/**
*  @brief Returns the intersection between two given arraylists
*  @param al1 array list 1
*  @param al2 array list 2
*  @param compare callback function for comparison
*  @param[in,out] dest pointer to the destination array list
*  @return pointer to the destination array list
*  
*  Returns the intersection (the common elements) between two given arraylists
*/
array_list_bs_t* array_list_bs_intersect(array_list_bs_t *al1, array_list_bs_t *al2, int (*compare)(const void *a, const void *b), array_list_bs_t *dest);

//------------------------------------------------------------------------------------

/**
*  @brief Returns the complementary arraylist between two given arraylists
*  @param al1 array list 1
*  @param al2 array list 2
*  @param compare callback function for comparison
*  @param[in,out] dest pointer to the destination array list
*  @return pointer to the destination array list
*  
*  Returns the complementary between two given arraylists. The complementary elements are those of 
*  list 2 that are not present in list 1
*/
array_list_bs_t* array_list_bs_complement(array_list_bs_t *al1, array_list_bs_t *al2,  int (*compare)(const void *a, const void *b), array_list_bs_t *dest);

//------------------------------------------------------------------------------------

/**
*  @brief Compare function for strings
*  @param a pointer to string
*  @param b pointer to string
*  @return 0: if equal, <>0: if not equal
*  
*  Compare function for strings
*/
int compare(const void *a, const void *b);

//------------------------------------------------------------------------------------

#endif /* array_list_bs_H */
