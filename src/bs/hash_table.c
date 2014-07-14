#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "hash_table.h"

hash_table_t *hash_create(unsigned int (*hash_f)(void *k),
                            int (*compare)(const void* a, const void* b),
                            void (*destroy_key)(void* a),
                            void (*destroy_data)(void* a),
                            int initial_size)
{
    hash_table_t *hash_table;

    if( (hash_table = (hash_table_t *)malloc(sizeof(hash_table_t))) == NULL)
        return NULL;

    memset(hash_table, 0, sizeof(hash_table));

    hash_table->max_colitions = 0;
    hash_table->size = initial_size;
    hash_table->used = 0;
    hash_table->hash_f = hash_f;
    hash_table->compare = compare;
    hash_table->destroy_key = destroy_key;
    hash_table->destroy_data = destroy_data;

    hash_table->table = (hash_data_t *)malloc(sizeof(hash_data_t)*initial_size);

    memset(hash_table->table, 0, sizeof(hash_data_t)*initial_size);

    return hash_table;
}

hash_node_t *hash_insert(hash_table_t *table, void *key, void *data)
{
    unsigned int num_key;
    hash_node_t *node;

    num_key = table->hash_f(key) % table->size;
    
    node = (hash_node_t *)malloc(sizeof(hash_node_t));


    node->key = key;
    node->next = table->table[num_key].list;
    node->data = data;

    table->table[num_key].list = node;
    table->table[num_key].count++;
    
    if(table->max_colitions < table->table[num_key].count)
        table->max_colitions = table->table[num_key].count;
    
    return node;
}

void *hash_find_data(hash_table_t *table, void *key)
{
    hash_node_t *node;

    node = hash_find_node(table, key);

    if(node == NULL)
        return NULL;

    return node->data;
}

hash_node_t *hash_find_node(hash_table_t *table, void *key)
{
    unsigned int num_key;
    hash_node_t *node;

    num_key = table->hash_f(key) % table->size;

    for(node = table->table[num_key].list;
            node != NULL && table->compare(node->key, key) != 0;
            node = node->next);

    return node;
}

void hash_rebuild(hash_table_t *table, int new_size)
{
    hash_data_t *new_table;
    hash_node_t *node;
    hash_node_t *node_next;
    hash_data_t *old_table;
    int old_size;
    int i;

    new_table = (hash_data_t *)malloc(sizeof(hash_data_t) * new_size);

    memset(new_table, 0, sizeof(hash_data_t) * new_size);

    old_table = table->table;
    old_size = table->size;
    table->table = new_table;
    table->size = new_size;
    table->max_colitions = 0;

    for(i=0; i<old_size; i++)
    {
        node = old_table[i].list;
        while(node != NULL)
        {
            hash_insert(table, node->key, node->data);
            node_next = node->next;
            free(node);
            node = node_next;
        }
    }

    free(old_table);
}

void hash_delete(hash_table_t *table, void *key)
{
    unsigned int num_key;
    hash_node_t *node;
    hash_node_t *prev_node;

    num_key = table->hash_f(key) % table->size;

    for(prev_node = node = table->table[num_key].list;
            node != NULL && table->compare(node->key, key) != 0;
            prev_node = node, node = node->next);
    
    if(node == NULL)
    {
        return;
    }
    else if(prev_node == node)
    {
        table->table[num_key].list = node->next;
    }
    else
    {
        prev_node->next = node->next;
    }

    table->destroy_key(node->key);
    table->destroy_data(node->data);
    free(node);

}

void hash_destroy(hash_table_t *table)
{
    hash_node_t *node;
    hash_node_t *next_node;
    int i;

    for(i=0; i<table->size; i++)
    {
        node = table->table[i].list;
        
        while(node != NULL)
        {
            next_node = node->next;
            table->destroy_key(node->key);
            table->destroy_data(node->data);
            free(node);
            node = next_node;
        }
    }

    free(table->table);
    free(table);
}

unsigned int jenkins_one_at_a_time_hash(void *vkey)
{
    unsigned hash = 0;
    size_t i;
    unsigned char *key;
    int key_len;

    key = (unsigned char *)vkey;
    
    key_len = strlen(key);

    for (i = 0; i < key_len; i++) {
        hash += key[i];
        hash += (hash << 10);
        hash ^= (hash >> 6);
    }
    hash += (hash << 3);
    hash ^= (hash >> 11);
    hash += (hash << 15);
    return hash;
}


int scmp(const void *a1, const void *a2)
{
    char *s1, *s2;
    s1 = (char *)a1;
    s2 = (char *)a2;

    if(strcmp(s1, s2) < 0) return -1;
    if(strcmp(s1, s2) == 0) return 0;

    return 1;
}

void nulldes(void *a)
{
    ;
}

void destr_key(void *s)
{
    free((char *)s);
}

void hash_init(hash_table_t *table, char *cad) {
  int i, j, k;
  int len = strlen(cad);
  char *nkey;
  char *data;
  nkey = (char *)malloc(3 * sizeof(char));
  data = (char *)malloc(sizeof(char));

  for (i = 0; i < len; i++) {
    for (j = 0; j < len; j++) {
      for (k = 0; k < len; k++) {
	nkey[0] = cad[i];
	nkey[1] = cad[j];
	nkey[2] = cad[k];
	if (cad[i] == 'C') {
	  if (cad[j] == 'G') {
	    data[0] = 'Z';
	    //data = ZONE_CpG;
	  }
	  else {
	    if (cad[k] == 'G') {
	      data[0] = 'X';
	      //data = ZONE_CHG;
	    }
	    else {
	      data[0] = 'H';
	      //data = ZONE_CHH;
	    }
	  }
	}
	else {
	  data[0] = '-';
	  //data = ZONE_OTHER;
	}
	hash_insert(table, nkey, &data);
      }
    }
  }
  free(nkey);
  free(data);
}
