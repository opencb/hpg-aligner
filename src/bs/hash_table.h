#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#define ZONE_CpG   1
#define ZONE_CHG   2
#define ZONE_CHH   3
#define ZONE_OTHER 0

typedef struct hash_node_s
{
    void* key;
    void* data;
    struct hash_node_s *next;
} hash_node_t;

typedef struct 
{
    hash_node_t *list;
    short count;
} hash_data_t;

typedef struct
{
    hash_data_t *table;
    int max_colitions;
    int size;
    int used;
    unsigned int (*hash_f)(void *k);
    int (*compare)(const void* a, const void* b); 
    void (*destroy_key)(void* a);
    void (*destroy_data)(void* a);
} hash_table_t;

hash_table_t *hash_create(unsigned int (hash_f)(void *k),
                            int (*Compare)(const void* a, const void* b),
                            void (*DestroyKey)(void* a),
                            void (*DestroyInfo)(void* a),
                            int initial_size);

hash_node_t *hash_insert(hash_table_t *table, void *key, void *data);
void hash_delete(hash_table_t *table, void *key);
hash_node_t *hash_find_node(hash_table_t *table, void *key);
void *hash_find_data(hash_table_t *table, void *key);
void hash_destroy(hash_table_t *table);
void hash_rebuild(hash_table_t *table, int new_size);

unsigned int jenkins_one_at_a_time_hash(void *vkey);
int scmp(const void *a1, const void *a2);
void nulldes(void *a);
void destr_key(void *s);

void hash_init(hash_table_t *table, char *cad);
