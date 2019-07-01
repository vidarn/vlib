#pragma once

struct HashMap;

// You want num_buckets to be ~ 1.5*the number of entries 
struct HashMap *hash_map_create(int num_buckets)
;
void hash_map_free(struct HashMap* hash_map)
;
void hash_map_insert(struct HashMap *hash_map, void *key, unsigned int key_size, void *value, unsigned int value_size)
;
void *hash_map_find(struct HashMap *hash_map, void *key, unsigned int key_size, unsigned int *value_size_out)
;

struct HashMapIterator *hash_map_iterator_create(struct HashMap *hash_map)
;
int hash_map_iterator_next(struct HashMapIterator *hash_map_iterator, void **key_out, int *key_size_out, 
	void **value_out, int *value_size_out)
;
void hash_map_iterator_free(struct HashMapIterator *hash_map_iterator)
;
