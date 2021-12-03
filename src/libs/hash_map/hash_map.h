#pragma once
#include <stdint.h>

struct HashMap;

uint32_t hash_murmur3_32(const uint8_t* key, size_t key_size, uint32_t seed)
;

// You want num_buckets to be ~ 1.5*the number of entries 
struct HashMap *hash_map_create(int num_buckets)
;
struct HashMap* hash_map_clone(struct HashMap* src)
;
void hash_map_free(struct HashMap* hash_map)
;
void *hash_map_insert(struct HashMap *hash_map, void *key, unsigned int key_size, void *value, unsigned int value_size)
;
void *hash_map_insert_and_grow(struct HashMap *hash_map, void *key, unsigned int key_size, void *value, unsigned int value_size)
;
void *hash_map_find(struct HashMap *hash_map, void *key, unsigned int key_size, unsigned int *value_size_out)
;
int hash_map_remove(struct HashMap *hash_map, void *key, unsigned int key_size)
;
void hash_map_clear(struct HashMap *hash_map)
;
int hash_map_num_keys(struct HashMap *hash_map)
;
int hash_map_num_buckets(struct HashMap* hash_map)
;
void hash_map_resize(struct HashMap* hash_map, int num_buckets)
;

struct HashMapIterator *hash_map_iterator_create(struct HashMap *hash_map)
;
int hash_map_iterator_next(struct HashMapIterator *hash_map_iterator, void **key_out, int *key_size_out, 
	void **value_out, int *value_size_out)
;
void hash_map_iterator_remove_current(struct HashMapIterator* hash_map_iterator)
;
void hash_map_iterator_free(struct HashMapIterator *hash_map_iterator)
;
