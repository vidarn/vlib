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
