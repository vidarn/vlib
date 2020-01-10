#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "hash_map.h"

struct HashMapBucket
{
	unsigned int  *key_sizes;
	unsigned int  *value_sizes;
	unsigned char *keys;
	unsigned char *values;
	unsigned int alloc, num;
	unsigned int key_alloc, key_size;
	unsigned int value_alloc, value_size;
};

struct HashMap
{
	struct HashMapBucket *buckets;
	int num_buckets;
	int num_keys;
};

struct HashMapIterator
{
	struct HashMap *hash_map;
	int next_bucket;
	int next_bucket_entry;
};

//NOTE(Vidar):From here: http://www.cse.yorku.ca/~oz/hash.html
unsigned long hash_djb2(unsigned char *key, unsigned int key_size)
{
	unsigned long hash = 5381;
	for(unsigned int i=0;i<key_size;i++){
		int c = *key++;
		hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
	}
	return hash;
}

//NOTE(Vidar):From here: https://en.wikipedia.org/wiki/MurmurHash
uint32_t hash_murmur3_32(const uint8_t* key, size_t len, uint32_t seed)
{
	uint32_t h = seed;
	if (len > 3) {
		size_t i = len >> 2;
		do {
			uint32_t k;
			memcpy(&k, key, sizeof(uint32_t));
			key += sizeof(uint32_t);
			k *= 0xcc9e2d51;
			k = (k << 15) | (k >> 17);
			k *= 0x1b873593;
			h ^= k;
			h = (h << 13) | (h >> 19);
			h = h * 5 + 0xe6546b64;
		} while (--i);
	}
	if (len & 3) {
		size_t i = len & 3;
		uint32_t k = 0;
		do {
			k <<= 8;
			k |= key[i - 1];
		} while (--i);
		k *= 0xcc9e2d51;
		k = (k << 15) | (k >> 17);
		k *= 0x1b873593;
		h ^= k;
	}
	h ^= len;
	h ^= h >> 16;
	h *= 0x85ebca6b;
	h ^= h >> 13;
	h *= 0xc2b2ae35;
	h ^= h >> 16;
	return h;
}

struct HashMap *hash_map_create(int num_buckets)
{
	struct HashMap *hash_map = calloc(1, sizeof(struct HashMap));
	hash_map->num_buckets = num_buckets;
	hash_map->buckets = calloc(num_buckets, sizeof(struct HashMapBucket));
	return hash_map;
}
void hash_map_free(struct HashMap* hash_map)
{
	int num_buckets = hash_map->num_buckets;
	for (int i = 0; i < num_buckets; i++) {
		struct HashMapBucket *bucket = hash_map->buckets + i;
		if (bucket->alloc > 0) {
			free(bucket->key_sizes);
			free(bucket->value_sizes);
			free(bucket->keys);
			free(bucket->values);
		}
	}
	free(hash_map->buckets);
	free(hash_map);
}

void hash_map_insert(struct HashMap *hash_map, void *key, unsigned int key_size, void *value, unsigned int value_size)
{
	unsigned long hash = hash_djb2(key, key_size);
	unsigned int bucket_i = hash % hash_map->num_buckets;
	struct HashMapBucket *bucket = hash_map->buckets + bucket_i;
	int index = bucket->num++;
	if (bucket->num > bucket->alloc) {
		bucket->alloc = bucket->alloc > 0 ? bucket->alloc << 1 : 4;
		bucket->key_sizes     = realloc(bucket->key_sizes,   bucket->alloc * sizeof(unsigned int));
		bucket->value_sizes   = realloc(bucket->value_sizes, bucket->alloc * sizeof(unsigned int));
	}

	bucket->key_sizes[index] = key_size;
	int key_offset = bucket->key_size;
	bucket->key_size += key_size;
	if (bucket->key_size > bucket->key_alloc) {
		bucket->key_alloc = bucket->key_size << 1;
		bucket->keys = realloc(bucket->keys, bucket->key_alloc);
	}
	memcpy(bucket->keys + key_offset, key, key_size);

	bucket->value_sizes[index] = value_size;
	int value_offset = bucket->value_size;
	bucket->value_size += value_size;
	if (bucket->value_size > bucket->value_alloc) {
		bucket->value_alloc = bucket->value_size << 1;
		bucket->values = realloc(bucket->values, bucket->value_alloc);
	}
	memcpy(bucket->values + value_offset, value, value_size);
	hash_map->num_keys++;
}

void hash_map_insert_and_grow(struct HashMap* hash_map, void* key, unsigned int key_size, void* value, unsigned int value_size)
{
	int num_keys = hash_map->num_keys;
	if (hash_map->num_buckets <= hash_map->num_keys) {
		hash_map_resize(hash_map, hash_map->num_buckets * 2);
	}
	hash_map_insert(hash_map, key, key_size, value, value_size);
}

void *hash_map_find(struct HashMap *hash_map, void *key, unsigned int key_size, unsigned int *value_size_out)
{
	unsigned long hash = hash_djb2(key, key_size);
	unsigned int bucket_i = hash % hash_map->num_buckets;
	struct HashMapBucket *bucket = hash_map->buckets + bucket_i;
	int num_entries = bucket->num;
	int key_offset = 0;
	int value_offset = 0;
	for (int i = 0; i < num_entries; i++) {
		unsigned int ks = bucket->key_sizes[i];
		unsigned int vs =  bucket->value_sizes[i];
		if (ks == key_size) {
			unsigned char *k = bucket->keys + key_offset;
			if (memcmp(key, k, key_size) == 0) {
				void *ret = bucket->values + value_offset;
				if (value_size_out) {
					*value_size_out = vs;
				}
				return ret;
			}
		}
		key_offset   += ks;
		value_offset += vs;
	}
	if (value_size_out) {
		*value_size_out = 0;
	}
	return 0;
}

static int _hash_map_bucket_remove_entry(struct HashMapBucket* bucket, int entry)
{
	int num_entries = bucket->num;
	if (entry != -1) {
		unsigned char *k_tmp = _alloca(bucket->key_size);
		unsigned char *v_tmp = _alloca(bucket->value_size);
		memcpy(k_tmp, bucket->keys, bucket->key_size);
		memcpy(v_tmp, bucket->values, bucket->value_size);
		int key_offset_tmp = 0;
		int value_offset_tmp = 0;
		int key_offset = 0;
		int value_offset = 0;
		for (int i = 0; i < num_entries; i++) {
			unsigned int ks = bucket->key_sizes[i];
			unsigned int vs = bucket->value_sizes[i];
			if (i != entry) {
				memcpy(bucket->keys + key_offset, k_tmp + key_offset_tmp, ks);
				memcpy(bucket->values + value_offset, v_tmp + value_offset_tmp, vs);
				key_offset += ks;
				value_offset += vs;
			}
			key_offset_tmp   += ks;
			value_offset_tmp += vs;
		}
		bucket->key_size -= bucket->key_sizes[entry];
		bucket->value_size -= bucket->value_sizes[entry];
		for (int i = entry+1; i < num_entries; i++) {
			bucket->key_sizes[i-1] = bucket->key_sizes[i];
			bucket->value_sizes[i-1] = bucket->value_sizes[i];
		}
		bucket->num--;
		return 1;
	}
	return 0;
}

int hash_map_remove(struct HashMap *hash_map, void *key, unsigned int key_size)
{
	unsigned long hash = hash_djb2(key, key_size);
	unsigned int bucket_i = hash % hash_map->num_buckets;
	struct HashMapBucket *bucket = hash_map->buckets + bucket_i;
	int num_entries = bucket->num;
	int key_offset = 0;
	int value_offset = 0;
	int entry = -1;
	for (int i = 0; i < num_entries; i++) {
		unsigned int ks = bucket->key_sizes[i];
		unsigned int vs =  bucket->value_sizes[i];
		if (ks == key_size) {
			unsigned char *k = bucket->keys + key_offset;
			if (memcmp(key, k, key_size) == 0) {
				entry = i;
			}
		}
		key_offset   += ks;
		value_offset += vs;
	}
	int ret = _hash_map_bucket_remove_entry(bucket, entry);
	if (ret) {
		hash_map->num_keys--;
	}
	return ret;
}

void hash_map_clear(struct HashMap *hash_map)
{
	for (int i_bucket = 0; i_bucket < hash_map->num_buckets; i_bucket++)
	{
		struct HashMapBucket *bucket = hash_map->buckets + i_bucket;
		bucket->num = 0;
		bucket->key_size = 0;
		bucket->value_size = 0;
	}
	hash_map->num_keys = 0;
}

int hash_map_num_keys(struct HashMap* hash_map)
{
	return hash_map->num_keys;
	/*
	int num = 0;
	for (int i_bucket = 0; i_bucket < hash_map->num_buckets; i_bucket++)
	{
		struct HashMapBucket* bucket = hash_map->buckets + i_bucket;
		num += bucket->num;
	}
	return num;
	*/
}

void hash_map_resize(struct HashMap* hash_map, int num_buckets)
{
	struct HashMapBucket *new_buckets = calloc(num_buckets, sizeof(struct HashMapBucket));
	struct HashMapBucket *old_buckets = hash_map->buckets;
	int old_num_buckets = hash_map->num_buckets;
	hash_map->buckets = new_buckets;
	hash_map->num_buckets = num_buckets;
	hash_map->num_keys = 0;

	for (int i_bucket = 0; i_bucket < old_num_buckets; i_bucket++) {
		struct HashMapBucket* bucket = old_buckets + i_bucket;
		int num_entries = bucket->num;
		unsigned char* key = bucket->keys;
		unsigned char* value = bucket->values;
		for (int i = 0; i < num_entries; i++) {
			uint32_t key_size = bucket->key_sizes[i];
			uint32_t value_size = bucket->value_sizes[i];
			hash_map_insert(hash_map, key, key_size, value, value_size);
			key += key_size;
			value += value_size;
		}
	}

	free(old_buckets);
}

struct HashMapIterator *hash_map_iterator_create(struct HashMap *hash_map)
{
	struct HashMapIterator *hash_map_iterator = calloc(1, sizeof(struct HashMapIterator));
	hash_map_iterator->hash_map = hash_map;
	hash_map_iterator->next_bucket = 0;
	hash_map_iterator->next_bucket_entry = 0;
	return hash_map_iterator;
}

int hash_map_iterator_next(struct HashMapIterator *hash_map_iterator, void **key_out, int *key_size_out, 
	void **value_out, int *value_size_out)
{
	struct HashMap *hash_map = hash_map_iterator->hash_map;
	int i_bucket = hash_map_iterator->next_bucket;
	int i_bucket_entry =  hash_map_iterator->next_bucket_entry;
	if (i_bucket >= hash_map->num_buckets) {
		return 0;
	}
	struct HashMapBucket *bucket = hash_map->buckets + i_bucket;
	while(i_bucket_entry >= bucket->num) {
		i_bucket_entry = 0;
		i_bucket++;
		if (i_bucket >= hash_map->num_buckets) {
			return 0;
		}
		bucket = hash_map->buckets + i_bucket;
	}

	int key_offset = 0;
	int value_offset = 0;
	for (int i = 0; i < i_bucket_entry; i++) {
		unsigned int ks = bucket->key_sizes[i];
		unsigned int vs =  bucket->value_sizes[i];
		key_offset   += ks;
		value_offset += vs;
	}

	if(key_size_out)   *key_size_out   = bucket->key_sizes[i_bucket_entry];
	if(value_size_out) *value_size_out = bucket->value_sizes[i_bucket_entry];
	if(value_out) *value_out = bucket->values + value_offset;
	if(key_out) *key_out   = bucket->keys   + key_offset;

	hash_map_iterator->next_bucket = i_bucket;
	hash_map_iterator->next_bucket_entry = i_bucket_entry+1;

	return 1;
}

void hash_map_iterator_remove_current(struct HashMapIterator* hash_map_iterator)
{
	int remove_bucket;
	int remove_entry;
	if(hash_map_iterator->next_bucket_entry==0){
		remove_bucket = hash_map_iterator->next_bucket - 1;
		remove_entry = hash_map_iterator->hash_map->buckets[remove_bucket].num-1;
	}
	else {
		hash_map_iterator->next_bucket_entry--;
		remove_bucket = hash_map_iterator->next_bucket;
		remove_entry = hash_map_iterator->next_bucket_entry;
	}
	if (_hash_map_bucket_remove_entry(hash_map_iterator->hash_map->buckets + remove_bucket, remove_entry)) {
		hash_map_iterator->hash_map->num_keys--;
	}
}

void hash_map_iterator_free(struct HashMapIterator *hash_map_iterator)
{
	free(hash_map_iterator);
}