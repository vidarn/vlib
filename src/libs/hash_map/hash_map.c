#include <stdlib.h>
#include <string.h>

struct HashMapBucket
{
	unsigned int alloc, num;
	unsigned int key_alloc, key_size;
	unsigned int value_alloc, value_size;
	unsigned int  *key_sizes;
	unsigned int  *value_sizes;
	unsigned char *keys;
	unsigned char *values;
};

struct HashMap
{
	int num_buckets;
	struct HashMapBucket *buckets;
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