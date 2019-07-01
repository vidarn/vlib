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

void hash_map_iterator_free(struct HashMapIterator *hash_map_iterator)
{
	free(hash_map_iterator);
}