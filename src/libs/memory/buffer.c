#include <stdlib.h>
#include <string.h>

struct Buffer
{
	unsigned char *mem;
	int len, alloc;
};

struct Buffer *buffer_create(size_t initial_size)
{
	struct Buffer *buffer = calloc(1, sizeof(struct Buffer));
	buffer->alloc = initial_size;
	buffer->mem = calloc(1, initial_size);
	return buffer;
}
void buffer_free(struct Buffer *buffer)
{
	free(buffer->mem);
	free(buffer);
}

void buffer_reset(struct Buffer *buffer)
{
	buffer->len = 0;
}

void buffer_add(void *ptr, size_t len, struct Buffer *buffer)
{
	int resize = 0;
	while (buffer->alloc < len + buffer->len) {
		buffer->alloc <<= 1;
		resize = 1;
	}
	if (resize) {
		buffer->mem = realloc(buffer->mem, buffer->alloc);
	}
	memcpy(buffer->mem + buffer->len, ptr, len);
	buffer->len += len;
}

void buffer_shrink(size_t len, struct Buffer *buffer)
{
	buffer->len -= len;
}

void *buffer_get(size_t len, struct Buffer *buffer)
{
	int resize = 0;
	while (buffer->alloc < len + buffer->len) {
		buffer->alloc <<= 1;
		resize = 1;
	}
	if (resize) {
		buffer->mem = realloc(buffer->mem, buffer->alloc);
	}
	void *ret = buffer->mem + buffer->len;
	buffer->len += len;
	return ret;
}

void *buffer_ptr(struct Buffer *buffer)
{
	return buffer->mem;
}

size_t buffer_len(struct Buffer *buffer)
{
	return buffer->len;
}

