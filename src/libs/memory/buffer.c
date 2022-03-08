#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

struct Buffer
{
	unsigned char *mem;
	size_t len, alloc;
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

void *buffer_hand_over_memory(struct Buffer *buffer)
{
	void *mem = buffer->mem;
	free(buffer);
	return mem;
}

void buffer_reset(struct Buffer *buffer)
{
	buffer->len = 0;
}

void buffer_add(const void *ptr, size_t len, struct Buffer *buffer)
{
	int resize = 0;
	while (buffer->alloc < len + buffer->len) {
		if (buffer->alloc == 0) buffer->alloc = 32;
		buffer->alloc <<= 1;
		resize = 1;
	}
	if (resize) {
		buffer->mem = realloc(buffer->mem, buffer->alloc);
	}
	memcpy(buffer->mem + buffer->len, ptr, len);
	buffer->len += len;
}

void buffer_add_str(const char *str, struct Buffer *buffer)
{
    while(buffer->len > 0 && buffer->mem[buffer->len-1] == 0){
        buffer->len--;
    }
    buffer_add(str, strlen(str)+1, buffer);
}

void buffer_add_str_va(struct Buffer *buffer, ...)
{
    va_list va;
    va_start(va, buffer);
    int i=0;
    char *c;
    do{
        c = va_arg(va, char*);
        if(c) buffer_add_str(c, buffer);
    }while(c);
    va_end(va);
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

size_t buffer_num(struct Buffer *buffer, size_t elem_size)
{
	return buffer->len/elem_size;
}

struct BufferSized
{
	unsigned char *mem;
	size_t len, alloc;
	size_t elem_size;
};

struct BufferSized *buffer_s_create(size_t initial_count, size_t elem_size)
{
	struct BufferSized *buffer = calloc(1, sizeof(struct BufferSized));
	buffer->alloc = initial_count * elem_size;
	buffer->mem = calloc(1, buffer->alloc);
	buffer->elem_size = elem_size;
	return buffer;
}

void* buffer_s_get(struct BufferSized* buffer)
{
	return buffer_get(buffer->elem_size, (struct Buffer*)buffer);
}

void buffer_s_add(void* ptr, struct BufferSized* buffer)
{
	buffer_add(ptr, buffer->elem_size, (struct Buffer*)buffer);
}

void* buffer_s_elem(int i, struct BufferSized* buffer)
{
	return buffer->mem + i * buffer->elem_size;
}

size_t buffer_s_num(struct BufferSized *buffer)
{
	return buffer->len/buffer->elem_size;
}
