#pragma once

struct Buffer
;

struct Buffer *buffer_create(size_t initial_size)
;
void buffer_free(struct Buffer *buffer)
;
void *buffer_hand_over_memory(struct Buffer *buffer)
;
void buffer_reset(struct Buffer *buffer)
;
void buffer_add(const void *ptr, size_t len, struct Buffer *buffer)
;
void buffer_add_str(const char *str, struct Buffer *buffer)
;
void buffer_add_str_va(struct Buffer *buffer, ...)
;
void buffer_shrink(size_t len, struct Buffer *buffer)
;
void *buffer_get(size_t len, struct Buffer *buffer)
;
void *buffer_ptr(struct Buffer *buffer)
;
size_t buffer_len(struct Buffer *buffer)
;
//TODO:Swap the order of these arguments!
size_t buffer_num(struct Buffer* buffer, size_t elem_size)
;

// If you know that all elements will have the same size, these functions might be easier to use
// You can send a BufferSized * to any of the functions above. However, be careful that you don't mess up the 
// assumption that the lenth is always a multiple of elem_size
struct BufferSized
;
struct BufferSized* buffer_s_create(size_t initial_count, size_t elem_size)
;
void* buffer_s_get(struct BufferSized* buffer)
;
void buffer_s_add(void* ptr, struct BufferSized* buffer)
;
void* buffer_s_elem(int i, struct BufferSized* buffer)
;
size_t buffer_s_num(struct BufferSized *buffer)
;
