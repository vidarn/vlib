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
void buffer_add(void *ptr, size_t len, struct Buffer *buffer)
;
void buffer_shrink(size_t len, struct Buffer *buffer)
;
void *buffer_get(size_t len, struct Buffer *buffer)
;
void *buffer_ptr(struct Buffer *buffer)
;
size_t buffer_len(struct Buffer *buffer)
;
