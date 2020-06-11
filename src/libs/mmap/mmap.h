#pragma once

struct MMapContext
;
struct MMapContext *mmap_file_read(const char *filename)
;
struct MMapContext *mmap_file_write(const char *filename, size_t file_size)
;
void mmap_close(struct MMapContext *context)
;
size_t mmap_size(struct MMapContext *context)
;
void *mmap_get_ptr_read(struct MMapContext *context)
;
void *mmap_get_ptr_write(struct MMapContext *context)
;
