#pragma once

struct MMapContext
;

struct MMapContext *mmap_file_read(const char *filename)
;

void mmap_close(struct MMapContext *context)
;

size_t mmap_size(struct MMapContext *context)
;

void *mmap_get_ptr(struct MMapContext *context)
;
