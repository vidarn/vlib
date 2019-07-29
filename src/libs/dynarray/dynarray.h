#pragma once

void *dynarray_create(int num, int size);
void dynarray_destroy(void* a);
void *dynarray_append(void* a,void* element);
int dynarray_length(void* a);
void dynarray_clear(void* a);
