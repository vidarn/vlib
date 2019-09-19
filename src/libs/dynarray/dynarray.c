#include "dynarray.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

static const int num_int_offsets = 3;
static const int offset = 3*sizeof(int);

static void * create(int num, int element_size)
{
    int *a;
	a = calloc(1,num*element_size+offset);
    a[0] = element_size;
    a[1] = num;
    a[2] = 0;
    return (void*)(a+num_int_offsets);
}

void * dynarray_create(int num, int element_size)
{
    return create(num, element_size,0);
}

void dynarray_destroy(void* a)
{
	free((int*)a - num_int_offsets);
}

void *dynarray_append(void* a, void* element)
{
    int* d= (int*)a - num_int_offsets;
    if(d[2] == d[1]){
        int *tmp;
		tmp = calloc(1,d[0]*d[1]*2 + offset);
        memcpy(tmp,d,d[0]*d[1]+offset);
		free(d);
        d = tmp;
        d[1] *= 2;
        printf("New dynarray size: %d\n",d[1]);
    }
    char* c = (char*)d + offset;
    char* e = (char*)element;
    for(int i=0;i<d[0];i++){
        c[d[0]*d[2]+i] = e[i];
    }
    d[2]++;
    return (void*)(d+num_int_offsets);
}

int dynarray_length(void* a)
{
    return *((int*)a-num_int_offsets+2);
}

void dynarray_clear(void* a)
{
    int* d= (int*)a - num_int_offsets;
    d[2] = 0;
}
