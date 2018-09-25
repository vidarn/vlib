#pragma once

void merge_sort(float *a, float *b, int n);
void merge_sort_auxillary(float *a, float *b, int *aux_a, int *aux_b, int n);
void merge_sort_custom_auxillary(void *a, void *b, int stride, int *aux_a, int *aux_b, int n, int compar(void *a, void *b, void *data), void *data);
