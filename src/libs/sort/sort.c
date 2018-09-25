#include <string.h>
static void merge(float *a, int i1, int i2, int i3, float *b)
{
    int i = i1;
    int j = i2;
    for(int k=i1; k < i3; k++){
        if(i < i2 && (j >= i3 || a[i] <= a[j])){
            b[k] = a[i];
            i++;
        }else{
            b[k] = a[j];
            j++;
        }
    }
}

//Bottom-up merge sort
//Sorts a using b as temporary storage
void merge_sort(float *a, float *b, int n)
{
    float *original_a = a;
    for(int width = 1; width < n; width *= 2){
        for(int i=0; i<n; i += width*2){
            int i1 = i;
            int i2 = i+width < n ? i+width : n;
            int i3 = i+2*width < n ? i+2*width : n;
            merge(a, i1, i2, i3, b);
        }
        float *tmp =a;
        a = b;
        b = tmp;
    }
    if(a != original_a){
        memcpy(original_a,a,n*sizeof(float));
    }
}

static void merge_auxillary(float *a, int i1, int i2, int i3, float *b,
    int *aux_a, int *aux_b)
{
    int i = i1;
    int j = i2;
    for(int k=i1; k < i3; k++){
        if(i < i2 && (j >= i3 || a[i] <= a[j])){
            b[k] = a[i];
            aux_b[k] = aux_a[i];
            i++;
        }else{
            b[k] = a[j];
            aux_b[k] = aux_a[j];
            j++;
        }
    }
}

void merge_sort_auxillary(float *a, float *b, int *aux_a, int *aux_b, int n)
{
    float *original_a = a;
    int *original_aux_a = aux_a;
    for(int width = 1; width < n; width *= 2){
        for(int i=0; i<n; i += width*2){
            int i1 = i;
            int i2 = i+width < n ? i+width : n;
            int i3 = i+2*width < n ? i+2*width : n;
            merge_auxillary(a, i1, i2, i3, b, aux_a, aux_b);
        }
        float *tmp =a;
        a = b;
        b = tmp;
        int *tmp_aux = aux_a;
        aux_a = aux_b;
        aux_b = tmp_aux;
    }
    if(a != original_a){
        memcpy(original_a,a,n*sizeof(float));
        memcpy(original_aux_a,aux_a,n*sizeof(int));
    }
}

static void merge_custom_auxillary(char *a, int i1, int i2, int i3, char *b,
    int *aux_a, int *aux_b, int stride, int compar(void *a, void *b, void *data), void *data)
{
    int i = i1;
    int j = i2;
    for(int k=i1; k < i3; k++){
        if(i < i2 && (j >= i3 || compar(a+i*stride,a+j*stride, data))){
            memcpy(b+k*stride,a+i*stride,stride);
            aux_b[k] = aux_a[i];
            i++;
        }else{
            memcpy(b+k*stride,a+j*stride,stride);
            aux_b[k] = aux_a[j];
            j++;
        }
    }
}

void merge_sort_custom_auxillary(void *a, void *b, int stride, int *aux_a,
    int *aux_b, int n, int compar(void *a, void *b, void *data), void *data)
{
    void *original_a = a;
    int *original_aux_a = aux_a;
    for(int width = 1; width < n; width *= 2){
        for(int i=0; i<n; i += width*2){
            int i1 = i;
            int i2 = i+width < n ? i+width : n;
            int i3 = i+2*width < n ? i+2*width : n;
            merge_custom_auxillary(a, i1, i2, i3, b, aux_a, aux_b, stride, compar, data);
        }
        void *tmp =a;
        a = b;
        b = tmp;
        int *tmp_aux = aux_a;
        aux_a = aux_b;
        aux_b = tmp_aux;
    }
    if(a != original_a){
        memcpy(original_a,a,n*stride);
        memcpy(original_aux_a,aux_a,n*sizeof(int));
    }

}

