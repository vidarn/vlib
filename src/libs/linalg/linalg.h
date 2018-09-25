#pragma once

// A 3x2 matrix look like this:
// [ 0  3 ]
// [ 1  4 ]
// [ 2  5 ]


struct Matrix{
    int n[2];
    float *m;
};

#define MAT(matrix,i,j) (matrix.m[i+j*matrix.n[0]])

struct Matrix matrix_new(int m, int n);
void matrix_free(struct Matrix A);
struct Matrix matrix_init(int n, int m, ...);
struct Matrix matrix_identity(int n);
void matrix_copy(struct Matrix src, struct Matrix dest);
void qr_decomp(struct Matrix A, struct Matrix q, struct Matrix r);
void matrix_multiply(struct Matrix A, struct Matrix B, struct Matrix result);
void matrix_transpose(struct Matrix A, struct Matrix result);
void back_subst(struct Matrix A, struct Matrix x, struct Matrix y);
void matrix_print(struct Matrix A);
struct Matrix linear_regression(struct Matrix X, struct Matrix y);
void matrix_eigen_2x2(struct Matrix A, struct Matrix v1, struct Matrix v2,
    float *lambda1, float *lambda2);
