#include "linalg.h"
#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>

struct Matrix matrix_new(int m, int n)
{
    struct Matrix ret = {m,n,0};
    ret.m = calloc(m*n, sizeof(float));
    return ret;
}

void matrix_free(struct Matrix A)
{
    free(A.m);
}

struct Matrix matrix_init(int n, int m, ...)
{
    va_list args;
    va_start(args, m);
    struct Matrix ret = matrix_new(n,m);
    for(int j=0;j<ret.n[1];j++){
        for(int i=0;i<ret.n[0];i++){
            MAT(ret,i,j) =va_arg(args,double);
        }
    }
    va_end(args);
    return ret;
}

struct Matrix matrix_identity(int n)
{
    struct Matrix ret = matrix_new(n, n);
    for(int i=0;i<n;i++){
        MAT(ret,i,i) = 1.f;
    }
    return ret;
}

void matrix_copy(struct Matrix src, struct Matrix dest)
{
    int n = src.n[0]; int m = src.n[1];
    assert(dest.n[0] == n); assert(dest.n[1] == m);
    memcpy(dest.m, src.m, n*m*sizeof(float));
}

//TODO(Vidar): Solve for non-square matrices
void qr_decomp(struct Matrix A, struct Matrix q, struct Matrix r)
{

    //Using Householder reflections, according to
    //https://en.wikipedia.org/wiki/QR_decomposition
    int m = A.n[0];
    int n = A.n[1];
    assert(m >= n);
    assert(q.n[0] == m); assert(q.n[1] == m);
    assert(r.n[0] == m); assert(r.n[1] == n);

    struct Matrix r_tmp = matrix_new(m, n);
    matrix_copy(A, r_tmp);
    struct Matrix A_tmp = matrix_new(m, n);
    matrix_copy(A, A_tmp);
    for(int i=0;i<m;i++){
        for(int j=0;j<m;j++){
            MAT(q, i, j) = i==j ? 1.f: 0.f;
        }
    }
    
    int num_steps = m-1 < n ? m-1 : n;
    for(int step =0; step < num_steps; step++){
        struct Matrix x = matrix_new(m, 1);
        float alpha = 0.f;
        for(int i=step;i<m;i++){
            float a = MAT(A_tmp,i,step);
            MAT(x, i, 0) = a;
            alpha += a*a;
        }
        alpha = sqrtf(alpha);
    
        struct Matrix u = matrix_new(m, 1);
        for(int i=0;i<m;i++){
            MAT(u,i,0) = MAT(x,i,0);
        }
        MAT(u,step,0) -= alpha;
    
        float u_magn = 0.f;
        for(int i=0;i<m;i++){
            float a = MAT(u,i,0);
            u_magn += a*a;
        }
        u_magn = sqrtf(u_magn);
    
        struct Matrix v = matrix_new(m, 1);
        for(int i=0;i<m;i++){
            MAT(v,i,0) = MAT(u,i,0)/u_magn;
        }
    
        struct Matrix v_transpose = matrix_new(1, m);
        matrix_transpose(v, v_transpose);
    
        struct Matrix vv = matrix_new(m, m);
        matrix_multiply(v, v_transpose, vv);
    
        struct Matrix Q_tmp = matrix_new(m, m);
        for(int i=0;i<m;i++){
            for(int j=0;j<m;j++){
                float ident = i==j ? 1.f : 0.f;
                MAT(Q_tmp,i,j) = ident - 2.f * MAT(vv,i,j);
            }
        }
    
        matrix_multiply(Q_tmp, r_tmp, r);
        matrix_copy(r, A_tmp);
        matrix_copy(r, r_tmp);
        
        struct Matrix tmp = matrix_new(m, m);
        matrix_copy(q, tmp);
        matrix_transpose(Q_tmp, Q_tmp);
        matrix_multiply(tmp, Q_tmp, q);

        matrix_free(tmp);
        matrix_free(Q_tmp);
        matrix_free(vv);
        matrix_free(v_transpose);
        matrix_free(v);
        matrix_free(u);
        matrix_free(x);
    
        /*
    #define PRINT(M) printf("  ---" #M "---\n"); matrix_print(M);
        PRINT(A)
        PRINT(x)
        PRINT(u)
        PRINT(v)
        PRINT(v_transpose)
        PRINT(vv)
        PRINT(Q_tmp)
        PRINT(r)
        PRINT(q)
    #undef PRINT
         */
    }
    matrix_free(r_tmp);
    matrix_free(A_tmp);

    /*
    for(int j=0;j<n;j++){
        float s = 0.f;
        for(int i=0;i<n;i++){
            float val = A.m[i + j*n];
            s += val*val;
        }
        float a = sqrtf(s);
        r.m[j + j*n] = a;
        for(int i=0;i<n;i++){
            q.m[i + j*n] = q.m[i + j*n] / a;
        }
        for(int k=j+1;k<n;k++){
            float s = 0.f;
            for(int i=0;i<n;i++){
                float val1 = q.m[i + j*n];
                float val2 = q.m[i + k*n];
                s += val1*val2;
            }
            r.m[j+k*n] = s;
            for(int i=0;i<n;i++){
                q.m[i+k*n] = q.m[i+k*n] - q.m[i+j*n] * s;
            }
        }
    }
     */
    //BOOKMARK(Vidar): Q seems to be incorrect...
    /*
     for j = 1 to p
     {
         define r[j,j] = sqrt( sum_i x[i,j]^2 )
    
         # r[j,j] is the norm of the jth column of X
    
         for i = 1 to n
         {
             x[i,j] = x[i,j] / r[j,j]
         }
    
         for k = j+1 to p
         {
             r[j,k] = sum_{i=1}^n x[i,j]x[i,k]
             for i = 1 to n
             {
                 x[i,k] = x[i,k] - x[i,j] r[j,k]
             }
         }
     }
     */
}

void matrix_multiply(struct Matrix A, struct Matrix B, struct Matrix result)
{
    // A - nxm
    // B - mxp
    // => result - nxp
    int n = A.n[0];
    int m = B.n[0];
    int p = B.n[1];
    assert(A.n[1] == m);
    assert(result.n[0] == n);
    assert(result.n[1] == p);
    for(int i = 0;i<n;i++){
        for(int j = 0;j<p;j++){
            result.m[i+j*n] = 0.f;
            for(int k = 0;k<m;k++){
                result.m[i+j*n] += A.m[i+k*n]*B.m[k+j*m];
            }
        }
    }
}

void matrix_transpose(struct Matrix A, struct Matrix result){
    int n = A.n[0];
    int m = A.n[1];
    assert(result.n[0] == m);
    assert(result.n[1] == n);
    for(int i=0;i<n;i++){
        for(int j=0;j<m;j++){
            MAT(result,j,i) = MAT(A,i,j);
        }
    }
}

//NOTE(Vidar):Finds x from Lx = y
void back_subst(struct Matrix A, struct Matrix x, struct Matrix y)
{
    // A - mxn
    // x - nx1
    // y - mx1
    // m >= n
    int m = A.n[0];
    int n = A.n[1];
    assert(x.n[0] == n); assert(x.n[1] == 1);
    assert(y.n[0] == m); assert(y.n[1] == 1);
    assert(m >= n);
    for(int i=n-1;i>=0;i--){
        float val = MAT(y,i,0);
        for(int j=n-1;j>i;j--){
            val -= MAT(x,j,0) * MAT(A,i,j);
            //val -= x.m[j]*A.m[j+i*n];
        }
        MAT(x,i,0) = val/MAT(A,i,i);
        //x.m[i] = val/A.m[i+n*i];
    }
}

void matrix_print(struct Matrix A)
{
    for(int i=0;i<A.n[0];i++){
        printf("[  ");
        for(int j=0;j<A.n[1];j++){
            printf("%4.4f  ", MAT(A,i,j));
        }
        printf("]\n");
    }
    printf("\n");
}


struct Matrix linear_regression(struct Matrix X, struct Matrix y)
{
    int m = X.n[0];
    int n = X.n[1];
    assert(y.n[0] == m); assert(y.n[1] == 1);
    struct Matrix beta = matrix_new(n, 1);
    struct Matrix y_tmp = matrix_new(m, 1);
    struct Matrix Q = matrix_new(m, m);
    struct Matrix Q_transpose = matrix_new(m, m);
    struct Matrix R = matrix_new(m, n);
    qr_decomp(X, Q, R);
    matrix_transpose(Q, Q_transpose);
    matrix_multiply(Q_transpose,y,y_tmp);
    back_subst(R, beta, y_tmp);
    /*
#define PRINT(M)\
printf("--- Matrix \"" #M "\" ---\n");\
matrix_print(M);
    PRINT(X)
    PRINT(Q)
    PRINT(R)
    PRINT(y)
    PRINT(y_tmp)
    PRINT(beta)
#undef PRINT
    */
    matrix_free(y_tmp);
    matrix_free(Q);
    matrix_free(Q_transpose);
    matrix_free(R);
    return beta;
}

void matrix_eigen_2x2(struct Matrix A, struct Matrix v1, struct Matrix v2,
                      float *lambda1, float *lambda2)
{
    assert(A.n[0] == 2); assert(A.n[1] == 2);
    assert(v1.n[0] == 2); assert(v1.n[1] == 1);
    assert(v2.n[0] == 2); assert(v2.n[1] == 1);
    float trace = MAT(A, 0, 0) + MAT(A,1,1);
    float det = MAT(A, 0, 0)*MAT(A,1,1) - MAT(A,0,1)*MAT(A, 1, 0);
    float gap = sqrtf(trace*trace - 4.f*det);
    float l1 = (trace + gap)*0.5f;
    float l2 = (trace - gap)*0.5f;
    
    MAT(v1,0,0) = MAT(A, 0, 0)-l2;
    MAT(v1,1,0) = MAT(A, 1, 0);
    if(MAT(v1, 0, 0) == 0.f &&MAT(v1, 1, 0) == 0.f){
        MAT(v1,0,0) = MAT(A, 0, 1);
        MAT(v1,1,0) = MAT(A, 1, 1)-l2;
    }
    
    MAT(v2,0,0) = MAT(A, 0, 0)-l1;
    MAT(v2,1,0) = MAT(A, 1, 0);
    if(MAT(v2, 0, 0) ==  0.f && MAT(v2, 1, 0) == 0.f){
        MAT(v2,0,0) = MAT(A, 0, 1);
        MAT(v2,1,0) = MAT(A, 1, 1)-l1;
    }

    *lambda1 = l1;
    *lambda2 = l2;
}
