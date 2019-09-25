#include "linalg.h"
#include "sort/sort.h"
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

struct Matrix3 multiply_matrix3(struct Matrix3 A, struct Matrix3 B)
{
	struct Matrix3 ret={0};
	for(int i=0;i<3;i++){
	    for(int j=0;j<3;j++){
		for(int k=0;k<3;k++){
		    ret.m[i*3+j]+=A.m[i*3+k]*B.m[k*3+j];
		}
	    }
	}
	return ret;
}

struct Matrix3 get_rotation_matrix3(float x,float y,float z)
{
    float cx = cosf(x);
    float sx = sinf(x);
    float cy = cosf(y);
    float sy = sinf(y);
    float cz = cosf(z);
    float sz = sinf(z);
    struct Matrix3 rot_x =
    {
        1.f,0.f,0.f,
        0.f, cx,-sx,
        0.f, sx, cx,
    };
    struct Matrix3 rot_y =
    {
        cy,0.f, sy,
        0.f,1.f,0.f,
        -sy,0.f, cy,
    };
    struct Matrix3 rot_z =
    {
        cz,-sz,0.f,
        sz, cz,0.f,
        0.f,0.f,1.f,
    };
    struct Matrix3 tmp = multiply_matrix3(rot_x, rot_y);
    return multiply_matrix3(tmp, rot_z);
}

struct Matrix3 invert_homogeneous_matrix3(struct Matrix3 mat)
{
	float r11 = mat.m[0];
	float r21 = mat.m[1];
	float r12 = mat.m[3];
	float r22 = mat.m[4];
	float t1 = mat.m[6];
	float t2 = mat.m[7];
	float s1_sq = r11*r11 + r21 * r21;
	float s2_sq = r12*r12 + r22 * r22;
	struct Matrix3 ret = {
		r11/s1_sq, r12/s1_sq, 0.f,
		r21/s2_sq, r22/s2_sq, 0.f,
		-(t1*r11/s1_sq + t2*r21/s2_sq), -(t1*r12/s1_sq + t2*r22/s2_sq), 1.f
	};
	return ret;
}

struct Matrix3 lu_decompose_matrix3(struct Matrix3 A)
{
    struct Matrix3 L = {0};
    struct Matrix3 U = {
        1.f, 0.f, 0.f,
        0.f, 1.f, 0.f,
        0.f, 0.f, 1.f,
    };
    for(int i=0;i<3;i++){
        M3(L,i,0) = M3(A,i,0);
    }
    float L00 = M3(L,0,0);
    for(int i=1;i<3;i++){
        M3(U,0,i) = M3(A,0,i)/L00;
    }
    
    for(int j=1;j<3-1;j++){
        for(int i=j;i<3;i++){
            float sum = 0.f;
            for(int k=0;k<j;k++){
                sum += M3(L,i,k)*M3(U,k,j);
            }
            M3(L,i,j) = M3(A,i,j) - sum;
        }
        
        for(int k=j+1;k<3;k++){
            float sum = 0.f;
            for(int i=0;i<j;i++){
                sum += M3(L,j,i)*M3(U,i,k);
            }
            float tmp = (M3(A,j,k) - sum);
            if(fabs(tmp) < 1e-15f){
                M3(U, j, k) = 0.f;
            }else{
                float inv_ljj = 1.f/M3(L,j,j);
                M3(U,j,k) = tmp*inv_ljj;
            }
        }
    }
    float sum = 0.f;
    for(int k=0;k<2;k++){
        sum += M3(L,2,k) * M3(U,k,2);
    }
    M3(L,2,2) = M3(A,2,2) - sum;
    
    struct Matrix3 ret = {0};
    for(int c=0;c<3;c++){
        for(int r=0;r<3;r++){
            M3(ret,r,c)  = r < c ? M3(U,r,c) : M3(L,r,c);
        }
    }
    return ret;
}

struct Matrix3 get_translation_matrix3(float x, float y)
{
    struct Matrix3 ret;
    float m[8]=
    {
        1.f, 0.f,0.f,
        0.f, 1.f, 0.f,
        0.f, 0.f
    };
    memcpy(ret.m,m,8*sizeof(float));
	//ret = *(struct Matrix3*)m;
	ret.m[6] = x;
	ret.m[7] = y;
	ret.m[8] = 1.f;
    return ret;
}

struct Matrix3 get_scale_matrix3(float s)
{

    struct Matrix3 ret;
    float m[4]=
    {
        s,0.f,0.f, 0.f
    };
    memcpy(ret.m,m,4*sizeof(float));
    memcpy(ret.m+4,m,4*sizeof(float));
	ret.m[8] = 1.f;
	return ret;
}

struct Matrix3 get_scale_matrix3_non_uniform(float sx, float sy)
{
    struct Matrix3 ret;
    float mx[4]=
    {
        sx,0.f,0.f, 0.f
    };
    float my[4]=
    {
        sy,0.f,0.f, 0.f
    };
    memcpy(ret.m,mx,4*sizeof(float));
    memcpy(ret.m+4,my,4*sizeof(float));
	ret.m[8] = 1.f;
	return ret;
}

struct Matrix3 get_identity_matrix3(void)
{
    struct Matrix3 ret;
    float m[8]=
    {
        1.f,0.f,0.f,
		0.f,1.f,0.f,0.f,
		0.f
    };
    memcpy(ret.m,m,8*sizeof(float));
	ret.m[8] = 1.f;
	return ret;
}

static float det23(struct Matrix3 mat, int i, int j, int k, int l){
    return M3(mat, i,j)*M3(mat,k,l) - M3(mat,i,l)*M3(mat,k,j);
}

struct Matrix3 invert_matrix3(struct Matrix3 mat)
{
    struct Matrix3 ret = {
         det23(mat, 1, 1, 2, 2), -det23(mat, 0, 1, 2, 2), det23(mat, 0, 1, 1, 2),
        -det23(mat, 1, 0, 2, 2), det23(mat, 0, 0, 2, 2), -det23(mat, 0, 0, 1, 2),
         det23(mat, 1, 0, 2, 1), -det23(mat, 0, 0, 2, 1), det23(mat, 0, 0, 1, 1),

    };
    float det = M3(mat,0,0)*det23(mat, 1, 1, 2, 2) - M3(mat,0,1)*det23(mat, 1,0,2,2) + M3(mat,0,2)*det23(mat, 1, 0, 2, 1);
    for(int c=0;c<3;c++){
        for(int r=0;r<3;r++){
            M3(ret,r,c) /= det;
        }
    }
    return ret;
}

struct Matrix3 transpose_matrix3(struct Matrix3 mat)
{
    struct Matrix3 ret = {
        M3(mat,0,0), M3(mat,1,0), M3(mat,2,0),
        M3(mat,0,1), M3(mat,1,1), M3(mat,2,1),
        M3(mat,0,2), M3(mat,1,2), M3(mat,2,2),
    };
    return ret;
}


struct Vec3 solve_LU_matrix3_vec3(struct Matrix3 lu, struct Vec3 y)
{
    struct Vec3 x = {0};

    //Solve Lz = y
    struct Vec3 z = {0};
    for(int i=0;i<3;i++){
        float val = y.m[i];
        for(int j=0;j<i;j++){
            val -= M3(lu,i,j)*z.m[j];
        }
        val /= M3(lu,i,i);
        z.m[i] = val;
    }

    //Solve Ux = z
    for(int i=2;i>=0;i--){
        float val = z.m[i];
        for(int j=(i+1);j<3;j++){
            val -= M3(lu,i,j)*x.m[j];
        }
        x.m[i] = val;
    }

    return x;
}

struct Matrix3 solve_LU_matrix3_matrix3(struct Matrix3 lu, struct Matrix3 Y)
{
    struct Matrix3 ret = {0};
    for(int i_dim = 0; i_dim < 3; i_dim ++){
        struct Vec3 y = {0};
        for(int i=0;i<3;i++){
            y.m[i] = M3(Y,i,i_dim);
        }
        struct Vec3 x = {0};
    
        //Solve Lz = y
        struct Vec3 z = {0};
        for(int i=0;i<3;i++){
            float val = y.m[i];
            for(int j=0;j<i;j++){
                val -= M3(lu,i,j)*z.m[j];
            }
            float old_val = val;
            float div = M3(lu,i,i);
            if(div < 1e-15f){
                val = 0.f;
            }else{
                if(fabsf(val) > 1e-15f){
                    val /= div;
                }
            }
            if(isnan(val) || isinf(val)){
                printf("NAAAAN! %f\n", old_val);
            }
            z.m[i] = val;
        }
    
        //Solve Ux = z
        for(int i=2;i>=0;i--){
            float val = z.m[i];
            for(int j=(i+1);j<3;j++){
                val -= M3(lu,i,j)*x.m[j];
            }
            if(isnan(val)){
                printf("NAAAAN!\n");
            }
            x.m[i] = val;
        }
        for(int i=0;i<3;i++){
            M3(ret,i,i_dim) = x.m[i];
        }
    }
    
    
    return ret;
}

void print_matrix3(struct Matrix3 m)
{
    printf("[\n");
    for(int r=0;r<3;r++){
        for(int c=0;c<3;c++){
            printf("%3.3f ", M3(m,r,c));
        }
        printf(";\n");
    }
    printf("]\n");
}

static const int off_diagonal_i[] = { 0, 0, 1 };
static const int off_diagonal_j[] = { 1, 2, 2 };
static const int off_diagonal_k[] = { 2, 1, 0 };
static const int num_off_diagonal = 3;

void jacobi_diagonalize_matrix3(struct Matrix3 A, struct Matrix3 *eigenvectors,
                                struct Vec3 *eigenvalues, int num_iterations)
{
    struct Matrix3 V = {
        1.f, 0.f, 0.f,
        0.f, 1.f, 0.f,
        0.f, 0.f, 1.f,
    };
    struct Matrix3 G = {0};
    for(int i_iteration=0;i_iteration<num_iterations;i_iteration++){
        int i = 0;
        int j = 0;
        int k = 0;
        float Aij = 0.f;
        for(int d=0;d<num_off_diagonal;d++){
            int tmp_i = off_diagonal_i[d];
            int tmp_j = off_diagonal_j[d];
            float val = fabsf(M3(A, tmp_i, tmp_j));
            if(val > Aij){
                Aij = val;
                i = tmp_i;
                j = tmp_j;
                k = off_diagonal_k[d];
            }
        }
        if(Aij < 1e-15f){
            break;
        }
        Aij = M3(A,i,j);
        float Ajj = M3(A,j,j);
        float Aii = M3(A,i,i);
        float theta = atan2f(2.f*Aij, Ajj - Aii)*0.5f;
        float s = sinf(theta);
        float c = cosf(theta);
        memset(G.m, 0, 9*sizeof(float));
        M3(G,k,k) = 1.f;
        M3(G,i,i) = c;
        M3(G,j,j) = c;
        M3(G,i,j) = -s;
        M3(G,j,i) = s;
        //TODO(Vidar):Calculate this directly instead...
        A = multiply_matrix3(G, multiply_matrix3(A, transpose_matrix3(G)));
        V = multiply_matrix3(V, transpose_matrix3(G));
    }
    int indices[3];
    float abs_eigen[3];
    for(int i=0;i<3;i++){
        abs_eigen[i] = fabsf(M3(A,i,i));
        indices[i] = i;
    }
    float tmp[3] = {0};
    int tmp_i[3] = {0};
    merge_sort_auxillary(abs_eigen, tmp, indices, tmp_i, 3);
    for(int i=0;i<3;i++){
        int i2 = indices[i];
        eigenvalues->m[i] = M3(A,i2,i2);
        for(int j=0;j<3;j++){
            M3((*eigenvectors),j,i) = M3(V,j,i2);
        }
    }
}

float det_matrix3(struct Matrix3 A)
{
    float det = M3(A,0,0)*det23(A, 1, 1, 2, 2) - M3(A,0,1)*det23(A, 1,0,2,2) + M3(A,0,2)*det23(A, 1, 0, 2, 1);
    return det;
}

struct Matrix4 matrix3_to_matrix4(struct Matrix3 m)
{
    struct Matrix4 ret ={
        M3(m,0,0), M3(m,0,1), M3(m,0,2), 0.f,
        M3(m,1,0), M3(m,1,1), M3(m,1,2), 0.f,
        M3(m,2,0), M3(m,2,1), M3(m,2,2), 0.f,
        0.f,       0.f,       0.f,       1.0f
    };
    return ret;
}

struct Matrix4 multiply_matrix4(struct Matrix4 a,struct Matrix4 b)
{
    struct Matrix4 ret={0};
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            for(int k=0;k<4;k++){
                ret.m[i*4+j]+=a.m[i*4+k]*b.m[k*4+j];
            }
        }
    }
    return ret;
}


struct Matrix4 get_trackball_matrix4(float theta, float phi)
{
    struct Matrix4 ret;
    float proj[4][4]=
    {

         {cosf(theta), -sinf(theta)*cosf(phi), sinf(theta)*sinf(phi), 0.0f},
         {0.0f, sinf(phi), cosf(phi), 0.0f},
         {-sinf(theta), -cosf(theta)*cosf(phi), cosf(theta)*sinf(phi), 0.0f},
         {0.0f, 0.0f, 0.0f, 1.0f},
         
        /*
         {cosf(theta),               0.0f,         -sinf(theta),              0.0f},
         {-sinf(theta)*cosf(phi),    sinf(phi),    -cosf(theta)*cosf(phi),    0.0f},
         {sinf(theta)*sinf(phi),     cosf(phi),    cosf(theta)*sinf(phi),     0.0f},
         {0.0f,                      0.0f,         0.0f,                      1.0f},
         
        {cosf(theta), -sinf(theta)*cosf(phi), sinf(theta)*sinf(phi), 0.0f},
        {sinf(theta), cosf(theta)*cosf(phi), -cosf(theta)*sinf(phi), 0.0f},
        {0.0f, sinf(phi), cosf(phi), 0.0f},
        {0.0f, 0.0f, 0.0f, 1.0f},
         */
    };
    memcpy(ret.m,proj,16*sizeof(float));
    return ret;
}

struct Matrix4 get_aspect_correction_matrix4(float w,float h)
{
    float aspect=w/h;
    struct Matrix4 ret={0};
    if(aspect<1.f){
        ret.m[0]=1.f;
        ret.m[5]=aspect;
    } else{
        ret.m[0]=1.f/aspect;
        ret.m[5]=1.f;
    }
    ret.m[10]=1.f;
    ret.m[15]=1.f;
    return ret;
}

struct Matrix4 get_orthographic_matrix4(float scale,float n,float f)
{
    struct Matrix4 ret;
    float proj[4][4]=
    {
        {1.f/scale,0.f,0.f,0.0f},
        {0.f, 1.f/scale, 0.f, 0.0f},
        {0.0f, 0.f, -2.f/(f-n), -(f+n)/(f-n)},
        {0.0f, 0.0f, 0.0f, 1.0f},
    };
    memcpy(ret.m,proj,16*sizeof(float));
    return ret;
}

struct Matrix4 get_perspective_matrix4(float fov,float n,float f)
{
    float top=tanf(fov*0.5f);
    struct Matrix4 ret;
    float proj[4][4]=
    {
        {1.f/top,0.f,0.f,0.0f},
        {0.f, 1.f/top, 0.f, 0.0f},
        {0.0f, 0.f, -(f+n)/(f-n), -1.f},
        {0.0f, 0.0f, -2.f*f*n/(f-n), 0.0f},
    };
    memcpy(ret.m,proj,16*sizeof(float));
    return ret;
}

struct Matrix4 get_translation_matrix4(float x,float y,float z)
{
    struct Matrix4 ret;
    float proj[4][4]=
    {
        {1.f,0.f,0.f, 0.f},
        {0.f, 1.f, 0.f, 0.f},
        {0.0f, 0.f, 1.f, 0.f},
        {x, y, z, 1.0f},
    };
    memcpy(ret.m,proj,16*sizeof(float));
    return ret;
}

struct Matrix4 get_rotation_matrix4(float x,float y,float z)
{
    float cx = cosf(x);
    float sx = sinf(x);
    float cy = cosf(y);
    float sy = sinf(y);
    float cz = cosf(z);
    float sz = sinf(z);
    struct Matrix4 rot_x =
    {
        1.f,0.f,0.f, 0.f,
        0.f, cx,-sx, 0.f,
        0.f, sx, cx, 0.f,
		0.f, 0.f, 0.f, 1.f
    };
    struct Matrix4 rot_y =
    {
        cy,0.f, sy, 0.f,
        0.f,1.f,0.f, 0.f,
        -sy,0.f, cy, 0.f,
		0.f, 0.f, 0.f, 1.f
    };
    struct Matrix4 rot_z =
    {
        cz,-sz,0.f, 0.f,
        sz, cz,0.f, 0.f,
        0.f,0.f,1.f, 0.f,
		0.f, 0.f, 0.f, 1.f
    };
    struct Matrix4 tmp = multiply_matrix4(rot_x, rot_y);
    return multiply_matrix4(tmp, rot_z);
}

struct Matrix4 get_scale_matrix4(float s)
{
    struct Matrix4 ret;
    float proj[4][4]=
    {
        {s  , 0.f,0.f, 0.f},
        {0.f, s  , 0.f, 0.f},
        {0.f, 0.f, s  , 0.f},
        {0.f, 0.f, 0.f, 1.0f},
    };
    memcpy(ret.m,proj,16*sizeof(float));
    return ret;
}

struct Matrix4 get_scale_matrix4_non_uniform(float x, float y, float z)
{
    struct Matrix4 ret;
    float proj[4][4]=
    {
        {x  , 0.f,0.f, 0.f},
        {0.f, y  , 0.f, 0.f},
        {0.f, 0.f, z  , 0.f},
        {0.f, 0.f, 0.f, 1.0f},
    };
    memcpy(ret.m,proj,16*sizeof(float));
    return ret;
}

struct Matrix4 get_identity_matrix4(void)
{
    struct Matrix4 ret;
    float proj[4][4] =
    {
        { 1.f, 0.f, 0.f, 0.0f },
        { 0.f, 1.f, 0.f, 0.0f },
        { 0.0f, 0.f, 1.f, 0.f },
        { 0.f , 0.f, 0.f, 1.0f },
    };
    memcpy(ret.m, proj, 16 * sizeof(float));
    return ret;
}

struct Matrix4 transpose_matrix4(struct Matrix4 m)
{
    struct Matrix4 ret = {
        M4(m,0,0), M4(m,1,0), M4(m,2,0), M4(m,3,0),
        M4(m,0,1), M4(m,1,1), M4(m,2,1), M4(m,3,1),
        M4(m,0,2), M4(m,1,2), M4(m,2,2), M4(m,3,2),
        M4(m,0,3), M4(m,1,3), M4(m,2,3), M4(m,3,3),
    };
    return ret;
}

struct Matrix4 invert_matrix4_noscale(struct Matrix4 m)
{
	//NOTE(Vidar): Transpose the rotation and negate the translation
    struct Matrix3 rot = {
         M4(m,0,0), M4(m,0,1), M4(m,0,2),
         M4(m,1,0), M4(m,1,1), M4(m,1,2),
         M4(m,2,0), M4(m,2,1), M4(m,2,2),
    };
	struct Vec3 trans = {
		 -M4(m,3,0), -M4(m,3,1), -M4(m,3,2)
	};
	trans = multiply_matrix3_vec3(rot, trans);
    struct Matrix4 ret = {
        M3(rot,0,0), M3(rot,1,0), M3(rot,2,0), 0.f,
        M3(rot,0,1), M3(rot,1,1), M3(rot,2,1), 0.f,
        M3(rot,0,2), M3(rot,1,2), M3(rot,2,2), 0.f,
		trans.x,     trans.y,     trans.z,     1.f
    };
	return ret;
}

//https://www.gamedev.net/resources/_/technical/math-and-physics/matrix-inversion-using-lu-decomposition-r3637
//TODO(Vidar):Permutation matrix...
struct Matrix4 lu_decompose_matrix4(struct Matrix4 A)
{
    struct Matrix4 L = {0};
    struct Matrix4 U = get_identity_matrix4();
    for(int i=0;i<4;i++){
        M4(L,i,0) = M4(A,i,0);
    }
    float L00 = M4(L,0,0);
    for(int i=1;i<4;i++){
        M4(U,0,i) = M4(A,0,i)/L00;
    }
    
    for(int j=1;j<4-1;j++){
        for(int i=j;i<4;i++){
            float sum = 0.f;
            for(int k=0;k<j;k++){
                sum += M4(L,i,k)*M4(U,k,j);
            }
            M4(L,i,j) = M4(A,i,j) - sum;
        }
        float inv_ljj = 1.f/M4(L,j,j);
        for(int k=j+1;k<4;k++){
            float sum = 0.f;
            for(int i=0;i<j;i++){
                sum += M4(L,j,i)*M4(U,i,k);
            }
            M4(U,j,k) = (M4(A,j,k) - sum)*inv_ljj;
        }
    }
    float sum = 0.f;
    for(int k=0;k<3;k++){
        sum += M4(L,3,k) * M4(U,k,3);
    }
    M4(L,3,3) = M4(A,3,3) - sum;
    
    struct Matrix4 ret = {0};
    for(int c=0;c<4;c++){
        for(int r=0;r<4;r++){
            M4(ret,r,c)  = r < c ? M4(U,r,c) : M4(L,r,c);
        }
    }
    return ret;
}

//Finds x in Ax = y where A = LU
struct Vec4 solve_LU_matrix4_vec4(struct Matrix4 lu, struct Vec4 y)
{
    struct Vec4 x = {0};
    
    //Solve Lz = y
    struct Vec4 z = {0};
    for(int i=0;i<4;i++){
        float val = y.m[i];
        for(int j=0;j<i;j++){
            val -= M4(lu,i,j)*z.m[j];
        }
        val /= M4(lu,i,i);
        z.m[i] = val;
    }
    
    //Solve Ux = z
    for(int i=3;i>=0;i--){
        float val = z.m[i];
        for(int j=(i+1);j<4;j++){
            val -= M4(lu,i,j)*x.m[j];
        }
        x.m[i] = val;
    }
    
    return x;
}

void print_matrix4(struct Matrix4 m)
{
    for(int r=0;r<4;r++){
        printf("[ ");
        for(int c=0;c<4;c++){
            printf("%3.3f ", M4(m,r,c));
        }
        printf("]\n");
    }
}

struct Vec3 vec3(float x, float y, float z)
{
    struct Vec3 ret = {x,y,z};
    return ret;
}

struct Vec3 normalize_vec3(struct Vec3 v)
{
    float m = sqrtf(v.x*v.x + v.y*v.y + v.z*v.z);
    v.x /= m; v.y /= m; v.z /= m;
    return v;
}

float magnitude_vec3(struct Vec3 v)
{
    return sqrtf(v.x*v.x + v.y*v.y + v.z*v.z);
}

struct Vec2 transform_vector(struct Vec2 p, struct Vec2 x_axis, struct Vec2 y_axis)
{
    struct Vec2 ret = {0};
    ret.x = x_axis.x*p.x + x_axis.y*p.y;
    ret.y = y_axis.x*p.x + y_axis.y*p.y;
    return ret;
}

struct Vec2 untransform_vector(struct Vec2 p, struct Vec2 x_axis, struct Vec2 y_axis)
{
    struct Vec2 ret = {0};
    ret.x = x_axis.x*p.x + y_axis.x*p.y;
    ret.y = x_axis.y*p.x + y_axis.y*p.y;
    return ret;
}

struct Vec3 transform_vec3(struct Vec3 p, struct Vec3 x_axis,
    struct Vec3 y_axis, struct Vec3 z_axis)
{
    struct Vec3 ret = {0};
    ret.x = x_axis.x*p.x + x_axis.y*p.y + x_axis.z*p.z;
    ret.y = y_axis.x*p.x + y_axis.y*p.y + y_axis.z*p.z;
    ret.z = z_axis.x*p.x + z_axis.y*p.y + z_axis.z*p.z;
    return ret;
}

struct Vec3 untransform_vec3(struct Vec3 p, struct Vec3 x_axis,
    struct Vec3 y_axis, struct Vec3 z_axis)
{
    struct Vec3 ret = {0};
    ret.x = x_axis.x*p.x + y_axis.x*p.y + z_axis.x*p.z;
    ret.y = x_axis.y*p.x + y_axis.y*p.y + z_axis.y*p.z;
    ret.z = x_axis.z*p.x + y_axis.z*p.y + z_axis.z*p.z;
    return ret;
}


struct Vec2 add_vec2(struct Vec2 a, struct Vec2 b)
{
    struct Vec2 ret = {a.x + b.x, a.y + b.y};
    return ret;
}

struct Vec2 sub_vec2(struct Vec2 a, struct Vec2 b)
{
    struct Vec2 ret = {a.x - b.x, a.y - b.y};
    return ret;
}

struct Vec2 scale_vec2(float s, struct Vec2 a)
{
    struct Vec2 ret = {a.x*s, a.y*s};
    return ret;
}

float dot_vec2(struct Vec2 a, struct Vec2 b)
{
    return a.x*b.x + a.y*b.y;
}

struct Vec3 multiply_vec3_matrix3(struct Vec3 v,struct Matrix3 m)
{
    struct Vec3 ret = {0};
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            ret.m[i]+=M3(m,j,i)*v.m[j];
        }
    }
    return ret;
}

struct Vec3 multiply_matrix3_vec3(struct Matrix3 m,struct Vec3 v)
{
    struct Vec3 ret = {0};
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            ret.m[i]+=M3(m,i,j)*v.m[j];
        }
    }
    return ret;
}

struct Vec3 multiply_vec3_matrix4(struct Vec3 v, struct Matrix4 m)
{
    struct Vec3 ret = {0};
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            ret.m[i]+=M4(m,j,i)*v.m[j];
        }
    }
    return ret;
}

struct Vec3 multiply_matrix4_vec3(struct Matrix4 m,struct Vec3 v)
{
    struct Vec3 ret = {0};
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            ret.m[i]+=M4(m,i,j)*v.m[j];
        }
    }
    return ret;
}

struct Vec3 multiply_matrix4_vec3_point(struct Matrix4 m,struct Vec3 v)
{
    struct Vec4 v_tmp = {v.x,v.y,v.z,1.f};
    struct Vec4 res = {0};
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            res.m[i]+=M4(m,i,j)*v_tmp.m[j];
        }
    }
    struct Vec3 ret = {0};
    for(int i=0;i<3;i++){
        ret.m[i] = res.m[i]/res.m[3];
    }
    return ret;
}

struct Vec3 add_vec3(struct Vec3 a, struct Vec3 b)
{
    struct Vec3 ret = {a.x+b.x, a.y+b.y, a.z+b.z};
    return ret;
}

struct Vec3 sub_vec3(struct Vec3 a, struct Vec3 b)
{
    struct Vec3 ret = {a.x-b.x, a.y-b.y, a.z-b.z};
    return ret;
}

struct Vec3 scale_vec3(float s, struct Vec3 a)
{
    struct Vec3 ret = {s*a.x, s*a.y, s*a.z};
    return ret;
}

float dot_vec3(struct Vec3 a, struct Vec3 b)
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

struct Vec3 cross_vec3(struct Vec3 a, struct Vec3 b)
{
    struct Vec3 ret = {
        a.y*b.z-a.z*b.y, -(a.x*b.z-a.z*b.x), a.x*b.y-a.y*b.x
    };
    return ret;
}

