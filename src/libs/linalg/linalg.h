#pragma once
#include "compiler_features/compiler_features.h"

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


// A 3x3 matrix look like this:
// [ 0  1  2]
// [ 3  4  5]
// [ 6  7  8]
//  NOTE(Vidar):Or should it be like this??
// [ 0  3  6]
// [ 1  4  7]
// [ 2  5  8]
struct ALIGNED_(32) Matrix3 {
    float m[9];
};
#define M3(mat,x,y) mat.m[x*3+y]
struct Vec3{
    union{
        float m[3];
        struct{
            float x,y,z;
        };
    };
};
struct Matrix4{
    float m[16];
};
#define M4(mat,x,y) mat.m[x*4+y]
struct Vec4{
    union{
        float m[4];
        struct{
        float x,y,z,w;
        };
    };
};

struct Matrix3 multiply_matrix3(struct Matrix3 a, struct Matrix3 b);
struct Matrix3 get_rotation_matrix3(float x,float y,float z);
struct Matrix3 invert_homogeneous_matrix3(struct Matrix3 mat);
struct Matrix3 lu_decompose_matrix3(struct Matrix3 mat);
struct Matrix3 get_translation_matrix3(float x,float y);
struct Matrix3 get_scale_matrix3(float s);
struct Matrix3 get_scale_matrix3_non_uniform(float sx, float sy);
struct Matrix3 get_identity_matrix3(void);
struct Matrix3 invert_matrix3(struct Matrix3 mat);
struct Matrix3 transpose_matrix3(struct Matrix3 mat);
struct Vec3 solve_LU_matrix3_vec3(struct Matrix3 lu, struct Vec3 y);
struct Matrix3 solve_LU_matrix3_matrix3(struct Matrix3 lu, struct Matrix3 y);
void print_matrix3(struct Matrix3);
void jacobi_diagonalize_matrix3(struct Matrix3 A, struct Matrix3 *eigenvectors,
    struct Vec3 *eigenvalues, int num_iterations);
float det_matrix3(struct Matrix3 A);

struct Matrix4 matrix3_to_matrix4(struct Matrix3);

struct Matrix4 multiply_matrix4(struct Matrix4 a,struct Matrix4 b);
struct Matrix4 get_trackball_matrix4(float theta, float phi);
struct Matrix4 get_aspect_correction_matrix4(float w,float h);
struct Matrix4 get_perspective_matrix4(float fov,float n,float f);
struct Matrix4 get_orthographic_matrix4(float scale,float n,float f);
struct Matrix4 get_translation_matrix4(float x,float y,float z);
struct Matrix4 get_rotation_matrix4(float x, float y, float z);
struct Matrix4 get_scale_matrix4(float s);
struct Matrix4 get_identity_matrix4(void);
struct Matrix4 transpose_matrix4(struct Matrix4);
struct Matrix4 invert_matrix4_noscale(struct Matrix4 mat);
struct Matrix4 lu_decompose_matrix4(struct Matrix4 mat);
struct Vec4 solve_LU_matrix4_vec4(struct Matrix4 lu, struct Vec4 y);
void print_matrix4(struct Matrix4);

struct Vec3 vec3(float x, float y, float z);
struct Vec3 normalize_vec3(struct Vec3);
float magnitude_vec3(struct Vec3 v);

struct Vec2{
    float x,y;
};
struct Vec2 transform_vector(struct Vec2 p, struct Vec2 x_axis, struct Vec2 y_axis);
struct Vec2 untransform_vector(struct Vec2 p, struct Vec2 x_axis, struct Vec2 y_axis);

struct Vec3 transform_vec3(struct Vec3 p, struct Vec3 x_axis, struct Vec3 y_axis,
    struct Vec3 z_axis);
struct Vec3 untransform_vec3(struct Vec3 p, struct Vec3 x_axis, struct Vec3 y_axis,
    struct Vec3 z_axis);


struct Vec2 add_vec2(struct Vec2 a, struct Vec2 b);
struct Vec2 sub_vec2(struct Vec2 a, struct Vec2 b);
struct Vec2 scale_vec2(float s, struct Vec2 a);
float dot_vec2(struct Vec2 a, struct Vec2 b);

struct Vec3 multiply_vec3_matrix3(struct Vec3 v,struct Matrix3 m);
struct Vec3 multiply_matrix3_vec3(struct Matrix3 m,struct Vec3 v);
struct Vec3 multiply_vec3_matrix4(struct Vec3 v, struct Matrix4 m);
struct Vec3 multiply_matrix4_vec3(struct Matrix4 m,struct Vec3 v);
struct Vec3 multiply_matrix4_vec3_point(struct Matrix4 m,struct Vec3 v);
struct Vec3 add_vec3(struct Vec3 a, struct Vec3 b);
struct Vec3 sub_vec3(struct Vec3 a, struct Vec3 b);
struct Vec3 scale_vec3(float s, struct Vec3 a);
float dot_vec3(struct Vec3 a, struct Vec3 b);
struct Vec3 cross_vec3(struct Vec3 a, struct Vec3 b);
