#pragma once
#include "flags.h"
#if RF_OUTPUT_SIZE < 32
    #define CART_MAX_PARAM_DIM 32
#else
    #define CART_MAX_PARAM_DIM RF_OUTPUT_SIZE
#endif
#include "rng/rng.h"

struct CARTSettings{
    int num_features_per_node;
    int min_interior_size;
    float (*impurity_func)(float *inp, int input_dim, float *trg, int target_dim,
       int stride, int num_patterns, float *parameters);
    int random_splits;
    int max_depth;
};
struct CARTNode{
    int predictor;
    int left_node_index;
    float split;
    float val[CART_MAX_PARAM_DIM];
};
struct CART{
    int num_nodes;
    int input_dim, output_dim;
    struct CARTNode *nodes;
};
struct CARTConstructionProgress{
    int num_nodes_created;
    int max_depth_reached;
    int current_depth;
    int num_splits;
    int created_nodes_max_depth;
    unsigned char *created_nodes;
    // created_nodes should either be 0 or  a pointer to a buffer with
    // $2^{\text{created\_nodes\_max\_depth}+1}$ elements 
};
float cart_anova_impurity(float *inp, int input_dim, float *trg, int target_dim,
    int stride, int num_patterns, float *parameters);
struct CART cart_construct(int num_patterns, const float *inputs_in, int input_dim,
    const float *targets_in, int target_dim,
    struct CARTSettings settings, struct CARTConstructionProgress *progress,
    struct RNGState *rng);
void cart_free(struct CART cart);
void cart_default_evaluate(float *input, int input_dim, float *output,
    int output_dim, float* parameters);
void cart_run(struct CART cart, float *input, float *output,
    void (*evaluate_func)(float *input, int input_dim, float *output,
        int output_dim, float* parameters)
    );
void cart_run_default_evaluate(struct CART cart, float *input, float *output);
void cart_print(struct CART cart);

//More cache-friendly representation used for evaluation
#include <stdint.h>
#define CART_LEAF_SENTINEL16 0x8000
#define CART_LEAF_SENTINEL32 0x80000000
#define CART_SPLIT_MAX 65535
struct CARTNodeForEvaluation
{
    union{
        struct{
            //The split location
            uint16_t split;
            //The feature, or highest bit set if leaf
            uint16_t feature;
        };
        //The leaf index;
        uint32_t leaf_index;
    };
};
struct CARTSubtreeForEvaluation
{
    struct CARTNodeForEvaluation nodes[15];
    int offset;
};
struct CARTForEvaluation
{
    struct CARTSubtreeForEvaluation *subtrees;
    float *leaf_values;
    int num_subtrees;
    int num_leaf_values;
    int leaf_value_dim;
};

struct CARTForEvaluation cart_convert_for_evaluation(struct CART cart);
void cart_for_evaluation_print(struct CARTForEvaluation cart);
void cart_evaluate(struct CARTForEvaluation cart, uint16_t *input, float *output);
