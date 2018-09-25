#pragma once
#include "cart.h"
struct RandomForestSettings{
    int num_trees;
    int training_size;
    struct CARTSettings cart_settings;
};
struct RandomForest{
    int num_trees;
    int input_dim, output_dim;
    int num_samples_per_tree;
    struct CART *trees;
    int **trees_used_samples;
};
struct RandomForestConstructionProgress{
    int num_trees_constructed;
    struct CARTConstructionProgress cart_progress;
};
struct RandomForest rf_construct(int num_patterns, const float *inputs_in,
    int input_dim, const float *targets_in, int target_dim,
    struct RandomForestSettings settings, struct RNGState *rng);

struct RandomForest rf_construct_parallel(int num_patterns,
    void (*get_pattern)(float *input, float *target, int i, void *data),
    void *get_pattern_data, int input_dim, int target_dim,
    struct RandomForestSettings settings, int num_threads,
    struct RandomForestConstructionProgress *progress,
    struct RNGState **rngs);

void rf_free(struct RandomForest rf);
void rf_run(struct RandomForest rf, float *input, float *output,
    void (*evaluate_func)(float *input, int input_dim, float *output,
        int output_dim, float* parameters)
    );
void rf_run_default_evaluate(struct RandomForest rf, float *input,
    float *output);

size_t rf_to_mem_size(struct RandomForest rf);
void rf_to_mem(struct RandomForest rf, unsigned char *mem);
struct RandomForest rf_from_mem(unsigned char *mem, int max_num_trees);

void rf_oob_error(struct RandomForest rf, int num_patterns,
    void (*get_pattern)(float *input, float *target, int i, void *data),
    void *get_pattern_data, struct RNGState *rng, float *errors,
    int num_eval);

//More cache-friendly representation used for evaluation
struct RandomForestForEvaluation{
    struct CARTForEvaluation *trees;
    int num_trees;
    int input_dim, output_dim;
};

struct RandomForestForEvaluation rf_convert_for_evaluation(
    struct RandomForest rf);
void rf_evaluate(struct RandomForestForEvaluation rf, uint16_t *input,
    float *output);
#include "rng/rng.h"
void rf_evaluate_subset(struct RandomForestForEvaluation rf, uint16_t *input,
    float *output, int num, struct RNGState *rng);

size_t rf_eval_to_mem_size(struct RandomForestForEvaluation rf);
void rf_eval_to_mem(struct RandomForestForEvaluation rf, unsigned char *mem);
struct RandomForestForEvaluation rf_eval_from_mem(unsigned char *mem,
    int max_num_trees);
