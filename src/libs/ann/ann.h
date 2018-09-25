#pragma once

enum ActivationFunction{
    ACTIVATION_FUNCTION_TANH,
    ACTIVATION_FUNCTION_RELU,
    ACTIVATION_FUNCTION_LINEAR,
};

struct ANNLayerDescription{
    int num_neurons;
    enum ActivationFunction activation_function;
};

struct ANNLayer{
    float *weights;
    float *neurons;
    float *local_field;
    float *local_field_deriv;
    float *weight_updates;
    float *weight_inertia;
    float *gradient_running_average;
    int num_neurons;
    int num_left_neurons;
    int num_right_neurons;
    int num_weights;
    enum ActivationFunction activation_function;
};

struct ANN{
    int num_layers;
    struct ANNLayer *layers;
};

struct ANN *ann_create(int num_layers,
    const struct ANNLayerDescription *descriptions);
void ann_free(struct ANN* ann);
void ann_train(struct ANN *ann, int num_iterations,
    int num_patterns, float *inputs, float *target_outputs, float *energy_result);
void ann_train_parallel(int num_threads, struct ANN *ann, int num_syncs,
    int num_patterns, float *inputs, float *target_outputs);
void ann_run(struct ANN *ann, float *input, float *output);
void ann_print_weights(struct ANN *ann);
void ann_copy_weights(struct ANN *dest, struct ANN *src);
float ann_get_rmse(struct ANN *ann, int num_points, float *inputs,
    float *target_outputs);
