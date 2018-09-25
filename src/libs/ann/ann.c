#include "ann.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "thread/thread.h"
#include "fp_except/fp_except.h"
#include <immintrin.h>

struct ANN *ann_create(int num_layers, const struct ANNLayerDescription *descriptions)
{
    struct ANN *ann = calloc(1,sizeof(struct ANN));
    ann->num_layers = num_layers;
    ann->layers = calloc(num_layers, sizeof(struct ANNLayer));
    for(int i=0;i<num_layers;i++){
        struct ANNLayer *layer = ann->layers + i;
        layer->num_neurons = descriptions[i].num_neurons;
        layer->neurons     = calloc(layer->num_neurons+1,sizeof(float));
        layer->neurons[layer->num_neurons] = 1.f; // Bias
        layer->activation_function = descriptions[i].activation_function;
        if(i<num_layers-1){
            layer->num_right_neurons = descriptions[i+1].num_neurons;
        }
        if(i>0){
            layer->num_left_neurons = descriptions[i-1].num_neurons+1;
            layer->local_field = calloc(layer->num_neurons,sizeof(float));
            layer->local_field_deriv = calloc(layer->num_neurons,sizeof(float));
            int num_weights = layer->num_neurons * (layer->num_left_neurons);
            layer->weights = calloc(num_weights,sizeof(float));
            layer->weight_updates = calloc(num_weights,sizeof(float));
            layer->weight_inertia = calloc(num_weights,sizeof(float));
            layer->gradient_running_average = calloc(num_weights,sizeof(float));
            layer->num_weights    = num_weights;
            //Weight initialization according to http://proceedings.mlr.press/v9/glorot10a/glorot10a.pdf eq[16]
            float weight_scale = sqrtf(6.f/(float)(layer->num_left_neurons +
                layer->num_neurons));
            for(int j=0;j<num_weights;j++){
                //TODO(Vidar):Better initial weights
                float r = (float)rand()/(float)RAND_MAX*2.f - 1.f;
                layer->weights[j] = r*weight_scale;
            }
        }
    }
    return ann;
}

void ann_free(struct ANN* ann)
{
    for(int i= 0;i<ann->num_layers;i++){
        free(ann->layers[i].neurons);
        free(ann->layers[i].weights);
    }
    free(ann->layers);
    free(ann);
}

struct ANN *ann_clone(struct ANN *ann)
{
    int num_layers = ann->num_layers;
    struct ANNLayerDescription *layers = calloc(num_layers,
        sizeof(struct ANNLayerDescription));
    for(int i =0;i<num_layers;i++){
        layers[i].num_neurons = ann->layers[i].num_neurons;
        layers[i].activation_function = ann->layers[i].activation_function;
    }
    struct ANN *ret = ann_create(num_layers,layers);
    free(layers);
    return ret;
}

//Assumes that the two networks have identical structure...
void ann_copy_weights(struct ANN *dest, struct ANN *src){
    for(int i=1;i<dest->num_layers;i++){
        struct ANNLayer *dst_layer = dest->layers + i;
        struct ANNLayer *src_layer = src->layers  + i;
        int num_weights = dst_layer->num_weights;
        for(int j=0;j<num_weights;j++){
            dst_layer->weights[j] = src_layer->weights[j];
        }
    }
}

static const float beta = 1.f;
static const float max_perf_inc = 1.04f;
static const float lr_dec = 0.7f;
static const float lr_inc = 1.05f;
static const float momentum = 0.2f;//1e-4f;
static const float weight_decay = 0.f;//.001f;
void ann_train(struct ANN *ann, int num_iterations, 
    int num_patterns, float *inputs, float *target_outputs, float *energy_result)
{

    //fp_except_enable();

    float eta = 0.001f;
    float forgetting_factor = 0.9f;
    int input_size = ann->layers[0].num_neurons;
    int output_size = ann->layers[ann->num_layers-1].num_neurons;


    float *tmp_input  = calloc(input_size,  sizeof(float));
    float *tmp_output = calloc(output_size, sizeof(float));

    int num_layers = ann->num_layers;
    int max_num_neurons = 0;
    for(int i_layer = num_layers-1; i_layer >= 0; i_layer--){
        int num_neurons = ann->layers[i_layer].num_neurons;
        max_num_neurons = num_neurons > max_num_neurons ? num_neurons :
            max_num_neurons;
    }

    float *curr_deltas = calloc(max_num_neurons,sizeof(float));
    float *prev_deltas = calloc(max_num_neurons,sizeof(float));

    int end_layer = num_layers-1;
    for(int i_pattern = 0; i_pattern < num_iterations; i_pattern++){
        int pattern = rand() % num_patterns;
        float *input = inputs + pattern*input_size;
        ann_run(ann,input,0);
        for(int i_layer = end_layer; i_layer > 0; i_layer--){
            struct ANNLayer *layer = ann->layers + i_layer;
            int num_neurons = layer->num_neurons;
            struct ANNLayer *left_layer = layer - 1;
            int num_left_neurons = layer->num_left_neurons;
            struct ANNLayer *right_layer = layer + 1;
            int num_right_neurons = layer->num_right_neurons;
            for(int i=0;i<num_neurons;i++){
                float g_prime = layer->local_field_deriv[i];

                float delta;
                if(i_layer == end_layer){
                    delta = (target_outputs[pattern*output_size + i]-
                        layer->neurons[i]) * g_prime;
                }else{
                    delta = 0.f;
                    for(int j=0;j<num_right_neurons;j++){
                        delta += prev_deltas[j] *
                            right_layer->weights[j*num_neurons + i];
                    }
                    delta *= g_prime;
                }

                int a = i*num_left_neurons;
                for(int j=0;j<num_left_neurons;j++)
                {
                    float w =  delta * left_layer->neurons[j];
                    layer->weight_updates[a + j] = w;
                }
                curr_deltas[i] = delta;
            }
            float *tmp = curr_deltas;
            curr_deltas = prev_deltas;
            prev_deltas = tmp;
        }

        for(int i_layer = end_layer; i_layer > 0; i_layer--){
            struct ANNLayer *layer = ann->layers + i_layer;
            int num_weights = layer->num_weights;

            for(int i=0;i<num_weights;i++){
                float grad = layer->weight_updates[i];
                float d = 0.f;
                if((0)){
                    //rmsprop according to http://www.cs.toronto.edu/~tijmen/csc321/slides/lecture_slides_lec6.pdf
                    layer->gradient_running_average[i] =
                        forgetting_factor*layer->gradient_running_average[i]
                        + (1.f-forgetting_factor)*grad*grad;
                    d = eta * (sqrtf(layer->gradient_running_average[i]) + 1e-8);
                }else{
                    d = eta;
                }
                layer->weights[i] += d * (grad 
                    + momentum * layer->weight_inertia[i]
                    - weight_decay * layer->weights[i]);
                layer->weight_inertia[i] = grad;
            }
        }
    }
    free(tmp_input);
    free(tmp_output);
    free(curr_deltas);
    free(prev_deltas);

    //fp_except_disable();
}

void ann_run(struct ANN *ann, float *input, float *output)
{
    int num_inputs = ann->layers[0].num_neurons;
    for(int i=0;i<num_inputs;i++){
        ann->layers[0].neurons[i] = input[i];
    }
    int num_layers = ann->num_layers;
    for(int i=1;i<num_layers;i++){

        struct ANNLayer *layer = ann->layers+i;
        struct ANNLayer *prev_layer = ann->layers+i-1;

        float *neurons = layer->neurons;
        int num_neurons = layer->num_neurons;

        float *prev_neurons = prev_layer->neurons;
        int prev_num_neurons = prev_layer->num_neurons;

        for(int j=0;j<num_neurons;j++){
            float *weights = layer->weights + j*layer->num_left_neurons;
            float b = 0.f;
            float b_prime = 0.f;
            for(int k=0;k<layer->num_left_neurons;k++){
                b += prev_neurons[k]*weights[k];
            }
            switch(layer->activation_function){
                case ACTIVATION_FUNCTION_TANH:
                    b = tanhf(beta*b);
                    b_prime = beta*(1.f - b*b);
                    break;
                case ACTIVATION_FUNCTION_RELU:
                    b = b > 0.f ? b : 0.f;
                    b_prime = b > 0.f ? 1.f : 0.f;
                    break;
                case ACTIVATION_FUNCTION_LINEAR:
                    b = b;
                    b_prime = 1.f;
                    break;
            }
            layer->local_field[j] = b;
            layer->local_field_deriv[j] = b_prime;
            neurons[j] = b;
        }
    }
    if(output){
        struct ANNLayer *output_layer = ann->layers+ann->num_layers-1;
        int num_outputs = output_layer->num_neurons;
        for(int i=0;i<num_outputs;i++){
            output[i] = output_layer->neurons[i];
        }
    }
}

static void ann_run_avx2(struct ANN *ann, float *input, float *output)
{
    __m256 _one = _mm256_set1_ps(1.f);
    int num_inputs = ann->layers[0].num_neurons;
    for(int i=0;i<num_inputs;i++){
        __m256 tmp = _mm256_load_ps(input+i*8);
        _mm256_store_ps(ann->layers[0].neurons+i*8, tmp);
    }
    int num_layers = ann->num_layers;
    for(int i=1;i<num_layers;i++){

        struct ANNLayer *layer = ann->layers+i;
        struct ANNLayer *prev_layer = ann->layers+i-1;

        float *neurons = layer->neurons;
        int num_neurons = layer->num_neurons;

        float *prev_neurons = prev_layer->neurons;

        for(int j=0;j<num_neurons;j++){
            float *weights = layer->weights + j*layer->num_left_neurons;
            __m256 b = _mm256_setzero_ps();
            __m256 b_prime = _mm256_setzero_ps();
            for(int k=0;k<layer->num_left_neurons;k++){
                __m256 p = _mm256_load_ps(prev_neurons+k*8);
                __m256 w = _mm256_set1_ps(weights[k]);
                b = _mm256_add_ps(b,_mm256_mul_ps(p,w));
                //b += prev_neurons[k]*weights[k];
            }
            switch(layer->activation_function){
                case ACTIVATION_FUNCTION_TANH:
                    {
                        /*
                        __m256 _beta = _mm256_set1_ps(beta);
                        b = _mm256_tanh_ps(_mm256_mul_ps(_beta,b));
                        b_prime = _mm256_mul_ps(_beta,
                            _mm256_sub_ps(_one, _mm256_mul_ps(b,b)));
                        */
                    }
                    break;
                case ACTIVATION_FUNCTION_RELU:
                    {
                        __m256 _zero = _mm256_setzero_ps();
                        b_prime =_mm256_and_ps(_one,
                            _mm256_cmp_ps(b,_zero,_CMP_GE_OS));
                        b = _mm256_max_ps(b,_zero);
                    }
                    break;
                case ACTIVATION_FUNCTION_LINEAR:
                    //b = b;
                    b_prime = _one;
                    break;
            }
            _mm256_store_ps(layer->local_field+j*8,b);
            _mm256_store_ps(layer->local_field_deriv+j*8,b_prime);
            _mm256_store_ps(neurons+j*8,b);
        }
    }
}

void ann_print_weights(struct ANN *ann)
{
    int num_layers = ann->num_layers;
    for(int i=1;i<num_layers;i++){
        int num_neurons = ann->layers[i].num_neurons;
        int num_left_neurons = ann->layers[i].num_left_neurons;
        printf("\nWeights between layer %d and %d\n", i-1, i);
        for(int j=0;j<num_neurons;j++){
            for(int k=0;k<num_left_neurons;k++){
                printf("%f ", ann->layers[i].weights[j*num_left_neurons + k]);
            }
            printf("\n");
        }
    }
}

struct ANNThreadData{
    int id;
    struct ANN *ann;
    int num_patterns;
    int num_iterations;
    int quit;
    float *inputs, *target_outputs;
    struct Semaphore *sync_start_semaphore;
    struct Semaphore *sync_done_semaphore;
};

static unsigned long ann_thread_func(void *data)
{
    struct ANNThreadData *thread_data = (struct ANNThreadData*)data;
    while(!thread_data->quit){
        ann_train(thread_data->ann, thread_data->num_iterations,
            thread_data->num_patterns, thread_data->inputs,
            thread_data->target_outputs,0);
        thread_semaphore_post(thread_data->sync_start_semaphore);
        thread_semaphore_wait(thread_data->sync_done_semaphore);
    }
    return 0;
}

void ann_train_parallel(int num_threads, struct ANN *ann, int num_syncs,
    int num_patterns, float *inputs, float *target_outputs)
{
    int input_size = ann->layers[0].num_neurons;
    int output_size = ann->layers[ann->num_layers-1].num_neurons;

    struct ANN *ann_working = ann_clone(ann);

    int num_validation_patterns = num_patterns/10;
    int num_training_patterns   = num_patterns - num_validation_patterns;

    char *sem_name_template = "sem_%c %d";
    size_t len = strlen(sem_name_template)+10;
    char *buffer = calloc(len,sizeof(char));
    float rmse_best = ann_get_rmse(ann,num_validation_patterns,
        inputs + num_training_patterns*input_size,
        target_outputs + num_training_patterns*output_size
        );
    printf("initial rmse: %f\n", rmse_best);
    rmse_best = FLT_MAX;
    struct ThreadHandle **threads = calloc(num_threads,
        sizeof(struct ThreadHandle*));
    struct ANNThreadData **thread_data = calloc(num_threads,
        sizeof(struct ANNThreadData*));
    for(int i=0;i<num_threads;i++){
        thread_data[i] = calloc(1,sizeof(struct ANNThreadData));
        thread_data[i]->id = i;
        thread_data[i]->num_patterns = num_training_patterns;
        thread_data[i]->inputs = inputs;
        thread_data[i]->quit = 0;
        thread_data[i]->target_outputs = target_outputs;
        thread_data[i]->ann = ann_clone(ann_working);
        thread_data[i]->num_iterations = 5e4;
        sprintf(buffer, sem_name_template, 's', i);
        thread_data[i]->sync_start_semaphore = thread_semaphore_create(0,buffer);
        sprintf(buffer, sem_name_template, 'w', i);
        thread_data[i]->sync_done_semaphore = thread_semaphore_create(0,buffer);
        threads[i] = thread_start(ann_thread_func, thread_data[i]);
    }
    free(buffer);

    for(int i_syncs = 0; i_syncs < num_syncs; i_syncs++){
        for(int i=0;i<num_threads;i++){
            //printf("Waiting for thread %d\n", i);
            thread_semaphore_wait(thread_data[i]->sync_start_semaphore);
        }

        float inv_num_threads = 1.f/(float)num_threads;
        int num_layers = ann_working->num_layers;
        for(int i=0;i<num_layers;i++){
            struct ANNLayer *layer = ann_working->layers + i;
            int num_weights = layer->num_weights;
            float *weights = layer->weights;
            for(int j=0;j<num_weights;j++){
                weights[j] = 0.f;
                //TODO(Vidar):move the loop over threads out...
                for(int k=0;k<num_threads;k++){
                    weights[j] += thread_data[k]->ann->layers[i].weights[j]
                        * inv_num_threads;
                }
            }
        }

        //printf("Updating weights\n");
        for(int i=0;i<num_threads;i++){
            ann_copy_weights(thread_data[i]->ann,ann_working);
            if(i_syncs == num_syncs-1){
                thread_data[i]->quit = 1;
            }
            thread_semaphore_post(thread_data[i]->sync_done_semaphore);
        }
        float rmse = ann_get_rmse(ann_working,num_validation_patterns,
                inputs + num_training_patterns*input_size,
                target_outputs + num_training_patterns*output_size
                );
        if(rmse < rmse_best){
            rmse_best = rmse;
            //printf("Updated weights, rmse: %f\n", rmse);
            ann_copy_weights(ann, ann_working);
        }
        if(i_syncs % 10 == 0){
            printf("rmse: %f\n", rmse);
        }
    }

    for(int i=0;i<num_threads;i++){
        thread_wait(threads[i]);
        thread_semaphore_destroy(thread_data[i]->sync_start_semaphore);
        thread_semaphore_destroy(thread_data[i]->sync_done_semaphore);
        free(thread_data[i]);
    }

    float rmse_final = ann_get_rmse(ann,num_validation_patterns,
        inputs + num_training_patterns*input_size,
        target_outputs + num_training_patterns*output_size
        );
    printf("final rmse: %f\n", rmse_final);

    free(thread_data);
    ann_free(ann_working);

    ann_print_weights(ann);

    printf("All done!\n");
}

float ann_get_rmse(struct ANN *ann, int num_points, float *inputs, 
    float *target_outputs)
{
    int input_size = ann->layers[0].num_neurons;
    int output_size = ann->layers[ann->num_layers-1].num_neurons;
    float ret = 0.f;
    struct ANNLayer* layer = ann->layers + ann->num_layers-1;
    for(int i=0;i<num_points;i++){
        ann_run(ann,inputs+i*input_size,0);
        for(int j=0;j<output_size;j++){
            float a = target_outputs[i*output_size + j]-
                layer->neurons[j];
            ret += a*a;
        }
    }
    ret /= (float)num_points;
    return sqrtf(ret);
}




        /*
        //Shuffle the patterns
        for(int i=num_patterns-1; i > 0;i--){
            int j = rand()%(i+1);
            {
                int len = input_size*sizeof(float);
                float * a = inputs + i*input_size;
                float * b = inputs + j*input_size;
                memcpy(tmp_input, a, len);
                memcpy(a, b, len);
                memcpy(b, tmp_input, len);
            }
            {
                int len = output_size*sizeof(float);
                float * a = target_outputs + i*output_size;
                float * b = target_outputs + j*output_size;
                memcpy(tmp_output, a, len);
                memcpy(a, b, len);
                memcpy(b, tmp_output, len);
            }
        }
        */




            //Check gradient
            /*
            if(0){
                printf("Numerical gradient:\n");
                printf("input:\n");
                for(int i=0;i<input_size;i++){
                    printf("%f, ", input[i]);
                }
                printf("\noutput:\n");
                for(int i=0;i<output_size;i++){
                    printf("%f, ", target_outputs[pattern*output_size + i]);
                }
                printf("\n");
                struct ANN *tmp_ann = ann_clone(ann);
                float epsilon = 1e-5f;
                struct ANNLayer *last_layer = tmp_ann->layers + end_layer;
                for(int i_layer = 1; i_layer < num_layers; i_layer++){
                    struct ANNLayer *layer = tmp_ann->layers + i_layer;
                    int num_weights = layer->num_weights;
                    printf("Layer %d:\n",i_layer);

                    int num_neurons = layer->num_neurons;
                    int num_left_neurons = layer->num_left_neurons;
                    for(int i=0;i<num_neurons;i++){
                        for(int j=0;j<num_left_neurons;j++){
                            float v_plus = 0.f;
                            float v_minus = 0.f;
                            int a= i*num_left_neurons + j;
                            float tmp = layer->weights[a];
                            layer->weights[a] = tmp + epsilon;
                            ann_run(tmp_ann,input,0);
                            for(int k=0;k<output_size;k++){
                                float f = target_outputs[pattern*output_size + k] - last_layer->neurons[k];
                                v_plus += f*f;
                            }
                            layer->weights[a] = tmp - epsilon;
                            ann_run(tmp_ann,input,0);
                            for(int k=0;k<output_size;k++){
                                float f = target_outputs[pattern*output_size + k] - last_layer->neurons[k];
                                v_minus += f*f;
                            }
                            layer->weights[a] = tmp;
                            printf("%f (%f), ",
                                0.5f*(v_plus - v_minus)/epsilon*0.5f,
                                -ann->layers[i_layer].weight_updates[a]
                            );
                        }
                        printf("\n");
                    }
                }
                ann_free(tmp_ann);
            }
            */
