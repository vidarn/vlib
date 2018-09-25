#include "rf.h"
#include <stdlib.h>
#include <stdio.h>
#ifdef _MSC_VER
#include <malloc.h>
#endif
#include "thread/thread.h"

struct RandomForestThreadParameters 
{
    int *trees_built;
    int input_dim, target_dim;
    void (*get_pattern)(float *input, float *target, int i, void *data);
    void *get_pattern_data;
    struct CART *trees;
    struct RandomForestSettings settings;
    int num_patterns, num_trees;
    struct RandomForestConstructionProgress *progress;
    struct RNGState *rng;
    int **used_samples;
};

static unsigned long rf_construct_thread(void *data)
{
    struct RandomForestThreadParameters p =
        *(struct RandomForestThreadParameters*)data;

    float *inputs  = calloc(p.settings.training_size, p.input_dim*sizeof(float));
    float *targets = calloc(p.settings.training_size,p.target_dim*sizeof(float));
    float *tmp_input  = alloca(p.input_dim*sizeof(float));
    float *tmp_target = alloca(p.target_dim*sizeof(float));
    int *index_shuffle = calloc(p.num_patterns, sizeof(int));
    for(int i=0;i<p.num_patterns;i++){
        index_shuffle[i] = i;
    }

    for(int i_tree=0;i_tree<p.num_trees;i_tree++){
        p.used_samples[i_tree] = calloc(p.settings.training_size, sizeof(int));
        rng_shuffle((uint8_t*)index_shuffle, p.num_patterns, sizeof(int), p.rng);
        for(int i_pattern=0;i_pattern<p.settings.training_size;i_pattern++){
            int p_index = index_shuffle[i_pattern];
            p.used_samples[i_tree][i_pattern] = p_index;
            p.get_pattern(tmp_input, tmp_target, p_index, p.get_pattern_data);
            for(int i_inp =0;i_inp<p.input_dim;i_inp++){
                inputs[i_pattern + i_inp*p.settings.training_size]
                    = tmp_input[i_inp];
            }
            for(int i_tar =0;i_tar<p.target_dim;i_tar++){
                targets[i_pattern + i_tar*p.settings.training_size]
                    = tmp_target[i_tar];
            }
        }
        p.trees[i_tree] = cart_construct(p.settings.training_size, inputs,
            p.input_dim, targets, p.target_dim, p.settings.cart_settings,
            &p.progress->cart_progress, p.rng);
#ifdef _MSC_VER
#else
        int num_trees_built = __sync_fetch_and_add(p.trees_built, 1);
        if(num_trees_built%20 == 19){
            printf("%d trees built\n", num_trees_built+1);
        }
#endif
        if(p.progress){
            p.progress->num_trees_constructed++;
            p.progress->cart_progress.max_depth_reached = 0;
            p.progress->cart_progress.num_nodes_created = 0;
        }
    }
    free(index_shuffle);

    return 0;
}

struct RandomForest rf_construct_parallel(int num_patterns,
    void (*get_pattern)(float *input, float *target, int i, void *data),
    void *get_pattern_data, int input_dim, int target_dim,
    struct RandomForestSettings settings, int num_threads,
    struct RandomForestConstructionProgress *progress, struct RNGState **rngs)
{
    int num_trees_built = 0;
    struct RandomForest rf = {0};
    rf.num_trees = settings.num_trees;
    rf.trees = calloc(settings.num_trees, sizeof(struct CART));
    rf.input_dim  =  input_dim;
    rf.output_dim = target_dim;
    rf.num_samples_per_tree = settings.training_size;
    rf.trees_used_samples = calloc(settings.num_trees, sizeof(int**));
    struct RandomForestThreadParameters *thread_params = 
        calloc(num_threads, sizeof(struct RandomForestThreadParameters));
    int trees_per_thread = settings.num_trees/num_threads;
    int extra_trees = settings.num_trees - trees_per_thread*num_threads;
    struct ThreadHandle **threads = calloc(num_threads,
        sizeof(struct ThreadHandle*));
    int num_trees_accum = 0;
    for(int i_thread=0;i_thread < num_threads;i_thread++){
        thread_params[i_thread].trees_built  = &num_trees_built;
        thread_params[i_thread].input_dim  = input_dim;
        thread_params[i_thread].target_dim = target_dim;
        thread_params[i_thread].get_pattern  = get_pattern;
        thread_params[i_thread].get_pattern_data  = get_pattern_data;
        thread_params[i_thread].trees = rf.trees + num_trees_accum;
        thread_params[i_thread].settings = settings;
        thread_params[i_thread].num_patterns = num_patterns;
        thread_params[i_thread].num_trees = trees_per_thread
            + (i_thread < extra_trees ? 1 : 0);
        thread_params[i_thread].progress = progress + i_thread;
        thread_params[i_thread].rng = rngs[i_thread];
        thread_params[i_thread].used_samples = rf.trees_used_samples + num_trees_accum;
        num_trees_accum += thread_params[i_thread].num_trees;
        threads[i_thread] =
            thread_start(rf_construct_thread, thread_params+i_thread);
    }
    for(int i_thread=0;i_thread < num_threads;i_thread++){
        thread_wait(threads[i_thread]);
    }
    free(thread_params);
    free(threads);
    return rf;
}


struct RandomForest rf_construct(int num_patterns, const float *inputs_in,
    int input_dim, const float *targets_in, int target_dim,
    struct RandomForestSettings settings, struct RNGState *rng)
{
    struct RandomForest rf = {0};
    rf.num_trees = settings.num_trees;
    rf.trees = calloc(settings.num_trees, sizeof(struct CART));
    rf.input_dim  =  input_dim;
    rf.output_dim = target_dim;
    float *inputs  = calloc(settings.training_size, input_dim*sizeof(float));
    float *targets = calloc(settings.training_size,target_dim*sizeof(float));
    for(int i_tree=0;i_tree<settings.num_trees;i_tree++){
        for(int i_pattern=0;i_pattern<settings.training_size;i_pattern++){
            int r = rng_below(rng, num_patterns);
            for(int i_dim=0;i_dim<input_dim;i_dim++){
                inputs[i_pattern*input_dim + i_dim]
                    = inputs_in[r*input_dim + i_dim];
            }
            for(int i_dim=0;i_dim<target_dim;i_dim++){
                targets[i_pattern*target_dim + i_dim]
                    = targets_in[r*target_dim + i_dim];
            }
        }
        rf.trees[i_tree] = cart_construct(settings.training_size, inputs,
            input_dim, targets, target_dim, settings.cart_settings, 0, rng);
    }
    return rf;
}

void rf_free(struct RandomForest rf)
{
    for(int i_tree=0;i_tree<rf.num_trees;i_tree++){
        cart_free(rf.trees[i_tree]);
    }
}

void rf_run(struct RandomForest rf, float *input, float *output,
    void (*evaluate_func)(float *input, int input_dim, float *output,
        int output_dim, float* parameters)
    )
{
    float *tmp_out = alloca(rf.output_dim*sizeof(float));
    for(int i_dim=0;i_dim < rf.output_dim;i_dim++){
        output[i_dim] = 0.f;
    }
    for(int i_tree=0;i_tree<rf.num_trees;i_tree++){
        cart_run(rf.trees[i_tree],input,tmp_out, evaluate_func);
        for(int i_dim=0;i_dim < rf.output_dim;i_dim++){
            output[i_dim] += tmp_out[i_dim];
        }
    }
    for(int i_dim=0;i_dim < rf.output_dim;i_dim++){
        output[i_dim] *= 1.f/(float)rf.num_trees;
    }
}

void rf_run_default_evaluate(struct RandomForest rf, float *input, float *output)
{
    float *tmp_out = alloca(rf.output_dim*sizeof(float));
    for(int i_dim=0;i_dim < rf.output_dim;i_dim++){
        output[i_dim] = 0.f;
    }
    for(int i_tree=0;i_tree<rf.num_trees;i_tree++){
        cart_run_default_evaluate(rf.trees[i_tree],input,tmp_out);
        for(int i_dim=0;i_dim < rf.output_dim;i_dim++){
            output[i_dim] += tmp_out[i_dim];
        }
    }
    for(int i_dim=0;i_dim < rf.output_dim;i_dim++){
        output[i_dim] *= 1.f/(float)rf.num_trees;
    }
}


size_t rf_to_mem_size(struct RandomForest rf)
{
    size_t s = 0;
    s += sizeof(int);
    s += sizeof(int);
    s += sizeof(int);
    s += sizeof(int);
    for(int i_tree=0;i_tree<rf.num_trees;i_tree++){
        struct CART tree = rf.trees[i_tree];
        s += sizeof(int);
        for(int i_node=0;i_node<tree.num_nodes;i_node++){
            s += sizeof(struct CARTNode);
        }
    }
    s += sizeof(int);
    s += rf.num_trees*rf.num_samples_per_tree*sizeof(int);
    return s;
}

void rf_to_mem(struct RandomForest rf, unsigned char *mem)
{
    *(int*)mem = 2; // Version
    mem += sizeof(int);
    *(int*)mem = rf.num_trees;
    mem += sizeof(int);
    *(int*)mem = rf.input_dim;
    mem += sizeof(int);
    *(int*)mem = rf.output_dim;
    mem += sizeof(int);
    for(int i_tree=0;i_tree<rf.num_trees;i_tree++){
        struct CART tree = rf.trees[i_tree];
        *(int*)mem = tree.num_nodes;
        mem += sizeof(int);
        for(int i_node=0;i_node<tree.num_nodes;i_node++){
            *(struct CARTNode*)mem = tree.nodes[i_node];
            mem += sizeof(struct CARTNode);
        }
    }
    *(int*)mem = rf.num_samples_per_tree;
    mem += sizeof(int);
    for(int i_tree=0;i_tree<rf.num_trees;i_tree++){
        size_t len = rf.num_samples_per_tree*sizeof(int);
        memcpy(mem,rf.trees_used_samples[i_tree],len);
        mem += len;
    }
}

struct RandomForest rf_from_mem(unsigned char *mem, int max_num_trees)
{
    int version = *(int*)mem; // Version
    mem += sizeof(int);

    struct RandomForest rf;

    rf.num_trees = *(int*)mem;
    mem += sizeof(int);
    rf.input_dim = *(int*)mem;
    mem += sizeof(int);
    rf.output_dim = *(int*)mem;
    mem += sizeof(int);
    
    if(rf.num_trees > max_num_trees && max_num_trees > 0){
        rf.num_trees = max_num_trees;
    }

    rf.trees = calloc(rf.num_trees,sizeof(struct CART));

    for(int i_tree=0;i_tree<rf.num_trees;i_tree++){
        struct CART *tree = rf.trees+i_tree;
        tree->num_nodes = *(int*)mem;
        mem += sizeof(int);

        tree->nodes = calloc(tree->num_nodes,sizeof(struct CARTNode));
        tree->input_dim = rf.input_dim;
        tree->output_dim = rf.output_dim;

        for(int i_node=0;i_node<tree->num_nodes;i_node++){
            tree->nodes[i_node] = *(struct CARTNode*)mem;
            mem += sizeof(struct CARTNode);
        }
    }
    if(version > 1){
        rf.num_samples_per_tree = *(int*)mem;
        mem += sizeof(int);
        rf.trees_used_samples = calloc(rf.num_trees,sizeof(int*));
        for(int i_tree=0;i_tree<rf.num_trees;i_tree++){
            size_t len = rf.num_samples_per_tree*sizeof(int);
            rf.trees_used_samples[i_tree] = calloc(1,len);
            memcpy(mem,rf.trees_used_samples[i_tree],len);
            mem += len;
        }
    }else{
        rf.num_samples_per_tree = 0;
        rf.trees_used_samples = 0;
    }
    return rf;
}

void rf_oob_error(struct RandomForest rf, int num_patterns,
    void (*get_pattern)(float *input, float *target, int i, void *data),
    void *get_pattern_data, struct RNGState *rng, float *errors,
    int num_eval)
{
    float *tmp_input  = alloca(rf.input_dim*sizeof(float));
    float *tmp_target = alloca(rf.output_dim*sizeof(float));
    float *tmp_output = alloca(rf.output_dim*sizeof(float));
    int **used_samples = alloca(rf.num_trees*sizeof(int*));
    for(int i_tree=0;i_tree<rf.num_trees;i_tree++){
        used_samples[i_tree] = calloc(num_patterns, sizeof(int));
        for(int i=0;i<rf.num_samples_per_tree;i++){
            int i_pattern = rf.trees_used_samples[i_tree][i];
            used_samples[i_tree][i_pattern] = 1;
        }
    }
    for(int i_pattern =0;i_pattern<num_patterns;i_pattern++){
        get_pattern(tmp_input, tmp_target, i_pattern, get_pattern_data);
        double err=0.f;
        int num_valid = 0;
        for(int i_tree=0;i_tree<rf.num_trees;i_tree++){
            int valid = !used_samples[i_tree][i_pattern];
            if(valid){
                cart_run_default_evaluate(rf.trees[i_tree],tmp_input,tmp_output);
                for(int d=0;d<rf.output_dim;d++){
                    double delta = tmp_output[d] - tmp_target[d];
                    err += delta*delta;
                }
                num_valid++;
            }
        }
        if(num_valid > 0){
            errors[i_pattern] = (float)(err/(double)(rf.output_dim*(double)num_valid));
        }
    }
    for(int i_tree=0;i_tree<rf.num_trees;i_tree++){
        free(used_samples[i_tree]);
    }
}


#include "mkl.h"

struct RandomForestForEvaluation rf_convert_for_evaluation(
    struct RandomForest rf)
{
    struct RandomForestForEvaluation ret = {0};
    ret.num_trees = rf.num_trees;
    ret.input_dim = rf.input_dim;
    ret.output_dim = rf.output_dim;
    ret.trees = mkl_calloc(ret.num_trees,sizeof(struct CARTForEvaluation), 64);
    for(int i=0;i<ret.num_trees;i++){
        ret.trees[i] = cart_convert_for_evaluation(rf.trees[i]);
    }
    return ret;
}

void rf_evaluate(struct RandomForestForEvaluation rf, uint16_t *input,
    float *output)
{
    float *tmp_out = alloca(rf.output_dim*sizeof(float));
    for(int i_dim=0;i_dim < rf.output_dim;i_dim++){
        output[i_dim] = 0.f;
    }
    for(int i_tree=0;i_tree<rf.num_trees;i_tree++){
        cart_evaluate(rf.trees[i_tree],input,tmp_out);
        for(int i_dim=0;i_dim < rf.output_dim;i_dim++){
            output[i_dim] += tmp_out[i_dim];
        }
    }
    for(int i_dim=0;i_dim < rf.output_dim;i_dim++){
        output[i_dim] *= 1.f/(float)rf.num_trees;
    }
}

void rf_evaluate_subset(struct RandomForestForEvaluation rf, uint16_t *input,
    float *output, int num, struct RNGState *rng)
{
    if(num > rf.num_trees) num = rf.num_trees;
    int *visited = alloca(num*sizeof(int));
    float *tmp_out = alloca(rf.output_dim*sizeof(float));
    for(int i_dim=0;i_dim < rf.output_dim;i_dim++){
        output[i_dim] = 0.f;
    }
    for(int i_tree=0;i_tree<num;i_tree++){
        int unique_found = 0;
        int tree_index;
        while(!unique_found){
            tree_index = rng_below(rng, rf.num_trees);
            unique_found = 1;
            for(int i=0;i<i_tree;i++){
                if(visited[i] == tree_index){
                    unique_found = 0;
                    break;
                }
            }
        }
        cart_evaluate(rf.trees[tree_index],input,tmp_out);
        for(int i_dim=0;i_dim < rf.output_dim;i_dim++){
            output[i_dim] += tmp_out[i_dim];
        }
    }
    float f = 1.f/(float)num;
    for(int i_dim=0;i_dim < rf.output_dim;i_dim++){
        output[i_dim] *= f;
    }
}

size_t rf_eval_to_mem_size(struct RandomForestForEvaluation rf)
{
    size_t s = 0;
    s += sizeof(int);
    s += sizeof(int);
    s += sizeof(int);
    s += sizeof(int);
    for(int i_tree=0;i_tree<rf.num_trees;i_tree++){
        struct CARTForEvaluation tree = rf.trees[i_tree];
        s += sizeof(int);
        s += sizeof(int);
        s += tree.num_subtrees*sizeof(struct CARTNode);
        s += tree.num_leaf_values*sizeof(float)*rf.output_dim;
    }
    return s;
}

void rf_eval_to_mem(struct RandomForestForEvaluation rf, unsigned char *mem)
{
    *(int*)mem = 1; // Version
    mem += sizeof(int);
    *(int*)mem = rf.num_trees;
    mem += sizeof(int);
    *(int*)mem = rf.input_dim;
    mem += sizeof(int);
    *(int*)mem = rf.output_dim;
    mem += sizeof(int);
    for(int i_tree=0;i_tree<rf.num_trees;i_tree++){
        struct CARTForEvaluation tree = rf.trees[i_tree];
        *(int*)mem = tree.num_subtrees;
        mem += sizeof(int);
        *(int*)mem = tree.num_leaf_values;
        mem += sizeof(int);
        for(int i_subtree=0;i_subtree<tree.num_subtrees;i_subtree++){
            *(struct CARTSubtreeForEvaluation*)mem = tree.subtrees[i_subtree];
            mem += sizeof(struct CARTSubtreeForEvaluation);
        }
        for(int i_leaf=0;i_leaf<tree.num_leaf_values*rf.output_dim;i_leaf++){
            *(float*)mem = tree.leaf_values[i_leaf];
            mem += sizeof(float);
        }
    }
}

struct RandomForestForEvaluation rf_eval_from_mem(unsigned char *mem,
    int max_num_trees)
{
    int version = *(int*)mem; // Version
    mem += sizeof(int);

    struct RandomForestForEvaluation rf;

    rf.num_trees = *(int*)mem;
    mem += sizeof(int);
    rf.input_dim = *(int*)mem;
    mem += sizeof(int);
    rf.output_dim = *(int*)mem;
    mem += sizeof(int);
    
    if(rf.num_trees > max_num_trees){
        rf.num_trees = max_num_trees;
    }

    rf.trees = calloc(rf.num_trees,sizeof(struct CARTForEvaluation));

    for(int i_tree=0;i_tree<rf.num_trees;i_tree++){
        struct CARTForEvaluation *tree = rf.trees+i_tree;
        tree->num_subtrees = *(int*)mem;
        mem += sizeof(int);
        tree->num_leaf_values = *(int*)mem;
        mem += sizeof(int);

        tree->subtrees = calloc(tree->num_subtrees,sizeof(struct CARTSubtreeForEvaluation));
        tree->leaf_values = calloc(tree->num_leaf_values,sizeof(float)*rf.output_dim);
        tree->leaf_value_dim = rf.output_dim;

        for(int i_subtree=0;i_subtree<tree->num_subtrees;i_subtree++){
            tree->subtrees[i_subtree] = *(struct CARTSubtreeForEvaluation*)mem;
            mem += sizeof(struct CARTSubtreeForEvaluation);
        }
        for(int i_leaf=0;i_leaf<tree->num_leaf_values*rf.output_dim;i_leaf++){
            tree->leaf_values[i_leaf] = *(float*)mem;
            mem += sizeof(float);
        }
    }
    return rf;
}
