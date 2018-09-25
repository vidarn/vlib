#include "cart.h"
#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <string.h>
#ifdef _WIN32
#include <malloc.h>
#endif

#include "sort/sort.h"

struct CARTConstructionNode{
    struct CARTNode node;
    int pattern_offset;
    int num_patterns;
    int next_in_queue;
    int depth;
    int index;
};

int comp_float(const void* _a, const void* _b)
{
    float a = *(float*)_a;
    float b = *(float*)_b;
    if(a < b) return -1;
    if(a > b) return 1;
    return 0;
}

float cart_anova_impurity(float *inp, int input_dim, float *trg, int target_dim,
   int stride, int num_patterns, float *parameters)
{
    //NOTE(Vidar):Handle the sum of a large number of objects by summing 2^6=64
    // at a time. Do this in 6 levels, which means that we can sum 2^36 items,
    // more than can be enumerated by the 32-bit num_patterns parameter
    float *accumulation     = alloca(6*sizeof(float));
    int   *accumulation_num = alloca(6*sizeof(int));

    float impurity = 0.f;
    for(int j_dim = 0;j_dim < target_dim; j_dim++){
        memset(accumulation,    0,6*sizeof(float));
        memset(accumulation_num,0,6*sizeof(int));
        for(int i_pattern = 0; i_pattern < num_patterns; i_pattern++) {
            accumulation[0] += trg[i_pattern +  j_dim*stride];
            accumulation_num[0]++;
            int accum_level = 0;
            while(accumulation_num[accum_level] == (1<<6)){
                accumulation_num[accum_level+1]++;
                accumulation_num[accum_level]=0;
                accumulation[accum_level+1] += accumulation[accum_level];
                accumulation[accum_level] = 0.f;
                accum_level++;
            }
        }
        parameters[j_dim] = 0.f;
        for(int i_accum = 0; i_accum<6; i_accum++){
            parameters[j_dim] += accumulation[i_accum];
        }
        parameters[j_dim] /= (float)num_patterns;
    }
    for(int j_dim = 0;j_dim < target_dim; j_dim++){
        memset(accumulation,    0,6*sizeof(float));
        memset(accumulation_num,0,6*sizeof(int));
        for(int i_pattern = 0; i_pattern < num_patterns; i_pattern++) {
            double d = trg[i_pattern + j_dim*stride] - parameters[j_dim];
            accumulation[0] += d*d;
            accumulation_num[0]++;
            int accum_level = 0;
            while(accumulation_num[accum_level] == (1<<6)){
                accumulation_num[accum_level+1]++;
                accumulation[accum_level+1] += accumulation[accum_level];
                accumulation[accum_level] = 0.f;
                accum_level++;
            }
        }
        for(int i_accum = 0; i_accum<6; i_accum++){
            impurity += accumulation[i_accum];
        }
    }
    return impurity;
}

struct CART cart_construct(int num_patterns, const float *inputs_in, int input_dim,
    const float *targets_in, int target_dim, struct CARTSettings settings,
    struct CARTConstructionProgress *progress,
    struct RNGState *rng)
{
    int num_features_per_node = settings.num_features_per_node;
    if(num_features_per_node > input_dim){
        num_features_per_node = input_dim;
    }
    int *used_dims = 0;
    if(num_features_per_node > 0){
        used_dims = calloc(input_dim, sizeof(int));
    }
    float *inputs_buffer  = calloc(num_patterns, input_dim *sizeof(float));
    float *targets_buffer = calloc(num_patterns, target_dim*sizeof(float));
    memcpy(inputs_buffer, inputs_in, num_patterns* input_dim*sizeof(float));
    memcpy(targets_buffer,targets_in,num_patterns*target_dim*sizeof(float));
    int   *index_buffer = calloc(num_patterns, sizeof(int));

    float *sort_value_storage = calloc(num_patterns, sizeof(float));
    int   *sort_index_storage = calloc(num_patterns, sizeof(int));
    
    int created_nodes_size = 0;
    if(progress->created_nodes){
        created_nodes_size = 1 << (progress->created_nodes_max_depth+1);
        memset(progress->created_nodes,0,created_nodes_size);
    }

    int alloc_nodes = 64;
    struct CARTConstructionNode *nodes =
        calloc(alloc_nodes,sizeof(struct CARTConstructionNode));
    int num_nodes = 1;
    nodes[0].num_patterns = num_patterns;
    nodes[0].pattern_offset = 0;
    nodes[0].next_in_queue = -1;
    nodes[0].depth = 0;
    nodes[0].index = 0;
    int current_node = 0;
    while(current_node >= 0){
        struct CARTConstructionNode *node = nodes + current_node;

        if(progress){
            progress->current_depth = node->depth;
        }

        int best_dim = 0;
        float best_split = 0.f;
        float best_impurity = FLT_MAX;
        float best_impurity_left = FLT_MAX;
        float best_impurity_right = FLT_MAX;
        float best_mean_left [CART_MAX_PARAM_DIM] = {0};
        float best_mean_right[CART_MAX_PARAM_DIM] = {0};
        int best_num_left  = 0;
        int best_num_right = 0;
        float *inp = inputs_buffer  + node->pattern_offset;
        float *trg = targets_buffer + node->pattern_offset;
        int n = num_features_per_node;
        if(n<=0){
            n = input_dim;
        }else{
            memset(used_dims,0,input_dim*sizeof(int));
        }
        for(int i = 0;i < n; i++)
        {
            int i_dim = i;
            if(num_features_per_node > 0){
                do{
                    i_dim = rng_below(rng, input_dim);
                }while(used_dims[i_dim]);
                used_dims[i_dim] = 1;
            }

            //Sort according to current dimension
            {
                for(int i_index=0;i_index<num_patterns;i_index++){
                    index_buffer[i_index] = i_index;
                }
                merge_sort_auxillary(inp + i_dim*num_patterns, sort_value_storage,
                    index_buffer, sort_index_storage, node->num_patterns);
                for(int j_dim = 0; j_dim < input_dim; j_dim++){
                    if(j_dim != i_dim){
                        float *tmp = sort_value_storage;
                        memcpy(tmp, inp+j_dim*num_patterns, node->num_patterns*sizeof(float));
                        for(int i_index=0;i_index<node->num_patterns;i_index++){
                            inp[i_index + j_dim*num_patterns] =
                                tmp[index_buffer[i_index]];
                        }
                    }
                }
                for(int j_dim = 0; j_dim < target_dim; j_dim++){
                    float *tmp = sort_value_storage;
                    memcpy(tmp, trg+j_dim*num_patterns, node->num_patterns*sizeof(float));
                    for(int i_index=0;i_index<node->num_patterns;i_index++){
                        trg[i_index + j_dim*num_patterns] =
                            tmp[index_buffer[i_index]];
                    }
                }
            }

            //printf("Considering dim %d\n", i_dim);
            
            int num_splits = node->num_patterns;
            if(settings.random_splits){
                //NOTE(Vidar):"Extremely randomized trees"
                num_splits = 8;
            }

            for(int i_split = 0; i_split < num_splits; i_split++){
                float split;
                int i_lower;
                if(settings.random_splits){
                    int split_pattern= rng_below(rng, node->num_patterns);
                    /*
                    printf("Random split at pattern %d of %d (%f)\n",
                        split_pattern,node->num_patterns,
                        (float)split_pattern/(float)node->num_patterns);
                    */
                    split = inp[split_pattern + i_dim*num_patterns];
                    i_lower = i_split-1;
                    while(i_lower < node->num_patterns-1 && inp[i_lower+1 + i_dim*num_patterns] < split){
                        i_lower++;
                    }
                }else{
                    split = inp[i_split + i_dim*num_patterns];
                    i_lower = i_split-1;
                    while(i_lower >= 0 && inp[i_lower + i_dim*num_patterns] == split){
                        i_lower--;
                    }
                }
                int num_left  = i_lower+1;
                int num_right = node->num_patterns - num_left;
                if(num_left > 0 && num_right > 0){
                    split = (split + inp[i_lower + i_dim*num_patterns])*0.5f;
                    //printf("Split at %f in dim %d:\nnum_left %d num_right %d\n",split,i_dim,num_left,num_right);
                    float impurity_left  = 0.f;
                    float impurity_right = 0.f;
                    float mean_left [CART_MAX_PARAM_DIM] = {0};
                    float mean_right[CART_MAX_PARAM_DIM] = {0};
                    impurity_left  = settings.impurity_func(inp,
                         input_dim,trg,target_dim, num_patterns,
                        num_left, mean_left);
                    impurity_right = settings.impurity_func(inp + num_left,
                        input_dim, trg + num_left, target_dim, num_patterns,
                        num_right, mean_right);
                    float impurity = impurity_left + impurity_right;
                    //printf("impurity %f\n",impurity);
                    if(impurity < best_impurity){
                        best_impurity = impurity;
                        best_impurity_left  = impurity_left;
                        best_impurity_right = impurity_right;
                        best_dim = i_dim;
                        best_split = split;
                        for(int j_dim = 0;j_dim < CART_MAX_PARAM_DIM; j_dim++){
                            best_mean_left [j_dim]  = mean_left [j_dim];
                            best_mean_right[j_dim] = mean_right[j_dim];
                        }
                        best_num_left  = num_left;
                        best_num_right = num_right;
                        //printf("    Currently best split is at %f in dim %d, impurity %f (%f,%f), num left: %d, num right: %d\n", best_split, best_dim, best_impurity, best_impurity_left, best_impurity_right, best_num_left, best_num_right);
                    }
                    if(progress){
                        progress->num_splits++;
                    }
                }
            }
        }
        //printf("Best split is at %f in dim %d\n", best_split, best_dim);
        int left_node_index = num_nodes;
        int right_node_index = num_nodes+1;
        node->node.predictor = best_dim;
        node->node.split = best_split;
        node->node.left_node_index = num_nodes;
        if(num_nodes+2 >= alloc_nodes){
            alloc_nodes *= 2;
            nodes = realloc(nodes,alloc_nodes*sizeof(struct CARTConstructionNode));
            node = nodes + current_node;
        }
        struct CARTConstructionNode *left_node  = nodes + left_node_index;
        struct CARTConstructionNode *right_node = nodes + right_node_index;
        num_nodes+=2;

        //Sort patterns according to split
        {
            /*printf("node->num_patterns = %d, split = %f, dim = %d\n", node->num_patterns, best_split, best_dim);
            for(int i_dim =0;i_dim<input_dim;i_dim++){
                for(int i = 0;i<node->num_patterns;i++){
                    printf("%f, ", inp[i*input_dim + i_dim]);
                }
                printf("\n");
            }
            */
            int li = 0;
            int ri = node->num_patterns-1;
            //TODO(Vidar): Use merge sort instead??
            while(li < ri){
                while(li < node->num_patterns && inp[li + best_dim*num_patterns] < best_split)
                {
                    li++;
                }
                while(ri >= 0 && inp[ri + best_dim*num_patterns] >= best_split)
                {
                    ri--;
                }
                if(li < ri){
                    //printf("li: %d, ri: %d\n", li, ri);
                    for(int i_dim =0;i_dim<input_dim;i_dim++){
                        float tmp = inp[ri + i_dim*num_patterns];
                        inp[ri + i_dim*num_patterns] = inp[li + i_dim*num_patterns];
                        inp[li + i_dim*num_patterns] = tmp;
                        /*
                        for(int i = 0;i<node->num_patterns;i++){
                            printf("%f, ", inp[i*input_dim + i_dim]);
                        }
                        printf("\n");
                        */
                    }
                    for(int i_dim =0;i_dim<target_dim;i_dim++){
                        float tmp = trg[ri + i_dim*num_patterns];
                        trg[ri + i_dim*num_patterns] = trg[li + i_dim*num_patterns];
                        trg[li + i_dim*num_patterns] = tmp;
                    }
                }
            }
        }

        for(int j_dim = 0;j_dim < CART_MAX_PARAM_DIM; j_dim++){
            left_node->node.val [j_dim] = best_mean_left [j_dim];
            right_node->node.val[j_dim] = best_mean_right[j_dim];
        }
        
        /*
        
                0
           1          2             *2 +1
        3   4     5     6           *2 +1
       7 8 9 10 11 12 13 14         *2 +1
         */
        
        if(node->index < created_nodes_size){
            progress->created_nodes[node->index] = 1; //TODO(Vidar):Use bits instead of entire bytes?
        }
        
        int depth = node->depth;

        left_node->pattern_offset = node->pattern_offset;
        left_node->num_patterns = best_num_left;
        left_node->depth = depth+1;
        left_node->index = node->index*2 + 1;

        right_node->pattern_offset = node->pattern_offset + best_num_left;
        right_node->num_patterns = best_num_right;
        right_node->depth = depth+1;
        right_node->index = node->index*2 + 2;

        int max_depth = settings.max_depth;
        if(max_depth == 0) max_depth = 9999;


        if(best_impurity_right > 1e-8f
            && best_num_right > settings.min_interior_size
            && depth < max_depth-1)
        {
            right_node->next_in_queue = node->next_in_queue;
            node->next_in_queue = right_node_index;
        }else{
            right_node->node.left_node_index = 0;
            right_node->next_in_queue = -1;
        }

        if(best_impurity_left > 1e-8f
            && best_num_left > settings.min_interior_size
            && depth < max_depth-1)
        {
            left_node->next_in_queue = node->next_in_queue;
            node->next_in_queue = left_node_index;
        }else{
            left_node->node.left_node_index = 0;
            left_node->next_in_queue = -1;
        }
        current_node = node->next_in_queue;
        
        if(progress){
            progress->num_nodes_created = num_nodes;
            if(depth > progress->max_depth_reached){
                progress->max_depth_reached = depth;
            }
            progress->num_splits = 0;
        }
    }

    free(sort_value_storage);
    free(sort_index_storage);
    free(index_buffer);
    free(inputs_buffer);
    free(targets_buffer);
    free(used_dims);

    struct CART ret = {0};
    ret.input_dim  =  input_dim;
    ret.output_dim = target_dim;
    ret.num_nodes = num_nodes;
    ret.nodes = calloc(num_nodes,sizeof(struct CARTNode));
    for(int i_node =0; i_node < num_nodes; i_node++){
        ret.nodes[i_node] = nodes[i_node].node;
    }
    free(nodes);
    //cart_print(ret);
    return ret;
}

void cart_free(struct CART cart)
{
    free(cart.nodes);
}

void cart_default_evaluate(float *input, int input_dim, float *output,
    int output_dim, float* parameters)
{
    for(int i_dim=0;i_dim<output_dim;i_dim++){
        output[i_dim] = parameters[i_dim];
    }
}

void cart_run(struct CART cart, float *input, float *output,
    void (*evaluate_func)(float *input, int input_dim, float *output,
        int output_dim, float* parameters)
    )
{
    //printf("input: %f\n", input[0]);
    int input_dim  = cart.input_dim;
    int output_dim = cart.output_dim;
    int done = 0;
    struct CARTNode *node = cart.nodes;
    while(!done){
        if(node->left_node_index == 0){
            evaluate_func(input, input_dim,  output, output_dim, node->val);
            /*
            for(int i_dim=0;i_dim<output_dim;i_dim++){
                output[i_dim] = node->val[i_dim];
            }
            */
            //printf("output: %f\n", output[0]);
            done = 1;
        }
        else{
            if(input[node->predictor] < node->split){
                //printf("%f < %f\n",input[node->predictor] , node->split);
                node = cart.nodes + node->left_node_index;
            }else{
                //printf("%f >= %f\n",input[node->predictor] , node->split);
                node = cart.nodes + node->left_node_index+1;
            }
        }
    }
}

void cart_run_default_evaluate(struct CART cart, float *input, float *output)
{
    int output_dim = cart.output_dim;
    int done = 0;
    struct CARTNode *node = cart.nodes;
    while(!done){
        if(node->left_node_index == 0){
            for(int i_dim=0;i_dim<output_dim;i_dim++){
                output[i_dim] = node->val[i_dim];
            }
            done = 1;
        }
        else{
            if(input[node->predictor] < node->split){
                node = cart.nodes + node->left_node_index;
            }else{
                node = cart.nodes + node->left_node_index+1;
            }
        }
    }
}

void cart_print(struct CART cart)
{
    printf("Printing tree:\n");
    int num_nodes = cart.num_nodes;
    for(int i_node =0; i_node < num_nodes; i_node++){
        struct CARTNode *node = cart.nodes + i_node;
        if(node->left_node_index == 0){
            printf("[%d] Params are {", i_node);
            for(int i_dim =0; i_dim < CART_MAX_PARAM_DIM; i_dim++){
                if(i_dim > 0){
                    printf(", ");
                }
                printf("%f",node->val[i_dim]);
            }
            printf("}\n");
        }else{
            printf("[%d] If feature %d < %f, go to node %d, otherwise go to %d\n",
                i_node, node->predictor, node->split, node->left_node_index,
                node->left_node_index+1);
        }
    }
}

#include "mkl.h"
#include <assert.h>
#define DIM_BITS 5

//BOOKMARK(Vidar):Finish this implementation, then fix the construction of the
// trees, make sure that they can handle a huge number of patterns
struct CARTForEvaluation cart_convert_for_evaluation(struct CART cart)
{
    static const int traversal_order[15*3] = {
        -1,-1,-1,
         0,-1,-1,
         1,-1,-1,
         0, 0,-1,
         0, 1,-1,
         1, 0,-1,
         1, 1,-1,
         0, 0, 0,
         0, 0, 1,
         0, 1, 0,
         0, 1, 1,
         1, 0, 0,
         1, 0, 1,
         1, 1, 0,
         1, 1, 1,
    };
    //int subtree_size = sizeof(struct CARTSubtreeForEvaluation);
    //printf("Creating CART tree for evaluation. Subtree size:%d\n", subtree_size);

    int alloc_subtrees = 128;
    int num_subtrees = 0;
    struct CARTSubtreeForEvaluation *subtrees = mkl_calloc(alloc_subtrees,
        sizeof(struct CARTSubtreeForEvaluation), 64);
    
    int alloc_leaf_values = 128;
    int num_leaf_values = 0;
    int leaf_value_dim = cart.output_dim;
    float *leaf_values
        = mkl_calloc(alloc_leaf_values,sizeof(float)*leaf_value_dim,64);

    struct QueueNode{
        int i;
        struct QueueNode *next;
    };

    struct CARTNode *nodes = cart.nodes;

    struct QueueNode *queue = calloc(1,sizeof(struct QueueNode));
    queue[0].i = 0;
    int queue_length = 1;

    while(queue!=0){
        struct QueueNode *tmp = queue;
        int subtree_root = tmp->i;
        queue = queue[0].next;
        queue_length--;
        free(tmp);

        int offset = queue_length;
        struct CARTNode *subtree_nodes[15] = {0};
        int num_children = 0;
        int children_bits = 0;
        if(subtree_root > -1){
            for(int i=0;i<15;i++){
                const int *t = traversal_order + i*3;
                struct CARTNode *node = nodes + subtree_root;
                int depth = 0;
                for(int d=0;d<3;d++){
                    if(t[d] != -1){
                        if(node){
                            if(node->left_node_index == 0){
                                node = 0;
                            }else{
                                node = nodes + node->left_node_index + t[d];
                                depth++;
                            }
                        }
                    }
                }
                subtree_nodes[i] = node;
                if(depth == 3 && node->left_node_index != 0){
                    for(int q=0;q<2;q++){
                        struct QueueNode *queue_node
                            = calloc(1,sizeof(struct QueueNode));
                        
                        queue_node->i = node->left_node_index+q;
                        num_children++;
                        if(queue == 0){
                            queue = queue_node;
                        }else{
                            struct QueueNode *queue_last = queue;
                            while(queue_last->next != 0){
                                queue_last = queue_last->next;
                            }
                            queue_last->next = queue_node;
                        }
                        queue_length++;
                    }
                    children_bits |= 1<<(i-7);
                }
            }
        }


        if(num_subtrees+1 > alloc_subtrees){
            alloc_subtrees *= 2;
            struct CARTSubtreeForEvaluation *tmp = mkl_calloc(alloc_subtrees,
              sizeof(struct CARTSubtreeForEvaluation), 64);
            memcpy(tmp,subtrees,num_subtrees
                *sizeof(struct CARTSubtreeForEvaluation));
            mkl_free(subtrees);
            subtrees = tmp;
        }
        struct CARTSubtreeForEvaluation *subtree = subtrees+num_subtrees;
        num_subtrees++;
        
        for(int i=0;i<15;i++){
            struct CARTNodeForEvaluation n = {0};
            struct CARTNode *node = subtree_nodes[i];
            if(node){
                if(node->left_node_index == 0){
                    n.leaf_index = num_leaf_values;
                    assert((n.feature & CART_LEAF_SENTINEL16) == 0);
                    n.feature |= CART_LEAF_SENTINEL16;
                    if(num_leaf_values+1 > alloc_leaf_values){
                        alloc_leaf_values *= 2;
                        float *tmp = mkl_calloc(alloc_leaf_values,
                            sizeof(float)*leaf_value_dim, 64);
                        memcpy(tmp, leaf_values,
                            num_leaf_values*sizeof(float)*leaf_value_dim);
                        mkl_free(leaf_values);
                        leaf_values = tmp;
                    }
                    memcpy(leaf_values+num_leaf_values*leaf_value_dim,
                        node->val,leaf_value_dim*sizeof(float));
                    num_leaf_values++;
                }else{
                    int dim = node->predictor;
                    double f = node->split;
                    if(f >  10.0) f =  10.0;
                    if(f < -10.0) f = -10.0;
                    n.split
                        = (uint16_t)((f+10.0)*((double)CART_SPLIT_MAX)/20.0);
                    n.feature = dim+1;
                }
            }
            subtree->nodes[i] = n;
        }
        subtree->offset = (offset<<8) | children_bits;
    }
    
    struct CARTForEvaluation ret = {0};
    ret.subtrees = subtrees;
    ret.num_subtrees = num_subtrees;
    ret.leaf_value_dim = leaf_value_dim;
    ret.num_leaf_values = num_leaf_values;
    ret.leaf_values = leaf_values;
    return ret;
}

int get_subtree_offset(int offset, int index, int right)
{
    int ret = offset>>8;
    index -= 7;
    for(int i=0;i<index;i++){
        ret += 2*((offset>>i) & 1);
    }
    return ret+right+1;
}

void cart_for_evaluation_print(struct CARTForEvaluation cart)
{
    int leaf_value_dim = cart.leaf_value_dim;
    for(int i=0;i<cart.num_subtrees;i++){
        printf("Subtree %d:\n", i);
        struct CARTSubtreeForEvaluation s = cart.subtrees[i];
        for(int j=0;j<15;j++){
            struct CARTNodeForEvaluation n = s.nodes[j];
            if(n.feature & CART_LEAF_SENTINEL16)
            {
                int leaf_index = n.leaf_index & (~CART_LEAF_SENTINEL32);
                float *l = cart.leaf_values + leaf_index*leaf_value_dim;
                printf("[Node %d]: Leaf %d (", j, n.leaf_index);
                for(int k=0;k<leaf_value_dim-1;k++){
                    printf("%f, ", l[k]);
                }
                printf("%f)\n", l[leaf_value_dim-1]);
            }else{
                if(n.feature == 0){
                    printf("[Node %d]: ---\n", j);
                }else{
                    double split = ((double)n.split)
                        /(double)CART_SPLIT_MAX*20.0 - 10.0;
                    if(j<7){
                        printf("[Node %d]: if feature %d < %f, go to Node %d, "
                            "otherwise go to Node %d\n", j, n.feature-1, split,
                            2*j+1, 2*j+2);
                    }else{
                        printf("[Node %d]: if feature %d < %f, go to Subtree %d, "
                            "otherwise go to Subtree %d\n", j, n.feature-1, split,
                            i + get_subtree_offset(s.offset,j,0),
                            i + get_subtree_offset(s.offset,j,1));
                    }
                }
            }
        }
    }
}

#include <mmintrin.h>
#include <xmmintrin.h>

void cart_evaluate(struct CARTForEvaluation cart, uint16_t *input, float *output)
{
    int leaf_value_dim = cart.leaf_value_dim;
    int done = 0;
    int leaf_value_index = 0;
    struct CARTSubtreeForEvaluation *subtree = cart.subtrees;
    while(!done){
        int i0 = 0;
        int feature0 = subtree->nodes[i0].feature;
        _mm_prefetch(subtree+(subtree->offset>>8)+1,_MM_HINT_T1);
        if((feature0 & CART_LEAF_SENTINEL16) == 0){
            int i1 = 1
                + (input[feature0-1] >= subtree->nodes[i0].split);
            int feature1 = subtree->nodes[i1].feature;
            if((feature1 & CART_LEAF_SENTINEL16) == 0){
                int i2 = i1*2+1
                    + (input[feature1-1] >= subtree->nodes[i1].split);
                int feature2 = subtree->nodes[i2].feature;
                if((feature2 & CART_LEAF_SENTINEL16) == 0){
                    int i3 = i2*2+1
                        + (input[feature2-1] >= subtree->nodes[i2].split);
                    int feature3 = subtree->nodes[i3].feature;
                    if((feature3 & CART_LEAF_SENTINEL16) == 0){
                        //TODO(Vidar):Inline this?
                        int subtree_offset = get_subtree_offset(subtree->offset,
                            i3,(input[feature3-1] >= subtree->nodes[i3].split));
                        subtree += subtree_offset;
                    }
                    else{
                        done = 1;
                        leaf_value_index = subtree->nodes[i3].leaf_index;
                    }
                }
                else{
                    done = 1;
                    leaf_value_index = subtree->nodes[i2].leaf_index;
                }
            }
            else{
                done = 1;
                leaf_value_index = subtree->nodes[i1].leaf_index;
            }
        }
        else{
            done = 1;
            leaf_value_index = subtree->nodes[i0].leaf_index;
        }
    }
    memcpy(output, cart.leaf_values
        + (leaf_value_index & (~CART_LEAF_SENTINEL32))*leaf_value_dim,
        leaf_value_dim*sizeof(float));
}
