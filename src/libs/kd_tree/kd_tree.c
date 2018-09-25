#include "kd_tree.h"
#include <stdlib.h>
#include <float.h>
#include <string.h>
#ifdef _WIN32
#include <malloc.h>
#endif
#include <stdio.h>
#include "sort/sort.h"

struct ConstructionSubset {
    float *points;
    int *indices;
    int node;
    int num_points;
};

int compare_points(void *a_in, void *b_in, void*data) {
    int axis = *(int*)data;
    float *a = (float*)a_in;
    float *b = (float*)b_in;
    return a[axis] <= b[axis];
}

struct KDTree kd_tree_construct(int num_points, float *points, int dim)
{
    float *sort_points_tmp    = calloc(num_points, sizeof(float)*dim);
    int *point_indices = calloc(num_points, sizeof(int));
    int *sort_indices_tmp    = calloc(num_points, sizeof(int));
    for (int i = 0; i < num_points; i++) {
        point_indices[i] = i;
    }

    struct KDTree ret = { 0 };
    ret.dim = dim;
    int alloc_nodes = 128;
    ret.nodes = calloc(alloc_nodes, sizeof(struct KDTreeNode));
    ret.num_nodes = 1;
    ret.nodes[0].axis = 0;
    ret.nodes[0].parent = -1;

    int queue_num = 0;
    int queue_alloc = 128;
    struct ConstructionSubset *queue = calloc(queue_alloc,sizeof(struct ConstructionSubset));

    queue[0].points = points;
    queue[0].indices = point_indices;
    queue[0].num_points = num_points;
    queue[0].node = 0;
    queue_num++;

    while (queue_num > 0){
        if (queue_alloc < queue_num + 2) {
            queue_alloc *= 2;
            queue = realloc(queue, queue_alloc * sizeof(struct ConstructionSubset));
        }
        if (alloc_nodes < ret.num_nodes + 2) {
            alloc_nodes *= 2;
            ret.nodes = realloc(ret.nodes, alloc_nodes * sizeof(struct KDTreeNode));
        }

        queue_num--;
        float *current_points = queue[queue_num].points;
        int *current_indices = queue[queue_num].indices;
        int num_current_points = queue[queue_num].num_points;
        int current_node_index = queue[queue_num].node;
        struct KDTreeNode *current_node = ret.nodes + current_node_index;
        int current_axis = current_node->axis;

        int median_index = num_current_points / 2;

        int num_left = median_index;
        int num_right = num_current_points - median_index -1;
        
        merge_sort_custom_auxillary(current_points, sort_points_tmp, dim*sizeof(float), current_indices,
            sort_indices_tmp, num_current_points, compare_points, &current_axis);

        current_node->id = current_indices[median_index];
        memcpy(current_node->point, current_points+median_index*dim, sizeof(float)*dim);


        int new_axis = current_axis + 1;
        if (new_axis >= dim) new_axis = 0;

        current_node->left_node = -1;
        current_node->right_node = -1;

        if (num_left > 0) {
            struct ConstructionSubset *left_subset = queue + queue_num;
            queue_num++;
            left_subset->points = current_points;
            left_subset->indices = current_indices;
            left_subset->num_points = num_left;
            struct KDTreeNode *new_node = ret.nodes + ret.num_nodes;
            left_subset->node = ret.num_nodes;
            current_node->left_node = ret.num_nodes;
            ret.num_nodes++;
            new_node->axis = new_axis;
            new_node->parent = current_node_index;
        }

        if (num_right > 0) {
            struct ConstructionSubset *right_subset = queue + queue_num;
            queue_num++;
            right_subset->points = current_points + (median_index+1)*dim;
            right_subset->indices = current_indices + median_index+1;
            right_subset->num_points = num_right;
            struct KDTreeNode *new_node = ret.nodes + ret.num_nodes;
            right_subset->node = ret.num_nodes;
            current_node->right_node = ret.num_nodes;
            ret.num_nodes++;
            new_node->axis = new_axis;
            new_node->parent = current_node_index;
        }
    }

    free(queue);
    return ret;
}

void kd_tree_nearest_recursive(struct KDTree *tree, int base_index, float *pos,
    float *best_dist, int *best_id, void(*print_func)(char *str))
{
    char buffer[512];
    int index = base_index;
    int last_index = index;
    while (index >= 0)
    {
        int other_index;
        last_index = index;
        struct KDTreeNode *node = tree->nodes + index;
        int axis = node->axis;
        if (pos[axis] < node->point[axis]) {
            index = node->left_node;
            other_index=node->right_node;
            if(print_func)                {
                print_func("Going left\n");
            }
        }
        else {
            index = node->right_node;
            other_index=node->left_node;
            if(print_func)                {
                print_func("Going right\n");
            }
        }
        if(index==-1&&other_index>=0){
            index=other_index;
        }
    }
    struct KDTreeNode *leaf_node = tree->nodes + last_index;
    float dist = 0.f;
    for (int d = 0; d < tree->dim; d++) {
        float diff = leaf_node->point[d] - pos[d];
        dist += diff*diff;
    }
    if(print_func)                {
        sprintf(buffer,"Checking node %d: dist %f\n",leaf_node->id,dist);
        print_func(buffer);
    }
    if (dist < *best_dist) {
        *best_dist = dist;
        *best_id = leaf_node->id;
    }
    index = leaf_node->parent;
    while (last_index != base_index) {
        struct KDTreeNode *node = tree->nodes + index;
        float dist = 0.f;
        for (int d = 0; d < tree->dim; d++) {
            float diff = node->point[d] - pos[d];
            dist += diff*diff;
        }
        if(print_func)                {
            sprintf(buffer,"Checking node %d: dist %f\n",node->id,dist);
            print_func(buffer);
        }
        if (dist < *best_dist) {
            *best_dist = dist;
            *best_id = node->id;
        }
        float dist_to_split=0.f;
        {
            float d=node->point[node->axis]-pos[node->axis];
            dist_to_split=d*d;
        }
        if(dist_to_split<*best_dist){
            if(print_func)                {
                sprintf(buffer, "Entering subtree at node %d:\n", node->id);
                print_func(buffer);
            }
            int other_index=node->left_node;
            if(other_index==last_index) other_index=node->right_node;
            if(other_index>=0){
                kd_tree_nearest_recursive(tree,other_index,pos,best_dist,best_id,print_func);
            }
        }
        last_index = index;
        index = node->parent;
    }
}

int kd_tree_nearest(struct KDTree *tree, float *pos, void(*print_func)(char *str))
{
    float best_dist = FLT_MAX;
    int best_id = 0;
    if((0)){ //Debug check...
        float best_dist=FLT_MAX;
        for(int i=0;i<tree->num_nodes;i++){
            struct KDTreeNode *node = tree->nodes + i;
            float dist = 0.f;
            for (int d = 0; d < tree->dim; d++) {
                float diff = node->point[d] - pos[d];
                dist += diff*diff;
            }
            if (dist < best_dist) {
                best_dist = dist;
                best_id = node->id;
            }
        }
        int best_id_tree = 0;
        float best_dist_tree=FLT_MAX;
        kd_tree_nearest_recursive(tree, 0, pos, &best_dist_tree, &best_id_tree, print_func);
        if(best_id!=best_id_tree){
            print_func("Does not match!!\n");
            kd_tree_print(tree,print_func);
        }
        print_func("\n\n");
    } else{
        kd_tree_nearest_recursive(tree, 0, pos, &best_dist, &best_id, 0);
    }
    return best_id;
}


void kd_tree_print(struct KDTree *tree, void(*print_func)(char *str)) {
    char buffer[512];
    for (int i = 0; i < tree->num_nodes; i++) {
        struct KDTreeNode *n = tree->nodes + i;
        sprintf(buffer, "[%d] axis: %d, split: %f, left: %d, right: %d, parent: %d, id: %d\n    (",
            i, n->axis, n->point[n->axis], n->left_node, n->right_node, n->parent, n->id);
        print_func(buffer);
        for (int d = 0; d < tree->dim-1; d++) {
            sprintf(buffer, "%f, ", n->point[d]);
            print_func(buffer);
        }
        sprintf(buffer, "%f)\n", n->point[tree->dim-1]);
        print_func(buffer);
    }
}
