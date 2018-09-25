#pragma once

#ifndef KD_MAX_DIM
#define KD_MAX_DIM 3
#endif // !KD_MAX_DIM


struct KDTree {
    struct KDTreeNode *nodes;
    int dim;
    int num_nodes;
};

struct KDTreeNode {
    float point[KD_MAX_DIM];
    int axis;
    int id;
    int left_node;
    int right_node;
    int parent;
};

struct KDTree kd_tree_construct(int num_points, float *points, int dim);
int kd_tree_nearest(struct KDTree *tree,float *pos,void(*print_func)(char *str));
void kd_tree_print(struct KDTree *tree, void(*print_func)(char *str));
