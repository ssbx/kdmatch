#ifndef _KDTREE_H_
#define _KDTREE_H_

typedef struct kd_node kd_node;
struct kd_node {
    double x[3];
	short int dim;
    struct kd_node *left, *right;
    void *data;
	double data_x[3];
};

struct kd_node *kd_new();
void kd_insert(struct kd_node *root, void *data, double* v);
double dist(kd_node*, kd_node*, int);
void swap(kd_node*, kd_node*);
kd_node* find_median(kd_node*, kd_node*, int);
kd_node* make_tree(kd_node*, int, int, int);
void nearest(kd_node*, kd_node*, int, int, kd_node**, double*,int*);

#endif /* SRC_KDTREE_H_ */
