#include "kdtree.h"

#include <stdlib.h>
#include <string.h>
#include <time.h>

#define DIM_X 0
#define DIM_Y 1
#define DIM_Z 2
#define DIM_NEXT(dim) (dim > 1) ? 0 : 1

void kd_insert(struct kd_node *node, void* data, double *v) {
}

/*
double dist(kd_node *a, kd_node *b, int dim) {
    double t, d = 0;
    while (dim--) {
        t = a->x[dim] - b->x[dim];
        d += t * t;
    }
    return d;
}

kd_node*
find_median(kd_node *t, kd_node *t2, int a)
{
	return NULL;
}
kd_node*
make_tree(kd_node *t, int len, int i, int dim) {
    kd_node *n;

    if (!len)
        return 0;

    if ((n = find_median(t, t + len, i))) {
        i = (i + 1) % dim;
        n->left = make_tree(t, n - t, i, dim);
        n->right = make_tree(n + 1, t + len - (n + 1), i, dim);
    }
    return n;
}

void nearest(kd_node *root, kd_node *nd, int i, int dim, kd_node **best,
        double *best_dist, int *visited) {
    double d, dx, dx2;

    if (!root)
        return;
    d = dist(root, nd, dim);
    dx = root->x[i] - nd->x[i];
    dx2 = dx * dx;

    (*visited)++;

    if (!*best || d < *best_dist) {
        *best_dist = d;
        *best = root;
    }

    // if chance of exact match is high 
    if (!*best_dist)
        return;

    if (++i >= dim)
        i = 0;

    nearest(dx > 0 ? root->left : root->right, nd, i, dim, best, best_dist,
            visited);
    if (dx2 >= *best_dist)
        return;
    nearest(dx > 0 ? root->right : root->left, nd, i, dim, best, best_dist,
            visited);
}
*/
