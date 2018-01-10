#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "catalog.h"
#include "kdtree.h"

enum token_type {
	LEAF_TOKEN,
	BRANCH_TOKEN
};

struct token {
	enum token_type type;
	struct token *up;
	struct token *left;
	struct token *right;
	unsigned char dim;
	Sample *sample;
};

static int sort_on_x(const void *a, const void *b) {
	Sample *sa = (Sample*) a;
	Sample *sb = (Sample*) b;

	/* TODO
use union to store a - b (double), and use unsigned char to have eather 0 or 1
char value. The use sort_return array.
	*/

	return sa->vector[0] > sb->vector[0] ? 1 : -1;
}

static int sort_on_y(const void *a, const void *b) {
	Sample *sa = (Sample*) a;
	Sample *sb = (Sample*) b;
	return sa->vector[1] > sb->vector[1] ? 1 : -1;
}

static int sort_on_z(const void *a, const void *b) {
	Sample *sa = (Sample*) a;
	Sample *sb = (Sample*) b;
	return sa->vector[2] > sb->vector[2] ? 1 : -1;
}

static int (*sorters[3]) (const void*, const void*) = {
	sort_on_x,
	sort_on_y,
	sort_on_z
};

#define DIM_X 0
#define DIM_Y 1
#define DIM_Z 2
static int dim_switch[3] = {1,2,0};
struct token *build_token(struct token *up, Sample *samples, int nsamples, int dim) {

	struct token *self = malloc(sizeof(struct token));
	self->up = up;
	if (nsamples == 1) {
		printf("end of node for sample %p\n", samples);
		self->type = LEAF_TOKEN;
		self->right = self->left = NULL;
		self->sample = samples;
	} else {
		self->type = BRANCH_TOKEN;
		self->dim = dim;
		qsort(samples, nsamples, sizeof(Sample), sorters[dim]);
		int half = nsamples / 2;
		self->left = build_token(self, samples, half, dim_switch[dim]);
		self->right = build_token(self, &samples[half], nsamples - half, dim_switch[dim]);
	}
	return self;
}

void free_token(struct token *node) {

	if (node->left)
		free_token(node->left);

	if (node->right)
		free_token(node->right);

	free(node);
}

static int sort_return[] = {-1, 1};
union sort_sp {
	double d;
	unsigned char c;
};

int main(int argc, char **argv) {
	int i;
	Field field;
	Set set;

	Catalog_open(argv[1], &field);

	int total_samples = 0;
	for (i=0; i<field.nsets; i++) {
		set = field.sets[i];
		total_samples += set.nsamples;
	}

	Sample *allsamples = malloc(sizeof(Sample) * total_samples);
	Sample *pos = allsamples;
	for (i=0; i<field.nsets; i++) {
		set = field.sets[i];
		memcpy(pos, set.samples, sizeof(Sample) * set.nsamples);
		pos += set.nsamples;
	}

	printf("total %i\n", total_samples);

	struct token *root = build_token(NULL, allsamples, total_samples, DIM_X);
	free_token(root);

	union sort_sp a;
	a.d = 0.2;
	printf("char %i %p\n", a.c, &a.c);
	printf("double %lf %p\n", a.d, &a.d);

  for (i = 0; i < 8; i++) {
      printf("%d", !!((a.c << i) & 0x80));
  }

	a.d = -0.2;
	printf("char %i %p\n", a.c, &a.c);
	printf("double %lf %p\n", a.d, &a.d);

  for (i = 0; i < 8; i++) {
      printf("%d", !!((a.c << i) & 0x80));
  }


	Catalog_freeField(&field);
	return 0;
}
