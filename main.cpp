// g++ -std=c++11 -O3 -march=native main.cpp -o main
// g++ -std=c++11 -O3 -march=native -fopenmp main.cpp -o main
#include <algorithm>
#include <stdio.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <iostream>
#include "utils.h"

#define VERB 0
#define THREADNUM 8
#define PVER 2

struct Node {
  double* x = NULL;  // array of points
  double* q = NULL;  // array of charges
  int*    I = NULL;  // array of indices
  int     n = 0;     // number of points 
  double  c = 0;     // center of cell
  double* w = NULL;  // weights
  int     p = 0;     // number of weights
  int     level = 0; // level in tree (root is 0)
  int     box   = 0; // box in level (0 to 2^l-1)
  struct Node* parent = NULL; // parent node
  struct Node* left   = NULL; // left child node
  struct Node* right  = NULL; // right child node
 
  Node(
    double* _x, double* _q, int* _I, 
    int _n, double _c, int _p, 
    int _level, int _box,
    Node* _parent
    ) {
    x = _x;
    q = _q;
    I = _I;
    n = _n;
    c = _c;
    p = _p;
    level = _level;
    box   = _box;
    parent = _parent;
    w = (double*) malloc(p * sizeof(double));
    for (int m = 0; m < p; m++) w[m] = 0;
  }
};

Node* build_tree(
  double* x, double* q, int* I, int n, 
  int max_pts, double a, double b, int p, 
  int level, int box,
  Node* parent
  ) {
  struct Node* tree = new Node(x, q, I, n, (a+b)/2, p, level, box, parent);

  // if there are more than max_pts points in the cell, keep subdividing
  if (n > max_pts) {
    // find the number of elements in the left half of this cell 
    // using a for loop here is naive and slow, could use bisection
    int nleft = 0;
    for (int i = 0; i < n; i++) {
      if (x[i] < (a+b)/2) {
        nleft = i+1;
      }
    }
    int nright = n - nleft;
    if (nleft > 0 && nright > 0) { // avoid subdividing empty cells
      tree->left  = build_tree(
        tree->x, tree->q, tree->I, nleft, 
        max_pts, a, (a+b)/2, p, 
        level+1, 2*box, tree
        );
      tree->right = build_tree(
        tree->x + nleft, tree->q + nleft, tree->I + nleft, nright, 
        max_pts, (a+b)/2, b, p, 
        level+1, 2*box+1, tree
        );
    }
  }
  return tree;
}

void compute_weights(Node* tree) {

	#pragma omp parallel num_threads(THREADNUM) if(PVER==1) //if(tree->n > 50 && tree->p > 4)
	{
	#pragma omp for
	for (int m = 0; m < tree->p; m++) {
    for (int i = 0; i < tree->n; i++) {
      // add weight term q_i a_m(x[i] - c) 
      // where a_m(x[i] - c) = (c - x[i])^m
      if (VERB)  printf("w[%i]: %.3f ", m, tree->w[m]);
      tree->w[m] += tree->q[i] * std::pow((tree->c - tree->x[i]), m);
      if (VERB)  printf(
        "-> %.3f (q[%i] = %.3f, x[%i] = %.3f, c = %.3f)\n", 
        tree->w[m], i, tree->q[i], i, tree->x[i], tree->c
        );
    }
  }
	}

  // traverse the tree
  if (tree->left) {
    compute_weights(tree->left); 
    compute_weights(tree->right);
  }
}

void add_near_field(double* u, Node* source, Node* target) {
  // evaluate near-field directly
	
	#pragma omp parallel num_threads(THREADNUM) if(PVER==1)
	{
	#pragma omp for
  for (int i = 0; i < target->n; i++) {
    for (int j = 0; j < source->n; j++) {
      if (std::abs(target->x[i] - source->x[j]) > 1e-16) {
        // add exact potential q_i / |x[i] - x[j]|
        if (VERB) printf("u[%i]: %.3f ", target->I[i], u[target->I[i]]);
        u[target->I[i]] += source->q[j] / std::abs(target->x[i] - source->x[j]);
        if (VERB) printf(
          "-> %.3f (q[%i] = %.3f, x[%i] = %.3f, x[%i] = %.3f)\n", 
          u[target->I[i]], j, source->q[j], i, target->x[i], j, source->x[j]
          );
      }
    }
  } 
	} 
	
}

void add_far_field(double* u, Node* source, Node* target) {
  // evaluate far-field using multipole expansions

	#pragma omp parallel num_threads(THREADNUM) if(PVER==1) //if(target->n > 50 && target->p > 4)
	{
	#pragma omp for schedule(dynamic)
  for (int i = 0; i < target->n; i++) {
    for (int m = 0; m < target->p; m++) {
      // add approximate potential term w_{i,m} * S_m(c - x[i]) 
      // where S_m(c - x[i]) = 1 / (|c - x[i]|*(c - x[i])^m)
      if (VERB) printf("u[%i]: %.3f ", target->I[i], u[target->I[i]]);
      u[target->I[i]] += source->w[m] / (
        std::abs(source->c - target->x[i])*std::pow(source->c - target->x[i], m)
        );
      if (VERB) printf(
        "-> %.3f (w[%i] = %.3f, x[%i] = %.3f, c = %.3f)\n", 
        u[target->I[i]], m, source->w[m], i, target->x[i], source->c
        );
    }
  } 
	} 
 	
}

void add_potential(double* u, Node* target, bool left) {
  // all the comments here are assuming left = 1
  // switch right and left in the comments for the case left = 0
  bool leaf = !(left ? target->left : target->right); // assumes some symmetry

  // add left near-field terms if we're on the lowest level
  if (leaf) {
    if ((left ? target->parent->right : target->parent->left) == target) {
      if (VERB) printf(
        "adding near-field terms for source %.4f and target %.4f\n", 
        (left ? target->parent->left : target->parent->right)->c, target->c
      );
      add_near_field(
        u, left ? target->parent->left : target->parent->right, target
        );
    }
  }

  Node* source = target->parent; // initialize pointer we'll search tree with
  int l = 1; // how many levels have we moved up in our search
  if (VERB) printf("search cell center is %.4f, l = %i\n", source->c, l);
 
  // ascend the tree until we are not the left child of our parent
  while (
    source->parent && 
    (left ? source->parent->left : source->parent->right) == source
    ) { 
    source = source->parent;
    l++;
    if (VERB) printf("ascending cell center is %.4f, l = %i\n", source->c, l);
  }

  if (source->parent) {
    // now we're the right child of our parent
    // so we cross the tree to be in the left branch
    source = left ? source->parent->left : source->parent->right;
    if (VERB) printf("going opposite branch cell center is %.4f, l = %i\n", source->c, l);

    // descend the tree to the right until we are one level above the base node
    for (; l > 1; l--) {
      source = left ? source->right : source->left;
      if (VERB) printf("descending cell center is %.4f, l = %i\n", source->c, l);
    }

    // add left child interaction as it will always be in the interaction list
    if (VERB) printf(
        "adding interaction list terms for source %.4f and target %.4f\n", 
        (left ? source->left : source->right)->c, target->c
      );
    add_far_field(u, left ? source->left : source->right, target);

    // if the base node was originally a right node, both children here are
    // in the interaction list, so add the right child interaction also
    if ((left ? target->parent->right : target->parent->left) == target) {
      if (VERB) printf(
        "adding interaction list terms for source %.4f and target %.4f\n", 
        (left ? source->right : source->left)->c, target->c
      );
      add_far_field(u, left ? source->right : source->left, target);
    } else if (leaf) {
      // we're on the lowest level and right child is in the near field 
      // so add it directly
      if (VERB) printf(
        "adding near-field terms for source %.4f and target %.4f\n", 
        (left ? source->right : source->left)->c, target->c
      );
      add_near_field(u, left ? source->right : source->left, target);
    }
  }
  if (VERB) printf("\n");
}

void compute_potential(double* u, Node* tree) {
  // the top level has no near field or interaction list
  if (tree->parent) {     
    if (VERB) printf("adding left interaction for cell with center %.4f\n", tree->c);
    add_potential(u, tree, 1);

    if (VERB) printf("adding right interaction for cell with center %.4f\n", tree->c);
    add_potential(u, tree, 0);
  }

  // traverse the tree
  if (tree->left) {
    // continue computing far-field terms by expansion
		#pragma omp parallel num_threads(THREADNUM) if((tree->level < std::log2(THREADNUM)) && PVER==2) //if(tree->level==0 && PVER==2)
		{
			#pragma omp sections
			{
				#pragma omp section
				{
		  		compute_potential(u, tree->left); 
				}
				#pragma omp section
				{
			  	compute_potential(u, tree->right);
				}
			}
		}
  } else {
    if (VERB) printf("adding same-cell terms for cell with center %.4f\n", tree->c);
    add_near_field(u, tree, tree);
  }
}

void barnes_hut(double* u, int n, Node* tree, int p) {
  // compute the weights and store them in the tree
  compute_weights(tree);

  // use the precomputed weights to compute the potential
  compute_potential(u, tree);
}

int main(int argc, char** argv) {
  // number of source / target points
  int n = read_option<int>("-n", argc, argv, "128");  
  // maximum number of points per neighborhood
  int max_pts = read_option<int>("-m", argc, argv, "8");
  // number of terms in multipole expansion
  int p = read_option<int>("-p", argc, argv, "1"); 
	
	// draw and sort a set of n uniform random points in [0,1]
	double* x = (double*) malloc(n * sizeof(double));
	for (int i = 0; i < n; i++) x[i] = ((double) i)/n; // ((double)rand()/RAND_MAX);
	std::sort(x, x+n);

	// draw a set of n uniform random charges in [-1, 1]
	double* q = (double*) malloc(n * sizeof(double));
	for (int i = 0; i < n; i++) q[i] = 2*((double)rand()/RAND_MAX) - 1;

	// make a simple index vector [0,...,n]
	int* I = (int*) malloc(n * sizeof(int));
	for (int i = 0; i < n; i++) I[i] = i;

	// build a binary tree subdiving our points
	if (VERB) printf("p = %i\n", p);
	Node* tree = build_tree(x, q, I, n, max_pts, 0, 1, p, 0, 0, NULL);

	// allocate the potential and initialize to zero
	double* u = (double*) malloc(n * sizeof(double));
	for (int i = 0; i < n; i++) u[i] = 0;

	// if version 2, turn on nesting
	if (PVER==2) omp_set_nested(true);

	// run barnes_hut
	Timer tt;
	tt.tic();
	barnes_hut(u, n, tree, p);
	double runtime = tt.toc();


	// display potential at source / target points
	//if (VERB) {
	printf("\nPotential:\n");
	for (int i = n-5; i < n; i++) printf("u(%.3f) = %.3f\n", x[i], u[i]);
	//}

	printf("\nn=%d, m=%d, p=%d, time=%f s\n\n", n, max_pts, p, runtime);

	free(x);
	free(q);
	free(I);
	free(u);
	
}
