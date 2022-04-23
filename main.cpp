#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <iostream>
#include "utils.h"

struct Node {
  double* x = NULL; // array of points
  double* q = NULL; // array of charges
  int*    I = NULL; // array of indices
  int     n = 0;    // number of points 
  double  c = 0;    // center of cell
  double* w = NULL; // weights
  int     p = 0;    // number of weights
  struct Node* parent = NULL; // parent node
  struct Node* left   = NULL; // left child node
  struct Node* right  = NULL; // right child node
 
  Node(
    double* _x, double* _q, int* _I, 
    int _n, double _c, int _p, Node* _parent
    ) {
    x = _x;
    q = _q;
    I = _I;
    n = _n;
    c = _c;
    p = _p;
    parent = _parent;
    w = (double*) malloc(p * sizeof(double));
    for (int m = 0; m < p; m++) w[m] = 0;
  }
};

Node* build_tree(
  double* x, double* q, int* I, int n, 
  int max_pts, double a, double b, int p,
  Node* parent
  ) {
  struct Node* tree = new Node(x, q, I, n, (a+b)/2, p, parent);

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
        max_pts, a, (a+b)/2, p, tree
        );
      tree->right = build_tree(
        tree->x + nleft, tree->q + nleft, tree->I + nleft, nright, 
        max_pts, (a+b)/2, b, p, tree
        );
      
    }
  }
  return tree;
}

void compute_weights(Node* tree) {
  for (int i = 0; i < tree->n; i++) {
    for (int m = 0; m < tree->p; m++) {
      // add weight term q_i a_m(x[i] - c) 
      // where a_m(x[i] - c) = (c - x[i])^m
      printf("w[m]: %.3f ", tree->w[m]);
      tree->w[m] += tree->q[i] * std::pow((tree->c - tree->x[i]), m);
      printf(
        "-> %.3f (q[i] = %.3f, x[i] = %.3f, c = %.3f)\n", 
        tree->w[m], tree->q[i], tree->x[i], tree->c
        );
    }
  }

  // traverse the tree
  if (tree->left) {
    compute_weights(tree->left); 
    compute_weights(tree->right);
  }
}

Node* get_interaction(Node* tree, bool left) {
  // all the comments here are assuming left = 1
  // switch right and left in the comments for the case left = 0
  Node* node = tree->parent; // initialize node pointer we'll search tree with
  int l = 1; // how many levels have we moved up in our search
  printf("search cell center is %.4f, l = %i\n", node->c, l);
  
  // ascend the tree until we are not the left child of our parent
  while (
    node->parent && (left ? node->parent->left : node->parent->right) == node
    ) { 
    node = node->parent;
    l++;
    printf("ascending cell center is %.4f, l = %i\n", node->c, l);
  }

  if (node->parent) {
    // now we're the right child of our parent
    // so we cross the tree to be in the left branch
    node = left ? node->parent->left : node->parent->right;
    printf("going opposite branch cell center is %.4f, l = %i\n", node->c, l);

    // descend the tree to the right until we are one level above the base node
    for (; l > 1; l--) {
      node = left ? node->right : node->left;
      printf("descending cell center is %.4f, l = %i\n", node->c, l);
    }

    // if the base node was originally a left node, the right child here
    // is in the near field, so give only the left child
    if ((left ? tree->parent->left : tree->parent->right) == tree) {
      if (left ? node->left : node->right) {
        node = left ? node->left : node->right;
        l--;
      } else {
        // we can't go down to the correct level to get out of the near-field
        // so we have no left interactions
        return nullptr;
      }
    }

    return node;
  } else { 
    // we're at the root node, which means we're the leftmost cell
    // at this level, so we have no left interactions
    return nullptr;
  }
  
}

void add_far_field(double* u, Node* tree, Node* node) {
  // evaluate far-field using multipole expansions
  for (int i = 0; i < tree->n; i++) {
    for (int m = 0; m < tree->p; m++) {
      // add approximate potential term w_{i,m} * S_m(c - x[i]) 
      // where S_m(c - x[i]) = 1 / (|c - x[i]|*(c - x[i])^m)
      printf("u[i]: %.3f ", u[tree->I[i]]);
      u[tree->I[i]] += node->w[m] / (
        std::abs(node->c - node->x[i])*std::pow(node->c - node->x[i], m)
        );
      printf(
        "-> %.3f (w[m] = %.3f, x[i] = %.3f, c = %.3f)\n", 
        u[tree->I[i]], node->w[m], node->x[i], node->c
        );
    }
  }
}

void compute_potential(double* u, Node* tree) {
  // the top two levels have empty interaction lists
  if (tree->parent && tree->parent->parent) { 
    Node* node = NULL;
    
    printf("computing left interaction for cell with center %.4f\n", tree->c);
    node = get_interaction(tree, 1);
    if (node) {
      printf("left interaction cell center is %.4f\n\n", node->c);
    } else {
      printf("no left interaction\n\n");
    }
    if (node) add_far_field(u, tree, node);

    printf("computing right interaction for cell with center %.4f\n", tree->c);
    node = get_interaction(tree, 0);
    if (node) {
      printf("right interaction cell center is %.4f\n\n", node->c);
    } else {
      printf("no right interaction\n\n");
    }
    if (node) add_far_field(u, tree, node);
  }

  // traverse the tree
  if (tree->left) {
    // continue computing far-field terms by expansion
    compute_potential(u, tree->left); 
    compute_potential(u, tree->right);
  } else {
    printf("adding near-field terms for cell with center %.4f\n\n", tree->c);
    // evaluate near-field directly
    for (int i = 0; i < tree->n; i++) {
      for (int j = 0; j < tree->n; j++) {
        if (i != j) {
          // add exact potential q_i / |x[i] - x[j]|
          u[tree->I[i]] += tree->q[i] / std::abs(tree->x[i] - tree->x[j]);
        }
      }
    }
  }
}

void barnes_hut(double* u, int n, Node* tree, int p) {
  // allocate temporary vectors to store source and target points during the computation
  double* sources = (double*) malloc(n * sizeof(double));
  double* targets = (double*) malloc(n * sizeof(double));

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
  int p = read_option<int>("-m", argc, argv, "1"); 

  // draw and sort a set of n uniform random points in [0,1]
  double* x = (double*) malloc(n * sizeof(double));
  for (int i = 0; i < n; i++) x[i] = ((double)rand()/RAND_MAX); // ((double)i)/n;
  std::sort(x, x+n);

  // draw a set of n uniform random charges in [-1, 1]
  double* q = (double*) malloc(n * sizeof(double));
  for (int i = 0; i < n; i++) q[i] = 1; // 2*((double)rand()/RAND_MAX) - 1;

  // make a simple index vector [0,...,n]
  int* I = (int*) malloc(n * sizeof(int));
  for (int i = 0; i < n; i++) I[i] = i;

  // build a binary tree subdiving our points
  Node* tree = build_tree(x, q, I, n, max_pts, 0, 1, p, NULL);

  // allocate the potential and initialize to zero
  double* u = (double*) malloc(n * sizeof(double));
  for (int i = 0; i < n; i++) u[i] = 0;

  // run barnes_hut
  barnes_hut(u, n, tree, p);

  // display potential at source / target points
  for (int i = 0; i < n; i++) printf("u(%.3f) = %.3f\n", x[i], u[i]);

  free(x);
  free(q);
  free(I);
  free(u);
}