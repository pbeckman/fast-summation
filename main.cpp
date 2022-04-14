#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <iostream>

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
      tree->w[m] += tree->q[i] * std::pow((tree->c - tree->x[i]), m);
    }
  }

  // traverse the tree
  if (tree->left) {
    compute_weights(tree->left); 
    compute_weights(tree->right);
  }
}

Node* get_left_interaction(Node* tree) {
  Node* node = tree->parent; // initialize node pointer we'll search tree with
  int l = 1; // how many levels have we moved up in our search
  
  // ascend the tree until we are not the left child of our parent
  while (node->parent->left == node) { 
    node = node->parent;
    l++;
  }

  // now we're the right child of our parent
  // so we cross the tree to be in the left branch
  node = node->parent->left;

  // descend the tree to the right until we are one level above the base node
  for (; l > 1; l--) node = node->right;

  // if the base node was originally a left node, the right child here
  // is in the near field, so give only the left child
  if (tree->parent->left == tree) {
    node = node->left;
  }

  return node;
}

Node* get_right_interaction(Node* tree) {
  Node* node = tree->parent; // initialize node pointer we'll search tree with
  int l = 1; // how many levels have we moved up in our search
  
  // ascend the tree until we are not the right child of our parent
  while (node->parent->right == node) { 
    node = node->parent;
    l++;
  }

  // now we're the left child of our parent
  // so we cross the tree to be in the right branch
  node = node->parent->right;

  // descend the tree to the left until we are one level above the base node
  for (; l > 1; l--) node = node->left;

  // if the base node was originally a right node, the left child here
  // is in the near field, so give only the right child
  if (tree->parent->right == tree) {
    node = node->right;
  }

  return node;
}

void add_far_field(double* u, Node* tree, Node* node) {
  // evaluate far-field using multipole expansions
  for (int i = 0; i < tree->n; i++) {
    for (int m = 0; m < tree->p; m++) {
      // add approximate potential term w_{i,m} * S_m(c - x[i]) 
      // where S_m(c - x[i]) = 1 / (|c - x[i]|*(c - x[i])^m)
      u[tree->I[i]] += node->w[m] * (
        std::abs(node->c - node->x[i])*std::pow(node->c - node->x[i], m)
        );
    }
  }
}

void compute_potential(double* u, Node* tree) {
  // the top two levels have empty interaction lists
  if (tree->parent && tree->parent->parent) { 
    Node* node = NULL;
    
    node = get_left_interaction(tree);
    add_far_field(u, tree, node);

    node = get_right_interaction(tree);
    add_far_field(u, tree, node);
  }

  // traverse the tree
  if (tree->left) {
    // continue computing far-field terms by expansion
    compute_potential(u, tree->left); 
    compute_potential(u, tree->right);
  } else {
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

int main() {
  int n = 100; // number of source / target points
  int max_pts = 10; // maximum number of points per neighborhood
  int p = 5; // number of terms in multipole expansion

  // draw and sort a set of n uniform random points in [0,1]
  double* x = (double*) malloc(n * sizeof(double));
  for (int i = 0; i < n; i++) x[i] = ((double)rand()/RAND_MAX); // ((double)i)/n;
  std::sort(x, x+n);

  // draw a set of n uniform random charges in [-1, 1]
  double* q = (double*) malloc(n * sizeof(double));
  for (int i = 0; i < n; i++) q[i] = 2*((double)rand()/RAND_MAX) - 1;

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