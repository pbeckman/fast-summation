#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <iostream>

struct Node {
  double* x = NULL;    // array of points
  double* q = NULL;    // array of charges
  int n = 0;          // number of points 
  double c = 0;        // center
  double* w = NULL;    // weights
  int p = 0;           // number of weights
  struct Node* left  = NULL; // left child node
  struct Node* right = NULL; // right child node
 
  Node(double* _x, double* _q, int _n, double _c, int _p) {
    x = _x;
    q = _q;
    n = _n;
    c = _c;
    p = _p;
    w = (double*) malloc(p * sizeof(double));
  }
};

Node* build_tree(
  double* x, double* q, int n, int max_pts,
  double a, double b, 
  int p
  ) {
  struct Node* tree = new Node(x, q, n, (a+b)/2, p);

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
        tree->x, tree->q, nleft, max_pts, a, (a+b)/2, p
        );
      tree->right = build_tree(
        tree->x + nleft, tree->q + nleft, nright, max_pts, (a+b)/2, b, p
        );
      
    }
  }
  return tree;
}

void compute_weights(Node* tree) {
  for (int j = 0; j < tree->n; j++) {
    for (int m = 0; m < tree->p; m++) {
      // add weight q_j a_m(x[j] - c) = q_j (c - x[j])^m
      tree->w[m] += tree->q[j] * pow((tree->c - tree->x[j]), m);
    }
  }

  // traverse the tree
  if (tree->left) {
    compute_weights(tree->left); 
    compute_weights(tree->right);
  }
}

void compute_potential(double* u, int n, Node* tree) {
  
}

void barnes_hut(double* u, int n, Node* tree, int p) {
  // allocate temporary vectors to store source and target points during the computation
  double* sources = (double*) malloc(n * sizeof(double));
  double* targets = (double*) malloc(n * sizeof(double));

  // compute the weights and store them in the tree
  compute_weights(tree);

  // use the precomputed weights to compute the potential
  compute_potential(u, n, tree);

  // // Evaluate the potential (top to bottom)
  // for (int l = 0; l < j; l++) {
  //   for (int k = 0; k < pow(2, l); k++) {
  //     // write relevant targets into temporary targets vector, where nt_lk gives the number of targets in T_{l,k}
  //     int nt_lk = get_targets(targets, tree, l, k); 
  //     for (int i = 0; i < nt_lk; i++) {
  //       double y = targets[i];
  //       for (int m = 0; m < p; m++) { // far-field
  //         // write relevant sources into temporary sources vector, where ns_lk gives the number of sources in the interaction list of T_{l,k}
  //         int ns_lk = get_interaction_list(sources, tree, l, k); 
  //         for (j = 0; j < ns_lk; i++) {
  //           double s = sources[j];
  //           // TODO: u(i) = u(i) + w(l,k,m) * Sm(xstar(l,s) - y);
  //         }
  //       }
  //       if (l == j) { // near-field
  //         // TODO: u(i) = u(i) + near-field(y(i), T(J,:));
  //       }
  //     }
  //   }
  // }

}

int main() {
  int n = 100; // number of source / target points
  int max_pts = 10; // maximum number of points per neighborhood
  int p = 5; // number of terms in multipole expansion

  // draw and sort a set of n uniform random points in [0,1]
  double* x = (double*) malloc(n * sizeof(double));
  for (int i = 0; i < n; i++) x[i] = ((double)rand()/RAND_MAX);
  std::sort(x, x+n);

  // draw a set of n uniform random charges in [-1, 1]
  double* q = (double*) malloc(n * sizeof(double));
  for (int i = 0; i < n; i++) q[i] = 2*((double)rand()/RAND_MAX) - 1;

  // build a binary tree subdiving our points
  Node* tree = build_tree(x, q, n, max_pts, 0, 1, p);

  // allocate the potential
  double* u = (double*) malloc(n * sizeof(double));

  // run barnes_hut
  barnes_hut(u, n, tree, p);

  // // display charges at source points
  // for (int i = 0; i < n; i++) printf("u(%.3f) = %.3f\n", x[i], q[i]);
}