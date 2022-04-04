#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <iostream>

struct Node {
  float* points = NULL;
  int n = 0;
  struct Node* left;
  struct Node* right;
 
  Node(float* _points, int _n) {
    points = _points;
    n = _n;
    left = NULL;
    right = NULL;
  }
};

Node* build_tree(
  float* points, int n, int max_pts,
  float a, float b
  ) {
  struct Node* tree = new Node(points, n);

  // if there are more than max_pts points in the cell, keep subdividing
  if (n > max_pts) {
    // find the number of elements in the left half of this cell 
    // using a for loop here is naive and slow, could use bisection
    int nleft = 0;
    for (int i = 0; i < n; i++) {
      if (points[i] < (a+b)/2) {
        nleft = i+1;
      }
    }
    int nright = n - nleft;
    if (nleft > 0 && nright > 0) { // avoid subdividing empty cells
      tree->left  = build_tree(
        tree->points, nleft, max_pts, a, (a+b)/2
        );
      tree->right = build_tree(
        tree->points + nleft, nright, max_pts, (a+b)/2, b
        );
    }
  }
  return tree;
}

int get_targets(float* targets, Node* tree, int l, int k) {
  return 0;
}

int get_interaction_list(float* sources, Node* tree, int l, int k) {
  return 0;
}

void barnes_hut(float* u, int n, Node* tree, int p) {
  // allocate temporary vectors to store source and target points during the computation
  float* sources = (float*) malloc(n * sizeof(float));
  float* targets = (float*) malloc(n * sizeof(float));

  int j = floor(log2(n));
  // Compute the weight (bottom to top)
  for (int l = j-1; l >= 0; l--) {
    for (int k = 0; k < pow(2, l); k++) {
      for (int m = 0; m < p; m++) {
        // TODO: compute w(L,k,m);
      }
    }
  }

  // Evaluate the potential (top to bottom)
  for (int l = 0; l < j; l++) {
    for (int k = 0; k < pow(2, l); k++) {
      // write relevant targets into temporary targets vector, where nt_lk gives the number of targets in T_{l,k}
      int nt_lk = get_targets(targets, tree, l, k); 
      for (int i = 0; i < nt_lk; i++) {
        float y = targets[i];
        for (int m = 0; m < p; m++) { // far-field
          // write relevant sources into temporary sources vector, where ns_lk gives the number of sources in the interaction list of T_{l,k}
          int ns_lk = get_interaction_list(sources, tree, l, k); 
          for (j = 0; j < ns_lk; i++) {
            float s = sources[j];
            // TODO: u(i) = u(i) + w(l,k,m) * Sm(xstar(l,s) - y);
          }
        }
        if (l == j) { // near-field
          // TODO: u(i) = u(i) + near-field(y(i), T(J,:));
        }
      }
    }
  }
}

int main() {
  int n = 100; // number of source / target points
  int max_pts = 10; // maximum number of points per neighborhood
  int p = 5; // number of terms in multipole expansion

  // draw and sort a set of n uniform random points in [0,1]
  float* points = (float*) malloc(n * sizeof(float));
  for (int i = 0; i < n; i++) points[i]  = ((double)rand()/RAND_MAX);
  std::sort(points, points+n);

  // build a binary tree subdiving our points
  Node* tree = build_tree(points, n, max_pts, 0, 1);

  // allocate the potential
  float* u = (float*) malloc(n * sizeof(float));

  // run barnes_hut
  barnes_hut(u, n, tree, p);

  // display solution at source / target points
  for (int i = 0; i < n; i++) printf("u(%.3f) = %.3f\n", points[i], u[i]);
}