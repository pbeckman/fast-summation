#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <omp.h>

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

void build_tree(Node* tree, float* points, int n) {
  return;
}

int get_targets(float* targets, int l, int k) {
  return 0;
}

int get_interaction_list(float* sources, int l, int k) {
  return 0;
}

void barnes_hut(Node* tree, int n, int p) {
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
      int nt_lk = get_targets(targets, l, k); 
      for (int i = 0; i < nt_lk; i++) {
        float y = targets[i];
        for (int m = 0; m < p; m++) { // far-field
          // write relevant sources into temporary sources vector, where ns_lk gives the number of sources in the interaction list of T_{l,k}
          int ns_lk = get_interaction_list(sources, l, k); 
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
  int n = 100; // number of points
  int p = 5; // number of terms in multipole expansion

  float* points = (float*) malloc(n * sizeof(float));

  struct Node* tree = new Node(NULL, 0);

  build_tree(tree, points, n);

  barnes_hut(tree, n, p);
}