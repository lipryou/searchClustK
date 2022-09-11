#ifndef INCLUDE_twomeans_h_
#define INCLUDE_twomeans_h_

#include "base.h"
#include "list.h"

#define MAXIT 50

extern list *cluster_array;

extern void two_means(double *X, int p, int *ck, int nk, double *cnt1, double *cnt2, list ic1, list ic2);
extern void init_cluster_array(int *C, int n, int Kmax);
extern void free_cluster_array(int Kmax);

#endif
