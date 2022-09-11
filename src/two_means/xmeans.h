#ifndef INCLUDE_xmeans_h_
#define INCLUDE_xmeans_h_

#include "base.h"
#include "list.h"

/* ========================== for test ================================= */
extern double estimate_sigma(double *X, int p, int *c, int c_len, double *mu, int df);
extern double loglik(int n, int p, int k, double sigma);
extern boolean testCluster_bic(double *X, int p, int *ck, int nk, double *mk, double *vec1, double *vec2, list ic1, list ic2, double *bic_p, double *bic_c);
/* ===================================================================== */

#endif
