#ifndef INCLUDE_mixt4_h_
#define INCLUDE_mixt4_h_

#include "base.h"

/* ----------------- for test ------------------------ */
extern double mml(double *props, double loglik, double dmover2, int n, int M);
extern void estimate_mixprop(int comp, double *props, double n_m, double dmover2, int n, int M);
extern void lshift_props(double *props, int comp, int M);
extern void lshift_weights(double *weights, int comp, int n, int M, int Mmax);

extern params_link kill_comp(int comp, params_link pl);
/* --------------------------------------------------- */

#endif //INCLUDE_mixt4_h_
