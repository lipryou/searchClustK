#ifndef INCLUDE_pgmeans_h_
#define INCLUDE_pgmeans_h_

#include "base.h"

/* -------------- for test ---------------------*/
extern void create_projection(double *P, int p);
extern void projection_mean(double *P, double *mus, double *prj_mus, int p, int M);
extern void projection_fullcov(double *P, double *covs, double *prj_sigmas, int p, int M);
extern void projection_diagcov(double *P, double *sigs, double *prj_sigmas, int p, int M);
extern void projection_x(double *P, double *data, double *prj_x, int n, int p);
extern void add_onecluster_props(double *props, double *new_props, int M);
extern void add_onecluster_mean(double *xi, double *mus, double *new_mus, int p, int M);
extern void add_onecluster_fullcov(double *covs, double *new_covs, int p, int M);
extern void add_onecluster_diagcov(double *sigs, double *new_sigs, int p, int M);
extern double pmixnorm(double x, int M, double *pi, double *mu, double *sigma);
/* -------------------------------------------  */

#endif //INCLUDE_pgmeans_h_
