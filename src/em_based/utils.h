#ifndef INCLUDE_pgmeans_h_
#define INCLUDE_pgmeans_h_

#include "base.h"

extern void copy_prop(double *new_props, double *props, int M);
extern void copy_mean(double *new_mus, double *mus, int p, int M);
extern void copy_fullcov(double *new_covs, double *covs, int p, int M);
extern void copy_diagcov(double *new_sigs, double *sigs, int p, int M);
extern int asc_double( const void *x, const void *y );

#endif //INCLUDE_pgmeans_h_
