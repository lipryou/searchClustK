#ifndef INCLUDE_gmeans_h_
#define INCLUDE_gmeans_h_

#include "base.h"
#include "list.h"

extern void projection(double *X, int p, int *ck, int nk, double *vec1, double *vec2, double *prj_x);
extern boolean testCluster_ADtest(double *X, int p, int *ck, int nk, double *vec1, double *vec2, double *d, double cv);
extern double ADstatic(double *z, int n);
extern void covariance(double *X, int p, int *ck, int nk, double *mk, double **cov);
extern double eigen_power(double **cov, int p, double *mk,double *vec1, double *vec2);
extern double first_prcomp(double *X, int p, int *ck, int nk, double *mk, double **cov, double *vec1, double *vec2);
extern double normalCDF(double value);

#endif
