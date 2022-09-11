#ifndef INCLUDE_em_gmm_h_
#define INCLUDE_em_gmm_h_

#include "base.h"

typedef struct params_list {
  int ind;
  double *mu; //p-dimensional mean vector
  double *cov; // p-dim covariance matrix
  double *sig; // p-dim variances
  double *chol; // p-dim Cholesky factor of the covariance matrix
  double *chol_diag; // Cholesky factor of the covariance matrix
  struct params_list *next;
} params_list;

typedef params_list * params_link;

extern void em_gmm(double *data, double *weights, double *mus, double *covs, double *props, int *pcov_type,
                          int *pn, int *pp, int *pM, double *pth, int *pcountf, double *logliks, int *pitmax, int *psilent);

extern params_link init_list_fullcov(double *mus, double *covs, int M, int p);
extern params_link init_list_diagcov(double *mus, double *sigs, int M, int p);
extern void free_list_fullcov(params_link params);
extern void free_list_diagcov(params_link params);

extern double mvnorm_fullcov(double *xi, params_link pl, int p);
extern double mvnorm_diagcov(double *xi, params_link pl, int p);

extern void estep(params_link params, double (*mvnorm)(double *, params_link, int),
           double *data, double *weights, double *props, int n, int p, int M);
extern double mixloglik(params_link params, double (*mvnorm)(double *, params_link, int),
                        double *data, double *props, int n, int p, int M);

extern double mstep_fullcov(params_link pl, double *data, double *weights, int n, int p);
extern double mstep_diagcov(params_link pl, double *data, double *weights, int n, int p);

#endif //INCLUDE_em_gmm_h_
