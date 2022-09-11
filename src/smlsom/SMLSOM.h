#ifndef INCLUDE_SMLSOM_h_
#define INCLUDE_SMLSOM_h_

#include "base.h"

extern void uSMLSOM(double *data, double *moments, char **dtype, int *adjmatrix,
                    double *alphas, double *pbeta, double *radii, double *changes,
                    int *pn, int *pp, int *pM, int *pK, int *prlen, int *pchgbyns,
                    int *classes, double *logliks, double *consts, int *ptau, int *pdebug, int *psilent);
extern void mvn_SMLSOM(double *data, double *mu1, double *mu2, int *adjmatrix,
                       double *alphas, double *pbeta, double *radii, double *changes,
                       int *pn, int *pp, int *pM, int *prlen, int *pchgbyns,
                       int *classes, double *logliks, int *ptau, int *pdebug, int *psilent);

#endif //INCLUDE_SMLSOM_h_
