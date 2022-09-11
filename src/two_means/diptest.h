#ifndef INCLUDE_diptest_h_
#define INCLUDE_diptest_h_

#include "base.h"

double testCluster_unimodality(double *D, int nk, int *ck, double a, double vthd, int debug);
double calculate_score(double *D, int nk, int *ck, int nn_len, int *nn, double *qdips, double vthd, int debug);

#endif //INCLUDE_diptest_h_
