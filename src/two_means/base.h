#ifndef INCLUDE_base_h_
#define INCLUDE_base_h_

#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <math.h>
#include <float.h>


#define RANDIN  GetRNGstate()
#define RANDOUT PutRNGstate()
#define UNIF unif_rand()

#define EPS 1e-4

#endif //INCLUDE_base_h_
