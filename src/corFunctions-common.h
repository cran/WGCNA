/* 
 * Common functions for fast calculations of correlations
 *
 */


#ifndef __corFunctions_common_h__
#define __corFunctions_common_h__

#define LDOUBLE 	long double

#include <R.h>
#include <Rinternals.h>

//#include "pivot.h"

enum { noWarning, warnZeroMAD };

double median(double * x, size_t n, int copy, int * err);

double quantile(double * x, size_t n, double q, int copy, int * err);
double quantile_noCopy(double * x, size_t n, double q);


void testMedian(double *x, int * n, double * res);

void prepareColBicor(double * col, size_t nr, double maxPOutliers, int fallback, int cosine, double * res, 
                     size_t * nNAentries, int * NAmed, volatile int *zeroMAD, double * aux, double * aux2);

void prepareColCor(double * x, size_t nr, int cosine, double * res, size_t * nNAentries, int * NAmean);

#endif
