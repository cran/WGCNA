/* 
 * Common functions for fast calculations of correlations
 *
 */


/*
 * General notes about handling missing data, zero MAD etc:
 * The idea is that bicor should switch to cor whenever it is feasible, it helps, and it is requested:
 * (1) if median is NA, the mean would be NA as well, so there's no point in switching to Pearson
 * (2) In the results, columns and rows corresponding to input with NA means/medians are NA'd out.
 * (3) The convention is that if zeroMAD is set to non-zero, it is the index of the column in which MAD is
 *     zero (plus one for C indexing)
 */ 

#ifndef __corFunctions_utils_h__

#define __corFunctions_utils_h__

enum { noWarning, warnZeroMAD };

double median(double * x, size_t n, int copy, int * err);
double quantile(double * x, size_t n, double q, int copy, int * err);
double quantile_noCopy(double * x, size_t n, double q);
void testMedian(double *x, int * n, double * res);
void testQuantile(double *x, int *n, double *q, double *res);

void prepareColBicor(double * col, size_t nr, double maxPOutliers, int fallback,
                     int cosine,
                     double * res, size_t * nNAentries, 
                     int * NAmed, volatile int * zeroMAD,
                     double * aux, double * aux2);
void prepareColCor(double * x, size_t nr, int cosine, double * res, size_t * nNAentries, int * NAmean);

void prepareColCor_weighted(double * x, double * weights, 
      size_t nr, int cosine, double * res, size_t * nNAentries, int * NAmean);

int basic2variableCorrelation(
   double *xx, double *yy,
   size_t nr,
   double *res,
   int cosineX, int cosineY);

int basic2variableCorrelation_weighted(
   double *xx, double *yy, 
   double *wx, double *wy,
   size_t nr,
   double *res,
   int cosineX, int cosineY);

void * threadPrepColBicor(void * par);
void * threadPrepColCor(void * par);
void * threadSymmetrize(void * par);
void * threadSlowCalcBicor(void * par);
void * threadSlowCalcCor(void * par);
void * threadNAing(void * par);
void * threadSlowCalcBicor2(void * par);
void * threadSlowCalcCor2(void * par);
void * threadPrepColCor_weighted(void * par);
void * threadSlowCalcCor_weighted(void * par);
void * threadSlowCalcCor2_weighted(void * par);

#endif
