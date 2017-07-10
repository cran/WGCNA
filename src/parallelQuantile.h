// "Parallel" quantile: for a list of numeric vectors or arrays, calculate a given quantile of each vector
// containing element 'i' of each component of the input list.
// NA's are treated as last.
//

#ifndef __parallelQuantile_h__

#define __parallelQuantile_h__

#include <Rcpp.h>
#include <R.h>
#include <iostream>
#include "array.h"

#include "corFunctions-common.h"

using namespace std;
using namespace Rcpp;

RcppExport SEXP parallelQuantile(SEXP data_s, SEXP prob_s);
RcppExport SEXP parallelMean(SEXP data_s, SEXP weight_s);
RcppExport SEXP parallelMin(SEXP data_s);
RcppExport SEXP minWhich_call(SEXP matrix_s, SEXP rowWise_s);
RcppExport SEXP quantileC_call(SEXP data_s, SEXP q_s);
RcppExport SEXP rowQuantileC_call(SEXP data_s, SEXP q_s)

void quantileC(double * data, int *nrow, int * ncol, double * q, double * res);
void rowQuantileC(double * data, int *nrow, int * ncol, double * q, double * res);

#endif


