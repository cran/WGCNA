// "Parallel" quantile: for a list of numeric vectors or arrays, calculate a given quantile of each vector
// containing element 'i' of each component of the input list.
// NA's are treated as last.
//

#ifndef __parallelQuantile_stdC_h__

#define __parallelQuantile_stdC_h__

SEXP parallelQuantile(SEXP data_s, SEXP prob_s);
SEXP parallelMean(SEXP data_s, SEXP weight_s);
SEXP parallelMin(SEXP data_s);
SEXP minWhich_call(SEXP matrix_s, SEXP rowWise_s);
SEXP quantileC_call(SEXP data_s, SEXP q_s);
SEXP rowQuantileC_call(SEXP data_s, SEXP q_s);

#endif


