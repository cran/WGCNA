// "Parallel" quantile: for a list of numeric vectors or arrays, calculate a given quantile of each vector
// containing element 'i' of each component of the input list.
// NA's are treated as last.

#include <Rcpp.h>
#include <R.h>
#include <iostream>
#include "array.h"

#include "corFunctions-common.h"

using namespace std;
using namespace Rcpp;

/*
 *
 * Main function parallelQuantile. I will assume that the R code made the necessary checks on the data being
 * non-empty and all components having the same length.
 *
 */

RcppExport SEXP parallelQuantile(SEXP data_s, SEXP prob_s)
{
  BEGIN_RCPP
  List data_lst = List(data_s);
  NumericVector prob_v = NumericVector(prob_s);
  double prob = prob_v[0];

  size_t nSets = data_lst.size();

  // cout << "nSets: " << nSets << endl;

  vector <NumericVector> data(nSets);
  data.clear();

  for (size_t i=0; i<nSets; i++) 
  {
    NumericVector v(data_lst[i]);
    // cout << "Set: " << i << ", length: " << v.length() << endl;
    data.push_back(v);
  }

  size_t nElements = data[0].size();

  // cout << "nElements: " << nElements << endl;

  NumericVector quantiles(nElements);

  double * vec = new double[nSets];

  int err = 0;

  for (size_t i=0; i<nElements; i++)
  {
    for (size_t set=0; set<nSets; set++)
      vec[set] = data[set][i];

    quantiles[i] = quantile_noCopy(vec, nSets, prob);
  }
  delete[] vec;
  // cout << "nElements: " << nElements << endl;

  quantiles.attr("dim") = data[0].attr("dim");
  return(quantiles);
 
  END_RCPP 
}




/*===============================================================================================
 *
 * parallel weighted mean
 *
 *===============================================================================================*/

RcppExport SEXP parallelMean(SEXP data_s, SEXP weight_s)
{
  BEGIN_RCPP
  List data_lst = List(data_s);
  NumericVector weight = NumericVector(weight_s);

  size_t nSets = data_lst.size();
  if (nSets!=weight.length())
     throw "Compiled parallelMean: Length of 'weights' must equal length of 'data'.";

  // cout << "nSets: " << nSets << endl;

  vector <NumericVector> data(nSets);
  data.clear();

  for (size_t i=0; i<nSets; i++) 
  {
    NumericVector v(data_lst[i]);
    // cout << "Set: " << i << ", length: " << v.length() << endl;
    data.push_back(v);
  }

  size_t nElements = data[0].size();

  // cout << "nElements: " << nElements << endl;

  NumericVector mean(nElements);

  for (size_t i=0; i<nElements; i++)
  {
    double sum = 0, count = 0;
    for (size_t set=0; set<nSets; set++) if (!ISNAN(data[set][i]) && !ISNAN(weight[set]))
    {
       sum += data[set][i] * weight[set];
       count += weight[set];
    }
    if (count==0)
       mean[i] = NA_REAL;
    else 
       mean[i] = sum/count;
  }
  // cout << "nElements: " << nElements << endl;

  mean.attr("dim") = data[0].attr("dim");
  return(mean);
 
  END_RCPP 
}

/*===============================================================================================
 *
 * parallel min and index of minimum
 *
 *===============================================================================================*/

RcppExport SEXP parallelMin(SEXP data_s)
{
  BEGIN_RCPP
  List data_lst = List(data_s);

  size_t nSets = data_lst.size();
  // cout << "nSets: " << nSets << endl;

  vector <NumericVector> data(nSets);
  data.clear();

  for (size_t i=0; i<nSets; i++) 
  {
    NumericVector v(data_lst[i]);
    // cout << "Set: " << i << ", length: " << v.length() << endl;
    data.push_back(v);
  }

  size_t nElements = data[0].size();

  // cout << "nElements: " << nElements << endl;

  NumericVector minv(nElements), which(nElements);

  for (size_t i=0; i<nElements; i++)
  {
    double min1 = NA_REAL, index1 = NA_REAL;
    for (size_t set=0; set<nSets; set++) if (!ISNAN(data[set][i]))
      if (ISNAN(min1) || (min1 > data[set][i])) 
      { 
        min1 = data[set][i]; 
        index1 = set; 
      } 

    minv[i] = min1;
    which[i] = index1+1;
  }
  // cout << "nElements: " << nElements << endl;

  minv.attr("dim") = data[0].attr("dim");
  which.attr("dim") = data[0].attr("dim");
  List out;
  out["min"] = minv;
  out["which"] = which;
  return(out);
 
  END_RCPP 
}

// if rowWise is non-zero, the min and which will be by rows, otherwise by columns.
//

RcppExport SEXP minWhich_call(SEXP matrix_s, SEXP rowWise_s)
{
  BEGIN_RCPP
  NumericMatrix matrix(matrix_s);
  size_t nrows = matrix.nrow(), ncols = matrix.ncol();
  IntegerVector rowWise(rowWise_s);

  size_t nouter, ninner, outerStride, innerStride;

  if (rowWise[0] == 0)
  {
    nouter = ncols;
    ninner = nrows;
    outerStride = nrows;
    innerStride = 1;
  } else {
    nouter = nrows;
    ninner = ncols;
    outerStride = 1;
    innerStride = nrows;
  }

  NumericVector min(nouter), which(nouter);

  for (size_t i=0; i<nouter; i++)
  {
    size_t ind = i * outerStride;
    double curmin = NA_REAL;
    double curwhich = NA_REAL;
    for (size_t j=0; j<ninner; j++)
    {
      if ( ISNAN(curmin) || (!ISNAN(matrix[ind]) && (matrix[ind] < curmin))) 
      { 
         curmin = matrix[ind]; 
         curwhich = (double) j+1; 
      }
      ind+= innerStride;
    }
    min[i] = curmin;
    which[i] = curwhich;
  }

  List out;
  out["min"] = min;
  out["which"] = which;
  return(out);
  END_RCPP
}


/* ==================================================================================================
 *
 * Older functions for column and row-wise quantiles, and their .Call-suitable wrappers.
 *
 ====================================================================================================*/


extern "C" {

void quantileC(double * data, int *nrow, int * ncol, double * q, double * res)
{
  try
  {
    int nr = *nrow, nc = *ncol;
    dArray d;

    d.wrap(data, nr, nc);

    if ((*q<0) || (*q>1)) 
      throw(Exception(string("quantileC: given quantile is out of range 0 to 1.")));

    dArray quant;

    quant.wrap(res, nc);

    d.colQuantile(*q, quant);

  } catch (Exception & err)
  {
    Rprintf("Error in (compiled code) quantileC: %s\n", err.what().c_str());
  }
}

void rowQuantileC(double * data, int *nrow, int * ncol, double * q, double * res)
{
  try
  {
    int nr = *nrow, nc = *ncol;
    dArray d;

    d.wrap(data, nr, nc);

    if ((*q<0) || (*q>1)) 
      throw(Exception(string("quantileC: given quantile is out of range 0 to 1.")));

    dArray quant;

    quant.wrap(res, nr);

    d.rowQuantile(*q, quant);

  } catch (Exception & err)
  {
    Rprintf("Error in (compiled code) quantileC: %s\n", err.what().c_str());
  }

}

} // extern "C"


RcppExport SEXP quantileC_call(SEXP data_s, SEXP q_s)
{
  BEGIN_RCPP

  NumericMatrix data(data_s);
  // The dimensions below are defined as int to preserve compatibility with the original .C functions.
  int nrows = data.nrow(), ncols = data.ncol();

  NumericVector q(q_s);
  NumericVector res(ncols);
  quantileC(&data[0], &nrows, &ncols, &q[0], &res[0]);

  return(res);
  END_RCPP
}

RcppExport SEXP rowQuantileC_call(SEXP data_s, SEXP q_s)
{
  BEGIN_RCPP

  NumericMatrix data(data_s);
  // The dimensions below are defined as int to preserve compatibility with the original .C functions.
  int nrows = data.nrow(), ncols = data.ncol();

  NumericVector q(q_s);
  NumericVector res(nrows);
  rowQuantileC(&data[0], &nrows, &ncols, &q[0], &res[0]);

  return(res);
  END_RCPP
}

