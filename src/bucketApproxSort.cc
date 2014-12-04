

/*******************************************************************************************************
 *
 * Approximate sort using the bucket method
 *
 * Input min and max value. Split the range into N intervals. Keep a vector of indices of elements in each
 * interval. Assign each element to an internal. Optionally sort each interval using qsort; return order.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <math.h>
#include "exceptions.h"

#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>

#include "pivot.h"

using namespace std;

typedef vector <int> interval_int;
typedef vector <double> interval_double;

template <typename indexType> void bucketOrder(
           double * x, size_t n,
           double xmin, double xmax, size_t nIntervals,
           int exact,
           indexType * order)
{

  // Rprintf("Entering bucketOrder. Size of data: %ld, nIntervals: %ld\n", n, nIntervals);
  typedef vector <indexType> interval;
  interval * intervalContent = new interval [nIntervals];

  for (size_t t=0; t<nIntervals; t++) intervalContent[t].clear();

  xmax = xmax + (xmax-xmin)/(5*nIntervals);

  for (size_t i=0; i<n; i++)
  {
    size_t index = (size_t) ( (x[i]-xmin)/(xmax-xmin) * nIntervals);
    intervalContent[index].push_back((indexType) i+1);
  }

  if (exact) Rprintf("Exact sorting not implemented yet.\n");

  size_t index = 0;
  for (size_t t=0; t<nIntervals; t++) 
    for (typename interval::iterator i=intervalContent[t].begin(); i!=intervalContent[t].end(); i++) 
    {
       order[index] = *i;
       index++;
    }

  delete [] intervalContent;

}

extern "C" {

SEXP bucketOrder_R(SEXP data, SEXP min, SEXP max, SEXP nIntervals, SEXP exact)
{
  R_xlen_t n = Rf_xlength(data);

  double * x = REAL(data);
  double xmin = REAL(min)[0];
  double xmax = REAL(max)[0];


  size_t nInt_c = (size_t ) INTEGER(nIntervals)[0];
  int exact_c = INTEGER(exact)[0];

  SEXP ans;
  if (n<(size_t) 0x80000000)
  {
    // Rprintf("..returning integer order.\n");
    PROTECT (ans = allocVector(INTSXP, n));
    bucketOrder <int> (x, (size_t) n, xmin, xmax, nInt_c, exact_c, INTEGER(ans));
  } else {
    // Rprintf("..returning floating point (double) order.\n");
    PROTECT (ans = allocVector(REALSXP, n));
    bucketOrder <double> (x, (size_t) n, xmin, xmax, nInt_c, exact_c, REAL(ans));
  }

  UNPROTECT(1);
  return ans;
}

} // #extern "C"

