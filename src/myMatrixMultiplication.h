
#ifndef __myMatrixMultiplication_h__
#define __myMatrixMultiplication_h__

#include <R.h>
#include <Rinternals.h>

/* For a symmatrix matrix A, calculate A'A
 * result[i,j] = sum_k A[k,i] A[k,j]
 */

void squareSymmetricMatrix(const double * A, const size_t ncol, double * result);

#endif

