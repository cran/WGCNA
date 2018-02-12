/* For a symmatrix matrix A, calculate A'A
 * result[i,j] = sum_k A[k,i] A[k,j]
 */

#include "myMatrixMultiplication.h"

void squareSymmetricMatrix(const double * A, const size_t ncol,
                           double * result)
{
  for (size_t i=0; i<ncol; i++) for (size_t j=i; j<ncol; j++)
  {
    long double sum=0;
    for (size_t k=0; k<ncol; k++) sum+= A[i*ncol + k] * A[j*ncol + k];
    result[i*ncol + j] = result[j*ncol + i] = (double) sum;
  }
}
