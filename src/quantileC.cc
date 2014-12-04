#include <iostream>

#include "array.h"
#include <R.h>


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


