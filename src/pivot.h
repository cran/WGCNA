#ifndef __pivot_h__

#define __pivot_h__

#include <R.h>
#include <Rinternals.h>

double vMax(double * v, size_t len);
double vMin(double * v, size_t len);
double pivot(double * v, size_t len, double target);

typedef struct
{
  double val;
  size_t index;
} orderStructure;

int compareOrderStructure(const orderStructure * os1, const orderStructure * os2);
void qorder_internal(double * x, size_t n, orderStructure * os);

SEXP qorder(SEXP data);

#endif


