

/* 

Copyright (C) 2008 Peter Langfelder; parts based on R by R Development team

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

*/

int uselessFunction1()
{
  int a = 3;
  return a;
}

#ifndef WITH_THREADS

// This file should only be included if threads are not available.



#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/libextern.h>
#define LDOUBLE 	long double

#include "pivot.h"
#include "corFunctions-common.h"

//===================================================================================================

// Pearson correlation of a matrix with itself.
// This one uses matrix multiplication in BLAS to speed up calculation when there are no NA's

//===================================================================================================


void cor1Fast(double * x, int * nrow, int * ncol, double * quick, double * result, int *nNA, int * err, 
          int * nThreads,
          int * verbose, int * indent)
{
  int nr = *nrow, nc = *ncol;

  char          spaces[2* *indent+1];

  for (int i=0; i<2* *indent; i++) spaces[i] = ' ';
  spaces[2* *indent] = '\0';

  int maxDiffNA = (int) (*quick * nr);

  *nNA = 0;

  // const double asymptCorr = 1.4826, qnorm75 = 0.6744898;

  double *xx, *yy;

  // Allocate space for various variables

  double * multMatX;
  int * nNAentries, *NAmean;

  // This matrix will hold preprocessed entries that can be simply multiplied together to get the
  // numerator

  if ( (multMatX = malloc(nc*nr * sizeof(double)))==NULL )
  {
    *err = 1;
    Rprintf("cor1: memmory allocation error. If possible, please decrease block size.\n");
    return;
  }

  // Number of NA entries in each column

  if ( (nNAentries = malloc(nc * sizeof(int)))==NULL )
  {
    free(multMatX);
    *err = 1;
    Rprintf("cor1: memmory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  // Flag indicating whether the mean of each column is NA

  if ( (NAmean = malloc(nc * sizeof(int)))==NULL )
  {
    free(nNAentries); free(multMatX);
    *err = 1;
    Rprintf("cor1: memmory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  // Calculate multMatX = (x - mean(x))/sqrt(sum(x - mean(x)))

  for (int i = 0; i < nc; i++)
  {
    prepareColCor(x + i*nr, nr, multMatX + i*nr, nNAentries + i, NAmean + i);
  //  if (NAmean[i]) Rprintf("have a NA mean in column %d.\n", i);
  }

  // Rprintf("NAmean values: ");
  // for (int i = 0; i < nc; i++) Rprintf("%d, ", NAmean[i]);
  // Rprintf("\n");

  // The main loop is actually a matrix multiplication

  double alpha = 1.0, beta = 0.0;
  dsyrk_("L", "T", ncol, nrow, & alpha, multMatX, nrow, & beta, result, ncol);

  // Here I need to recalculate results that have NA's in them.

  int nSlow = 0;
  double vx=0, vy=0;

  // Rprintf("nNAentries values: ");
  // for (int i = 0; i < nc; i++) Rprintf("%d, ", nNAentries[i]);
  // Rprintf("\n");
  //
  for (int i=0; i<nc; i++) if ( NAmean[i] == 0 )
  {
     for (int j=i+1; j<nc; j++) if ( (NAmean[j]==0) && 
                                     ((nNAentries[i] > maxDiffNA) || ( nNAentries[j] > maxDiffNA)) )
     {
        // Rprintf("Recalculating row %d and column %d, column size %d\n", i, j, nr);
        xx = x + i * nr; yy = x + j * nr;
        LDOUBLE sumxy = 0, sumx = 0, sumy = 0, sumxs = 0, sumys = 0;
        int count = 0;

        for (int k=0; k<nr; k++)
        {
           vx = *xx; vy = *yy;
           if (!ISNAN(vx) && !ISNAN(vy))
           {
             count ++;
             sumxy += vx * vy;
             sumx += vx;
             sumy += vy;
             sumxs += vx*vx;
             sumys += vy*vy;
           }
           xx++; yy++;
        }
        if (count==0)
        {
            result[i*nc + j] = NA_REAL;
            (*nNA)++;
        } else {
            result[i*nc + j] = (double) ( (sumxy - sumx * sumy/count)/
                                sqrtl( (sumxs - sumx*sumx/count) * (sumys - sumy*sumy/count) ) );
        }
        // result[j*nc + i] = result[i*nc + j];
        nSlow++;

     }
  }

  if (*verbose) Rprintf("%s Fraction of slow calculations: %f\n", spaces, 
                         ( (double) nSlow*2 ) / (nc*(nc-1)) );
  
  // Symmetrize the result

  for (int i=0; i<nc; i++) if (NAmean[i] == 0)
  {
    double * resx = result + i*nc + i;
    for (int j=i; j<nc; j++) 
    {
      if (NAmean[j] == 0) // Start at j=i to check for values greater than 1.
      {
        if (!ISNAN(*resx))
        {
           if (*resx > 1.0) *resx = 1.0;
           if (*resx < -1.0) *resx = -1.0;
        }
        result[j*nc + i] = *resx;
      }
      resx ++;
    }
  } else {
     for (int j=0; j<nc; j++)
     {
        result[i*nc + j] = NA_REAL;
        result[j*nc + i] = NA_REAL;
     }
  }

  free(NAmean);
  free(nNAentries);
  free(multMatX);
}

//===================================================================================================

// Two-variable pearson correlation. Basically the same as cor1, just must calculate the whole matrix.

//===================================================================================================

void corFast(double * x, int * nrow, int * ncolx, double * y, int * ncoly,
           double * quick, double * result, int *nNA, int * err,
           int * nThreads,
           int * verbose, int * indent)
{
  int nr = *nrow, ncx = *ncolx, ncy = *ncoly;

  char          spaces[2* *indent+1];
  for (int i=0; i<2* *indent; i++) spaces[i] = ' ';
  spaces[2* *indent] = '\0';

  int maxDiffNA = (int) (*quick * nr);

  *nNA = 0;

  // Allocate space for auxiliary variables, NA counts and flags whether the corresponding means or
  // medians are NA

  double * multMatX, * multMatY;	
  int * nNAentriesX, * nNAentriesY;
  int *NAmeanX, *NAmeanY;

  if ( (multMatX = malloc(ncx*nr * sizeof(double)))==NULL )
  {
    *err = 1;
    Rprintf("cor: memmory allocation error. If possible, please decrease block size.\n");
    return;
  }

  if ( (multMatY = malloc(ncy*nr * sizeof(double)))==NULL )
  {
    free(multMatX);
    *err = 1;
    Rprintf("cor: memmory allocation error. If possible, please decrease block size.\n");
    return;
  }

  if ( (nNAentriesX = malloc(ncx * sizeof(int)))==NULL )
  {
    free(multMatY); free(multMatX);
    *err = 1;
    Rprintf("cor: memmory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  if ( (nNAentriesY = malloc(ncy * sizeof(int)))==NULL )
  {
    free(nNAentriesX); free(multMatY); free(multMatX);
    *err = 1;
    Rprintf("cor: memmory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  if ( (NAmeanX = malloc(ncx * sizeof(int)))==NULL )
  {
    free(nNAentriesY); free(nNAentriesX); free(multMatY); free(multMatX);
    *err = 1;
    Rprintf("cor: memmory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  if ( (NAmeanY = malloc(ncy * sizeof(int)))==NULL )
  {
    free(NAmeanX); free(nNAentriesY); free(nNAentriesX); free(multMatY); free(multMatX);
    *err = 1;
    Rprintf("cor: memmory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  // Calculate multMatX = (x - median(x))/(9.0 * qnorm75 * mad(x)); median and mad are columnwise.
  // or multMatX = x-mean(x) if robustX is 0.

  int NAmed = 0;
  for (int i = 0; i < ncx; i++)
  {
    prepareColCor(x + i*nr,  nr, multMatX + i*nr, nNAentriesX + i, &NAmed);
    // if (NAmed) Rprintf("median is NA in column %d of x\n", i);
    NAmeanX[i] = NAmed;
  }

  // Rprintf("NAmeanX: %d\n", NAmeanX[0])

  for (int j = 0; j < ncy; j++)
  {
    prepareColCor(y + j*nr,  nr, multMatY + j*nr, nNAentriesY + j, &NAmed);
    // if (NAmed) Rprintf("median is NA in column %d of y\n", j);
    NAmeanY[j] = NAmed;
  }
  // Rprintf("NAmeanY: %d\n", NAmeanY[0])


  double alpha = 1.0, beta = 1.0;
  dgemm_("T", "N", ncolx, ncoly, nrow, & alpha, multMatX, nrow, multMatY, nrow, & beta, result, ncolx);

  int nSlow = 0;
  double *xx, *yy, vx, vy;

  double * resx = result;
  for (int j=0; j<ncy; j++) for (int i=0; i<ncx; i++) 
  {
     if ((NAmeanX[i]>0) || (NAmeanY[j] > 0))
     {
       *resx = NA_REAL;
     } else 
     {
        if ( (nNAentriesX[i] > maxDiffNA) || ( nNAentriesY[j] > maxDiffNA))
        {
           xx = x + i * nr; yy = y + j * nr;
           LDOUBLE sumxy = 0, sumx = 0, sumy = 0, sumxs = 0, sumys = 0;
           int count = 0;
      
           for (int k=0; k<nr; k++)
           {
                 vx = *xx; vy = *yy;
                 if (!ISNAN(vx) && !ISNAN(vy))
                 {
                   count ++;
                   sumxy += vx * vy;
                   sumx += vx;
                   sumy += vy;
                   sumxs += vx*vx;
                   sumys += vy*vy;
                 }
                 xx++; yy++;
           }
           if (count==0)
           {
               *resx = NA_REAL;
               (*nNA)++;
           } else {
               *resx = (double) ( (sumxy - sumx * sumy/count)/
                                   sqrtl( (sumxs - sumx*sumx/count) * (sumys - sumy*sumy/count) ) );
           }
           nSlow++;
        }
        if ( !ISNAN(*resx) )
        {
          if (*resx > 1.0) *resx = 1.0;
          if (*resx < -1.0) *resx = -1.0;
        }
     }
     resx ++;
  }

  if (*verbose) Rprintf("%s Fraction of slow calculations: %f\n", spaces,
                         ( (double) nSlow) / (ncx*ncy) );
 
  free(NAmeanY);
  free(NAmeanX);
  free(nNAentriesY);
  free(nNAentriesX);
  free(multMatY);
  free(multMatX);
}

//===================================================================================================

// bicorrelation of a matrix with itself.
// This one uses matrix multiplication in BLAS to speed up calculation when there are no NA's

//===================================================================================================


void bicor1Fast(double * x, int * nrow, int * ncol, double * maxPOutliers, 
            double * quick, int * fallback, double * result, int *nNA, int * err, 
            int * nThreads,
            int * verbose, int * indent)
{
  int nr = *nrow, nc = *ncol;

  char          spaces[2* *indent+1];

  for (int i=0; i<2* *indent; i++) spaces[i] = ' ';
  spaces[2* *indent] = '\0';

  int maxDiffNA = (int) (*quick * nr);

  *nNA = 0;

  // const double asymptCorr = 1.4826, qnorm75 = 0.6744898;

  double *xx, *yy, *xxx, *yyy, *xx2, *yy2;

  if ( (xx=malloc(6*nr * sizeof(double)))==NULL)
  {
    *err = 1;
    Rprintf("bicor1: memmory allocation error. The needed block is very small... suspicious.\n");
    return;
  }

  yy = xx + nr;

  xxx = xx + 2*nr;
  yyy = yy + 2*nr;

  xx2 = xx + 4*nr;
  yy2 = yy + 4*nr;

  // Allocate space for various variables

  double * multMat;
  int * nNAentries, *NAmed;

  // This matrix will hold preprocessed entries that can be simply multiplied together to get the
  // numerator

  if ( (multMat = malloc(nc*nr * sizeof(double)))==NULL )
  {
    free(xx);
    *err = 1;
    Rprintf("bicor1: memmory allocation error. If possible, please decrease block size.\n");
    return;
  }

  // Number of NA entries in each column

  if ( (nNAentries = malloc(nc * sizeof(int)))==NULL )
  {
    free(multMat); free(xx);
    *err = 1;
    Rprintf("bicor1: memmory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  // Flag indicating whether the median of each column is NA

  if ( (NAmed = malloc(nc * sizeof(int)))==NULL )
  {
    free(nNAentries); free(multMat); free(xx);
    *err = 1;
    Rprintf("bicor1: memmory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  // Calculate multMat = (x - median(x))/(9.0 * qnorm75 * mad(x)); 
  // median and mad are columnwise. multMat is also
  // normalized so that sum of each column is 1.

  int zeroMAD;
  int pearson = 0;
  if (*fallback < 3)
  {
    for (int i = 0; i < nc; i++)
      prepareColBicor(x + i*nr, nr, *maxPOutliers, *fallback, multMat + i*nr, 
                      nNAentries + i, NAmed + i, &zeroMAD, xx, xxx);
  } else {
    int doFallback = 0;
    for (int i = 0; i < nc; i++)
    {
      prepareColBicor(x + i*nr, nr, *maxPOutliers, *fallback, multMat + i*nr, 
                      nNAentries + i, NAmed + i, &zeroMAD, xx, xxx);
      if (zeroMAD > 0) { doFallback = i+1; i = nc; }
    }
    if (doFallback)
    {
      if (*verbose > 0)
        Rprintf("Warning in bicor: zero-MAD column %d in variable x.\n   %s", doFallback, 
                "         Switching to Pearson correlation for variable x.\n\n");
      pearson = 1;
      for (int i = 0; i < nc; i++)
        prepareColCor(x + i*nr, nr, multMat + i*nr, nNAentries + i, NAmed + i);
    }
  }


  // The main loop is actually a matrix multiplication

  double alpha = 1.0, beta = 0.0;
  dsyrk_("L", "T", ncol, nrow, & alpha, multMat, nrow, & beta, result, ncol);


  // Here I need to recalculate results that have NA's in them.

  int nSlow = 0;

  int temp = 0;
  int fbx = *fallback;
  if (pearson) fbx = 4;

  // Rprintf("bicor(x): fbx = %d\n", fbx);
  // Rprintf("..starting remedial calculations. maxDiffNA = %d\n", maxDiffNA);
  // Rprintf("..nNAentries:\n");
  // for (int i=0; i<nc; i++) Rprintf(" %d,", nNAentries[i]);
  // Rprintf("\n");

  for (int i=0; i<nc; i++) if ( NAmed[i] == 0) 
  {
     for (int j=i+1; j<nc; j++) if ( (NAmed[j] == 0) && ((nNAentries[i] > maxDiffNA) 
                                                   || ( nNAentries[j] > maxDiffNA)))
     {
        memcpy((void *)xx, (void *)(x + i*nr), nr * sizeof(double));
        memcpy((void *)yy, (void *)(x + j*nr), nr * sizeof(double));

        int nNAx = 0, nNAy = 0;    
        for (int k=0; k<nr; k++)
        {
           if (ISNAN(xx[k])) yy[k] = NA_REAL;
           if (ISNAN(yy[k])) xx[k] = NA_REAL;
           if (ISNAN(xx[k])) nNAx++;
           if (ISNAN(yy[k])) nNAy++;
        }

        int NAx = 0, NAy = 0;
        double *xx3, *yy3;

        // Rprintf("nNAx: %d, nNAy: %d\n", nNAx, nNAy);

        if ((nNAx - nNAentries[i] > maxDiffNA) || (nNAy-nNAentries[j] > maxDiffNA))
        {
            // Rprintf("Recalculating row %d and column %d\n", i, j);
            // must recalculate the auxiliary variables for both columns
            if (nNAx - nNAentries[i] > maxDiffNA)
               {
                  prepareColBicor(xx, nr, *maxPOutliers, fbx, xxx, &temp, &NAx, &zeroMAD, xx2, yy2);
                  xx3 = xxx;
               }
               else
                  xx3 = multMat + i * nr;
            if (nNAy-nNAentries[j] > maxDiffNA)
               {    
                  prepareColBicor(yy, nr, *maxPOutliers, fbx, yyy, &temp, &NAy, &zeroMAD, xx2, yy2);
                  yy3 = yyy;
               }
               else
                  yy3 = multMat + j * nr;
            if (NAx + NAy==0)
            {
               //LDOUBLE sumxy = 0;
               double  sumxy = 0;
               int count = 0;
               for (int k=0; k<nr; k++)
               {
                 double vx = *xx3, vy = *yy3;
                 // Rprintf("i: %d, j: %d, k: %d: vx: %e, vy: %e\n", i,j,k,vx,vy);
                 if (!ISNAN(vx) && !ISNAN(vy))
                 {
                   sumxy += vx * vy; 
                   count++;
                 }
                 xx3++; yy3++;
               }
               if (count==0) 
               {
                  result[i*nc + j] = NA_REAL; 
                  (*nNA)++;
               } else {
                  result[i*nc + j] = (double) sumxy;
               }
             } else {
                result[i*nc + j] = NA_REAL; 
                (*nNA)++;
             }
             nSlow++;
        }
     }
  }

  // Symmetrize the result

  for (int i=0; i<nc; i++) if (NAmed[i] == 0)
  {
    double * resx = result + i*nc + i;
    for (int j=i; j<nc; j++) 
    {
      if (NAmed[j] == 0) // Start at j=i to check for values greater than 1.
      {
        if (!ISNAN(*resx))
        {
           if (*resx > 1.0) *resx = 1.0;
           if (*resx < -1.0) *resx = -1.0;
        }
        result[j*nc + i] = *resx;
      }
      resx ++;
    }
  } else
  {
     for (int j=0; j<nc; j++)
     {
        result[i*nc + j] = NA_REAL;
        result[j*nc + i] = NA_REAL;
     }
  }

  if (*verbose) Rprintf("%s Fraction of slow calculations: %f\n", spaces, 
                         ( (double) nSlow*2 ) / (nc*(nc-1)) );
  
  free(NAmed);
  free(nNAentries);
  free(multMat);
  free(xx);
}

//===================================================================================================

// Two-variable bicorrelation. Basically the same as bicor1, just must calculate the whole matrix.
// If robustX,Y is zero, the corresponding variable will be treated as in pearson correlation.

//===================================================================================================

void bicorFast(double * x, int * nrow, int * ncolx, double * y, int * ncoly,
           int * robustX, int * robustY, double * maxPOutliers, 
           double * quick, int * fallback, double * result, int *nNA, int * err, 
           int * nThreads,
           int * verbose, int * indent)
{
  int nr = *nrow, ncx = *ncolx, ncy = *ncoly;

  char          spaces[2* *indent+1];
  for (int i=0; i<2* *indent; i++) spaces[i] = ' ';
  spaces[2* *indent] = '\0';

  int maxDiffNA = (int) (*quick * nr);

  *nNA = 0;

  // const double asymptCorr = 1.4826, qnorm75 = 0.6744898;

  double *xx, *yy, *xxx, *yyy, *xx2, *yy2;

  if ( (xx=malloc(6*nr * sizeof(double)))==NULL)
  {
    *err = 1;
    Rprintf("bicor: memmory allocation error. The needed block is very small... suspicious.\n");
    return;
  }

  yy = xx + nr;

  xxx = xx + 2*nr;
  yyy = yy + 2*nr;

  xx2 = xx + 4*nr;
  yy2 = yy + 4*nr;

  // Allocate space for auxiliary variables, NA counts and flags whether the corresponding means or
  // medians are NA

  double * multMatX, * multMatY;	
  int * nNAentriesX, * nNAentriesY;
  int *NAmeanX, *NAmeanY;

  if ( (multMatX = malloc(ncx*nr * sizeof(double)))==NULL )
  {
    free(xx);
    *err = 1;
    Rprintf("bicor: memmory allocation error. If possible, please decrease block size.\n");
    return;
  }

  if ( (multMatY = malloc(ncy*nr * sizeof(double)))==NULL )
  {
    free(xx); free(multMatX);
    *err = 1;
    Rprintf("bicor: memmory allocation error. If possible, please decrease block size.\n");
    return;
  }

  if ( (nNAentriesX = malloc(ncx * sizeof(int)))==NULL )
  {
    free(multMatY); free(multMatX); free(xx);
    *err = 1;
    Rprintf("bicor: memmory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  if ( (nNAentriesY = malloc(ncy * sizeof(int)))==NULL )
  {
    free(nNAentriesX); free(multMatY); free(multMatX); free(xx);
    *err = 1;
    Rprintf("bicor: memmory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  if ( (NAmeanX = malloc(ncx * sizeof(int)))==NULL )
  {
    free(nNAentriesY); free(nNAentriesX); free(multMatY); free(multMatX); free(xx);
    *err = 1;
    Rprintf("bicor: memmory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  if ( (NAmeanY = malloc(ncy * sizeof(int)))==NULL )
  {
    free(NAmeanX); free(nNAentriesY); free(nNAentriesX); free(multMatY); free(multMatX); free(xx);
    *err = 1;
    Rprintf("bicor: memmory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  // Calculate multMatX = (x - median(x))/(9.0 * qnorm75 * mad(x)); median and mad are columnwise.
  // or multMatX = x-mean(x) if robustX is 0.

  int zeroMAD;
  int pearsonX = 0, pearsonY = 0;
  if (*robustX)
  {
    if (*fallback < 3)
    {
      for (int i = 0; i < ncx; i++)
        prepareColBicor(x + i*nr, nr, *maxPOutliers, *fallback, multMatX + i*nr, nNAentriesX + i, NAmeanX + i, 
                        &zeroMAD, xx, xxx);
    } else {
      int doFallback = 0;
      for (int i = 0; i < ncx; i++)
      {
        prepareColBicor(x + i*nr, nr, *maxPOutliers, *fallback, multMatX + i*nr, nNAentriesX + i, NAmeanX + i, 
                        &zeroMAD, xx, xxx);
        if (zeroMAD > 0) { doFallback = i+1; i=ncx; }
      }
      if (doFallback > 0)
      {
        if (*verbose > 0)
          Rprintf("Warning in bicor: zero-MAD column %d in variable x.\n   %s", doFallback, 
                  "         Switching to Pearson correlation for variable x.\n\n");
        pearsonX = 1;
        *robustX = 0;
        for (int i = 0; i < ncx; i++)
          prepareColCor(x + i*nr, nr, multMatX + i*nr, nNAentriesX + i, NAmeanX + i);
      }
    }   
  } else {
    for (int i = 0; i < ncx; i++)
      prepareColCor(x + i*nr,  nr, multMatX + i*nr, nNAentriesX + i, NAmeanX + i);
  }

  if (*robustY)
  {
    if (*fallback < 3)
    {
      for (int i = 0; i < ncy; i++)
        prepareColBicor(y + i*nr, nr, *maxPOutliers, *fallback, multMatY + i*nr, nNAentriesY + i, NAmeanY + i, 
                        &zeroMAD, xx, xxx);
    } else {
      int doFallback = 0;
      for (int i = 0; i < ncy; i++)
      {
        prepareColBicor(y + i*nr, nr, *maxPOutliers, *fallback, multMatY + i*nr, nNAentriesY + i, NAmeanY + i, 
                        &zeroMAD, xx, xxx);
        if (zeroMAD > 0) { doFallback = i+1; i = ncy;}
      }
      if (doFallback)
      {
        if (*verbose > 0)
          Rprintf("Warning in bicor: zero-MAD column %d in variable y.\n   %s", doFallback, 
                  "         Switching to Pearson correlation for variable y.\n\n");
        pearsonY = 1;
        *robustY = 0;
        for (int i = 0; i < ncy; i++)
          prepareColCor(y + i*nr, nr, multMatY + i*nr, nNAentriesY + i, NAmeanY + i);
      }
    }   
  } else {
    for (int i = 0; i < ncy; i++)
      prepareColCor(y + i*nr,  nr, multMatY + i*nr, nNAentriesY + i, NAmeanY + i);
  }

  // Rprintf("NAmeanX: %d\n", NAmeanX[0])

  // This is older code, kept here for reference.
  // for (int j = 0; j < ncy; j++)
  // {
    // if (*robustY)
      // prepareColBicor(y + j*nr, nr, *maxPOutliers, multMatY + j*nr, nNAentriesY + j, &NAmed, yy, yyy);
    // else
      // prepareColCor(y + j*nr,  nr, multMatY + j*nr, nNAentriesY + j, &NAmed);
    // // if (NAmed) Rprintf("median is NA in column %d of y\n", j);
    // NAmeanY[j] = NAmed;
  // }
  // Rprintf("NAmeanY: %d\n", NAmeanY[0])

  // main multiplication.

  double alpha = 1.0, beta = 1.0;
  dgemm_("T", "N", ncolx, ncoly, nrow, & alpha, multMatX, nrow, multMatY, nrow, & beta, result, ncolx);

  int nSlow = 0;
  double * xx3, *yy3;

  int fbx = *fallback;
  int fby = *fallback;

  if (!robustX || pearsonX) fbx = 4;
  if (!robustY || pearsonY) fby = 4;

  double * resx = result;
  for (int j=0; j<ncy; j++) for (int i=0; i<ncx; i++) 
  {
     if ((NAmeanX[i]>0) || (NAmeanY[j] > 0))
     {
       *resx = NA_REAL;
     } else 
     {  
        if ( (nNAentriesX[i] > maxDiffNA) || ( nNAentriesY[j] > maxDiffNA))
        {
           memcpy((void *)xx, (void *)(x + i*nr), nr * sizeof(double));
           memcpy((void *)yy, (void *)(y + j*nr), nr * sizeof(double));
   
           int nNAx = 0, nNAy = 0;    
           for (int k=0; k<nr; k++)
           {
              if (ISNAN(xx[k])) yy[k] = NA_REAL;
              if (ISNAN(yy[k])) xx[k] = NA_REAL;
              if (ISNAN(xx[k])) nNAx++;
              if (ISNAN(yy[k])) nNAy++;
           }
           int NAx = 0, NAy = 0;
   
           if ((nNAx - nNAentriesX[i] > maxDiffNA) || (nNAy-nNAentriesY[j] > maxDiffNA))
           {
               // Rprintf("Recalculating row %d and column %d, column size %d\n", i, j, nr);
               // must recalculate the auxiliary variables for both columns
               int temp = 0;
               if (nNAx - nNAentriesX[i] > maxDiffNA)
               {
                  // Rprintf("...Recalculating row... \n");
                  // if (robustX && (pearsonX==0))
                     prepareColBicor(xx, nr, *maxPOutliers, fbx, xxx, &temp, &NAx, &zeroMAD, xx2, yy2);
                  // else
                    //  prepareColCor(xx, nr, xxx, &temp, &NAx);
                  xx3 = xxx;
               } else
                  xx3 = multMatX + i * nr;
               if (nNAy-nNAentriesY[j] > maxDiffNA)
               {
                  // Rprintf("...Recalculating column... \n");
                  // if (robustY && (pearsonY==0))
                     prepareColBicor(yy, nr, *maxPOutliers, fby, yyy, &temp, &NAy, &zeroMAD, yy2, xx2);
                  // else
                     // prepareColCor(yy, nr, yyy, &temp, &NAy);
                  yy3 = yyy;
               } else
                  yy3 = multMatY + j * nr;
               if (NAx + NAy==0)
               {
                  LDOUBLE sumxy = 0;
                  int count = 0;
                  for (int k=0; k<nr; k++)
                  {
                    double vx = *(xx3 + k), vy = *(yy3 +  k);
                    // Rprintf("i: %d, j: %d, k: %d: vx: %e, vy: %e\n", i,j,k,vx,vy);
                    if (!ISNAN(vx) && !ISNAN(vy))
                    {
                      sumxy += vx * vy; 
                      count++;
                    }
                  }
                  if (count==0) 
                  {
                     *resx = NA_REAL; 
                     (*nNA)++;
                  } else {
                     *resx = (double) sumxy;
                  }
                } else {
                   *resx = NA_REAL; 
                   (*nNA)++;
                }
                nSlow++;
           }
        }
        if ( !ISNAN(*resx))
        { 
            if (*resx > 1.0) *resx = 1.0;
            if (*resx < -1.0) *resx = -1.0;
        }
     }
     resx++;
  }

  if (*verbose) Rprintf("%s Fraction of slow calculations: %f\n", spaces,
                        ( (double) nSlow) / (ncx*ncy) );

  free(NAmeanY);
  free(NAmeanX);
  free(nNAentriesY);
  free(nNAentriesX);
  free(multMatY);
  free(multMatX);
  free(xx);
}


#endif
