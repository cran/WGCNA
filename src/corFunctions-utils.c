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

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include <R.h>
#include <R_ext/BLAS.h>
#include <R_ext/libextern.h>

#include "pivot.h"
#include "conditionalThreading.h"
#include "corFunctions-typeDefs.h"
#include "corFunctions-utils.h"

#define RefUX	0.5

/*===================================================================================
 *
 * median
 *
 * ==================================================================================*/

// Here I first put all NAs to the end, then call the pivot function to find the median of the remaining
// (finite) entries.

double median(double * x, size_t n, int copy, int * err)
{
  double * xx, res;
  if (copy)
  {
    if ( (xx=(double *) malloc(n * sizeof(double)))==NULL ) 
    {
      Rprintf("Memory allocation error in median(). Could not allocate %d kB.\n", 
              (int) (n * sizeof(double) / 1024 + 1));
      *err = 1;
      return NA_REAL;
    }
    memcpy((void *)xx, (void *)x, n * sizeof(double));
  } else xx = x;

    
  *err = 0;
  // Put all NA's at the end.
  size_t bound = n;
  for (size_t i=n; i>0; ) 
  {
    i--;
    if (ISNAN(xx[i]))
    {
       bound--;
       xx[i] = xx[bound];
       xx[bound] = NA_REAL;
    }
  }

  // Rprintf("Median: n: %d, bound: %d\n", n, bound);
  // Any non-NA's left?

  if (bound==0)
    res = NA_REAL;
  else 
  // yes, return the appropriate pivot. 
    res = pivot(xx, bound, ( 1.0 * (bound-1))/2);

  if (copy) free(xx);

  return res;

}


/*===================================================================================
 *
 * quantile
 *
 * ==================================================================================*/

// Here I first put all NAs to the end, then call the pivot function to find the appropriate
// quantile of the remaining (finite) entries.

// q is the quantile: 1/2 will give exactly the median above.

double quantile(double * x, size_t n, double q, int copy, int * err)
{
  double * xx;
  double res;

  if (copy)
  {
    if ( (xx=(double *) malloc(n * sizeof(double)))==NULL ) 
    {
      Rprintf("Memory allocation error in quantile(). Could not allocate %d kB.\n", 
              (int) (n * sizeof(double) / 1024 + 1));
      *err = 1;
      return NA_REAL;
    }
    memcpy((void *)xx, (void *)x, n * sizeof(double));
  } else xx = x;

    
  *err = 0;
  // Put all NA's at the end.
  size_t bound = n;
  for (size_t i=n; i>0; ) 
  {
    i--;
    if (ISNAN(xx[i]))
    {
       bound--;
       xx[i] = xx[bound];
       xx[bound] = NA_REAL;
    }
  }

  // Rprintf("Quantile: q: %f, n: %d, bound: %d\n", q, n, bound);
  // Any non-NA's left?

  if (bound==0)
    res = NA_REAL;
  else
  // yes, return the appropriate pivot. 
    res = pivot(xx, bound, ( 1.0 * (bound-1))*q);

  if (copy) free(xx);

  return res;

}

double quantile_noCopy(double * x, size_t n, double q)
{
  double res;
  // Put all NA's at the end.
  size_t bound = n;
  for (size_t i=n; i>0; ) 
  {
    i--;
    if (ISNAN(x[i]))
    {
       bound--;
       x[i] = x[bound];
       x[bound] = NA_REAL;
    }
  }

  // Rprintf("Quantile: q: %f, n: %d, bound: %d\n", q, n, bound);
  // Any non-NA's left?

  if (bound==0)
    res = NA_REAL;
  else
  // yes, return the appropriate pivot. 
    res = pivot(x, bound, ( 1.0 * (bound-1))*q);

  return res;

}

/*==========================================================================================
 *
 * testMedian
 *
 * =========================================================================================*/


void testMedian(double *x, int * n, double * res)
{
  int err;
  *res = median(x, (size_t) *n, 0, &err);
} 

/*==========================================================================================
 *
 * testQuantile
 *
 * =========================================================================================*/


void testQuantile(double *x, int *n, double *q, double *res)
{
  int err;
  *res = quantile(x, (size_t) *n, *q, 0, &err);
} 


/*==========================================================================================
 *
 * prepareColBicor
 *
 * =========================================================================================*/

// prepareColBicor: calculate median and mad of x and put 
// (1-u^2)^2 * (x - median(x))/(9.0 * qnorm75 * mad(x))/ appropriate normalization 
// into res. 
// res must have enough space allocated to hold the result; 
// aux and aux2 each must also have enough space to hold a copy of x.

// maxQoutliers is the maximum allowed proportion of outliers on either side of the median.

// fallback: 1: none, 2: individual, 3: all, 4: force Pearson calculation. 4 is necessary for remedial
// calculations.

// In this case: Pearson pre-calculation entails normalizing columns by mean and variance. 

void prepareColBicor(double * col, size_t nr, double maxPOutliers, int fallback,
                     int cosine,
                     double * res, size_t * nNAentries, 
                     int * NAmed, volatile int * zeroMAD,
                     double * aux, double * aux2)
{
  // const double asymptCorr = 1.4826, qnorm75 = 0.6744898;
  // Note to self: asymptCorr * qnorm75 is very close to 1 and should equal 1 theoretically. Should
  // probably leave them out completely. 

  if (fallback==4)
  {
    prepareColCor(col, nr, cosine, res, nNAentries, NAmed);
    return;
  }

  int err = 0;

  // Calculate the median of col

  memcpy((void *)res, (void *)col, nr * sizeof(double));
  double med = median(res, nr, 0, &err);

  // Create a conditional copy of the median
  double medX;
  if (cosine) medX = 0; else medX = med;

  *zeroMAD = 0;
  // calculate absolute deviations from the median

  if (ISNAN(med))
  {
    *NAmed = 1;
    for (size_t k=0; k<nr; k++) *(res + k) = 0;
  } else {
    *NAmed = 0;
    *nNAentries = 0;
    for (size_t k=0; k<nr; k++)
      if (ISNAN(col[k]))
      {
        (*nNAentries)++;
        res[k] = NA_REAL;
        aux[k] = NA_REAL;
      } else {
        res[k] = col[k] - medX;
        aux[k] = fabs(col[k] - med);
      }

    // calculate mad, i.e. median absolute deviation
    double mad = median(aux, nr, 0, &err);

    // If mad is zero, value of fallback decides what is it we will do.
    if (mad==0)
    {
       *zeroMAD = 1;
       switch (fallback)
       {
          case 1: 
          {
             // Return after zeoring out results and setting the NAmed flag
             for (size_t k=0; k<nr; k++) res[k] = 0;
             *NAmed = 1;
             return;
          }
          case 2: 
             // Switch to Pearson correlation and return
             // Rprintf("mad is zero in a column. Switching to Pearson for this column.\n");
             prepareColCor(col, nr, cosine, res, nNAentries, NAmed);
          case 3: 
             // Do nothing: the setting of *zeroMAD above is enough.
             return;
       }
    } 

    // We now re-use aux to store a copy of the weights ux. To calculate them, first get (x-med)/(9*mad).

    // Rprintf("median: %6.4f, mad: %6.4f, cosine: %d\n", med, mad, cosine);

    double denom = 9.0 * mad;
    for (size_t k=0; k<nr; k++)
      if (!ISNAN(col[k]))
        aux[k] = (col[k] - med) / denom;
      else
        aux[k] = NA_REAL;

    // Get the low and high quantiles  of ux
    memcpy((void *)aux2, (void *)aux, nr * sizeof(double));
    double lowQ = quantile(aux2, nr, maxPOutliers, 0, &err);

    memcpy((void *)aux2, (void *)aux, nr * sizeof(double));
    double hiQ = quantile(aux2, nr, 1-maxPOutliers, 0, &err);

    // Rprintf("prepareColBicor: lowQ=%f, hiQ = %f\n", lowQ, hiQ);
    // If the low quantile is below -1, rescale the aux (that serve as ux below)
    // such that the low quantile will fall at -1; similarly for the high quantile

    if (lowQ > -RefUX) lowQ = -RefUX;
    if (hiQ < RefUX) hiQ = RefUX;
    lowQ = fabs(lowQ);

    for (size_t k=0; k<nr; k++) if (!ISNAN(aux[k]))
    {
      if (aux[k] < 0)
        aux[k] = aux[k] * RefUX / lowQ;
      else
        aux[k] = aux[k] * RefUX / hiQ;
    }

    // Calculate the (1-ux^2)^2 * (x-median(x)) 

    LDOUBLE sum = 0;
    for (size_t k=0; k<nr; k++)
      if (!ISNAN(res[k]))
      {
        double ux = aux[k];
        if (fabs(ux) > 1) ux = 1;  // sign of ux doesn't matter.
        ux = 1-ux*ux;
        res[k] *= ux*ux ;
        sum += res[k]*res[k];
      } else
        res[k] = 0;
    sum = sqrtl(sum);
    if (sum==0)
    {
       for (size_t k=0; k<nr; k++)
          res[k] = 0;
       *NAmed = 1;
    } else {
       for (size_t k=0; k<nr; k++)
          res[k] = res[k] / sum;
    }
  }
}

/*======================================================================================
 *
 * prepareColCor
 *
 * =====================================================================================*/

// Used for Pearson correlation fallback
// and when bicor is called with robustsX or robustY = 0
// if cosine is not zero, the cosine correlation will be calculated.

void prepareColCor(double * x, size_t nr, int cosine, double * res, size_t * nNAentries, int * NAmean)
{
  *nNAentries = 0;
  size_t count = 0;
  LDOUBLE mean = 0, sum = 0;
  for (size_t k = 0; k < nr; k++)
    if (!ISNAN(x[k]))
    {
      mean += x[k];
      sum += ((LDOUBLE) x[k])*( (LDOUBLE) x[k]);
      count ++;
    }
  if (count > 0)
  {
    *NAmean = 0;
    *nNAentries = nr-count;
    if (cosine) mean = 0; else mean = mean/count;
    sum = sqrtl(sum - count * mean*mean);
    if (sum > 0)
    {
      // Rprintf("sum: %Le\n", sum);
       for (size_t k=0; k<nr; k++)
         if (!ISNAN(x[k]))
            res[k] = (x[k] - mean)/sum;
         else
            res[k] = 0;
    } else {
       // Rprintf("prepareColCor: have zero variance.\n");
       *NAmean = 1;
       for (size_t k=0; k<nr; k++) res[k] = 0;
    }
  } else {
    *NAmean = 1;
    *nNAentries = nr;
    for (size_t k=0; k<nr; k++)
       res[k] = 0;
  }
}

/*======================================================================================
 *
 * prepareColCor
 *
 * =====================================================================================*/

// Used for Pearson correlation fallback
// and when bicor is called with robustsX or robustY = 0
// if cosine is not zero, the cosine correlation will be calculated.

void prepareColCor_weighted(double * x, double * weights, 
      size_t nr, int cosine, double * res, size_t * nNAentries, int * NAmean)
{
  *nNAentries = 0;
  size_t count = 0;
  LDOUBLE mean = 0, wsum = 0, wsumSq = 0, sumSq = 0, sumxwSq = 0;
  for (size_t k = 0; k < nr; k++)
    if (!ISNAN(x[k]) && !ISNAN(weights[k]))
    {
      wsum += weights[k];
      mean += x[k] * weights[k];
      sumSq += ((LDOUBLE) x[k])* x[k] * weights[k] * weights[k];  
      sumxwSq += ( (LDOUBLE) x[k]) * weights[k] * weights[k];
      wsumSq += ( (LDOUBLE) weights[k]) * weights[k];
      count ++;
    }
  if (count > 0)
  {
    *NAmean = 0;
    *nNAentries = nr-count;
    if (cosine) mean = 0; else mean = mean/wsum;
    sumSq = sqrtl(sumSq - 2*mean * sumxwSq  + mean*mean * wsumSq);
    //Rprintf("\nprepareColCor_weighted: \n");
    //Rprintf("  mean: %5.3Lf, sumSq: %5.3Lf\n", mean, sumSq);
    //Rprintf("  x: "); RprintV(x, nr);
    //Rprintf("  weights: "); RprintV(weights, nr);
    if ((wsum > 0) && (sumSq > 0))
    {
      // Rprintf("sum: %Le\n", sum);
       for (size_t k=0; k<nr; k++)
         if (!ISNAN(x[k]))
            res[k] = weights[k] * (x[k] - mean)/sumSq;
         else
            res[k] = 0;
    } else {
       // Rprintf("prepareColCor: have zero variance.\n");
       *NAmean = 1;
       for (size_t k=0; k<nr; k++) res[k] = 0;
    }
    //Rprintf("res: "); RprintV(res, nr);
  } else {
    *NAmean = 1;
    *nNAentries = nr;
    for (size_t k=0; k<nr; k++) res[k] = 0;
  }
}

/*===================================================================================================
 *
 * Threaded "slow" calculations for weighted pearson correlation
 *
 *===================================================================================================
*/

// basic function that calculates the (unweighted) correlation of two columns.

// The input pointers must point to the start of the rows in x, y
// The res pointer must point to the appropriate component of the output.
// The return value is 1 if the result is undefined (NA) and 0 if it is valid.

int basic2variableCorrelation(
   double *xx, double *yy,
   size_t nr,
   double *res,
   int cosineX, int cosineY)
{
  LDOUBLE sumxy = 0, sumx = 0, sumy = 0, sumxs = 0, sumys = 0;
  size_t count = 0;
  double vx, vy;
  for (size_t k=0; k<nr; k++)
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
      *res = NA_REAL; 
      return 1;
  } else {
      if (cosineX) sumx = 0;
      if (cosineY) sumy = 0;
      LDOUBLE varx = sumxs - sumx*sumx/count,
              vary = sumys - sumy*sumy/count;
      if (varx==0 || vary==0)
      {
         *res = NA_REAL;
         return 1;
      } else
         *res = (double) ( (sumxy - sumx * sumy/count)/ sqrtl( varx * vary));
  }
  return 0;
}

// basic function that calculates the weighted correlation of two columns with two different sets of weights.
// The input pointers must point to the start of the rows in x, y, weights for x and weights for y.
//
// The res pointer must point to the appropriate component of the output.
//
// The return value is 1 if the result is undefined (NA) and 0 if it is valid.

int basic2variableCorrelation_weighted(
   double *xx, double *yy, 
   double *wx, double *wy,
   size_t nr,
   // output
   double *res,
   // options
   int cosineX, int cosineY)
{
  LDOUBLE 
     sum_wx_x = 0, sum_w_x = 0, sum_wSq_x = 0, sum_xwSq_x = 0, sum_Sq_x = 0,
     sum_wx_y = 0, sum_w_y = 0, sum_wSq_y = 0, sum_xwSq_y = 0, sum_Sq_y = 0,
     sum_wwxy = 0, sum_wwy = 0, sum_wwx = 0, sum_ww = 0;

  size_t count = 0;

  for (size_t k=0; k<nr; k++)
  {
     double vx = *xx, vy = *yy, vwx = *wx, vwy = *wy;
     if (!ISNAN(vx) && !ISNAN(vy) && !ISNAN(vwx) && !ISNAN(vwy))
     {
       double ww = vwx * vwy;
       count ++;
       sum_wx_x += vx * vwx;
       sum_Sq_x += ((LDOUBLE) vx) * vx * vwx * vwx;
       sum_xwSq_x += ( (LDOUBLE) vx) * vwx * vwx;
       sum_w_x += vwx;
       sum_wSq_x += ( (LDOUBLE) vwx) * vwx;

       sum_wx_y += vy * vwy;
       sum_Sq_y += ((LDOUBLE) vy) * vy * vwy * vwy;
       sum_xwSq_y += ( (LDOUBLE) vy) * vwy * vwy;
       sum_w_y += vwy;
       sum_wSq_y += ( (LDOUBLE) vwy) * vwy;

       sum_wwxy += ( (LDOUBLE) vx) * vy * ww;
       sum_wwx += ( (LDOUBLE) vx) * ww;
       sum_wwy += ( (LDOUBLE) vy) * ww;
       sum_ww += ww;
     }
     xx++; yy++;
     wx++; wy++;
  }
  // Rprintf("Count of included observations: %d\n", count);
  if (count==0) 
  {
    *res = NA_REAL; 
    return 1;
  } else {
    double 
      mean_x = cosineX ? 0 : sum_wx_x/sum_w_x,
      mean_y = cosineY ? 0 : sum_wx_y/sum_w_y,
      varx = sum_Sq_x - 2*mean_x * sum_xwSq_x + mean_x * mean_x * sum_wSq_x,
      vary = sum_Sq_y - 2*mean_y * sum_xwSq_y + mean_y * mean_y * sum_wSq_y;

    if (varx==0 || vary==0)
    { 
      *res = NA_REAL;
      return 1;
    } else 
       *res = (double) ( (sum_wwxy - mean_x*sum_wwy - mean_y * sum_wwx + sum_ww *  mean_x * mean_y)/
                     sqrt( varx * vary));

  }
  return 0;
}

/*=========================================================================================================
 *
 *
 * Threaded functions
 *
 *
 *=========================================================================================================*/


/*======================================================================================
 *
 * prepareColBicor
 *
 * =====================================================================================*/

void * threadPrepColBicor(void * par)
{
  colPrepThreadData volatile * td = (colPrepThreadData *) par;
  cor1ThreadData volatile * x = td->x;

  // Rprintf("Preparing columns: nr = %d, nc = %d\n", x->nr, x->nc);
  while (td->pc->i < td->pc->n)
  {
      // Grab the next column that needs to be done
      pthread_mutex_lock_c( td->lock, x->threaded);
      if (td->pc->i < td->pc->n)
      {
         size_t col = td->pc->i;
         // Rprintf("...working on column %d in thread %d\n", col, td->x->id);
         td->pc->i++;
         pthread_mutex_unlock_c( td->lock, x->threaded );
 
         prepareColBicor(x->x + col * x->nr, 
                         x->nr, 
                         x->maxPOutliers, 
                         x->fallback, 
                         x->cosine,
                         x->multMat + col * x->nr,
                         x->nNAentries + col,
                         x->NAme + col,
                         &(x->zeroMAD),
                         x->aux,
                         x->aux + x->nr);
         // if (x->zeroMAD > 0) { Rprintf("threadPrepColBicor: mad was zero in column %d.\n", col); }
         if (x->zeroMAD > 0) *(x->warn) = warnZeroMAD;
         if ( (x->zeroMAD > 0) && (x->fallback==3)) 
         { 
           pthread_mutex_lock_c( td->lock, x->threaded );
           // Rprintf("threadPrepColBicor: Moving counter from %d %d to end at %d in thread %d.\n", 
                   // col, td->pc->i, td->pc->n, x->id);
           x->zeroMAD = col+1; td->pc->i = td->pc->n; 
           pthread_mutex_unlock_c( td->lock, x->threaded );
         }
      } else 
         pthread_mutex_unlock_c( td->lock, x->threaded );
  }
  return NULL;
} 
      

/*======================================================================================
 *
 * prepareColCor
 *
 * =====================================================================================*/

// Used for the fast calculation of Pearson correlation
// and when bicor is called with robustsX or robustY = 0


void * threadPrepColCor(void * par)
{
  colPrepThreadData volatile * td = (colPrepThreadData *) par;
  cor1ThreadData volatile * x = td->x;
  //Rprintf("threadPrepColCor: starting in thread %d: counter.i = %d, counter.n = %d, nc = %d.\n", 
  //         td->x->id, td->pc->i, td->pc->n, td->x->nc);
  while (td->pc->i < td->pc->n)
  {
      // Grab the next column that needs to be done
      pthread_mutex_lock_c( td->lock, x->threaded );
      int col = td->pc->i;
      if (col < td->x->nc)
      {
         td->pc->i++;
    //     Rprintf("threadPrepColCor: preparing column %d in thread %d.\n", col, td->x->id);
         pthread_mutex_unlock_c( td->lock, x->threaded );
 
         prepareColCor(x->x + col * x->nr, 
                       x->nr, 
                       x->cosine,
                       x->multMat + col * x->nr,
                       x->nNAentries + col,
                       x->NAme + col);
      } else 
         pthread_mutex_unlock_c( td->lock, x->threaded );
  }
  return NULL;
} 
      
/*===================================================================================================
 *
 * Threaded symmetrization and NA'ing out of rows and columns with NA means
 *
 *===================================================================================================
*/

void * threadSymmetrize(void * par)
{
  symmThreadData * td = (symmThreadData *) par;
  cor1ThreadData * x = td->x;

  size_t nc = x->nc;
  double * result = x->result;
  int * NAmean = x->NAme;
  size_t col = 0;
  while ( (col = td->pc->i) < nc)
  {
      // Symmetrize the column
      // point counter to the next column
      td->pc->i = col+1;
      // and update the matrix. Start at j=col to check for values greater than 1.
      if (NAmean[col] == 0)
      {
        double * resx = result + col*nc + col;
        // Rprintf("Symmetrizing row %d to the same column.\n", col);
        for (size_t j=col; j<nc; j++) 
        {
          // Rprintf("..symmetrizing element %d...\n", j);
          if (NAmean[j] == 0)
          {
            if (!ISNAN(*resx))
            {
               if (*resx > 1.0) *resx = 1.0;
               if (*resx < -1.0) *resx = -1.0;
            }
            result[j*nc + col] = *resx;
          }
          resx ++;
        }
      } else {
        // Rprintf("NA-ing out column and row %d\n", col);
        for (size_t j=0; j<nc; j++)
        {
           result[col*nc + j] = NA_REAL;
           result[j*nc + col] = NA_REAL;
        }
      }
  }
  return NULL;
} 
      
/*===================================================================================================
 *
 * Threaded "slow" calculations for bicor
 *
 *===================================================================================================
*/
// This can actually be relatively slow, since the search for calculations that need to be done is not
// parallel, so one thread may have to traverse the whole matrix. I can imagine parallelizing even that
// part, but for now leave it as is as this will at best be a minuscule improvement.

void * threadSlowCalcBicor(void * par)
{
  slowCalcThreadData * td = (slowCalcThreadData *) par;
  size_t * nSlow = td->nSlow;
  size_t * nNA = td->nNA;
  double * x = td->x->x;
  double * multMat = td->x->multMat;
  double * result = td->x->result;
  int fbx = td->x->fallback;
  int cosine = td->x->cosine;
  size_t nc = td->x->nc, nc1 = nc-1, nr = td->x->nr;
  int * NAmean = td->x->NAme;
  size_t * nNAentries = td->x->nNAentries;
  progressCounter * pci = td->pci, * pcj = td->pcj;

  double maxPOutliers = td->x->maxPOutliers;

  double * xx = td->x->aux, * yy = xx + nr; 
  double * xxx = xx + 2*nr, * yyy = xx + 3*nr;
  double * xx2 = xx + 4*nr, * yy2 = xx + 5*nr;

  size_t maxDiffNA = (size_t) (td->x->quick * nr);

  if (fbx==3) fbx = 2; // For these calculations can't go back and redo everything

 
  // Rprintf("Checking %d rows and %d columns\n", nc1, nc);
  // Rprintf("starting at %d and %d\n", pci->i, pcj->i);
  while (pci->i < nc1)
  {
     pthread_mutex_lock_c( td->lock, td->x->threaded );
     size_t i = pci->i, ii = i;
     size_t j = pcj->i, jj = j;
     do
     {
       i = ii;
       j = jj;
       jj++;
       if (jj==nc) 
       {
         ii++;
         jj = ii+1;
       }
     } while ((i<nc1) && (j<nc) && 
              ((NAmean[i] > 0) || (NAmean[j] > 0) ||
               ( (nNAentries[i] <= maxDiffNA) && ( nNAentries[j] <= maxDiffNA))));
     pci->i = ii;
     pcj->i = jj;
     pthread_mutex_unlock_c( td->lock, td->x->threaded );
 
     if ((i < nc1) && (j < nc) )
     {
        // Rprintf("Recalculating row %d and column %d, column size %d\n", i, j, nr);
        memcpy((void *)xx, (void *)(x + i*nr), nr * sizeof(double));
        memcpy((void *)yy, (void *)(x + j*nr), nr * sizeof(double));

        size_t nNAx = 0, nNAy = 0;    
        for (size_t k=0; k<nr; k++)
        {
           if (ISNAN(xx[k])) yy[k] = NA_REAL;
           if (ISNAN(yy[k])) xx[k] = NA_REAL;
           if (ISNAN(xx[k])) nNAx++;
           if (ISNAN(yy[k])) nNAy++;
        }
        int NAx = 0, NAy = 0;

        if ((nNAx - nNAentries[i] > maxDiffNA) || (nNAy-nNAentries[j] > maxDiffNA))
        {
            // must recalculate the auxiliary variables for both columns
            size_t temp = 0;
            int zeroMAD = 0;
            if (nNAx - nNAentries[i] > maxDiffNA)
               {
                  prepareColBicor(xx, nr, maxPOutliers, fbx, cosine, xxx, &temp, &NAx, &zeroMAD, xx2, yy2);
                  if (zeroMAD) *(td->x->warn) = warnZeroMAD;
               }
               else
                  memcpy((void *) xxx, (void *) (multMat + i * nr),  nr * sizeof(double));
            if (nNAy-nNAentries[j] > maxDiffNA)
               {
                  prepareColBicor(yy, nr, maxPOutliers, fbx, cosine, yyy, &temp, &NAy, &zeroMAD, xx2, yy2);
                  if (zeroMAD) *(td->x->warn) = warnZeroMAD;
               }
               else
                  memcpy((void *) yyy, (void *) (multMat + j * nr),  nr * sizeof(double));
            if (NAx + NAy==0)
            {
               LDOUBLE sumxy = 0;
               size_t count = 0;
               for (size_t k=0; k<nr; k++)
               {
                 double vx = *(xxx + k), vy = *(yyy +  k);
                 // Rprintf("i: %d, j: %d, k: %d: vx: %e, vy: %e\n", i,j,k,vx,vy);
                 if (!ISNAN(vx) && !ISNAN(vy))
                 {
                   sumxy += vx * vy; 
                   count++;
                 }
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
             // result[j*nc + i] = result[i*nc + j];
             (*nSlow)++;
        }
     }
  }
  return NULL;
}

/*===================================================================================================
 *
 * Threaded "slow" calculations for pearson correlation
 *
 *===================================================================================================
*/
// This can actually be relatively slow, since the search for calculations that need to be done is not
// parallel, so one thread may have to traverse the whole matrix. I can imagine parallelizing even that
// part, but for now leave it as is as this will at best be a minuscule improvement.

void * threadSlowCalcCor(void * par)
{
  slowCalcThreadData * td = (slowCalcThreadData *) par;
  size_t * nSlow = td->nSlow;
  size_t * nNA = td->nNA;
  double * x = td->x->x;
  double * result = td->x->result;
  size_t nc = td->x->nc, nc1 = nc-1, nr = td->x->nr;
  int cosine = td->x->cosine;
  int * NAmean = td->x->NAme;
  size_t * nNAentries = td->x->nNAentries;
  progressCounter * pci = td->pci, * pcj = td->pcj;

  double *xx, *yy;
  double vx, vy;

  size_t maxDiffNA = (size_t) (td->x->quick * nr);

  // Rprintf("quick:%f\n", td->x->quick);


  // Rprintf("Checking %d rows and %d columns\n", nc1, nc);
  // Rprintf("starting at %d and %d\n", pci->i, pcj->i);
  while (pci->i < nc1)
  {
     pthread_mutex_lock_c( td->lock, td->x->threaded );
     size_t i = pci->i, ii = i;
     size_t j = pcj->i, jj = j;
     do
     {
       i = ii;
       j = jj;
       jj++;
       if (jj==nc) 
       {
         ii++;
         jj = ii+1;
       }
     } while ((i<nc1) && (j<nc) && 
               ((NAmean[i] > 0) || (NAmean[j] > 0) ||
                ( (nNAentries[i] <= maxDiffNA) && ( nNAentries[j] <= maxDiffNA))));

     pci->i = ii;
     pcj->i = jj;
     pthread_mutex_unlock_c( td->lock, td->x->threaded );
 
     if ((i < nc1) && (j < nc))
     {
        // Rprintf("Recalculating column %d and row %d, column size %d\n", i, j, nr);
        *nNA += basic2variableCorrelation(
                       x + i * nr, x + j * nr,
                       nr, result + i*nc + j,
                       cosine, cosine);
        (*nSlow)++;
     }
  }
  return NULL;
}


/*===================================================================================================
 *
 * Threaded NA-ing
 *
 *===================================================================================================
*/

void * threadNAing(void * par)
{
  NA2ThreadData * td = (NA2ThreadData *) par;

  double * result = td->x->x->result;
  size_t ncx = td->x->x->nc;
  int * NAmedX = td->x->x->NAme;

  size_t ncy = td->x->y->nc;
  int * NAmedY = td->x->y->NAme;

  progressCounter * pci = td->pci;
  progressCounter * pcj = td->pcj;

  // Go row by row

  size_t row = 0, col = 0;

  while  ((row = pci->i) < ncx)
  {
      pci->i = row + 1;
      if (NAmedX[row])
      {
         // Rprintf("NA-ing out column and row %d\n", col);
         for (size_t j=0; j<ncy; j++)
              result[row + j * ncx] = NA_REAL;
      } 
  }

  // ... and column by column

  while ( (col = pcj->i) < ncy)
  {
     pcj->i = col + 1;
     if (NAmedY[col])
     {
        // Rprintf("NA-ing out column and row %d\n", col);
        for (size_t i=0; i<ncx; i++)
             result[i + col * ncx] = NA_REAL;
     } else {
        double *resx = result + col*ncx;
        for (size_t i=0; i<ncx; i++) 
        {
           if (!ISNAN(*resx))
           {
             if (*resx > 1.0) *resx = 1.0;
             if (*resx < -1.0) *resx = -1.0;
           }
           resx++;
        }
     }
  }

  return NULL;
} 


/*===================================================================================================
 *
 * Threaded "slow" calculations for bicor(x,y)
 *
 *===================================================================================================
*/
// This can actually be relatively slow, since the search for calculations that need to be done is not
// parallel, so one thread may have to traverse the whole matrix. I can imagine parallelizing even that
// part, but for now leave it as is as this will at best be a minuscule improvement.

void * threadSlowCalcBicor2(void * par)
{
  slowCalc2ThreadData * td = (slowCalc2ThreadData *) par;
  size_t * nSlow = td->nSlow;
  size_t * nNA = td->nNA;

  double * x = td->x->x->x;
  double * multMatX = td->x->x->multMat;
  double * result = td->x->x->result;
  size_t ncx = td->x->x->nc, nr = td->x->x->nr;
  int * NAmeanX = td->x->x->NAme;
  size_t * nNAentriesX = td->x->x->nNAentries;
  int robustX = td->x->x->robust;
  int fbx = td->x->x->fallback;
  int cosineX = td->x->x->cosine;

  double * y = td->x->y->x;
  double * multMatY = td->x->y->multMat;
  size_t ncy = td->x->y->nc;
  int * NAmeanY = td->x->y->NAme;
  size_t * nNAentriesY = td->x->y->nNAentries;
  int robustY = td->x->y->robust;
  int fby = td->x->y->fallback;
  int cosineY = td->x->y->cosine;

  double maxPOutliers = td->x->x->maxPOutliers;

  progressCounter * pci = td->pci, * pcj = td->pcj;

  double * xx = td->x->x->aux;
  double * xxx = xx + nr;
  double * xx2 = xx + 2*nr;

  double * yy = td->x->y->aux;
  double * yyy = yy + nr;
  double * yy2 = yy + 2*nr;

  double * xx3, *yy3;

  int maxDiffNA = (int) (td->x->x->quick * nr);

  if (fbx==3) fbx = 2;
  if (fby==3) fby = 2;

  if (!robustX) fbx = 4;
  if (!robustY) fby = 4;

  // Rprintf("Remedial calculation thread #%d: starting at %d and %d\n", td->x->x->id, 
  //            pci->i, pcj->i);
  //

  while (pci->i < ncx)
  {
     pthread_mutex_lock_c( td->lock, td->x->x->threaded );
     size_t i = pci->i, ii = i;
     size_t j = pcj->i, jj = j;
     do
     {
       i = ii;
       j = jj;
       jj++;
       if (jj==ncy) 
       {
         ii++;
         jj = 0;
       }
     } while ((i<ncx) && (j<ncy) && 
              ((NAmeanX[i] > 0) || (NAmeanY[j] > 0) || 
                ( (nNAentriesX[i] <= maxDiffNA) && ( nNAentriesY[j] <= maxDiffNA))));
     pci->i = ii;
     pcj->i = jj;
     pthread_mutex_unlock_c( td->lock, td->x->x->threaded );
 
     if ((i < ncx) && (j < ncy))
     {
        memcpy((void *)xx, (void *)(x + i*nr), nr * sizeof(double));
        memcpy((void *)yy, (void *)(y + j*nr), nr * sizeof(double));

        size_t nNAx = 0, nNAy = 0;    
        for (size_t k=0; k<nr; k++)
        {
           if (ISNAN(xx[k])) yy[k] = NA_REAL;
           if (ISNAN(yy[k])) xx[k] = NA_REAL;
           if (ISNAN(xx[k])) nNAx++;
           if (ISNAN(yy[k])) nNAy++;
        }
        int NAx = 0, NAy = 0;

        if ((nNAx - nNAentriesX[i] > maxDiffNA) || (nNAy-nNAentriesY[j] > maxDiffNA))
        {
            // Rprintf("Recalculating row %d and column %d, column size %d in thread %d\n", i, j, nr,
            //         td->x->x->id);
            // must recalculate the auxiliary variables for both columns
            size_t temp = 0;
            int zeroMAD = 0;
            if (nNAx - nNAentriesX[i] > maxDiffNA)
            {
               // Rprintf("...Recalculating row... \n");
               //if (robustX && (fbx!=4))
                  prepareColBicor(xx, nr, maxPOutliers, fbx, cosineX, xxx, &temp, &NAx, &zeroMAD, xx2, yy2);
                  if (zeroMAD) *(td->x->x->warn) = warnZeroMAD;
               //else
               //   prepareColCor(xx, nr, xxx, &temp, &NAx);
               xx3 = xxx;
            } else
               xx3 = multMatX + i * nr;
            if (nNAy-nNAentriesY[j] > maxDiffNA)
            {
               // Rprintf("...Recalculating column... \n");
               //if (robustY && (fby!=4))
                  prepareColBicor(yy, nr, maxPOutliers, fby, cosineY, yyy, &temp, &NAy, &zeroMAD, xx2, yy2);
                  if (zeroMAD) *(td->x->y->warn) = warnZeroMAD;
               //else
               //   prepareColCor(yy, nr, yyy, &temp, &NAy);
               yy3 = yyy;
            } else
               yy3 = multMatY + j * nr;
            if (NAx + NAy==0)
            {
               // LDOUBLE sumxy = 0;
               double sumxy = 0;
               size_t count = 0;
               for (size_t k=0; k<nr; k++)
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
                  result[i + j*ncx] = NA_REAL; 
                  (*nNA)++;
               } else {
                  result[i + j*ncx] = (double) sumxy;
                  // Rprintf("Recalculated row %d and column %d, column size %d in thread %d: result = %12.6f\n", 
                  //         i, j, nr, td->x->x->id,  result[i + j*ncx]);
               }
            } else {
                result[i + j*ncx] = NA_REAL; 
                (*nNA)++;
            }
            (*nSlow)++;
        }
     }
  }
  return NULL;
}

/*===================================================================================================
 *
 * Threaded "slow" calculations for pearson correlation of 2 variables.
 *
 *===================================================================================================
*/
// This can actually be relatively slow, since the search for calculations that need to be done is not
// parallel, so one thread may have to traverse the whole matrix. I can imagine parallelizing even that
// part, but for now leave it as is as this will at best be a minuscule improvement.

void * threadSlowCalcCor2(void * par)
{
  slowCalc2ThreadData * td = (slowCalc2ThreadData *) par;
  size_t * nSlow = td->nSlow;
  size_t * nNA = td->nNA;

  double * x = td->x->x->x;
//  double * multMatX = td->x->x->multMat;
  double * result = td->x->x->result;
  size_t ncx = td->x->x->nc, nr = td->x->x->nr;
  int * NAmeanX = td->x->x->NAme;
  size_t * nNAentriesX = td->x->x->nNAentries;
  int cosineX = td->x->x->cosine;

  double * y = td->x->y->x;
//  double * multMatY = td->x->y->multMat;
  size_t ncy = td->x->y->nc;
  int * NAmeanY = td->x->y->NAme;
  size_t * nNAentriesY = td->x->y->nNAentries;
  int cosineY = td->x->y->cosine;

  size_t maxDiffNA = (size_t) (td->x->x->quick * nr);

  progressCounter * pci = td->pci, * pcj = td->pcj;
  double * xx, * yy;

  // Rprintf("Will tolerate %d additional NAs\n",  maxDiffNA);
  // Rprintf("Checking %d rows and %d columns\n", nc1, nc);
  // Rprintf("starting at %d and %d\n", pci->i, pcj->i);
  //

  while (pci->i < ncx)
  {
     pthread_mutex_lock_c( td->lock, td->x->x->threaded );
     size_t i = pci->i, ii = i;
     size_t j = pcj->i, jj = j;
     do
     {
       i = ii;
       j = jj;
       jj++;
       if (jj==ncy) 
       {
         ii++;
         jj = 0;
       }
     } while ((i<ncx) && (j<ncy) && 
              ((NAmeanX[i] > 0) || (NAmeanY[j] > 0) || 
                ( (nNAentriesX[i] <= maxDiffNA) && ( nNAentriesY[j] <= maxDiffNA))));
     pci->i = ii;
     pcj->i = jj;
     pthread_mutex_unlock_c( td->lock, td->x->x->threaded );
 
     if ((i < ncx) && (j < ncy))
     {
        // Rprintf("Recalculating row %d and column %d, column size %d; cosineX: %d, cosineY: %d\n", 
        //         i, j, nr, cosineX, cosineY);
        *nNA += basic2variableCorrelation(
                       x + i * nr, y + j * nr,
                       nr, result + i + j*ncx,
                       cosineX, cosineY);
        (*nSlow)++;
     }
  }
  return NULL;
}

/*======================================================================================
 *
 * threaded prepareColCor_weighted
 *
 * =====================================================================================*/

// Used for the fast calculation of Pearson correlation
// and when bicor is called with robustsX or robustY = 0


void * threadPrepColCor_weighted(void * par)
{
  colPrepThreadData volatile * td = (colPrepThreadData *) par;
  cor1ThreadData volatile * x = td->x;
  //Rprintf("threadPrepColCor: starting in thread %d: counter.i = %d, counter.n = %d, nc = %d.\n", 
  //         td->x->id, td->pc->i, td->pc->n, td->x->nc);
  while (td->pc->i < td->pc->n)
  {
      // Grab the next column that needs to be done
      pthread_mutex_lock_c( td->lock, x->threaded );
      int col = td->pc->i;
      if (col < td->x->nc)
      {
         td->pc->i++;
    //     Rprintf("threadPrepColCor: preparing column %d in thread %d.\n", col, td->x->id);
         pthread_mutex_unlock_c( td->lock, x->threaded );
 
         prepareColCor_weighted(x->x + col * x->nr, 
                       x->weights + col * x->nr,
                       x->nr, 
                       x->cosine,
                       x->multMat + col * x->nr,
                       x->nNAentries + col,
                       x->NAme + col);
      } else 
         pthread_mutex_unlock_c( td->lock, x->threaded );
  }
  return NULL;
} 
      
void * threadSlowCalcCor_weighted(void * par)
{
  slowCalcThreadData * td = (slowCalcThreadData *) par;
  size_t * nSlow = td->nSlow;
  size_t * nNA = td->nNA;
  double * x = td->x->x;
  double * weights = td->x->weights;
  double * result = td->x->result;
  size_t nc = td->x->nc, nc1 = nc-1, nr = td->x->nr;
  int cosine = td->x->cosine;
  int * NAmean = td->x->NAme;
  size_t * nNAentries = td->x->nNAentries;
  progressCounter * pci = td->pci, * pcj = td->pcj;

  size_t maxDiffNA = (size_t) (td->x->quick * nr);

  // Rprintf("quick:%f\n", td->x->quick);


  // Rprintf("Checking %d rows and %d columns\n", nc1, nc);
  // Rprintf("starting at %d and %d\n", pci->i, pcj->i);
  while (pci->i < nc1)
  {
     pthread_mutex_lock_c( td->lock, td->x->threaded );
     size_t i = pci->i, ii = i;
     size_t j = pcj->i, jj = j;
     do
     {
       i = ii;
       j = jj;
       jj++;
       if (jj==nc) 
       {
         ii++;
         jj = ii+1;
       }
     } while ((i<nc1) && (j<nc) && 
               ((NAmean[i] > 0) || (NAmean[j] > 0) ||
                ( (nNAentries[i] <= maxDiffNA) && ( nNAentries[j] <= maxDiffNA))));

     pci->i = ii;
     pcj->i = jj;
     pthread_mutex_unlock_c( td->lock, td->x->threaded );
 
     if ((i < nc1) && (j < nc))
     {
        // Rprintf("Recalculating column %d and row %d, column size %d\n", i, j, nr);
        *nNA += basic2variableCorrelation_weighted(x + i * nr, x + j * nr,
                        weights + i * nr, weights + j * nr,
                        nr, result + i*nc + j,
                        cosine, cosine);
        (*nSlow)++;
     }
  }
  return NULL;
}

/*===================================================================================================
 *
 * Threaded "slow" calculations for weighted pearson correlation of 2 variables.
 *
 *===================================================================================================
*/
// The search for calculations that need to be done is not
// parallel, so one thread may have to traverse the whole matrix. I can imagine parallelizing even that
// part, but for now leave it as is as this will at best be a minuscule improvement.

void * threadSlowCalcCor2_weighted(void * par)
{
 
  slowCalc2ThreadData * td = (slowCalc2ThreadData *) par;
  size_t * nSlow = td->nSlow;
  size_t * nNA = td->nNA;

  double * x = td->x->x->x;
  double * weights_x = td->x->x->weights;
//  double * multMatX = td->x->x->multMat;
  double * result = td->x->x->result;
  size_t ncx = td->x->x->nc, nr = td->x->x->nr;
  int * NAmeanX = td->x->x->NAme;
  size_t * nNAentriesX = td->x->x->nNAentries;
  int cosineX = td->x->x->cosine;

  double * y = td->x->y->x;
  double * weights_y = td->x->y->weights;
//  double * multMatY = td->x->y->multMat;
  size_t ncy = td->x->y->nc;
  int * NAmeanY = td->x->y->NAme;
  size_t * nNAentriesY = td->x->y->nNAentries;
  int cosineY = td->x->y->cosine;

  size_t maxDiffNA = (size_t) (td->x->x->quick * nr);

  progressCounter * pci = td->pci, * pcj = td->pcj;

  double * xx, * yy;
  double vx = 0, vy = 0;

  while (pci->i < ncx)
  {
     pthread_mutex_lock_c( td->lock, td->x->x->threaded );
     size_t i = pci->i, ii = i;
     size_t j = pcj->i, jj = j;
     do
     {
       i = ii;
       j = jj;
       jj++;
       if (jj==ncy) 
       {
         ii++;
         jj = 0;
       }
     } while ((i<ncx) && (j<ncy) && 
              ((NAmeanX[i] > 0) || (NAmeanY[j] > 0) || 
                ( (nNAentriesX[i] <= maxDiffNA) && ( nNAentriesY[j] <= maxDiffNA))));
     pci->i = ii;
     pcj->i = jj;
     pthread_mutex_unlock_c( td->lock, td->x->x->threaded );
 
     if ((i < ncx) && (j < ncy))
     {
        // Rprintf("Recalculating row %d and column %d, column size %d; cosineX: %d, cosineY: %d\n", 
        //         i, j, nr, cosineX, cosineY);
        *nNA += basic2variableCorrelation_weighted(x + i * nr, y + j * nr,
                        weights_x + i * nr, weights_y + j * nr,
                        nr, result + i + j * ncx,
                        cosineX, cosineY);
        (*nSlow)++;
     }
  }
  return NULL;
}



