/* 
Calculation of unweighted Pearson and biweght midcorrelation.

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



Some notes on handling of zero MAD:
(.) in the threaded calculations, each columns has its own NAmed, but the zeroMAD flag is one flag per thread.
    Thus, it should be zeroed out before the threaded calculation starts and checked at the end.

*/


#include "corFunctions.h"
#include "conditionalThreading.h"
#include <stdio.h>
//#include <stdlib.h>

#include <sys/time.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/BLAS.h>
#include <R_ext/libextern.h>

#include "pivot.h"

#include "corFunctions-typeDefs.h"
#include "corFunctions-utils.h"

/*========================================================================
 *
 * Short test code to see whether parallel code can be incorporated into R
 *
 * =======================================================================
 */

int nProcessors()
{
#ifdef WITH_THREADS
#ifdef _SC_NPROCESSORS_CONF
  long nProcessorsOnline = sysconf(_SC_NPROCESSORS_ONLN);
#else
  long nProcessorsOnline = 2;
#endif
#else
  long nProcessorsOnline = 1;
#endif
  return (int) nProcessorsOnline;
}

// Function to calculate suitable number of threads to use.

int useNThreads(size_t n, int nThreadsRequested)
{
#ifdef WITH_THREADS
  int nt = nThreadsRequested;
  if ((nt < 1) || (nt > MxThreads))
  {
    nt = nProcessors();
    if (nt >MxThreads) nt = MxThreads;
  }
  if (n < nt * minSizeForThreading) nt = (n/minSizeForThreading) + 1;
  return nt;
#else
  // Silence "unused argument" warning
  n = n+1;
  return 1;
#endif
}
  
//===================================================================================================

// Pearson correlation of a matrix with itself.
// This one uses matrix multiplication in BLAS to speed up calculation when there are no NA's
// and uses threading to speed up the rest of the calculation.

//===================================================================================================


// C-level correlation calculation

void cor1Fast(double * x, int * nrow, int * ncol, 
          double * weights, double * quick, 
          int * cosine, 
          double * result, int *nNA, int * err, 
          int * nThreads,
          int * verbose, int * indent)
{
  size_t nr = (size_t) *nrow, nc = (size_t) *ncol;

  char          spaces[2* *indent+1];

  for (int i=0; i<2* *indent; i++) spaces[i] = ' ';
  spaces[2* *indent] = '\0';

  *err = 0;

  size_t nNA_ext = 0;

  // Allocate space for various variables

  double * multMat;
  size_t * nNAentries;
  int *NAmean;

  // This matrix will hold preprocessed entries that can be simply multiplied together to get the
  // numerator

  if ( (multMat = (double *) malloc(nc*nr * sizeof(double)))==NULL )
  {
    *err = 1;
    Rprintf("cor1: memory allocation error. If possible, please decrease block size.\n");
    return;
  }

  // Number of NA entries in each column

  if ( (nNAentries = (size_t *) malloc(nc * sizeof(size_t)))==NULL )
  {
    free(multMat);
    *err = 1;
    Rprintf("cor1: memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  // Flag indicating whether the mean of each column is NA

  if ( (NAmean = (int *) malloc(nc * sizeof(int)))==NULL )
  {
    free(nNAentries); free(multMat); 
    *err = 1;
    Rprintf("cor1: memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  // Decide how many threads to use
  int nt = useNThreads( nc*nc, *nThreads);
  
  if (*verbose)
  {
    if (nt > 1)
      Rprintf("%s..will use %d parallel threads.\n", spaces, nt);
    else
      Rprintf("%s..will not use multithreading.\n", spaces, nt);
  }

  // double * aux[MxThreads];

  // for (int t=0; t < nt; t++)
  // {
     // if ( (aux[t] = (double *) malloc(6*nr * sizeof(double)))==NULL)
     // {
       // *err = 1;
       // Rprintf("cor1: memory allocation error. The needed block is very small... suspicious.\n");
       // for (int tt = t-1; tt>=0; tt--) free(aux[tt]);
       // free(NAmean); free(nNAentries); free(multMat);
       // return;
     // }
  // }

  // Put the general data of the correlation calculation into a structure that can be passed on to
  // threads.

  cor1ThreadData thrdInfo[MxThreads];
  for (int t = 0; t < nt; t++)
  {
     thrdInfo[t].x = x;
     thrdInfo[t].weights = weights;
     thrdInfo[t].nr = nr;
     thrdInfo[t].nc = nc;
     thrdInfo[t].multMat = multMat;
     thrdInfo[t].result = result;
     thrdInfo[t].nNAentries = nNAentries;
     thrdInfo[t].NAme = NAmean;
     thrdInfo[t].quick = *quick;
     thrdInfo[t].cosine = *cosine;
     thrdInfo[t].id = t;
     thrdInfo[t].threaded = (nt > 1);
  }

  // Column preparation (calculation of the matrix to be multiplied) in a threaded form.

  colPrepThreadData  cptd[MxThreads];
  pthread_t  thr[MxThreads];
  int       status[MxThreads];

  pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
  progressCounter pc;

  pc.i = 0;
  pc.n = nc;

  // Rprintf("Preparing columns...\n");
  for (int t=0; t<nt; t++)
  {
    cptd[t].x = &thrdInfo[t];
    cptd[t].pc = &pc;
    cptd[t].lock = &mutex1;
    status[t] = pthread_create_c(&thr[t], 
                    NULL, 
                    weights==NULL ? threadPrepColCor : threadPrepColCor_weighted,
                    (void *) &cptd[t], 
                    thrdInfo[t].threaded);
    if (status[t]!=0)
    {
      Rprintf("Error in cor(x): thread %d could not be started successfully. Error code: %d.\n%s",
              t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
      *err = 2;
    }
  }

  for (int t=0; t<nt; t++)
      if (status[t]==0) pthread_join_c(thr[t], NULL, thrdInfo[t].threaded);

  // Rprintf("done...\n");
  // Rprintf("NAmean:");
  // for (int i=0; i<nc; i++) Rprintf(" %d,", NAmean[i]);
  // Rprintf("\n");

  // The main loop is actually a matrix multiplication

  double alpha = 1.0, beta = 0.0;
  F77_NAME(dsyrk)("L", "T", ncol, nrow, & alpha, multMat, nrow, & beta, result, ncol);

  size_t nSlow = 0;

  // Rprintf("nNAentries values: ");
  // for (int i = 0; i < nc; i++) Rprintf("%d, ", nNAentries[i]);
  // Rprintf("\n");
  //
  //
  
  if (*quick < 1.0)
  {
      // Parallelized slow calculations
      slowCalcThreadData  sctd[MxThreads];
      progressCounter pci, pcj;
      pthread_mutex_t mutexSC = PTHREAD_MUTEX_INITIALIZER;

      pthread_t  thr3[MxThreads];

      pci.i = 0;
      pci.n = nc;
      pcj.i = 1;
      pcj.n = nc;

      // Rprintf("slow calculations... nt=%d\n", nt);
      for (int t=0; t<nt; t++)
      {
        sctd[t].x = &thrdInfo[t];
        sctd[t].pci = &pci;
        sctd[t].pcj = &pcj;
        sctd[t].nSlow = &nSlow;
        sctd[t].nNA = &nNA_ext;
        sctd[t].lock = &mutexSC;
        status[t] = pthread_create_c(&thr3[t], NULL, 
                weights==NULL? threadSlowCalcCor : threadSlowCalcCor_weighted, 
                (void *) &sctd[t], thrdInfo[t].threaded);
        if (status[t]!=0)
        {
          Rprintf("Error in cor(x): thread %d could not be started successfully. Error code: %d.\n%s",
                  t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
        }
      }
    
      for (int t=0; t<nt; t++)
         if (status[t]==0) pthread_join_c(thr3[t], NULL, thrdInfo[t].threaded);
    
      // Rprintf("done...\n");
    
      if (*verbose) Rprintf("%s Fraction of slow calculations: %f\n", spaces, 
                             ( (double) nSlow*2 ) / (nc*(nc-1)) );
  }
      
  // Symmetrize the result and set all rows and columns with NA means to zero
  
  symmThreadData  std[MxThreads];
  // reset the progress counter
  pc.i = 0;
  pc.n = nc;

  pthread_t  thr2[MxThreads];

  // Rprintf("symmetrizing... nt=%d\n", nt);
  for (int t=0; t<nt; t++)
  {
    std[t].x = &thrdInfo[t];
    std[t].pc = &pc;
    status[t] = pthread_create_c(&thr2[t], NULL, threadSymmetrize, (void *) &std[t], thrdInfo[t].threaded);
    if (status[t]!=0)
    {
      Rprintf("Error in cor(x): thread %d could not be started successfully. Error code: %d.\n%s",
              t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
    }
  }

  for (int t=0; t<nt; t++)
     if (status[t]==0) pthread_join_c(thr2[t], NULL, thrdInfo[t].threaded);

  // Rprintf("done... nt=%d\n", nt);
  // Here I need to recalculate results that have NA's in them.

  // for (int t=nt-1; t >= 0; t--) free(aux[t]);

  //Rprintf("End of cor1Fast (1): err = %d\n", *err);
  //*nNA = 1234;
  //Rprintf("End of cor1Fast (2): err = %d\n", *err);
  *nNA = (int) nNA_ext;
  //Rprintf("End of cor1Fast (3): err = %d\n", *err);
  free(NAmean);
  free(nNAentries);
  free(multMat);
}


//===================================================================================================

// bicorrelation of a matrix with itself.
// This one uses matrix multiplication in BLAS to speed up calculation when there are no NA's
// and is threaded to speed up the rest of the calculation.

//===================================================================================================


void bicor1Fast(double * x, int * nrow, int * ncol, double * maxPOutliers, 
            double * quick, int * fallback, int * cosine,
            double * result, int *nNA, int * err, 
            int * warn,
            int * nThreads,
            int * verbose, int * indent)
{
  size_t nr = *nrow, nc = *ncol;

  char          spaces[2* *indent+1];

  for (int i=0; i<2* *indent; i++) spaces[i] = ' ';
  spaces[2* *indent] = '\0';

  *nNA = 0;
  *warn = noWarning;
  *err = 0;

  size_t nNA_ext = 0;

  // Allocate space for various variables

  double * multMat;
  size_t * nNAentries;
  int  *NAmed;

  if ( (multMat = (double *) malloc(nc*nr * sizeof(double)))==NULL )
  {
    *err = 1;
    Rprintf("cor1: memory allocation error. If possible, please decrease block size.\n");
    return;
  }

  // Number of NA entries in each column

  if ( (nNAentries = (size_t *) malloc(nc * sizeof(size_t)))==NULL )
  {
    free(multMat);
    *err = 1;
    Rprintf("cor1: memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  // Flag indicating whether the mean of each column is NA

  if ( (NAmed = (int *) malloc(nc * sizeof(int)))==NULL )
  {
    free(nNAentries); free(multMat); 
    *err = 1;
    Rprintf("cor1: memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  // Decide how many threads to use
  int nt = useNThreads( nc*nc, *nThreads);
  
  if (*verbose)
  {
    if (nt > 1)
      Rprintf("%s..will use %d parallel threads.\n", spaces, nt);
    else
      Rprintf("%s..will not use multithreading.\n", spaces, nt);
  }

  double * aux[MxThreads];

  for (int t=0; t < nt; t++)
  {
     if ( (aux[t] = (double *) malloc(6*nr * sizeof(double)))==NULL)
     {
       *err = 1;
       Rprintf("cor1: memory allocation error. The needed block is very small... suspicious.\n");
       for (int tt = t-1; tt>=0; tt--) free(aux[tt]);
       free(NAmed); free(nNAentries); free(multMat);
       return;
     }
  }

  // Put the general data of the correlation calculation into a structure that can be passed on to
  // threads.

  cor1ThreadData thrdInfo[MxThreads];
  for (int t = 0; t < nt; t++)
  {
     thrdInfo[t].x = x;
     thrdInfo[t].weights = NULL;
     thrdInfo[t].nr = nr;
     thrdInfo[t].nc = nc;
     thrdInfo[t].multMat = multMat;
     thrdInfo[t].result = result;
     thrdInfo[t].nNAentries = nNAentries;
     thrdInfo[t].NAme = NAmed;
     thrdInfo[t].zeroMAD = 0; 
     thrdInfo[t].warn = warn;   // point the pointer 
     thrdInfo[t].aux = aux[t];
     thrdInfo[t].robust = 0;
     thrdInfo[t].fallback = *fallback;
     thrdInfo[t].quick = *quick;
     thrdInfo[t].cosine = *cosine;
     thrdInfo[t].maxPOutliers = *maxPOutliers;
     thrdInfo[t].id = t;
     thrdInfo[t].threaded = (nt > 1);
  }

  // Column preparation (calculation of the matrix to be multiplied) in a threaded form.

  colPrepThreadData  cptd[MxThreads];
  pthread_t  thr[MxThreads];
  int       status[MxThreads];

  pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
  progressCounter pc;

  pc.i = 0;
  pc.n = nc;

  // Rprintf("Preparing columns...\n");
  for (int t=0; t<nt; t++)
  {
    cptd[t].x = &thrdInfo[t];
    cptd[t].pc = &pc;
    cptd[t].lock = &mutex1;
    status[t] = pthread_create_c(&thr[t], NULL, threadPrepColBicor, (void *) &cptd[t], thrdInfo[t].threaded);
    if (status[t]!=0)
    {
      Rprintf("Error in bicor(x): thread %d could not be started successfully. Error code: %d.\n%s",
              t, status[t], "WARNING: RETURNED RESULTS WILL BE INCORRECT.");
          *err = 2;
    }
             
  }

  for (int t=0; t<nt; t++)
      if (status[t]==0) pthread_join_c(thr[t], NULL, thrdInfo[t].threaded);

  int pearson = 0;

  if (*fallback==3)
  {
    for (int t=0; t<nt; t++) if (thrdInfo[t].zeroMAD > 0)
    { 
      pearson = 1;
      if (*verbose)
        Rprintf("Warning in bicor(x): Thread %d (of %d) reported zero MAD in column %d. %s",
                t, nt, thrdInfo[t].zeroMAD, "Switching to Pearson correlation.\n");
    }
    if (pearson==1) // Re-do all column preparations using Pearson preparation.
    { 
      // Set fallback to 4 for slow calculations below.
      for (int t = 0; t < nt; t++) thrdInfo[t].fallback = 4;

      pthread_mutex_t mutex2 = PTHREAD_MUTEX_INITIALIZER;
      pc.i = 0;
      pc.n = nc;

      for (int t=0; t<nt; t++)
      {
        cptd[t].lock = &mutex2;
        status[t] = pthread_create_c(&thr[t], NULL, threadPrepColCor, (void *) &cptd[t], thrdInfo[t].threaded);
        if (status[t]!=0)
        {
          Rprintf("Error in bicor(x): thread %d could not be started successfully. Error code: %d.\n%s",
                  t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
        }
      }
      for (int t=0; t<nt; t++) 
         if (status[t]==0) pthread_join_c(thr[t], NULL, thrdInfo[t].threaded);
    }
  }

  // Rprintf("done...\n");
  // Rprintf("NAmed:");
  // for (int i=0; i<nc; i++) Rprintf(" %d,", NAmed[i]);
  // Rprintf("\n");

  // The main loop is actually a matrix multiplication

  double alpha = 1.0, beta = 0.0;
  // Rprintf("alpha: %f\n", alpha);
  F77_NAME(dsyrk)("L", "T", ncol, nrow, & alpha, multMat, nrow, & beta, result, ncol);

  // Here I need to recalculate results that have NA's in them.

  size_t nSlow = 0;

  // Rprintf("nNAentries values: ");
  // for (int i = 0; i < nc; i++) Rprintf("%d, ", nNAentries[i]);
  // Rprintf("\n");
  //
  //
  
  if (*quick < 1)
  {
      // Parallelized slow calculations
      slowCalcThreadData  sctd[MxThreads];
      progressCounter pci, pcj;
      pthread_mutex_t mutexSC = PTHREAD_MUTEX_INITIALIZER;

      pthread_t  thr3[MxThreads];

      pci.i = 0;
      pci.n = nc;
      pcj.i = 1;
      pcj.n = nc;

      // Rprintf("slow calculations... nt=%d\n", nt);
      for (int t=0; t<nt; t++)
      {
        sctd[t].x = &thrdInfo[t];
        sctd[t].pci = &pci;
        sctd[t].pcj = &pcj;
        sctd[t].nSlow = &nSlow;
        sctd[t].nNA = &nNA_ext;
        sctd[t].lock = &mutexSC;
        status[t] = pthread_create_c(&thr3[t], NULL, threadSlowCalcBicor, (void *) &sctd[t], 
                    thrdInfo[t].threaded);
        if (status[t]!=0)
        {
          Rprintf("Error in bicor(x): thread %d could not be started successfully. Error code: %d.\n%s",
                  t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
        }
      }
    
      for (int t=0; t<nt; t++)
         if (status[t]==0) pthread_join_c(thr3[t], NULL, thrdInfo[t].threaded);
    
      // Rprintf("done...\n");
    
      if (*verbose) Rprintf("%s Fraction of slow calculations: %f\n", spaces, 
                             ( (double) nSlow*2 ) / (nc*(nc-1)) );
  }
      
  // Symmetrize the result and set all rows and columns with NA means to zero
  
  symmThreadData  std[MxThreads];
  // reset the progress counter
  pc.i = 0;
  pc.n = nc;

  pthread_t  thr2[MxThreads];

  // Rprintf("symmetrizing... nt=%d\n", nt);
  for (int t=0; t<nt; t++)
  {
    std[t].x = &thrdInfo[t];
    std[t].pc = &pc;
    status[t] = pthread_create_c(&thr2[t], NULL, threadSymmetrize, (void *) &std[t], thrdInfo[t].threaded);
    if (status[t]!=0)
    {
       Rprintf("Error in bicor(x): thread %d could not be started successfully. Error code: %d.\n%s",
               t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
    }
  }

  for (int t=0; t<nt; t++)
      if (status[t]==0) pthread_join_c(thr2[t], NULL, thrdInfo[t].threaded);

  for (int t=nt-1; t >= 0; t--) free(aux[t]); 

  *nNA = (int) nNA_ext;

  free(NAmed);
  free(nNAentries);
  free(multMat);
}

//===================================================================================================
//
// Two-variable bicorrelation. Basically the same as bicor1, just must calculate the whole matrix.
// If robustX,Y is zero, the corresponding variable will be treated as in pearson correlation.
//
//===================================================================================================

void bicorFast(double * x, int * nrow, int * ncolx, double * y, int * ncoly,
           int * robustX, int * robustY, double *maxPOutliers, 
           double * quick, int * fallback,
           int * cosineX, int * cosineY, 
           double * result, int *nNA, int * err,
           int * warnX, int * warnY,
           int * nThreads,
           int * verbose, int * indent)
{
  size_t nr = *nrow, ncx = *ncolx, ncy = *ncoly;

  char          spaces[2* *indent+1];
  for (int i=0; i<2* *indent; i++) spaces[i] = ' ';
  spaces[2* *indent] = '\0';

  *warnX = noWarning;
  *warnY = noWarning;
  *err = 0;

  size_t nNA_ext = 0;

  double * multMatX, * multMatY;
  size_t * nNAentriesX, * nNAentriesY;
  int *NAmedX, *NAmedY;


  // Rprintf("nr: %d, ncx: %d, ncy: %d\n", nr, ncx, ncy);
  // Rprintf("robustX: %d, robustY: %d, cosineX: %d, cosineY: %d\n", *robustX, *robustY, *cosineX, *cosineY);
  // Rprintf("quick: %12.6f, maxPOutliers: %12.6f\n", *quick, *maxPOutliers);

  // Rprintf("Last few entries of x:\n");
  // for (int i = nr-2; i<nr; i++)
  // {
  //   for (int j = ncx-3; j<ncx; j++)
  //     Rprintf("%12.6f ", x[j*nr + i]);
  //   Rprintf("\n");
  // }

  // Rprintf("Last few entries of y:\n");
  // for (int i = nr-2; i<nr; i++)
  // {
  //   for (int j = ncy-3; j<ncy; j++)
  //     Rprintf("%12.6f ", y[j*nr + i]);
  //   Rprintf("\n");
  // }

  if ( (multMatX = (double *) malloc(ncx*nr * sizeof(double)))==NULL )
  if ( (multMatX = (double *) malloc(ncx*nr * sizeof(double)))==NULL )
  {
    *err = 1;
    Rprintf("bicor: memory allocation error. If possible, please decrease block size.\n");
    return;
  }

  if ( (multMatY = (double *) malloc(ncy*nr * sizeof(double)))==NULL )
  {
    free(multMatX);
    *err = 1;
    Rprintf("bicor: memory allocation error. If possible, please decrease block size.\n");
    return;
  }

  if ( (nNAentriesX = (size_t *) malloc(ncx * sizeof(size_t)))==NULL )
  {
    free(multMatY); free(multMatX);
    *err = 1;
    Rprintf("bicor: memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  if ( (nNAentriesY = (size_t *) malloc(ncy * sizeof(size_t)))==NULL )
  {
    free(nNAentriesX); free(multMatY); free(multMatX);
    *err = 1;
    Rprintf("bicor: memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  if ( (NAmedX = (int *) malloc(ncx * sizeof(int)))==NULL )
  {
    free(nNAentriesY); free(nNAentriesX); free(multMatY); free(multMatX);
    *err = 1;
    Rprintf("bicor: memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  if ( (NAmedY = (int *) malloc(ncy * sizeof(int)))==NULL )
  {
    free(NAmedX); free(nNAentriesY); free(nNAentriesX); free(multMatY); free(multMatX);
    *err = 1;
    Rprintf("bicor: memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  // Decide how many threads to use
  int nt = useNThreads( ncx* ncy, *nThreads);
  
  if (*verbose)
  {
    if (nt > 1)
      Rprintf("%s..will use %d parallel threads.\n", spaces, nt);
    else
      Rprintf("%s..will not use multithreading.\n", spaces, nt);
  }

  double * aux[MxThreads];

  for (int t=0; t < nt; t++)
  {
     if ( (aux[t] = (double *) malloc(6*nr * sizeof(double)))==NULL)
     {
       *err = 1;
       Rprintf("cor1: memory allocation error. The needed block is very small... suspicious.\n");
       for (int tt = t-1; tt>=0; tt--) free(aux[tt]);
       free(NAmedY); free(NAmedX); free(nNAentriesY); free(nNAentriesX); free(multMatY); free(multMatX);
       return;
     }
  }

  cor1ThreadData thrdInfoX[MxThreads];
  cor1ThreadData thrdInfoY[MxThreads];
  cor2ThreadData thrdInfo[MxThreads];

  for (int t = 0; t < nt; t++)
  {
     thrdInfoX[t].x = x;
     thrdInfoX[t].weights = NULL;
     thrdInfoX[t].nr = nr;
     thrdInfoX[t].nc = ncx;
     thrdInfoX[t].multMat = multMatX;
     thrdInfoX[t].result = result;
     thrdInfoX[t].nNAentries = nNAentriesX;
     thrdInfoX[t].NAme = NAmedX;
     thrdInfoX[t].zeroMAD = 0;
     thrdInfoX[t].aux = aux[t];
     thrdInfoX[t].robust = *robustX;
     thrdInfoX[t].fallback = *fallback;
     thrdInfoX[t].maxPOutliers = *maxPOutliers;
     thrdInfoX[t].quick = *quick;
     thrdInfoX[t].cosine = *cosineX;
     thrdInfoX[t].warn = warnX;
     thrdInfoX[t].id = t;
     thrdInfoX[t].threaded = (nt > 1);

   
     thrdInfoY[t].x = y;
     thrdInfoY[t].weights = NULL;
     thrdInfoY[t].nr = nr;
     thrdInfoY[t].nc = ncy;
     thrdInfoY[t].multMat = multMatY;
     thrdInfoY[t].result = result;
     thrdInfoY[t].nNAentries = nNAentriesY;
     thrdInfoY[t].NAme = NAmedY;
     thrdInfoY[t].zeroMAD = 0;
     thrdInfoY[t].aux = aux[t] + 3 * nr;
     thrdInfoY[t].robust = *robustY;
     thrdInfoY[t].fallback = *fallback;
     thrdInfoY[t].maxPOutliers = *maxPOutliers;
     thrdInfoY[t].quick = *quick;
     thrdInfoY[t].cosine = *cosineY;
     thrdInfoY[t].warn = warnY;
     thrdInfoY[t].id = t;
     thrdInfoY[t].threaded = (nt > 1);

     thrdInfo[t].x = thrdInfoX + t;
     thrdInfo[t].y = thrdInfoY + t;
  }

  // Prepare the multMat columns in X and Y

  // Rprintf(" ..preparing columns in x\n");

  colPrepThreadData  cptd[MxThreads];
  pthread_t  thr[MxThreads];
  int       status[MxThreads];

  progressCounter pcX, pcY;
  int pearsonX = 0, pearsonY = 0;

  // Prepare columns in X
 
  pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;

  pcX.i = 0;
  pcX.n = ncx;
  for (int t=0; t<nt; t++)
  {
    cptd[t].x = &thrdInfoX[t];
    cptd[t].pc = &pcX;
    cptd[t].lock = &mutex1;
    if (* robustX)
         status[t] = pthread_create_c(&thr[t], NULL, threadPrepColBicor, (void *) &cptd[t],
                                      thrdInfoX[t].threaded);
       else
         status[t] = pthread_create_c(&thr[t], NULL, threadPrepColCor, (void *) &cptd[t],
                                      thrdInfoX[t].threaded);
    if (status[t]!=0)
    {
       Rprintf("Error in bicor(x,y): thread %d could not be started successfully. Error code: %d.\n%s",
               t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
    }
  }
  for (int t=0; t<nt; t++)
      if (status[t]==0) pthread_join_c(thr[t], NULL, thrdInfoX[t].threaded);

  // If the fallback method is to re-do everything in Pearson, check whether any columns had zero MAD.
  if (*fallback==3)
  {
    for (int t=0; t<nt; t++) if (thrdInfoX[t].zeroMAD > 0)
    { 
      pearsonX = 1;
      if (*verbose)
        Rprintf("Warning in bicor(x, y): thread %d of %d reported zero MAD in column %d of x. %s", 
                t, nt, thrdInfoX[t].zeroMAD, "Switching to Pearson calculation for x.\n");
    }
    if (pearsonX==1) // Re-do all column preparations 
    { 
      for (int t = 0; t < nt; t++) thrdInfoX[t].fallback = 4;

      pthread_mutex_t mutex2 = PTHREAD_MUTEX_INITIALIZER;
      pcX.i = 0;
      pcX.n = ncx;

      for (int t=0; t<nt; t++)
      {
        cptd[t].lock = &mutex2;
        status[t] = pthread_create_c(&thr[t], NULL, threadPrepColCor, (void *) &cptd[t], 
                                     thrdInfoX[t].threaded);
        if (status[t]!=0)
        {
           Rprintf("Error in bicor(x,y): thread %d could not be started successfully. Error code: %d.\n%s",
                   t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
        }
      }
      for (int t=0; t<nt; t++) if (status[t]==0) pthread_join_c(thr[t], NULL, thrdInfoX[t].threaded);
    }
  }


  // Prepare columns in Y
 
  // Rprintf(" ..preparing columns in y\n");
  pthread_mutex_t mutex1Y = PTHREAD_MUTEX_INITIALIZER;

  pcY.i = 0;
  pcY.n = ncy;
  for (int t=0; t<nt; t++)
  {
    cptd[t].x = &thrdInfoY[t];
    cptd[t].pc = &pcY;
    cptd[t].lock = &mutex1Y;
    if (* robustY)
         status[t] = pthread_create_c(&thr[t], NULL, threadPrepColBicor, (void *) &cptd[t], 
                                      thrdInfoX[t].threaded);
       else
         status[t] = pthread_create_c(&thr[t], NULL, threadPrepColCor, (void *) &cptd[t],
                                     thrdInfoX[t].threaded);
    if (status[t]!=0)
    {
       Rprintf("Error in bicor(x,y): thread %d could not be started successfully. Error code: %d.\n%s",
               t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
    }
  }

  for (int t=0; t<nt; t++)
    if (status[t]==0)  pthread_join_c(thr[t], NULL, thrdInfoX[t].threaded);

  // If the fallback method is to re-do everything in Pearson, check whether any columns had zero MAD.
  if (*fallback==3)
  {
    for (int t=0; t<nt; t++) if (thrdInfoY[t].zeroMAD > 0)
    { 
      pearsonY = 1;
      if (*verbose)
        Rprintf("Warning in bicor(x, y): thread %d of %d reported zero MAD in column %d of y. %s", 
                t, nt, thrdInfoY[t].zeroMAD, "Switching to Pearson calculation for y.\n");
    }
    if (pearsonY==1) // Re-do all column preparations 
    { 
      for (int t = 0; t < nt; t++) thrdInfoY[t].fallback = 4;

      pthread_mutex_t mutex2Y = PTHREAD_MUTEX_INITIALIZER;
      pcY.i = 0;
      pcY.n = ncy;

      for (int t=0; t<nt; t++)
      {
      //  Rprintf("Starting pearson re-calculation in thread %d of %d.\n", t, nt);
        cptd[t].lock = &mutex2Y;
        status[t] = pthread_create_c(&thr[t], NULL, threadPrepColCor, (void *) &cptd[t], thrdInfoX[t].threaded);
        if (status[t]!=0)
        {
           Rprintf("Error in bicor(x,y): thread %d could not be started successfully. Error code: %d.\n%s",
                   t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
        }
      }
      for (int t=0; t<nt; t++) if (status[t]==0) pthread_join_c(thr[t], NULL, thrdInfoX[t].threaded);
    }
  }

   // Rprintf("nNAentriesX:");
   // for (int i=0; i<ncx; i++) Rprintf(" %d,", nNAentriesX[i]);
   // Rprintf("\n");
   // Rprintf("nNAentriesY:");
   // for (int i=0; i<ncx; i++) Rprintf(" %d,", nNAentriesY[i]);
   // Rprintf("\n");
 

  // The main calculation: matrix multiplication
  
  double alpha = 1.0, beta = 0.0;
  F77_NAME(dgemm)("T", "N", ncolx, ncoly, nrow, & alpha, multMatX, nrow, multMatY, nrow, & beta, result, ncolx);

  // Rprintf("matrix multiplication result:\n");
  // for (int i=0; i<ncx; i++)
  // {
  //   for (int j=0; j<ncy; j++) Rprintf(" %12.6f ", result[i + ncx*j]);
  //   Rprintf("\n");
 //  }
/*
  Rprintf("Last few entries of result just after multiplication:\n");
  for (int i = 0; i<ncx; i++)
  {
    for (int j = 0; j<ncy; j++)
      Rprintf("%12.6f ", result[j*ncx + i]);
    Rprintf("\n");
  }

  Rprintf("multMatX:\n");
  for (int i=0; i<nr; i++)
  {
    for (int j=0; j<ncx; j++) Rprintf(" %12.6f ", multMatX[i + nr*j]);
    Rprintf("\n");
  }
 
  Rprintf("multMatY:\n");
  for (int i=0; i<nr; i++)
  {
    for (int j=0; j<ncy; j++) Rprintf(" %12.6f ", multMatY[i + nr*j]);
    Rprintf("\n");
  }

*/

  // Remedial calculations

  size_t nSlow = 0;
  if (*quick < 1.0)
  {
      slowCalc2ThreadData  sctd[MxThreads];
      pthread_mutex_t mutexSC = PTHREAD_MUTEX_INITIALIZER;

      pthread_t  thr3[MxThreads];

      pcX.i = 0; pcY.i = 0;

      // Rprintf("slow calculations... nt=%d\n", nt);
      for (int t=0; t<nt; t++)
      {
        sctd[t].x = &thrdInfo[t];
        sctd[t].pci = &pcX;
        sctd[t].pcj = &pcY;
        sctd[t].nSlow = &nSlow;
        sctd[t].nNA = &nNA_ext;
        sctd[t].lock = &mutexSC;
        sctd[t].quick = *quick;
        status[t] = pthread_create_c(&thr3[t], NULL, threadSlowCalcBicor2, (void *) &sctd[t], 
                                     thrdInfoX[t].threaded);
        if (status[t]!=0)
        {
           Rprintf("Error in bicor(x,y): thread %d could not be started successfully. Error code: %d.\n%s",
                   t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
        }
      }

      for (int t=0; t<nt; t++)
        if (status[t]==0)  pthread_join_c(thr3[t], NULL, thrdInfoX[t].threaded);

      if (*verbose) Rprintf("%s Fraction of slow calculations: %f\n", spaces,
                             ( (double) nSlow) / (ncx*ncy) );
  }

  // NA out all rows and columns that need it and check for values outside of [-1, 1]

  NA2ThreadData  natd[MxThreads];
  // reset the progress counter
  pcX.i = 0;
  pcY.i = 0;

  pthread_t  thr2[MxThreads];

  for (int t=0; t<nt; t++)
  {
    natd[t].x = &thrdInfo[t];
    natd[t].pci = &pcX;
    natd[t].pcj = &pcY;
    status[t] = pthread_create_c(&thr2[t], NULL, threadNAing, (void *) &natd[t], thrdInfoX[t].threaded);
    if (status[t]!=0)
    {
       Rprintf("Error in bicor(x,y): thread %d could not be started successfully. Error code: %d.\n%s",
               t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
    }
  }

  for (int t=0; t<nt; t++)
     if (status[t]==0) pthread_join_c(thr2[t], NULL, thrdInfoX[t].threaded);

  *nNA = (int) nNA_ext;

  
  // Rprintf("Last few entries of result:\n");
  // for (int i = ncx-4; i<ncx; i++)
  // {
  //  for (int j = ncy-4; j<ncy; j++)
  //    Rprintf("%12.6f ", result[j*ncx + i]);
  //  Rprintf("\n");
  //}

  // Clean up

  for (int t=nt-1; t >= 0; t--) free(aux[t]);
  free(NAmedY);
  free(NAmedX);
  free(nNAentriesY);
  free(nNAentriesX);
  free(multMatY);
  free(multMatX);
}

/*======================================================================================================
 *
 * corFast: fast correlation of 2 matrices
 *
 *======================================================================================================

One important note: if weights_x is not NULL, weights_y is also assumed to be valid.

*/

void corFast(double * x, int * nrow, int * ncolx, double * y, int * ncoly,
           double * weights_x, double * weights_y,
           double * quick, 
           int * cosineX, int * cosineY, 
           double * result, int *nNA, int * err,
           int * nThreads,
           int * verbose, int * indent)
{
  size_t nr = *nrow, ncx = *ncolx, ncy = *ncoly;

  char          spaces[2* *indent+1];
  for (int i=0; i<2* *indent; i++) spaces[i] = ' ';
  spaces[2* *indent] = '\0';

  size_t nNA_ext = 0;
  *err = 0;

  double * multMatX, * multMatY;
  size_t * nNAentriesX, * nNAentriesY;
  int *NAmeanX, *NAmeanY;

  if ( (weights_x == NULL) != (weights_y == NULL))
  {
    *err = 2;
    error("corFast: weights_x and weights_y must both be either NULL or non-NULL.\n");
  }

  if ( (multMatX = (double *) malloc(ncx*nr * sizeof(double)))==NULL )
  {
    *err = 1;
    Rprintf("cor(x,y): memory allocation error. If possible, please decrease block size.\n");
    return;
  }

  if ( (multMatY = (double *) malloc(ncy*nr * sizeof(double)))==NULL )
  {
    free(multMatX);
    *err = 1;
    Rprintf("cor(x,y): memory allocation error. If possible, please decrease block size.\n");
    return;
  }

  if ( (nNAentriesX = (size_t *) malloc(ncx * sizeof(size_t)))==NULL )
  {
    free(multMatY); free(multMatX);
    *err = 1;
    Rprintf("cor(x,y): memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  if ( (nNAentriesY = (size_t *) malloc(ncy * sizeof(size_t)))==NULL )
  {
    free(nNAentriesX); free(multMatY); free(multMatX);
    *err = 1;
    Rprintf("cor(x,y): memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  if ( (NAmeanX = (int *) malloc(ncx * sizeof(int)))==NULL )
  {
    free(nNAentriesY); free(nNAentriesX); free(multMatY); free(multMatX);
    *err = 1;
    Rprintf("cor(x,y): memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  if ( (NAmeanY = (int *) malloc(ncy * sizeof(int)))==NULL )
  {
    free(NAmeanX); free(nNAentriesY); free(nNAentriesX); free(multMatY); free(multMatX);
    *err = 1;
    Rprintf("cor(x,y): memory allocation error. The needed block is relatively small... suspicious.\n");
    return;
  }

  // Decide how many threads to use
  int nt = useNThreads( ncx* ncy, *nThreads);
  
  if (*verbose)
  {
    if (nt > 1)
      Rprintf("%s..will use %d parallel threads.\n", spaces, nt);
    else
      Rprintf("%s..will not use multithreading.\n", spaces, nt);
  }

  cor1ThreadData thrdInfoX[MxThreads];
  cor1ThreadData thrdInfoY[MxThreads];
  cor2ThreadData thrdInfo[MxThreads];

  for (int t = 0; t < nt; t++)
  {
     thrdInfoX[t].x = x;
     thrdInfoX[t].weights = weights_x;
     thrdInfoX[t].nr = nr;
     thrdInfoX[t].nc = ncx;
     thrdInfoX[t].multMat = multMatX;
     thrdInfoX[t].result = result;
     thrdInfoX[t].nNAentries = nNAentriesX;
     thrdInfoX[t].NAme = NAmeanX;
     thrdInfoX[t].quick = *quick;
     thrdInfoX[t].cosine = *cosineX;
     thrdInfoX[t].maxPOutliers = 1;
     thrdInfoX[t].id = t;
     thrdInfoX[t].threaded = (nt > 1);
   
     thrdInfoY[t].x = y;
     thrdInfoY[t].weights = weights_y;
     thrdInfoY[t].nr = nr;
     thrdInfoY[t].nc = ncy;
     thrdInfoY[t].multMat = multMatY;
     thrdInfoY[t].result = result;
     thrdInfoY[t].nNAentries = nNAentriesY;
     thrdInfoY[t].NAme = NAmeanY;
     thrdInfoY[t].quick = *quick;
     thrdInfoY[t].cosine = *cosineY;
     thrdInfoY[t].maxPOutliers = 1;
     thrdInfoY[t].id = t;
     thrdInfoY[t].threaded = (nt > 1);

     thrdInfo[t].x = thrdInfoX + t;
     thrdInfo[t].y = thrdInfoY + t;
  }

  // Prepare the multMat columns in X and Y

  colPrepThreadData  cptd[MxThreads];
  pthread_t  thr[MxThreads];
  int       status[MxThreads];

  progressCounter pcX, pcY;

  // Prepare columns in X
 
  pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;

  pcX.i = 0;
  pcX.n = ncx;
  for (int t=0; t<nt; t++)
  {
    cptd[t].x = &thrdInfoX[t];
    cptd[t].pc = &pcX;
    cptd[t].lock = &mutex1;
    status[t] = pthread_create_c(&thr[t], NULL, 
                   weights_x==NULL ? threadPrepColCor : threadPrepColCor_weighted, 
                   (void *) &cptd[t], thrdInfoX[t].threaded);
    if (status[t]!=0)
    {
       Rprintf("Error in cor(x,y): thread %d could not be started successfully. Error code: %d.\n%s",
               t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
    }
  }
  for (int t=0; t<nt; t++)
    if (status[t]==0) pthread_join_c(thr[t], NULL, thrdInfoX[t].threaded);

  // Prepare columns in Y
 
  pthread_mutex_t mutex1Y = PTHREAD_MUTEX_INITIALIZER;

  pcY.i = 0;
  pcY.n = ncy;
  for (int t=0; t<nt; t++)
  {
    cptd[t].x = &thrdInfoY[t];
    cptd[t].pc = &pcY;
    cptd[t].lock = &mutex1Y;
    status[t] = pthread_create_c(&thr[t], NULL, 
           weights_y==NULL ? threadPrepColCor: threadPrepColCor_weighted, 
           (void *) &cptd[t], thrdInfoX[t].threaded);
    if (status[t]!=0)
    {
       Rprintf("Error in cor(x,y): thread %d could not be started successfully. Error code: %d.\n%s",
               t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
    }
  }

  for (int t=0; t<nt; t++)
    if (status[t]==0) pthread_join_c(thr[t], NULL, thrdInfoX[t].threaded);

  //Rprintf("multMatX:\n");
  //for (int i=0; i<nr; i++)
  //{
    //for (int j=0; j<ncx; j++) Rprintf(" %12.6f ", multMatX[i + nr*j]);
    //Rprintf("\n");
  //}
 
  //Rprintf("multMatY:\n");
  //for (int i=0; i<nr; i++)
  //{
    //for (int j=0; j<ncy; j++) Rprintf(" %12.6f ", multMatY[i + nr*j]);
    //Rprintf("\n");
  //}

  // Rprintf("nNAentriesX:");
  // for (int i=0; i<ncx; i++) Rprintf(" %d,", nNAentriesX[i]);
  // Rprintf("\n");
  // Rprintf("nNAentriesY:");
  // for (int i=0; i<ncx; i++) Rprintf(" %d,", nNAentriesY[i]);
  // Rprintf("\n");
 

  // The main calculation: matrix multiplication
  
  double alpha = 1.0, beta = 0.0;
  F77_NAME(dgemm)("T", "N", ncolx, ncoly, nrow, & alpha, multMatX, nrow, multMatY, nrow, & beta, result, ncolx);

  // Remedial calculations

  size_t nSlow = 0;
  if (*quick < 1.0)
  {
      slowCalc2ThreadData  sctd[MxThreads];
      pthread_mutex_t mutexSC = PTHREAD_MUTEX_INITIALIZER;

      pthread_t  thr3[MxThreads];

      pcX.i = 0; pcY.i = 0;

      // Rprintf("slow calculations... nt=%d\n", nt);
      for (int t=0; t<nt; t++)
      {
        sctd[t].x = &thrdInfo[t];
        sctd[t].pci = &pcX;
        sctd[t].pcj = &pcY;
        sctd[t].nSlow = &nSlow;
        sctd[t].nNA = &nNA_ext;
        sctd[t].lock = &mutexSC;
        sctd[t].quick = *quick;
        status[t] = pthread_create_c(&thr3[t], NULL, 
                 weights_x==NULL ? threadSlowCalcCor2 : threadSlowCalcCor2_weighted,
                 (void *) &sctd[t], thrdInfoX[t].threaded);
        if (status[t]!=0)
        {
           Rprintf("Error in cor(x,y): thread %d could not be started successfully. Error code: %d.\n%s",
                   t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
        }
      }

      for (int t=0; t<nt; t++)
          if (status[t]==0) pthread_join_c(thr3[t], NULL, thrdInfoX[t].threaded);

      if (*verbose) Rprintf("%s Fraction of slow calculations: %f\n", spaces,
                             ( (double) nSlow) / (ncx*ncy) );
  }

  // NA out all rows and columns that need it and check for values outside of [-1, 1]
  //
  NA2ThreadData  natd[MxThreads];
  // reset the progress counters
  pcX.i = 0;
  pcY.i = 0;

  pthread_t  thr2[MxThreads];

  for (int t=0; t<nt; t++)
  {
    natd[t].x = &thrdInfo[t];
    natd[t].pci = &pcX;
    natd[t].pcj = &pcY;
    status[t] = pthread_create_c(&thr2[t], NULL, threadNAing, (void *) &natd[t], thrdInfoX[t].threaded);
    if (status[t]!=0)
    {
       Rprintf("Error in cor(x,y): thread %d could not be started successfully. Error code: %d.\n%s",
               t, status[t], "*** WARNING: RETURNED RESULTS WILL BE INCORRECT. ***");
          *err = 2;
    }
  }

  for (int t=0; t<nt; t++)
    if (status[t]==0)  pthread_join_c(thr2[t], NULL, thrdInfoX[t].threaded);

  *nNA = (int) nNA_ext;

  // clean up and return

  free(NAmeanY);
  free(NAmeanX);
  free(nNAentriesY);
  free(nNAentriesX);
  free(multMatY);
  free(multMatX);
}

//===================================================================================================
//
// Two-variable Pearson correlation. 
//
//===================================================================================================

SEXP bicor2_call(SEXP x_s, SEXP y_s,
                 SEXP robustX_s, SEXP robustY_s,
                 SEXP maxPOutliers_s, SEXP quick_s, 
                 SEXP fallback_s, 
                 SEXP cosineX_s, SEXP cosineY_s,
                 SEXP nNA_s, SEXP err_s, 
                 SEXP warnX_s, SEXP warnY_s,
                 SEXP nThreads_s, SEXP verbose_s, SEXP indent_s)
{
  SEXP dimX, dimY, cor_s; 

  int nr, ncx, ncy;
  int *cosineX, *cosineY;
  int *err, *nThreads, *verbose, *indent, *fallback;
  int *warnX, *warnY, *robustX, *robustY;
  int *nNA;

  double *x, *y, *corMat, *quick, *maxPOutliers;

  /* Get dimensions of 'x'. */
  PROTECT(dimX = getAttrib(x_s, R_DimSymbol));
  nr = INTEGER(dimX)[0];
  ncx = INTEGER(dimX)[1];
  // Rprintf("Matrix x dimensions: %d %d\n", nr, ncx);
  /* Get dimensions of 'y'. */
  PROTECT(dimY = getAttrib(y_s, R_DimSymbol));
  ncy = INTEGER(dimY)[1];
  // Rprintf("Matrix y dimensions: %d %d\n", INTEGER(dimY)[0], ncy);

  x = REAL(x_s);
  y = REAL(y_s);

  // Rprintf("First three elements of x: %f %f %f\n", x[0], x[1], x[2]);

  quick = REAL(quick_s);
  maxPOutliers = REAL(maxPOutliers_s);
  cosineX = INTEGER(cosineX_s);
  cosineY = INTEGER(cosineY_s);
  robustX = INTEGER(robustX_s);
  robustY = INTEGER(robustY_s);
  nThreads = INTEGER(nThreads_s);
  verbose = INTEGER(verbose_s);
  indent = INTEGER(indent_s);
  fallback = INTEGER(fallback_s);

  // Allocate space for the result
  PROTECT(cor_s = allocMatrix(REALSXP, ncx, ncy));
  // PROTECT(nNA_s = allocVector(REALSXP, 1));
  // PROTECT(err_s = allocVector(REALSXP, 1));

  corMat = REAL(cor_s);
  nNA = INTEGER(nNA_s);
  err = INTEGER(err_s);
  warnX = INTEGER(warnX_s);
  warnY = INTEGER(warnY_s);

  bicorFast(x, &nr, &ncx, y, &ncy,
           robustX, robustY,
           maxPOutliers, quick, fallback, 
           cosineX, cosineY, 
           corMat, nNA, err,
           warnX, warnY, 
           nThreads, verbose, indent);

  // Rprintf("Done...\n");
  UNPROTECT(3);
  return cor_s;
} 

SEXP corFast_call(SEXP x_s, SEXP y_s,
                 SEXP weights_x_s, SEXP weights_y_s,
                 SEXP quick_s, 
                 SEXP cosineX_s, SEXP cosineY_s,
                 SEXP nNA_s, SEXP err_s, 
                 SEXP nThreads_s, SEXP verbose_s, SEXP indent_s)
{
  SEXP dimX, dimY, cor_s; 

  int nr, ncx, ncy;
  int *cosineX, *cosineY;
  int *err, *nThreads, *verbose, *indent;
  int *nNA;

  double *x, *y, *weights_x, *weights_y, *corMat, *quick;

  /* Get dimensions of 'x'. */
  PROTECT(dimX = getAttrib(x_s, R_DimSymbol));
  nr = INTEGER(dimX)[0];
  ncx = INTEGER(dimX)[1];
  // Rprintf("Matrix dimensions: %d %d\n", nr, nc);
  /* Get dimensions of 'y'. */
  PROTECT(dimY = getAttrib(y_s, R_DimSymbol));
  ncy = INTEGER(dimY)[1];

  x = REAL(x_s);
  y = REAL(y_s);

  if (isNull(weights_x_s)) weights_x = NULL; else weights_x = REAL(weights_x_s);
  if (isNull(weights_y_s)) weights_y = NULL; else weights_y = REAL(weights_y_s);

  // Rprintf("First three elements of x: %f %f %f\n", x[0], x[1], x[2]);

  quick = REAL(quick_s);
  cosineX = INTEGER(cosineX_s);
  cosineY = INTEGER(cosineY_s);
  nThreads = INTEGER(nThreads_s);
  verbose = INTEGER(verbose_s);
  indent = INTEGER(indent_s);

  // Allocate space for the result
  PROTECT(cor_s = allocMatrix(REALSXP, ncx, ncy));
  // PROTECT(nNA_s = allocVector(REALSXP, 1));
  // PROTECT(err_s = allocVector(REALSXP, 1));

  corMat = REAL(cor_s);
  nNA = INTEGER(nNA_s);
  err = INTEGER(err_s);

  // Rprintf("Calling cor1Fast...\n");
  corFast(x, &nr, &ncx, y, &ncy,
          weights_x, weights_y,
          quick, 
          cosineX, cosineY, 
          corMat, nNA, err,
          nThreads, verbose, indent);

  // Rprintf("Done...\n");
  UNPROTECT(3);
  return cor_s;
} 

// Re-write cor1Fast as a function that can be called using .Call
// Since I don't know how to create and fill lists in C code, I will for now return the nNA and err results
// via supplied arguments. Not ideal but will do.

SEXP cor1Fast_call(SEXP x_s, SEXP weights_s, SEXP quick_s, SEXP cosine_s,
                   SEXP nNA_s, SEXP err_s,
                   SEXP nThreads_s, SEXP verbose_s, SEXP indent_s)
{
  SEXP dim, cor_s; 
  // SEXP out, nNA_s, err_s;

  int nr, nc;
  int *cosine, *err, *nThreads, *verbose, *indent;
  int *nNA;

  double *x, *weights, *corMat, *quick;

  /* Get dimensions of 'x'. */
  PROTECT(dim = getAttrib(x_s, R_DimSymbol));
  nr = INTEGER(dim)[0];
  nc = INTEGER(dim)[1];
  // Rprintf("Matrix dimensions: %d %d\n", nr, nc);

  x = REAL(x_s);

  if (isNull(weights_s)) weights = NULL; else weights = REAL(weights_s);

  // Rprintf("First three elements of x: %f %f %f\n", x[0], x[1], x[2]);

  quick = REAL(quick_s);
  cosine = INTEGER(cosine_s);
  nThreads = INTEGER(nThreads_s);
  verbose = INTEGER(verbose_s);
  indent = INTEGER(indent_s);

  // Allocate space for the result
  PROTECT(cor_s = allocMatrix(REALSXP, nc, nc));
  // PROTECT(nNA_s = allocVector(REALSXP, 1));
  // PROTECT(err_s = allocVector(REALSXP, 1));

  corMat = REAL(cor_s);
  nNA = INTEGER(nNA_s);
  err = INTEGER(err_s);

  // Rprintf("Difference of nNA and err pointers: %d\n", (int) (err - nNA));

  // Rprintf("Calling cor1Fast...\n");
  cor1Fast(x, &nr, &nc, weights, quick, cosine, 
           corMat, nNA, err,
           nThreads, verbose, indent);

  // Rprintf("Done...\n");
  UNPROTECT(2);
  return cor_s;
} 

SEXP bicor1_call(SEXP x_s, 
                 SEXP maxPOutliers_s, SEXP quick_s, 
                 SEXP fallback_s, SEXP cosine_s,
                 SEXP nNA_s, SEXP err_s, SEXP warn_s,
                 SEXP nThreads_s, SEXP verbose_s, SEXP indent_s)
{
  SEXP dim, cor_s; 
  // SEXP out, nNA_s, err_s;

  int nr, nc;
  int *cosine, *err, *nThreads, *verbose, *indent, *fallback;
  int *nNA, *warn;

  double *x, *corMat, *quick, *maxPOutliers;

  /* Get dimensions of 'x'. */
  PROTECT(dim = getAttrib(x_s, R_DimSymbol));
  nr = INTEGER(dim)[0];
  nc = INTEGER(dim)[1];
  // Rprintf("Matrix dimensions: %d %d\n", nr, nc);

  x = REAL(x_s);

  // Rprintf("First three elements of x: %f %f %f\n", x[0], x[1], x[2]);

  quick = REAL(quick_s);
  maxPOutliers = REAL(maxPOutliers_s);
  cosine = INTEGER(cosine_s);
  nThreads = INTEGER(nThreads_s);
  verbose = INTEGER(verbose_s);
  indent = INTEGER(indent_s);
  fallback = INTEGER(fallback_s);

  // Allocate space for the result
  PROTECT(cor_s = allocMatrix(REALSXP, nc, nc));
  // PROTECT(nNA_s = allocVector(REALSXP, 1));
  // PROTECT(err_s = allocVector(REALSXP, 1));

  corMat = REAL(cor_s);
  nNA = INTEGER(nNA_s);
  err = INTEGER(err_s);
  warn = INTEGER(warn_s);

  // Rprintf("Calling cor1Fast...\n");
  bicor1Fast(x, &nr, &nc, 
           maxPOutliers, quick, 
           fallback, cosine, 
           corMat, nNA, err,
           warn, 
           nThreads, verbose, indent);

  // Rprintf("Done...\n");
  UNPROTECT(2);
  return cor_s;
} 

