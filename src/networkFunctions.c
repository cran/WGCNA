/*
 *
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

#include "networkFunctions.h"

/* =============================================================================================
 *
 * tomSimilarity
 *
 * =============================================================================================

     * expr: the expression matrix, array of nSamples * nGenes doubles
             (NA's allowed)
     corType: CorTypePearson, CorTypeBicor, CorTypeSpearman
     tomType: TomTypeNone, TomTypeUnsigned, TomTypeSigned

*/

enum { CorTypePearson = 0, CorTypeBicor = 1, CorTypeSpearman = 2 };
enum { TomTypeNone = 0, TomTypeUnsigned = 1, TomTypeSigned = 2, TomTypeSignedNowick = 3,
       TomTypeUnsigned2 = 4, TomTypeSigned2 = 5, TomTypeSignedNowick2 = 6 };
enum { TomDenomMin = 0, TomDenomMean = 1 };
enum { AdjTypeUnsigned = 0, AdjTypeSigned = 1, AdjTypeHybrid = 2, AdjTypeUnsignedKeepSign = 3 };

#define MxStr	200
typedef char	plString[MxStr];
const plString 	AdjErrors[] = {"No error. Just a placeholder.",
                               "Standard deviation of some genes is zero.",
                               "Unrecognized correlation type.",
                               "Unrecognized adjacency type."};

void tomSimilarityFromAdj(double * adj, int * nGenes,
                   int * tomType, int * denomType,
                   int * suppressTOMForZeroAdj, int * suppressNegativeTOM,
                   int * useInternalMatrixAlgebra, 
                   double * tom,
                   int * verbose, int * indent);

//===============================================================================================

// adjacency 

//===============================================================================================

void adjacency(double * expr, double * weights, int nSamples, int nGenes, int corType, int adjType, double power, 
               double maxPOutliers, double quick, int fallback, int cosine, 
               int replaceMissing, 
               double * adj, int * errCode, int *warn, int * nThreads, int verbose, int indent)
{
  size_t nElems = ((size_t) nGenes) * ((size_t ) nGenes);
  int nNA = 0;
  double replacementValue = 0;

  // Rprintf("Received nGenes: %d, nSamples: %d\n", nGenes, nSamples);

  // Rprintf("adjacency: adjType: %d\n", adjType);
  // Rprintf("adjacency: replaceMissing: %d\n", replaceMissing);
  int err = 0;
  switch (corType)
  {
     case CorTypePearson :
        // Rprintf("Calling cor_pairwise1...");
        cor1Fast(expr, &nSamples, &nGenes, weights, &quick, &cosine, adj, &nNA, &err, nThreads, &verbose, &indent);
        // Rprintf("..done.\n");
        if ((nNA > 0) && (!replaceMissing))
        {
          * errCode = 1;
          return;
        }
        break;
     case CorTypeBicor :
        // Rprintf("Calling bicor1...");
        bicor1Fast(expr, &nSamples, &nGenes, &maxPOutliers, &quick, &fallback, &cosine, adj, &nNA, &err,
                   warn, nThreads, &verbose, &indent);
        // Rprintf("..done.\n");
        if ((nNA > 0) && (!replaceMissing))
        {
          // Rprintf("nNA: %d\n", nNA);
          * errCode = 1;
          return;
        }
        if (err>0)
        {
          // Rprintf("bicor1 returned err: %d\n", err);
          * errCode = 3;
          return;
        }
        break;
     default : 
        * errCode = 2;
        return;
  }

  if ((*errCode==1) && replaceMissing) 
  {
    Rprintf("Replacing missing adjacency values.\n");
    *errCode = 0;
    if (adjType==AdjTypeSigned) replacementValue = -1;
    for (size_t i=0; i < nElems; i++)
       if (ISNAN(adj[i])) adj[i] = replacementValue;
    
  }

  // Rprintf("ADJ 1\n");
  

  switch (adjType) 
  {
    case AdjTypeUnsigned :
       for (size_t i=0; i < nElems; i++)
          adj[i] = pow(fabs(adj[i]), power);
       break;
    case AdjTypeUnsignedKeepSign : 
       for (size_t i=0; i < nElems; i++)
         adj[i] = (signbit(adj[i])? -1: 1) * pow(fabs(adj[i]), power);
       break;
    case AdjTypeSigned :
       for (size_t i=0; i < nElems; i++)
          adj[i] = pow((1+adj[i])/2, power);
       break;
    case AdjTypeHybrid :
       for (size_t i=0; i < nElems; i++)
          adj[i] = adj[i] > 0 ? pow(adj[i], power) : 0;
       break;
    default : 
       * errCode = 3;
  }
       
}
    
void testAdjacency(double * expr, double * weights, int * nSamples, int * nGenes, int * corType, int * adjType, 
                   double * power, double * maxPOutliers, double * quick, int * fallback, int * cosine, 
                   double * adj, 
                   int * errCode, int * warn, int * nThreads)
{
 adjacency(expr, weights, * nSamples, * nGenes, * corType, * adjType, * power, 
           *maxPOutliers, *quick, *fallback, *cosine, 0,
           adj, errCode, warn, nThreads, 1, 0);
}


/********************************************************************************************
 *
 * checkAvailableMemory
 *
 ********************************************************************************************/

size_t checkAvailableMemory(void)
{
  size_t guess;
  if ( sizeof (size_t)==4 ) 
     guess = 16384;  // 2^14
  else
     guess = 131072;  // power of 2 nearest to 100k

  int tooLarge = 1;
  double * pt;
  while ( (tooLarge) && (guess > 1000))
  {
     // Rprintf("trying matrix of size %d\n", guess);
     tooLarge = ( (pt=malloc(guess*guess*sizeof(double))) == NULL );
     if (tooLarge) guess = (guess * 3) / 4;
     // Rprintf("next size will be %d\n", guess);
  }

  if (!tooLarge) free(pt);   

  // Rprintf("Returning %d.\n", guess * guess);

  return guess*guess;
}

  
//====================================================================================================

// TOM similarity from adjacency

//====================================================================================================

void tomSimilarityFromAdj(double * adj, int * nGenes, 
                   int * tomType, int * denomType, 
                   int * suppressTOMForZeroAdj, 
                   int * suppressNegativeTOM,
                   int * useInternalMatrixAlgebra, 
                   double * tom, 
                   int * verbose, int * indent)
{
  size_t		ng = (size_t) *nGenes;
  // int 		err = 0;

  char		spaces[2* *indent+1];

  for (int i=0; i<2* *indent; i++) spaces[i] = ' ';
  spaces[2* *indent] = '\0';

  double * conn;
  if ( (conn = malloc(ng * sizeof(double))) == NULL)
    error("Memmory allocation error (connectivity)");

  if (*verbose > 0) Rprintf("%s..connectivity..\n", spaces);
  double *tom2 = tom;
  for (size_t gene = 0; gene < ng; gene++)
  {
    double * adj2 = (adj + gene * ng);
    // set diagonal to 1.
    *(adj2 + gene) = 1;
    // Calculate connectivity
    double sum = 0.0;
    for (size_t g2 = 0; g2 < ng; g2++)
    {
      sum += fabs(adj2[g2]);
      *tom2 = 0;
      tom2++;
    }
    conn[gene] = sum;
  }

  if (*useInternalMatrixAlgebra > 0)
  {
    if (*verbose > 0) Rprintf("%s..matrix multiplication (custom code)..\n", spaces);
    squareSymmetricMatrix(adj, *nGenes, tom);
  } else {
    if (*verbose > 0) Rprintf("%s..matrix multiplication (system BLAS)..\n", spaces);
    double alpha = 1.0, beta = 0.0;
    F77_NAME(dsyrk)("L", "N", nGenes, nGenes, & alpha, adj, nGenes, 
           & beta, tom, nGenes FCONE FCONE);
  }

  if (*verbose > 0) Rprintf("%s..normalization..\n", spaces);
  // Rprintf("Using denomType %d\n", *denomType);
  tom2 = tom;
  double * adj2;
  size_t ng1 = ng-1;
  size_t nAbove1 = 0, nSuppressed = 0;
  if (*suppressTOMForZeroAdj)
  {
    Rprintf("%s..will suppress TOM for pairs of nodes with zero adjacency.\n", spaces);
  }
  int form = *tomType > TomTypeSignedNowick;

  switch (* tomType)
  {
    case TomTypeUnsigned:
    case TomTypeUnsigned2:
      for (size_t j=0; j< ng1; j++)
      {
        tom2 = tom + (ng+1)*j + 1;  
        adj2 = adj + (ng+1)*j + 1;
        for (size_t i=j+1; i< ng; i++)
        {
          double den1;
          if ((* denomType) == TomDenomMin)
             den1 = fmin(conn[i], conn[j]);
          else
             den1 = (conn[i] + conn[j])/2;
          double den = den1 - * adj2;
          if (form > 0)
          {
            double r;
            if (den <= 1) r = 0; else r = (*tom2 - *adj2 * 2)/(den-1);
            *tom2 = (*adj2 + r)/2;
          } else {
            if (den==0)
              *tom2 = 0;
            else
              *tom2 = ( *tom2 - *adj2) / den ;
          }
          *(tom + ng*i + j) = *tom2;
          if (*tom2 > 1) nAbove1++;
          tom2++;
          adj2++;
        }
      }
      break;
    case TomTypeSigned:
    case TomTypeSigned2:
      for (size_t j=0; j < ng1; j++)
      {
        tom2 = tom + (ng+1)*j + 1;  
        adj2 = adj + (ng+1)*j + 1;
        for (size_t i=j+1; i< ng; i++)
        {
          if ((*suppressTOMForZeroAdj == 0) || (*adj2 > 0))
          {
            double den1;
            if ((* denomType) == TomDenomMin)
               den1 = fmin(conn[i], conn[j]);
            else
               den1 = (conn[i] + conn[j])/2;
            double den = den1 - fabs(*adj2);
            if (form > 0)
            {
              double r;
              if (den <= 1) r = 0; else r = (*tom2 - *adj2 * 2)/(den-1);
              *tom2 = fabs(*adj2 + r)/2;
            } else {
              if (den==0)
                *tom2 = 0;
              else
                *tom2 = fabs( *tom2 - *adj2) / den ;
            }
            *(tom + ng*i + j) = *tom2;
            if (*tom2 > 1) 
            {
              Rprintf("TOM greater than 1: actual value: %f, i: %d, j: %d\n", *tom2, i, j);
              nAbove1++;
            }
          } else {
            *tom2 = 0;
            *(tom + ng*i + j) = 0;
            nSuppressed++;
          }
          tom2++;
          adj2++;
        }
      }
      break;
    case TomTypeSignedNowick: 
    case TomTypeSignedNowick2:  // Differs from the above only in one missing fabs and potential suppression of negative values
      // Rprintf("Calculating Nowick-type TOM. SuppressNegativeTOM: %d\n", *suppressNegativeTOM);
      for (size_t j=0; j < ng1; j++)
      {
        tom2 = tom + (ng+1)*j + 1;
        adj2 = adj + (ng+1)*j + 1;
        for (size_t i=j+1; i< ng; i++)
        {
          if ((*suppressTOMForZeroAdj == 0) || (*adj2 > 0))
          {
            double den1;
            if ((* denomType) == TomDenomMin)
               den1 = fmin(conn[i], conn[j]);
            else
               den1 = (conn[i] + conn[j])/2;
            double den = den1 - fabs(*adj2);
            if (form > 0)
            {
              double r;
              if (den <= 1) r = 0; else r = (*tom2 - *adj2 * 2)/(den-1);
              *tom2 = (*adj2 + r)/2;
            } else {
              if (den==0)
                *tom2 = 0;
              else
                *tom2 = ( *tom2 - *adj2) / den ;
            }
            *(tom + ng*i + j) = *tom2;
            if (fabs(*tom2) > 1)
            {
              Rprintf("TOM greater than 1: actual value: %f, i: %d, j: %d\n", *tom2, i, j);
              nAbove1++;
            }
          } else {
            *tom2 = 0;
            *(tom + ng*i + j) = 0;
            nSuppressed++;
          }
          if (*suppressNegativeTOM && (*tom2 < 0)) { *tom2 = 0; *(tom + ng*i + j) = 0; }
          tom2++;
          adj2++;
        }
      }
      break;
  }

  if (nSuppressed > 0)
    Rprintf("%s.. %lu TOM elements were set to zero because of zero adjacencies.\n", spaces, 
             (unsigned long) nSuppressed); 
  
  // Set the diagonal of tom to 1
  for (size_t i=0; i<ng; i++)
    *(tom + ng*i + i) = 1;

  if (nAbove1 > 0)
    Rprintf("problem: %d TOM entries are larger than 1.\n", nAbove1);

  free(conn);
  if (*verbose > 0) Rprintf("%s..done.\n", spaces);
}



//===========================================================================================
//
// tomSimilarity
//
//===========================================================================================
//

void tomSimilarity(double * expr, double * weights, int * nSamples, int * nGenes, 
                   int * corType, 
                   int * adjType,
                   double * power, 
                   int * tomType, 
                   int * denomType,
                   double * maxPOutliers, 
                   double * quick,
                   int * fallback,
                   int * cosine,
                   int * replaceMissing,
                   int * suppressTOMForZeroAdj,
                   int * suppressNegativeTOM,
                   int * useInternalMatrixAlgebra, 
                   double * tom, 
                   int * warn,
                   int * nThreads,
                   int * verbose, int * indent)
{
  // Rprintf("Starting tomSimilarity...\n");
  double 	* adj;

  int		ng = *nGenes, ns = *nSamples;
  size_t	matSize = ( (size_t)ng) * ( (size_t) ng);

  int 		err = 0;

  char		spaces[2* *indent+1];

  for (int i=0; i<2* *indent; i++) spaces[i] = ' ';
  spaces[2* *indent] = '\0';

/*
  int size=4000;
  int success = 0;
  double * pt;
  while ( (pt=malloc(size*size*sizeof(double)))!=NULL )
  {
    size+=1000;
    success = 1;
    free(pt);
  }

  size -= 1000;
  if ((*verbose > 0) && success) 
     Rprintf("%sRough guide to maximum array size: about %d x %d array of doubles..\n", 
            spaces, size, size);
*/

  if (*verbose > 0) Rprintf("%sTOM calculation: adjacency..\n", spaces);
  if (* tomType==TomTypeNone)  // just calculate adjacency.
  {
    adjacency(expr, weights, ns, ng, *corType, *adjType, *power, *maxPOutliers, *quick, *fallback, *cosine, 
              *replaceMissing, tom, &err, warn, nThreads, *verbose, *indent);
    if (*verbose > 0) Rprintf("\n");
    if (err) error(AdjErrors[err]);
    return;
  }
    
  if ((adj = malloc(matSize * sizeof(double))) == NULL)
    error("Memmory allocation error.");

  if ((* tomType == TomTypeSigned) && (* adjType == AdjTypeUnsigned))
    * adjType = AdjTypeUnsignedKeepSign;
    
  if (*tomType == TomTypeSignedNowick) *adjType = AdjTypeUnsignedKeepSign;

  if ((* tomType == TomTypeUnsigned) && (* adjType == AdjTypeUnsignedKeepSign))
    * adjType = AdjTypeUnsigned;

  adjacency(expr, weights, ns, ng, * corType, * adjType, * power, * maxPOutliers, * quick, *fallback, *cosine, 
            *replaceMissing, adj, & err, warn, nThreads, *verbose, *indent);

  // Rprintf("TOM 1\n");
  if (err) 
  {
     Rprintf("TOM: exit because 'adjacency' reported an error.\n");
     free(adj);
     error(AdjErrors[err]);
  } else {
    tomSimilarityFromAdj(adj, nGenes, tomType, denomType, suppressTOMForZeroAdj, suppressNegativeTOM,
                         useInternalMatrixAlgebra, 
                         tom, verbose, indent);
    free(adj);
  }
}

/*======================================================================================================

  Function returning the column-wise minimum and minimum index. For easier integration with R, the index
will also be stored as a double. NA's are ignored.

========================================================================================================*/

void minWhichMin(double * matrix, int * nRows, int * nColumns, double * min, double * whichMin)
{
  int nrows = *nRows, ncols = *nColumns;

  for (size_t i=0; i<ncols; i++)
  {
    double * col = matrix + i*nrows;
    double curmin = *col;
    double curwhich = 0;
    for (size_t j=1; j<nrows; j++)
    {
      col++;
      if ( ISNAN(curmin) || (!ISNAN(*col) && (*col < curmin))) { curmin = *col; curwhich = (double) j; }
    }
    min[i] = curmin;
    whichMin[i] = curwhich;
  }
}

void minWhichMin_row(double * matrix, int * nRows, int * nColumns, double * min, double * whichMin)
{
  int nrows = *nRows, ncols = *nColumns;

  for (size_t i=0; i<nrows; i++)
  {
    double * val = matrix + i;
    double curmin = *val;
    double curwhich = 0;
    for (size_t j=1; j<ncols; j++)
    {
      val+=nrows;
      if ( ISNAN(curmin) || (!ISNAN(*val) && (*val < curmin))) { curmin = *val; curwhich = (double) j; }
    }
    min[i] = curmin;
    whichMin[i] = curwhich;
  }
}


/*

  Function returning the column-wise mean. NAs are ignored.

*/

void mean(double * matrix, int * nRows, int * nColumns, double * mean)
{
  int nrows = *nRows, ncols = *nColumns;

  for (size_t i=0; i<ncols; i++)
  {
    double * col = matrix + i*nrows;
    double sum = 0;
    size_t count = 0;
    for (size_t j=1; j<nrows; j++)
    {
      col++;
      if (!ISNAN(*col))
      {
        sum += *col;
        count++;
      }
    }
    if (count > 0)
      mean[i] = sum/count;
    else
      mean[i] = NA_REAL;
  }
}


// Version to be called via .Call
//

SEXP tomSimilarityFromAdj_call(SEXP adj_s, 
                        SEXP tomType_s, SEXP denomType_s,
                        SEXP suppressTOMForZeroAdj_s,
                        SEXP suppressNegativeTOM_s,
                        SEXP useInternalMatrixAlgebra_s,
                        SEXP verbose_s, SEXP indent_s)
{
  // Rprintf("Step 1\n");
  SEXP dim, tom_s;

  int *nGenes, *verbose, *indent;
  int *tomType, *denomType, *suppressTOMForZeroAdj, *suppressNegativeTOM, *useInternalMatrixAlgebra;

  double *adj, *tom, *maxPOutliers;
  // Rprintf("Step 2\n");

  /* Get dimensions of 'expr'. */
  PROTECT(dim = getAttrib(adj_s, R_DimSymbol));
  nGenes = INTEGER(dim);
  if (*nGenes!= *(nGenes+1)) 
  {
    UNPROTECT(1);
    error("Input adjacency is not symmetric.");
  }
  // nGenes = INTEGER(dim)+1;
  adj = REAL(adj_s);

  // Rprintf("Step 3\n");
  tomType = INTEGER(tomType_s);
  denomType = INTEGER(denomType_s);
  suppressTOMForZeroAdj = INTEGER(suppressTOMForZeroAdj_s);
  suppressNegativeTOM = INTEGER(suppressNegativeTOM_s);
  useInternalMatrixAlgebra = INTEGER(useInternalMatrixAlgebra_s);
  verbose = INTEGER(verbose_s);
  indent = INTEGER(indent_s);

  // Rprintf("Step 4\n");
  PROTECT(tom_s = allocMatrix(REALSXP, *nGenes, *nGenes));   
  tom = REAL(tom_s);

  // Rprintf("Calling tomSimilarity...\n");
  tomSimilarityFromAdj(adj, nGenes,
                tomType, denomType, suppressTOMForZeroAdj, suppressNegativeTOM,
                useInternalMatrixAlgebra, tom, 
                verbose, indent);

  // Rprintf("Returned from tomSimilarity...\n");
  UNPROTECT(2);

  return tom_s;
}

// Version to be called via .Call
//

SEXP tomSimilarity_call(SEXP expr_s, 
                        SEXP weights_s,
                        SEXP corType_s, SEXP adjType_s, SEXP power_s,
                        SEXP tomType_s, SEXP denomType_s,
                        SEXP maxPOutliers_s, SEXP quick_s,  
                        SEXP fallback_s, SEXP cosine_s, 
                        SEXP replaceMissing_s,
                        SEXP suppressTOMForZeroAdj_s,
                        SEXP suppressNegativeTOM_s,
                        SEXP useInternalMatrixAlgebra_s,
                        SEXP warn_s, // This is an "output" variable
                        SEXP nThreads_s, SEXP verbose_s, SEXP indent_s)
{
  // Rprintf("Step 1\n");
  SEXP dim, tom_s;

  int *nSamples, *nGenes, *fallback, *cosine, *warn, *nThreads, *verbose, *indent;
  int *corType, *adjType, *tomType, *denomType, *replaceMissing, *suppressTOMForZeroAdj, *suppressNegativeTOM,
      *useInternalMatrixAlgebra;

  double *expr, *weights, *power, *quick, *tom, *maxPOutliers;
  // Rprintf("Step 2\n");

  /* Get dimensions of 'expr'. */
  PROTECT(dim = getAttrib(expr_s, R_DimSymbol));
  nSamples = INTEGER(dim);
  nGenes = INTEGER(dim)+1;
  expr = REAL(expr_s);

  weights = isNull(weights_s)? NULL : REAL(weights_s);

  // Rprintf("Step 3\n");
  corType = INTEGER(corType_s);
  adjType = INTEGER(adjType_s);
  tomType = INTEGER(tomType_s);
  denomType = INTEGER(denomType_s);
  fallback = INTEGER(fallback_s);
  cosine = INTEGER(cosine_s);
  replaceMissing = INTEGER(replaceMissing_s);
  suppressTOMForZeroAdj = INTEGER(suppressTOMForZeroAdj_s);
  suppressNegativeTOM = INTEGER(suppressNegativeTOM_s);
  useInternalMatrixAlgebra = INTEGER(useInternalMatrixAlgebra_s);
  warn = INTEGER(warn_s);
  nThreads = INTEGER(nThreads_s);
  verbose = INTEGER(verbose_s);
  indent = INTEGER(indent_s);

  power = REAL(power_s);
  quick = REAL(quick_s);
  maxPOutliers = REAL(maxPOutliers_s);

  // Rprintf("Step 4\n");
  PROTECT(tom_s = allocMatrix(REALSXP, *nGenes, *nGenes));   
  tom = REAL(tom_s);

  // Rprintf("Calling tomSimilarity...\n");
  tomSimilarity(expr, weights, nSamples, nGenes,
                corType, adjType, power,
                tomType, denomType, 
                maxPOutliers, quick, fallback, cosine,
                replaceMissing,
                suppressTOMForZeroAdj, suppressNegativeTOM, 
                useInternalMatrixAlgebra,
                tom, warn, nThreads, verbose, indent);

  // Rprintf("Returned from tomSimilarity...\n");
  UNPROTECT(2);

  return tom_s;
}

void checkAvailableMemoryForR(double * size)
{
  *size = 1.0 *  checkAvailableMemory() ;
}

/* =============================================================================================
 *
 *  Register native routines here.
 *
 * =============================================================================================*/

void R_init_WGCNA(DllInfo * info)
{
  static const R_CallMethodDef callMethods[]  = {
    {"tomSimilarity_call", (DL_FUNC) &tomSimilarity_call, 19},
    {"tomSimilarityFromAdj_call", (DL_FUNC) &tomSimilarityFromAdj_call, 8},
    {"cor1Fast_call", (DL_FUNC) &cor1Fast_call, 9},
    {"bicor1_call", (DL_FUNC) &bicor1_call, 11},
    {"bicor2_call", (DL_FUNC) &bicor2_call, 16},
    {"corFast_call", (DL_FUNC) &corFast_call, 12},
    {"parallelQuantile", (DL_FUNC) &parallelQuantile, 2},
    {"parallelMean", (DL_FUNC) &parallelMean, 2},
    {"parallelMin", (DL_FUNC) &parallelMin, 1},
    {"minWhich_call", (DL_FUNC) &minWhich_call, 2},
    {"quantileC_call", (DL_FUNC) &quantileC_call, 2},
    {"rowQuantileC_call", (DL_FUNC) &rowQuantileC_call, 2},
    {"qorder", (DL_FUNC) &qorder, 1},
    {NULL, NULL, 0}
  };

  static R_NativePrimitiveArgType checkAvailableMemoryForR_t[] = { REALSXP };
  static const R_CMethodDef CMethods[] = {
    {"checkAvailableMemoryForR", (DL_FUNC) &checkAvailableMemoryForR, 1, checkAvailableMemoryForR_t},
    {NULL, NULL, 0}
  };
  R_registerRoutines(info, CMethods, callMethods, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

