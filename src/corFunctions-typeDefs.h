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

#ifndef __corFunctions_internal_h__

#define __corFunctions_internal_h__

#define LDOUBLE 	long double

typedef struct 
{
   volatile size_t i, n;
}  progressCounter;

/* For each parallel operation will presumably need a separate structure to hold its
 * information, but can define a common structure holding the general information that is needed to
 * calculate correlation. Can keep two versions, one for calculating cor(x), one for cor(x,y).
 * Each specific thread-task specific struct can contain a pointer to the general structure.
 */

// General information for a [bi]cor(x) calculation

typedef struct
{
   double * x, * weights;
   size_t nr, nc; 
   double * multMat, * result;
   double * aux;
   size_t *nNAentries;
   int *NAme;
   int zeroMAD;
   int * warn;
   double maxPOutliers;
   double quick;
   int robust, fallback;
   int cosine;
   int id;
   int threaded; 	// This flag will be used to indicate whether the calculation really is threaded. 
			// For small problems it doesn't make sense to use threading.
}  cor1ThreadData;

// General information for a [bi]cor(x,y) calculation

typedef struct
{
   cor1ThreadData * x, * y;
}  cor2ThreadData;

// Information for column preparation

typedef struct
{
   cor1ThreadData * x;
   progressCounter * pc;
   pthread_mutex_t * lock;
}  colPrepThreadData;

// Information for symmetrization

typedef struct
{
   cor1ThreadData * x;
   progressCounter * pc;
}  symmThreadData;

// Information for threaded slow calculations for cor1

typedef struct
{
   cor1ThreadData * x;
   progressCounter * pci, * pcj;
   size_t * nSlow, * nNA;
   pthread_mutex_t * lock;
}  slowCalcThreadData;

/*==============================================================================================
 *
 * Threaded 2-variable versions of the correlation functions
 *
 *==============================================================================================
*/


typedef struct
{
   cor2ThreadData * x;
   progressCounter * pci, *pcj;
   size_t * nSlow, * nNA;
   pthread_mutex_t * lock;
   double quick;
}  slowCalc2ThreadData;

// Data for NAing out appropriate rows and columns

typedef struct
{
   cor2ThreadData * x; 
   progressCounter * pci, *pcj;
}  NA2ThreadData;

#endif
