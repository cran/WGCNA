
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

#ifndef __corFunctions_h__

#define __corFunctions_h__

#define minSizeForThreading	100

#include <R.h>
#include <Rinternals.h>



void cor1Fast(double * x, int * nrow, int * ncol, double * quick, 
          int * cosine,
          double * result, int *nNA, int * err, 
          int * nThreads,
          int * verbose, int * indent);

void bicor1Fast(double * x, int * nrow, int * ncol, 
            double * maxPOutliers, double * quick, 
            int * fallback, int * cosine,
            double * result, int *nNA, int * err, 
            int * warn,
            int * nThreads,
            int * verbose, int * indent);

void bicorFast(double * x, int * nrow, int * ncolx, double * y, int * ncoly,
           int * robustX, int * robustY, 
           double * maxPOutliers, double * quick, 
           int * fallback, 
           int * cosineX, int * cosineY, 
           double * result, int *nNA, int * err,
           int * warnX, int * warnY,
           int * nThreads,
           int * verbose, int * indent);

void corFast(double * x, int * nrow, int * ncolx, double * y, int * ncoly,
           double * quick,
           int * cosineX, int * cosineY,
           double * result, int *nNA, int * err,
           int * nThreads,
           int * verbose, int * indent);

SEXP cor1Fast_call(SEXP x_s, SEXP quick_s, SEXP cosine_s,
                   SEXP nNA_s, SEXP err_s,
                   SEXP nThreads_s, SEXP verbose_s, SEXP indent_s);

SEXP corFast_call(SEXP x_s, SEXP y_s,
                 SEXP quick_s,
                 SEXP cosineX_s, SEXP cosineY_s,
                 SEXP nNA_s, SEXP err_s,
                 SEXP nThreads_s, SEXP verbose_s, SEXP indent_s);

SEXP bicor1_call(SEXP x_s,
                 SEXP maxPOutliers_s, SEXP quick_s,
                 SEXP fallback_s, SEXP cosine_s,
                 SEXP nNA_s, SEXP err_s, SEXP warn_s,
                 SEXP nThreads_s, SEXP verbose_s, SEXP indent_s);

SEXP bicor2_call(SEXP x_s, SEXP y_s,
                 SEXP robustX_s, SEXP robustY_s,
                 SEXP maxPOutliers_s, SEXP quick_s,
                 SEXP fallback_s,
                 SEXP cosineX_s, SEXP cosineY_s,
                 SEXP nNA_s, SEXP err_s,
                 SEXP warnX_s, SEXP warnY_s,
                 SEXP nThreads_s, SEXP verbose_s, SEXP indent_s);


#endif
