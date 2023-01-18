
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

#ifndef __networkFunctions_h__

#define __networkFunctions_h__

#include <stdio.h>
#include <stdlib.h>

#include <sys/time.h>

#include <R.h>
#include <Rinternals.h>

#include <R_ext/BLAS.h>
#include <R_ext/libextern.h>
#include <R_ext/Rdynload.h>

#define LDOUBLE         long double

#include "parallelQuantile_stdC.h"
#include "pivot_declarations.h"
#include "corFunctions-utils.h"
#include "corFunctions.h"
#include "myMatrixMultiplication.h"


size_t checkAvailableMemory(void);

#endif
