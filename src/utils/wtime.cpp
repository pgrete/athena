//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================

// Athena headers
#include <time.h>
#ifndef DOS
#include <sys/time.h>
#endif
//======================================================================================
//! \file wtime.cpp
//======================================================================================



void wtime(double *t)
{
   /* a generic timer */
   static int sec = -1;
   struct timeval tv;
   gettimeofday(&tv, (struct timezone *)0);
   if (sec < 0) sec = tv.tv_sec;
   *t = (tv.tv_sec - sec) + 1.0e-6*tv.tv_usec;
}

