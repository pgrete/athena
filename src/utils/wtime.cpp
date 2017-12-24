//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

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

