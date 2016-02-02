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
#include "../athena.hpp"


//======================================================================================
//! \file Gauleg.cpp
//======================================================================================



//--------------------------------------------------------------------------------------
//! \fn  void Gauleg(int n, Real x1, Real x2,  AthenaArray<Real> &x,
//          AthenaArray<Real> &w);
// \brief gauss-legendre weight routine from numerical recipes

void Gauleg(int n, Real x1, Real x2,  AthenaArray<Real> &x,
            AthenaArray<Real> &w)
{
  
  Real eps = 3.0e-14;
  Real xm, xl, z, z1;
  Real p1, p2, p3, pp;
  int m;
  
  m = (n + 1) / 2;
  xm = 0.5 * (x2 + x1);
  xl = 0.5 * (x2 - x1);
  
  for(int i=1; i<=m; ++i) {
    z = cos(PI * ((Real)i - 0.25) / ((Real)n + 0.5));
    do {
      p1=1.0;
      p2=0.0;
      for(int j=1; j<=n; ++j) {
        p3 = p2;
        p2 = p1;
        p1 = ((2.0 * (Real)j - 1.0) * z * p2 - ((Real)j - 1.0) * p3) / (Real)j;
      }
      pp = (Real)n * (z * p1 - p2) / (z * z - 1.0);
      z1 = z;
      z = z1 - p1 / pp;
    }  while(fabs(z - z1) > eps);
    x(i-1) = xm - xl * z;
    x(n-i) = xm + xl * z;
    w(i-1) = 2.0 * xl / ((1.0 - z * z) * pp * pp);
    w(n-i) = w(i-1);
  }
  
}