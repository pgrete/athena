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
//! \file matrixmult.cpp
//======================================================================================

//--------------------------------------------------------------------------------------
//! \fn  void MatrixMult(int m, int n, AthenaArray<Real> &a,
//    AthenaArray<Real> &b, AthenaArray<Real> &c)
// \brief perform a(m,n)*b(n)=c(m)

void MatrixMult(int m, int n, AthenaArray<Real> &a,
                AthenaArray<Real> &b, AthenaArray<Real> &c)
{
  
  for(int i=0; i<m; ++i) {
    c(i) = 0.0;
    for(int j=0; j<n; ++j){
      Real& ap = a(i,j);
      Real& bp = b(j);
      c(i) += ap * bp;
    }}
}