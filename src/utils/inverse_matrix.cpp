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
#include "./utils.hpp"


//======================================================================================
//! \file inverse_matrix.cpp
//======================================================================================


//--------------------------------------------------------------------------------------
//! \fn  void  InverseMatrix(int n, AthenaArray<Real> &a, AthenaArray<Real> &b)
// \brief Inverse matrix solver
// a: input matrix; n: matrix size, b: return matrix
// Note: the input matrix will be DESTROYED

void InverseMatrix(int n, AthenaArray<Real> &a, AthenaArray<Real> &b)
{
  AthenaArray<int> indx;
  AthenaArray<Real> col;
  Real d;
  
  indx.NewAthenaArray(n);
  col.NewAthenaArray(n);
  
  
  
  Ludcmp_nr(n,a,indx,&d);
  
  for (int j=0; j<n; ++j) {
    for (int i=0; i<n; ++i) col(i)=0.0;
    col(j)=1.0;
    Lubksb_nr(n, a, indx, col);
    for (int i=0; i<n; i++)    b(i,j) = col(i);
  }
  
  
  indx.DeleteAthenaArray();
  col.DeleteAthenaArray();
  
  return;
  
}



