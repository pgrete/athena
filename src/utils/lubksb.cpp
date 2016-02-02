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
//! \file lubksb.cpp
//======================================================================================



//--------------------------------------------------------------------------------------
//! \fn  void  Lubksb_nr(int n, AthenaArray<Real> &a, AthenaArray<int> &indx,
//      AthenaArray<Real> &b)
// \brief Backward substitution (from numerical recipies)
// a is the input matrix done with LU decomposition, n is the matrix size
// indx id the history of row permutation
// b is the vector on the right (AX=b), and is returned with the solution
void Lubksb_nr(int n, AthenaArray<Real> &a, AthenaArray<int> &indx,
               AthenaArray<Real> &b)
{
  int ii=-1,ip;
  Real sum;
  // Solve L*y=b
  for(int i=0;i<n;++i) {
    ip=indx(i);
    sum=b(ip);
    b(ip)=b(i);
    if (ii>=0)
      for (int j=ii;j<=i-1;++j) sum -= a(i,j)*b(j);
    else if (sum) ii=i;
    b(i)=sum;
  }
  // Solve U*x=y
  for(int i=n-1;i>=0;--i) {
    sum=b(i);
    for (int j=i+1;j<n;++j){
      Real& ap = a(i,j);
      Real& bp = b(j);
      sum -= ap * bp;
    }
    b(i)=sum/a(i,i);
  }
  
}
