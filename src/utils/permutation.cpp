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
//! \file permutation.cpp
//======================================================================================

//--------------------------------------------------------------------------------------
//! \fn int Permutation(int i, int j, int k, int np, AthenaArray<int> &pl)
// \brief permutate the array element

int Permutation(int i, int j, int k, int np, AthenaArray<int> &pl)
{
  int ip=-1;
  
  for(int l=0; l<np; l++) {
    // check each permutation in the table
    for(int m=0; m<3; m++)
      if(i == pl(l,m))
        for(int n=0; n<3; n++)
          if(n != m)
            if(j == pl(l,n))
              for(int o=0;o<3;o++)
                if((o != m) && (o != n))
                  if(k == pl(l,o))
                    ip = l;
  }
  
  return ip;
}
