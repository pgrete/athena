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
//! \file radiation.cpp
//  \brief implementation of functions in class Radiation
//======================================================================================

// Athena++ headers
#include "radiation.hpp"
#include "../mesh.hpp"

Radiation::Radiation(MeshBlock *pmb, ParameterInput *pin)
{
  // read in the parameters
  int six_ray_flag = pin->GetOrAddInteger("radiation","six_ray_flag",0);
	nfreq = pin->GetOrAddInteger("radiation","n_frequency",1);
  
  pmy_block = pmb;
  
	if (six_ray_flag) {
		nang = 6;
	} else {
		nang = 1;
	}
  
  n_fre_ang = nang * nfreq;
  
  
  // allocate arrays
  int ncells1 = pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmy_block->block_size.nx2 > 1) ncells2 = pmy_block->block_size.nx2 + 2*(NGHOST);
  if (pmy_block->block_size.nx3 > 1) ncells3 = pmy_block->block_size.nx3 + 2*(NGHOST);
  // store frequency and angles as [nfre][ang]
  ir.NewAthenaArray(ncells3, ncells2, ncells1, n_fre_ang);
}

// destructor

Radiation::~Radiation()
{
  ir.DeleteAthenaArray();
}
