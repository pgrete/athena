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
//! \file const.cpp
//  \brief implementation of radiation integrators: constant radiation
//======================================================================================


// Athena++ headers
#include "../radiation.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"

// Class header
#include "rad_integrators.hpp"

RadIntegrator::RadIntegrator(Radiation *prad, ParameterInput *pin)
{
  pmy_rad = prad;
  rad_G0_ = pin->GetReal("problem", "G0");
#ifdef INCLUDE_CHEMISTRY
  MeshBlock* pmy_block = prad->pmy_block;
  pmy_chemnet = pmy_block->pspec->pchemnet;
  n_cols_ang = pmy_rad->nfreq * pmy_chemnet->n_cols_;
  //allocate array for column density
  int ncells1 = pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmy_block->block_size.nx2 > 1) ncells2 = pmy_block->block_size.nx2 + 2*(NGHOST);
  if (pmy_block->block_size.nx3 > 1) ncells3 = pmy_block->block_size.nx3 + 2*(NGHOST);
  col_tot.NewAthenaArray(ncells3, ncells2, ncells1, n_cols_ang);
#endif
}

RadIntegrator::~RadIntegrator() {
#ifdef INCLUDE_CHEMISTRY
  col_tot.DeleteAthenaArray();
#endif 
}

#ifdef INCLUDE_CHEMISTRY
void RadIntegrator::GetColMB() {}
void RadIntegrator::UpdateRadiation() {}
#endif

