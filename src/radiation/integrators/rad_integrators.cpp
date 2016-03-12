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
//! \file rad_integrators.cpp
//  \brief implementation of radiation integrators
//======================================================================================


// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh.hpp"
#include "../radiation.hpp"
#include "rad_integrators.hpp"



RadIntegrator::RadIntegrator(Radiation *prad, ParameterInput *pin)
{

  pmy_rad = prad;
  
      // factor to separate the diffusion and advection part
  taufact_ = pin->GetOrAddInteger("radiation","taucell",5);
  compton_flag_=pin->GetOrAddInteger("radiation","Compton",0);


  int nthreads = prad->pmy_block->pmy_mesh->GetNumMeshThreads();
  int ncells1 = prad->pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1;
  int ncells3 = 1;

  flx_.NewAthenaArray(nthreads,ncells1,prad->n_fre_ang);
  vel_.NewAthenaArray(nthreads,ncells1,prad->n_fre_ang);
  flx2_.NewAthenaArray(nthreads,ncells1,prad->n_fre_ang);
  vel2_.NewAthenaArray(nthreads,ncells1,prad->n_fre_ang);
  
  x1face_area_.NewAthenaArray(nthreads,ncells1+1);
  if(prad->pmy_block->block_size.nx2 > 1) {
    x2face_area_.NewAthenaArray(nthreads,ncells1);
    x2face_area_p1_.NewAthenaArray(nthreads,ncells1);
    ncells2 = prad->pmy_block->block_size.nx2 + 2*(NGHOST);
  }
  if(prad->pmy_block->block_size.nx3 > 1) {
    x3face_area_.NewAthenaArray(nthreads,ncells1);
    x3face_area_p1_.NewAthenaArray(nthreads,ncells1);
    ncells3 = prad->pmy_block->block_size.nx3 + 2*(NGHOST);
  }
  cell_volume_.NewAthenaArray(nthreads,ncells1);
  
  temp_i1_.NewAthenaArray(nthreads,ncells3,ncells2,ncells1,prad->n_fre_ang);
  temp_i2_.NewAthenaArray(nthreads,ncells3,ncells2,ncells1,prad->n_fre_ang);
  
  vncsigma_.NewAthenaArray(prad->nang);
  vncsigma2_.NewAthenaArray(prad->nang);
  wmu_cm_.NewAthenaArray(prad->nang);
  tran_coef_.NewAthenaArray(prad->nang);
  cm_to_lab_.NewAthenaArray(prad->nang);
  ir_cm_.NewAthenaArray(prad->n_fre_ang);

}

// destructor

RadIntegrator::~RadIntegrator()
{
  flx_.DeleteAthenaArray();
  vel_.DeleteAthenaArray();
  flx2_.DeleteAthenaArray();
  vel2_.DeleteAthenaArray();
  
  x1face_area_.DeleteAthenaArray();
  if(pmy_rad->pmy_block->block_size.nx2 > 1) {
    x2face_area_.DeleteAthenaArray();
    x2face_area_p1_.DeleteAthenaArray();
  }
  if(pmy_rad->pmy_block->block_size.nx3 > 1) {
    x3face_area_.DeleteAthenaArray();
    x3face_area_p1_.DeleteAthenaArray();
  }
  cell_volume_.DeleteAthenaArray();
  temp_i1_.DeleteAthenaArray();
  temp_i2_.DeleteAthenaArray();
  
  vncsigma_.DeleteAthenaArray();
  vncsigma2_.DeleteAthenaArray();
  wmu_cm_.DeleteAthenaArray();
  tran_coef_.DeleteAthenaArray();
  cm_to_lab_.DeleteAthenaArray();
  ir_cm_.DeleteAthenaArray();
  
}

