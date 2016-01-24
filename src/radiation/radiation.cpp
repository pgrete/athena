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
#include "../athena.hpp"
#include "../athena_arrays.hpp" 
#include "radiation.hpp"

// constructor, initializes data structures and parameters

Radiation::Radiation(MeshBlock *pmb, ParameterInput *pin)
{
  // read in the parameters
  int nmu = pin->GetInteger("radiation","nmu");
  int angle_flag = pin->GetOrAddInteger("radiation","angle_flag",0);
  prat = pin->GetReal("radiation","Prat");
  crat = pin->GetReal("radiation","Crat");
  reduced_c  = crat * pin->GetOrAddReal("radiation","reduced_factor",1.0);
  int nfreq = pin->GetOrAddInteger("radiation","n_frequency",1);
  
  // calculate noct based on dimension
  int ndim = 1;
  if(pmb->block_size.nx2 > 1) ndim = 2;
  if(pmb->block_size.nx3 > 1) ndim = 3;
  
  int n1z = pmy_block->block_size.nx1 + 2*(NGHOST);
  int n2z = 1;
  int n3z = 1;
  
  
  
  int n_ang, noct;
  // total calculate total number of angles based on dimensions
  if(ndim == 1){
    n_ang = nmu;
    noct = 2;
  }else if(ndim == 2){
    noct = 4;
    n2z = pmy_block->block_size.nx2 + 2*(NGHOST);
    if(angle_flag == 0){
      n_ang = nmu * (nmu + 1)/2;
    }else if(angle_flag == 10){
      n_ang = nmu;
    }
  }else if(ndim == 3){
    noct = 8;
    n2z = pmy_block->block_size.nx2 + 2*(NGHOST);
    n3z = pmy_block->block_size.nx3 + 2*(NGHOST);
    if(angle_flag == 0){
      n_ang = nmu * (nmu + 1)/2;
    }else if(angle_flag == 10){
      n_ang = nmu * nmu/2;
    }
  }// end 3D
  
  nang = n_ang * noct;
  
  // allocate arrays
  ir.NewAthenaArray(nfreq,n3z,n2z,n1z,nang);
  ir1.NewAthenaArray(nfreq,n3z,n2z,n1z,nang);
  
  sigma_s.NewAthenaArray(nfreq,n3z,n2z,n1z);
  sigma_a.NewAthenaArray(nfreq,n3z,n2z,n1z);
  
  mu.NewAthenaArray(3,n3z,n2z,n1z,nang);
  wmu.NewAthenaArray(nang);

  
  

}

// destructor

Radiation::~Radiation()
{
  ir.DeleteAthenaArray();
  ir1.DeleteAthenaArray();
  sigma_s.DeleteAthenaArray();
  sigma_a.DeleteAthenaArray();
  mu.DeleteAthenaArray();
  wmu.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
// \!fn 
// \brief
