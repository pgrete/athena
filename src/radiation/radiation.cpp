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
#include "../parameter_input.hpp"
#include "../mesh.hpp"
#include "integrators/rad_integrators.hpp"

// constructor, initializes data structures and parameters

// The default opacity function.
// Do nothing. Keep the opacity as the initial value
inline void DefaultOpacity(MeshBlock *pmb, AthenaArray<Real> &prim)
{
  
}


Radiation::Radiation(MeshBlock *pmb, ParameterInput *pin)
{
  // read in the parameters
  int nmu = pin->GetInteger("radiation","nmu");
  int angle_flag = pin->GetOrAddInteger("radiation","angle_flag",0);
  prat = pin->GetReal("radiation","Prat");
  crat = pin->GetReal("radiation","Crat");
  reduced_c  = crat * pin->GetOrAddReal("radiation","reduced_factor",1.0);
  nfreq = pin->GetOrAddInteger("radiation","n_frequency",1);
  
  pmy_block = pmb;
  
  // calculate noct based on dimension
  int ndim = 1;
  if(pmb->block_size.nx2 > 1) ndim = 2;
  if(pmb->block_size.nx3 > 1) ndim = 3;
  
  int n1z = pmb->block_size.nx1 + 2*(NGHOST);
  int n2z = 1;
  int n3z = 1;
  
  
  
  int n_ang; // number of angles per octant and number of octant
  // total calculate total number of angles based on dimensions
  if(ndim == 1){
    n_ang = nmu;
    noct = 2;
  }else if(ndim == 2){
    noct = 4;
    n2z = pmb->block_size.nx2 + 2*(NGHOST);
    if(angle_flag == 0){
      n_ang = nmu * (nmu + 1)/2;
    }else if(angle_flag == 10){
      n_ang = nmu;
    }
  }else if(ndim == 3){
    noct = 8;
    n2z = pmb->block_size.nx2 + 2*(NGHOST);
    n3z = pmb->block_size.nx3 + 2*(NGHOST);
    if(angle_flag == 0){
      n_ang = nmu * (nmu + 1)/2;
    }else if(angle_flag == 10){
      n_ang = nmu * nmu/2;
    }
  }// end 3D
  
  nang = n_ang * noct;
  
  n_fre_ang = nang * nfreq;
  
  
  // allocate arrays
  // store frequency and angles as [nfre][ang]
  ir.NewAthenaArray(n3z,n2z,n1z,n_fre_ang);
  ir1.NewAthenaArray(n3z,n2z,n1z,n_fre_ang);
  
  rad_mom.NewAthenaArray(13,n3z,n2z,n1z);
  sigma_s.NewAthenaArray(n3z,n2z,n1z,nfreq);
  sigma_a.NewAthenaArray(n3z,n2z,n1z,nfreq);
  
  grey_sigma_s.NewAthenaArray(n3z,n2z,n1z);
  grey_sigma_a.NewAthenaArray(n3z,n2z,n1z);
  
  mu.NewAthenaArray(3,n3z,n2z,n1z,nang);
  wmu.NewAthenaArray(nang);

  wfreq.NewAthenaArray(nfreq);
  
  AngularGrid(angle_flag, nmu);
  
  //allocate memory to store the flux
  flux[x1face].NewAthenaArray(n3z,n2z,n1z,n_fre_ang);
  if(n2z > 1) flux[x2face].NewAthenaArray(n3z,n2z,n1z,n_fre_ang);
  if(n3z > 1) flux[x3face].NewAthenaArray(n3z,n2z,n1z,n_fre_ang);
  
  // Initialize the frequency weight
  FrequencyGrid();
  
  // set a default opacity function
  UpdateOpacity = DefaultOpacity;
  
  pradintegrator = new RadIntegrator(this, pin);



}

// destructor

Radiation::~Radiation()
{
  ir.DeleteAthenaArray();
  ir1.DeleteAthenaArray();
  rad_mom.DeleteAthenaArray();
  sigma_s.DeleteAthenaArray();
  sigma_a.DeleteAthenaArray();
  grey_sigma_s.DeleteAthenaArray();
  grey_sigma_a.DeleteAthenaArray();
  mu.DeleteAthenaArray();
  wmu.DeleteAthenaArray();
  wfreq.DeleteAthenaArray();
  
  flux[x1face].DeleteAthenaArray();
  if(pmy_block->block_size.nx2 > 1) flux[x2face].DeleteAthenaArray();
  if(pmy_block->block_size.nx3 > 1) flux[x3face].DeleteAthenaArray();
  
  delete pradintegrator;
  
}


//Enrol the function to update opacity

void Radiation::EnrollOpacityFunction(Opacity_t MyOpacityFunction)
{
  UpdateOpacity = MyOpacityFunction;
  
}

// calculate the frequency integrated moments of the radiation field
// including the ghost zones
void Radiation::CalculateMoment()
{
  Real er, frx, fry, frz, prxx, pryy, przz, prxy, prxz, pryz;
  int n1z = pmy_block->block_size.nx1 + 2*(NGHOST);
  int n2z = pmy_block->block_size.nx2;
  int n3z = pmy_block->block_size.nx3;
  if(n2z > 1) n2z += (2*(NGHOST));
  if(n3z > 1) n3z += (2*(NGHOST));
  
  Real *weight = &(wmu(0));
  
  AthenaArray<Real> i_mom;
  
  i_mom.InitWithShallowCopy(rad_mom);

  
  // reset the moment arrays to be zero
  // There are 13 3D arrays
  for(int n=0; n<13; ++n)
    for(int k=0; k<n3z; ++k)
      for(int j=0; j<n2z; ++j)
#pragma simd
        for(int i=0; i<n1z; ++i){
          i_mom(n,k,j,i) = 0.0;
        }
  
  
  for(int k=0; k<n3z; ++k){
    for(int j=0; j<n2z; ++j){
      for(int i=0; i<n1z; ++i){
        for(int ifr=0; ifr<nfreq; ++ifr){
          er=0.0; frx=0.0; fry=0.0; frz=0.0;
          prxx=0.0; pryy=0.0; przz=0.0; prxy=0.0;
          prxz=0.0; pryz=0.0;
          Real *intensity = &(ir(k,j,i,ifr*nang));
          Real *cosx = &(mu(0,k,j,i,0));
          Real *cosy = &(mu(1,k,j,i,0));
          Real *cosz = &(mu(2,k,j,i,0));
#pragma simd
          for(int n=0; n<nang; ++n){
            er   += weight[n] * intensity[n];
            frx  += weight[n] * intensity[n] * cosx[n];
            fry  += weight[n] * intensity[n] * cosy[n];
            frz  += weight[n] * intensity[n] * cosz[n];
            prxx += weight[n] * intensity[n] * cosx[n] * cosx[n];
            pryy += weight[n] * intensity[n] * cosy[n] * cosy[n];
            przz += weight[n] * intensity[n] * cosz[n] * cosz[n];
            prxy += weight[n] * intensity[n] * cosx[n] * cosy[n];
            prxz += weight[n] * intensity[n] * cosx[n] * cosz[n];
            pryz += weight[n] * intensity[n] * cosy[n] * cosz[n];
          }
          //multiply the frequency weight
          er *= wfreq(ifr);
          frx *= wfreq(ifr);
          fry *= wfreq(ifr);
          frz *= wfreq(ifr);
          prxx *= wfreq(ifr);
          pryy *= wfreq(ifr);
          przz *= wfreq(ifr);
          prxy *= wfreq(ifr);
          prxz *= wfreq(ifr);
          pryz *= wfreq(ifr);
          

          
          //assign the moments
          i_mom(IER,k,j,i) += er;
          i_mom(IFR1,k,j,i) += frx;
          i_mom(IFR2,k,j,i) += fry;
          i_mom(IFR3,k,j,i) += frz;
          i_mom(IPR11,k,j,i) += prxx;
          i_mom(IPR12,k,j,i) += prxy;
          i_mom(IPR13,k,j,i) += prxz;
          i_mom(IPR21,k,j,i) += prxy;
          i_mom(IPR22,k,j,i) += pryy;
          i_mom(IPR23,k,j,i) += pryz;
          i_mom(IPR31,k,j,i) += prxz;
          i_mom(IPR32,k,j,i) += pryz;
          i_mom(IPR33,k,j,i) += przz;
          
        }// End frequency loop
        // Now calculate frequency inetgrated opacity
        Real sum_sigma_s=0.0, sum_sigma_a = 0.0;
        Real *sigmas=&(sigma_s(k,j,i,0));
        Real *sigmaa=&(sigma_a(k,j,i,0));
#pragma simd
        for(int ifr=0; ifr<nfreq; ++ifr){
          sum_sigma_s += sigmas[ifr] * wfreq(ifr);
          sum_sigma_a += sigmaa[ifr] * wfreq(ifr);
        }
        grey_sigma_s(k,j,i) = sum_sigma_s;
        grey_sigma_a(k,j,i) = sum_sigma_a;
      }
    }
    
  }
  
  
}


