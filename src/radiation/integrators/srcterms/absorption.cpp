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
//! \file absorption.cpp
//  \brief  Add absorption source terms
//======================================================================================

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../radiation.hpp"
#include "../../../hydro/hydro.hpp"
#include "../../../hydro/eos/eos.hpp"
#include "../../../mesh.hpp"
#include "../../../coordinates/coordinates.hpp"
#include "../../../utils/utils.hpp"

// this class header
#include "../rad_integrators.hpp"

//--------------------------------------------------------------------------------------
//! \fn RadIntegrator::Absorption()
//  \brief 

// wmu_cm is the weight in the co-moving frame
// wmu_cm=wmu * 1/(1-vdotn/Crat)^2 / Lorz^2
// tran_coef is (1-vdotn/Crat)*Lorz
// rho is gas density
// tgas is gas temperature
// This function only update the absorption term in each cell

void RadIntegrator::Absorption(const AthenaArray<Real> &wmu_cm,
          const AthenaArray<Real> &tran_coef, Real *sigma_a,
          Real *sigma_ae, Real dt, Real rho, Real &tgas,
          AthenaArray<Real> &ir_cm)
{

  Real& prat = pmy_rad->prat;
  Real ct = dt * pmy_rad->reduced_c;
  int& nang=pmy_rad->nang;
  int& nfreq=pmy_rad->nfreq;
  Real gamma = pmy_rad->pmy_block->phydro->peos->GetGamma();
  
  // Temporary array
  AthenaArray<Real> vncsigma, vncsigma2;
  
  vncsigma.InitWithShallowCopy(vncsigma_);
  vncsigma2.InitWithShallowCopy(vncsigma2_);
  
  int badcell=0;
  
  
  Real coef[5];
  for (int i=0; i<5; ++i)
    coef[i] = 0.0;
  
  Real tgasnew = tgas;
  
  for(int ifr=0; ifr<nfreq; ++ifr){
    
    Real suma1=0.0, suma2=0.0;
    Real jr_cm=0.0;
    
    Real dtcsigmat = ct * sigma_a[ifr];
    Real dtcsigmae = ct * sigma_ae[ifr];
    
#pragma simd
    for(int n=0; n<nang; n++){
       vncsigma(n) = 1.0/(1.0 + dtcsigmae * tran_coef(n));
       vncsigma2(n) = tran_coef(n) * vncsigma(n);
       Real ir_weight = ir_cm(n+nang*ifr) * wmu_cm(n);
       jr_cm += ir_weight;
       suma1 += (dtcsigmat * wmu_cm(n) * vncsigma2(n));
       suma2 += (ir_weight * vncsigma(n));
    }
    // Now solve the equation
    // rho dT/gamma-1=-Prat c(sigma T^4 - sigma(a1 T^4 + a2))
    // make sure jr_cm is positive
    jr_cm = std::max(jr_cm, TINY_NUMBER);
    Real tr = sqrt(sqrt(jr_cm));
    // In general, it should be sigmat(T^4 - sigmae E/sigmat)
    if(dtcsigmat > TINY_NUMBER)
      tr *= sqrt(sqrt(dtcsigmae/dtcsigmat));
    
    // No need to do this if already in thermal equilibrium
    if(fabs(tr - tgas) > 1.e-12){
      coef[4] = prat * (dtcsigmat - dtcsigmae * suma1);
      coef[1] = rho /(gamma - 1.0);
      coef[0] = -rho * tgas/(gamma-1.0) - dtcsigmae * prat * suma2;
      if(coef[4] > 0.0 && coef[1] > 0.0 && coef[0] > 0.0){
        badcell = 1;
        tgasnew = tgas;
      }else if(coef[4]/coef[1] > 1.e6){
        tgasnew = pow(-coef[0]/coef[4],0.25);
      }else if(fabs(coef[4]) < 1.e-13){
        tgasnew = -coef[0]/coef[1];
      } else {
      
        int flag = ExactPolynomial(coef[4], coef[1], coef[0], tgasnew);
        if(flag == -1){
         // take the larger value to avoid negative root
         tgasnew = std::max(tgas, tr);

         Laguer(coef, 4, tgasnew);
        }
      }
      
        
    }
    
    // even if tr=told, there can be change for intensity, making them isotropic
    
    if(!badcell){
    
      Real emission = tgasnew * tgasnew * tgasnew * tgasnew;
    
    // Update the co-moving frame specific intensity
#pragma simd
      for(int n=0; n<nang; n++){
        ir_cm(n+nang*ifr) = ir_cm(n+nang*ifr) * vncsigma(n)
                        + dtcsigmat * emission * vncsigma2(n);
      }
    }
    
  }// End Frequency
  
  // Update gas temperature
  tgas = tgasnew;


  return;
}
