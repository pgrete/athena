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

void RadIntegrator::Scattering(const AthenaArray<Real> &wmu_cm,
          const AthenaArray<Real> &tran_coef, Real *sigma_s,
          Real dt, Real rho, Real *tgas, AthenaArray<Real> &ir_cm)
{

  Real& prat = pmy_rad->prat;
  Real ct = dt * pmy_rad->reduced_c;
  int& nang=pmy_rad->nang;
  int& nfreq=pmy_rad->nfreq;
  Real gamma = pmy_rad->pmy_block->phydro->peos->GetGamma();
  Real telectron = 1.0/pmy_rad->telectron;
  
  // Temporary array
  AthenaArray<Real> vncsigma2;
  
  vncsigma2.InitWithShallowCopy(vncsigma2_);
  
  
  
  Real tgasnew = (*tgas);
  
  for(int ifr=0; ifr<nfreq; ++ifr){
    
    Real suma1=0.0, suma2=0.0;
    
    Real dtcsigma = ct * sigma_s[ifr];
    
#pragma simd
    for(int n=0; n<nang; n++){
       Real vncsigma = 1.0/(1.0 + dtcsigma * tran_coef(n));
       vncsigma2(n) = tran_coef(n) * dtcsigma * vncsigma;
       suma1 += (wmu_cm(n) * vncsigma);
       suma2 += (ir_cm(n+nang*ifr) * wmu_cm(n) * vncsigma);
    }
    // Now solve the equation
    // rho dT/gamma-1=-Prat c(sigma T^4 - sigma(a1 T^4 + a2))
    // make sure jr_cm is positive
    Real jr_cm = suma2/suma1;
    
    // Add Simple Compton scattering
    // for compton scattering, call a function to calculate the mean energy change
    // in operator spliting way, using the partially update jrcom
    Real compte = 0.0;
    if(compton_flag_ > 0){
    
     // Calculate the sum \int gamma (1-vdotn/c) dw_0 4 dt csigma_s /T_e
      suma1 = 0.0; suma2 = 0.0;
#pragma simd
      for(int n=0; n<nang; n++){
         suma1 += tran_coef(n) * wmu_cm(n) * 4.0 * dtcsigma * telectron;
      }
      suma2 = 4.0 * prat * dtcsigma*(gamma-1.0)*telectron/rho;
      compte = Compton(suma1,suma2,(*tgas),jr_cm,&tgasnew);
      jr_cm += compte * suma1;
      // update co-moving specific intensity
#pragma simd
      for(int n=0; n<nang; n++){
        ir_cm(n+nang*ifr) += (tran_coef(n) * 4.0 * dtcsigma * compte * telectron);
      }
      
    }// End Compton
    
    
    // Update the co-moving frame specific intensity
#pragma simd
    for(int n=0; n<nang; n++){
      ir_cm(n+nang*ifr) +=  (jr_cm - ir_cm(n+nang*ifr)) * vncsigma2(n);
    }
    
    
  }// End Frequency
  
  // Update gas temperature
  (*tgas) = tgasnew;


  return;
}


// Add simple isotropic energy exchange due to Compton process
//
Real RadIntegrator::Compton(const Real suma1, const Real suma2,
  Real tgas, Real ercom, Real *tgas_new)
{

  ercom = std::max(ercom, TINY_NUMBER);
  Real tr = sqrt(sqrt(ercom));
  
  Real compte = 0.0;
  
  if(fabs(tr - tgas) > 1.e-12){
    Real coef8 = suma2/suma1;
    Real coef5 = 1.0;
    Real coef4 = 1.0/suma1 - suma2 * ercom/suma1 - tgas;
    Real tconst = -ercom/suma1;
    
    Real tmin = std::min(tr, tgas);
    Real tmax = std::max(tr, tgas);
    
    tmin *= 0.5;
    tmax *= 1.5;
        
    Real trnew = Rtsafe(Tcompton, tmin, tmax, 1.e-10, coef8, coef5, coef4, tconst);
      
    Real jrnew = trnew * trnew * trnew * trnew;
    
    Real tgasnew = 1.0/suma1 - ercom/(suma1 * jrnew) + trnew;

    *tgas_new = tgasnew;
    
    compte = (tgasnew - trnew) * jrnew;
    
  }// end if tr and tgas are different
  else{
    *tgas_new = tgas;
  }
 
  return compte;
}


