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
          Real dt, Real rho, Real &tgas, AthenaArray<Real> &ir_cm)
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
  
  // Polynomial coefficients for Compton
  Real coef[5];
  for (int i=0; i<5; ++i)
  coef[i] = 0.0;

  
  Real tgasnew = tgas;
  
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
    Real jr_cm = std::max(suma2/suma1,TINY_NUMBER);
    

    
    // Add Simple Compton scattering using the partically updated jr and ir

    if((compton_flag_ > 0) && (dtcsigma > TINY_NUMBER)){
    
     // Calculate the sum \int gamma (1-vdotn/c) dw_0 4 dt csigma_s /T_e
      suma1 = 0.0; suma2 = 0.0;
#pragma simd
      for(int n=0; n<nang; n++){
         suma1 += tran_coef(n) * wmu_cm(n) * 4.0 * dtcsigma * telectron;
      }
      suma2 = 4.0 * prat * dtcsigma*(gamma-1.0)*telectron/rho;
      
      Real tr = sqrt(sqrt(jr_cm));
      
      if(fabs(tr - tgas) > 1.e-12){
        coef[4] = (1.0 + suma2* jr_cm)/(suma1 * jr_cm);
        coef[1] = 1.0;
        coef[0] = -(1.0+suma2*jr_cm)/suma1-tgas;
        
        Real trnew;
        
        if(coef[4] > 1.e8){
           trnew = pow(-coef[0]/coef[4],0.25);
        }else{
        
          int flag = ExactPolynomial(coef[4], coef[1], coef[0], trnew);
          if(flag == -1){
          
            trnew = std::max(tgas,tr);
            Laguer(coef, 4, trnew);
          }// end flag
        }
        Real jrnew = trnew * trnew * trnew * trnew;
    
        tgasnew = jrnew/(suma1*jr_cm) - 1.0/suma1 + trnew;
    
        jr_cm = jrnew;
      
      }
      
    }// End Compton
    
    // Update the co-moving frame specific intensity
#pragma simd
    for(int n=0; n<nang; n++){
      ir_cm(n+nang*ifr) +=  (jr_cm - ir_cm(n+nang*ifr)) * vncsigma2(n);
    }
    
    
  }// End Frequency
  
  // Update gas temperature
  tgas = tgasnew;


  return;
}



