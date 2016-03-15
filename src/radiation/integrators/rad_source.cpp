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
//! \file rad_source.cpp
//  \brief Add radiation source terms to both radiation and gas
//======================================================================================


// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh.hpp"
#include "../radiation.hpp"
#include "../../coordinates/coordinates.hpp" //


// class header
#include "rad_integrators.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif


void RadIntegrator::AddSourceTerms(MeshBlock *pmb, AthenaArray<Real> &u,
        AthenaArray<Real> &w, AthenaArray<Real> &ir, const int step)
{
  Radiation *prad=pmb->prad;
  Coordinates *pco = pmb->pcoord;
  
  Real& prat = prad->prat;
  Real invcrat = 1.0/prad->crat;
  
  AthenaArray<Real> wmu_cm, tran_coef, ir_cm, cm_to_lab;
  
  Real *sigma_at, *sigma_aer, *sigma_s;
  Real *lab_ir;
  
  int& nang =prad->nang;
  int& nfreq=prad->nfreq;
  
  
  // Get the temporary array
  wmu_cm.InitWithShallowCopy(wmu_cm_);
  tran_coef.InitWithShallowCopy(tran_coef_);
  ir_cm.InitWithShallowCopy(ir_cm_);
  cm_to_lab.InitWithShallowCopy(cm_to_lab_);
  

  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  Real dt;
  if (step == 1) {
    dt = 0.5*(pmb->pmy_mesh->dt);
  } else {
    dt = (pmb->pmy_mesh->dt);
  }
  
  
  for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){
      for(int i=is; i<=ie; ++i){
         // First, get the estimated velocity
         Real vx = w(IVX,k,j,i);
         Real vy = w(IVY,k,j,i);
         Real vz = w(IVZ,k,j,i);
        
         Real rho = w(IDN,k,j,i);
         Real tgas= w(IEN,k,j,i)/rho;
        
         // Do not use the velocity directly
         // use the predicted velocity based on moment equation
         // This is to avoid angular truncation error in the implicit method
         PredictVel(k,j,i, 0.5 * dt, rho, &vx, &vy, &vz);
        
         Real vel = vx * vx + vy * vy + vz * vz;
        
         Real ratio = sqrt(vel) * invcrat;
         // Limit the velocity to be smaller than the speed of light
         if(ratio > prad->vmax){
           Real factor = prad->vmax/ratio;
           vx *= factor;
           vy *= factor;
           vz *= factor;
           
           vel *= (factor*factor);
         }
        
        
         Real lorzsq = 1.0/(1.0 - vel  * invcrat * invcrat);
        

        
         sigma_at=&(prad->sigma_a(k,j,i,0));
         sigma_aer=&(prad->sigma_ae(k,j,i,0));
         sigma_s=&(prad->sigma_s(k,j,i,0));

         // Prepare the transformation coefficients
         Real numsum = 0.0;

#pragma simd
         for(int n=0; n<nang; ++n){
            Real vdotn = vx * prad->mu(0,k,j,i,n) + vy * prad->mu(1,k,j,i,n)
                        + vz * prad->mu(2,k,j,i,n);
            Real vnc = 1.0 - vdotn * invcrat;
            tran_coef(n) = sqrt(lorzsq) * vnc;
            wmu_cm(n) = prad->wmu(n)/(tran_coef(n) * tran_coef(n));
            numsum += wmu_cm(n);
            cm_to_lab(n) = tran_coef(n)*tran_coef(n)*tran_coef(n)*tran_coef(n);
           
         }
           // Normalize weight in co-moving frame to make sure the sum is one
         numsum = 1.0/numsum;
#pragma simd
         for(int n=0; n<nang; ++n){
            wmu_cm(n) *= numsum;
         }
        
         lab_ir=&(ir(k,j,i,0));
        
         for(int ifr=0; ifr<nfreq; ++ifr){
#pragma simd
           for(int n=0; n<nang; ++n){
             ir_cm(n+ifr*nang) = lab_ir[n+ifr*nang] * cm_to_lab(n);
           }
         }// End frequency
        
         // Add absorption opacity source
         Absorption(wmu_cm,tran_coef, sigma_at,sigma_aer, dt, rho, &tgas, ir_cm);
        
         // Add scattering opacity source
         Scattering(wmu_cm,tran_coef, sigma_s, dt, rho, &tgas, ir_cm);
        
       
        //After all the source terms are applied
        // Calculate energy and momentum source terms

        
        // first, calculate Er and Fr in lab frame before the step
        Real er0 = 0.0;
        Real frx0 = 0.0;
        Real fry0 = 0.0;
        Real frz0 = 0.0;
        
        for(int ifr=0; ifr<nfreq; ++ifr){
#pragma simd
           for(int n=0; n<nang; ++n){
               Real ir_weight = lab_ir[n+ifr*prad->nang] * prad->wmu(n);
               er0 += ir_weight;
               frx0 += ir_weight * prad->mu(0,k,j,i,n);
               fry0 += ir_weight * prad->mu(1,k,j,i,n);
               frz0 += ir_weight * prad->mu(2,k,j,i,n);
           }
           er0 *= prad->wfreq(ifr);
           frx0 *= prad->wfreq(ifr);
           fry0 *= prad->wfreq(ifr);
           frz0 *= prad->wfreq(ifr);
        }
        
        
       // now update the lab frame intensity
        for(int ifr=0; ifr<nfreq; ++ifr){
#pragma simd
          for(int n=0; n<nang; ++n){
              lab_ir[n+ifr*nang] = std::max(ir_cm(n+ifr*nang)/cm_to_lab(n),
                                   TINY_NUMBER);
          }
        }

       // now calculate the new moments
               // first, calculate Er and Fr in lab frame before the step
        Real er = 0.0;
        Real frx = 0.0;
        Real fry = 0.0;
        Real frz = 0.0;
        
        for(int ifr=0; ifr<nfreq; ++ifr){
#pragma simd
           for(int n=0; n<nang; ++n){
               Real ir_weight = lab_ir[n+ifr*prad->nang] * prad->wmu(n);
               er += ir_weight;
               frx += ir_weight * prad->mu(0,k,j,i,n);
               fry += ir_weight * prad->mu(1,k,j,i,n);
               frz += ir_weight * prad->mu(2,k,j,i,n);
           }
           er *= prad->wfreq(ifr);
           frx *= prad->wfreq(ifr);
           fry *= prad->wfreq(ifr);
           frz *= prad->wfreq(ifr);
        }
        
        
        // Now apply the radiation source terms to gas with energy and
        // momentum conservation
        u(IM1,k,j,i) += (-prat*(frx- frx0) * invcrat);
        u(IM2,k,j,i) += (-prat*(fry- fry0) * invcrat);
        u(IM3,k,j,i) += (-prat*(frz- frz0) * invcrat);
        u(IEN,k,j,i) += (-prat*(er - er0 ));
        
        //limit the velocity by speed of light
        vx = u(IM1,k,j,i)/u(IDN,k,j,i);
        vy = u(IM2,k,j,i)/u(IDN,k,j,i);
        vz = u(IM3,k,j,i)/u(IDN,k,j,i);
        vel = vx*vx+vy*vy+vz*vz;
        vel = sqrt(vel);
        ratio = vel * invcrat;
        if(ratio > prad->vmax){
          Real factor = prad->vmax/ratio;
          u(IM1,k,j,i) *= factor;
          u(IM2,k,j,i) *= factor;
          u(IM3,k,j,i) *= factor;
        
        }
        
      }// end i
    }// end j
  }// end k
    


  
}


void RadIntegrator::PredictVel(int k, int j, int i, Real dt, Real rho, Real *vx,
                              Real *vy, Real *vz)
{
    Radiation *prad = pmy_rad;
  
    Real &prat = prad->prat;
    Real invcrat = 1.0/prad->crat;
    Real ct = dt * prad->reduced_c;
  
    Real &er =prad->rad_mom(IER,k,j,i);
  
    Real &fr1=prad->rad_mom(IFR1,k,j,i);
    Real &fr2=prad->rad_mom(IFR1,k,j,i);
    Real &fr3=prad->rad_mom(IFR1,k,j,i);
  
    Real &pr11=prad->rad_mom(IPR11,k,j,i);
    Real &pr12=prad->rad_mom(IPR12,k,j,i);
    Real &pr13=prad->rad_mom(IPR13,k,j,i);
    Real &pr22=prad->rad_mom(IPR22,k,j,i);
    Real &pr23=prad->rad_mom(IPR23,k,j,i);
    Real &pr33=prad->rad_mom(IPR33,k,j,i);
  
    Real dtcsigma = ct * (prad->grey_sigma_s(k,j,i) + prad->grey_sigma_a(k,j,i));
  
    Real vx0 = (*vx);
    Real vy0 = (*vy);
    Real vz0 = (*vz);
  
    Real m0x = prat * fr1 * invcrat + rho * vx0;
    Real m0y = prat * fr2 * invcrat + rho * vy0;
    Real m0z = prat * fr3 * invcrat + rho * vz0;
  
    Real vx11 = rho * (1.0 + dtcsigma) + prat * dtcsigma * (er + pr11) * invcrat;
    Real vy11 = rho * (1.0 + dtcsigma) + prat * dtcsigma * (er + pr22) * invcrat;
    Real vz11 = rho * (1.0 + dtcsigma) + prat * dtcsigma * (er + pr33) * invcrat;
    Real vx12 = dtcsigma * prat * pr12 * invcrat;
    Real vx13 = dtcsigma * prat * pr13 * invcrat;
    Real vy12 = dtcsigma * prat * pr23 * invcrat;
  
    Real rhs1 = rho * vx0 + dtcsigma * m0x;
    Real rhs2 = rho * vy0 + dtcsigma * m0y;
    Real rhs3 = rho * vz0 + dtcsigma * m0z;
  
    Real factor = vx11 * vy11 * vz11 - vy11 * vx13 * vx13 + 2.0 * vx12 * vx13 * vy12
              - vx11 * vy12 * vy12 - vx12 * vx12 * vz11;
    factor = 1.0/factor;

    (*vx) = factor*(rhs3*(vx12*vy12 - vx13*vy11) + rhs2*(vy12*vx13
                - vx12*vz11) + rhs1*(vy11*vz11 - vy12*vy12));
          
    (*vy) = factor*(rhs3*(vx12*vx13 - vx11*vy12) + rhs2*(vx11*vz11
                - vx13*vx13) + rhs1*(vx13*vy12 - vx12*vz11));
          
    (*vz) = factor*(rhs3*(vx11*vy11 - vx12*vx12) + rhs2*(vx12*vx13
                - vx11*vy12) + rhs1*(vx12*vy12 - vx13*vy11));
  
}








