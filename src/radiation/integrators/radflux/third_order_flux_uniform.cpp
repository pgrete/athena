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
//! \file second_order_flux.cpp
//  \brief  piecewise linear reconstruction
//======================================================================================
#include <stdexcept>  // runtime_error
#include <iostream>
#include <sstream>

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../radiation.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../coordinates/coordinates.hpp"

// this class header
#include "../rad_integrators.hpp"

//--------------------------------------------------------------------------------------
//! \fn RadIntegrator::ThirdOrderFluxX1()
//  \brief 


void RadIntegrator::ThirdOrderFluxX1Uniform(Coordinates *pco, const int kl, const int ku, 
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const AthenaArray<Real> &vel, AthenaArray<Real> &flx)
{
 // CS08 constant used in second derivative limiter, >1 , independent of h
  const Real C2 = 1.25;
  Real qa,qb,qc,qd,qe,rho;

  // 1D scratch arrays
  int ncells1 = iu-il + 2*(NGHOST);
  int tot_ang=pmy_rad->n_fre_ang;
  AthenaArray<Real> dph, qplus, qminus, dqf_plus, dqf_minus, d2qf, d2qc;
  dph.NewAthenaArray(ncells1,tot_ang);
  qplus.NewAthenaArray(ncells1,tot_ang);
  qminus.NewAthenaArray(ncells1,tot_ang);
  dqf_plus.NewAthenaArray(ncells1,tot_ang);
  dqf_minus.NewAthenaArray(ncells1,tot_ang);
  d2qf.NewAthenaArray(ncells1,tot_ang);
  d2qc.NewAthenaArray(ncells1,tot_ang);

  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {

//--- Step 1. ----------------------------------------------------------------------------
// Reconstruct interface averages <a>_{i-1/2} using PPM for uniform mesh (CW eq 1.6)
//#pragma simd
    for (int i=il-1; i<=(iu+1); ++i) {
      for(int n=0; n<tot_ang;++n){
        dph(i,n) = (7.0*(q(k,j,i-1,n)+q(k,j,i,n)) - (q(k,j,i+1,n)+q(k,j,i-2,n)))/12.0;
        d2qc(i,n) = q(k,j,i-1,n) - 2.0*q(k,j,i,n) + q(k,j,i+1,n); //(CD eq 85a) (no 1/2)
      }
    }
    for(int n=0; n<tot_ang;++n)
      d2qc(il-2,n) = q(k,j,il-3,n) - 2.0*q(k,j,il-2,n) + q(k,j,il-1,n);

    // Limit interpolated interface states as in CD section 4.3.1
    // #pragma simd // poor vectorization efficiency
    for (int i=il-1; i<=(iu+1); ++i) {
      for(int n=0; n<tot_ang;++n){
        qa = dph(i,n) - q(k,j,i-1,n); // (CD eq 84a)
        qb = q(k,j,i,n) - dph(i,n);   // (CD eq 84b)
        if (qa*qb < 0.0) { // Local extrema detected at i-1/2 face
          qa = 3.0*(q(k,j,i-1,n) - 2.0*dph(i,n) + q(k,j,i,n));  // (CD eq 85b)
          qb = d2qc(i-1,n); // (CD eq 85a) (no 1/2)
          qc = d2qc(i,n);   // (CD eq 85c) (no 1/2)
          qd = 0.0;
          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)){
            qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
          }
          dph(i,n) = 0.5*(q(k,j,i-1,n)+q(k,j,i,n)) - qd/6.0;
        }
      }
    }
//#pragma simd
    for (int i=il-1; i<=iu; ++i) {
      for(int n=0; n<tot_ang;++n){
        qminus(i,n) = dph(i,n);
        qplus(i,n) =  dph(i+1,n);
      }
    }

//--- Step 2. ----------------------------------------------------------------------------
// Compute cell-centered difference stencils (MC section 2.4.1)
//#pragma simd
    for (int i=il-1; i<=iu; ++i) {
      for(int n=0; n<tot_ang;++n){
        dqf_minus(i,n) = q(k,j,i,n) - qminus(i,n); // (CS eq 25)
        dqf_plus(i,n)  = qplus(i,n) - q(k,j,i,n);
        d2qf(i,n) = 6.0*(dph(i,n) - 2.0*q(k,j,i,n) + dph(i+1,n)); // a6 coefficient * -2
      }
    }

//--- Step 3. ----------------------------------------------------------------------------
// Apply CS limiters to parabolic interpolant
    // #pragma simd // poor vectorization efficiency
    for (int i=il-1; i<=iu; ++i) {
      for(int n=0; n<tot_ang;++n){
        qa = dqf_minus(i,n)*dqf_plus(i,n);
        qb = (q(k,j,i+1,n) - q(k,j,i,n))*(q(k,j,i,n) - q(k,j,i-1,n));

        // Check for local extrema
        if (qa <= 0.0 || qb <= 0.0 ) {
          // Check if extrema is smooth
          qa = d2qc(i-1,n);
          qb = d2qc(i,n);
          qc = d2qc(i+1,n);
          qd = d2qf(i,n);
          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)) {
            // Extrema is smooth
            qe = SIGN(qd)* std::min(std::min(C2*fabs(qa),C2*fabs(qb)),
                                    std::min(C2*fabs(qc),fabs(qd))); // (CS eq 22)
          } else {
            // Extrema is at interface adjacent to a discontinuity: flatten derivative
            qe =0.0;
          }

          // Check if 2nd derivative is close to roundoff error
          qa = std::max(fabs(q(k,j,i-1,n)),fabs(q(k,j,i-2,n)));
          qb = std::max(std::max(fabs(q(k,j,i,n)),fabs(q(k,j,i+1,n))),
                          fabs(q(k,j,i+2,n)));

          if (fabs(qd) <= (1.0e-12)*std::max(qa,qb)) {
            // 2nd derivative magnitude is too small: flatten
            rho = 0.0;
          } else {
            // Limiter is not sensitive to roundoff. Use limited ratio (MC eq 27)
            rho = qe/qd;
          }

          // Check if relative change in limited 2nd deriv is > roundoff
          if (rho <= (1.0 - (1.0e-12))) {
            // Limit smooth extrema
            qminus(i,n) = q(k,j,i,n) - rho*dqf_minus(i,n); // (CS eq 23)
            qplus(i,n) = q(k,j,i,n) + rho*dqf_plus(i,n);
          }

        // No extrema detected
        } else {
          // Overshoot i-1/2,R / i,(-) state
          if (fabs(dqf_minus(i,n)) >= 2.0*fabs(dqf_plus(i,n))) {
            qminus(i,n) = q(k,j,i,n) - 2.0*dqf_plus(i,n);
          }
          // Overshoot i+1/2,L / i,(+) state
          if (fabs(dqf_plus(i,n)) >= 2.0*fabs(dqf_minus(i,n))) {
            qplus(i,n) = q(k,j,i,n) + 2.0*dqf_minus(i,n);
          }
        }
      }
    }

//--- Step 4. ----------------------------------------------------------------------------
// Convert limited cell-centered values to interface-centered L/R Riemann states
// both L/R values defined over [il,iu]
//#pragma simd
    for (int i=il; i<=iu; ++i) {
      for(int n=0; n<tot_ang;++n){
        if(vel(k,j,i,n) > 0.0){
          flx(k,j,i,n) += qplus(i-1,n) * vel(k,j,i,n);
        }else{
          flx(k,j,i,n) += qminus(i,n) * vel(k,j,i,n);
        }
      }
    }

  }}

  dph.DeleteAthenaArray();
  qplus.DeleteAthenaArray();
  qminus.DeleteAthenaArray();
  dqf_plus.DeleteAthenaArray();
  dqf_minus.DeleteAthenaArray();
  d2qf.DeleteAthenaArray();
  d2qc.DeleteAthenaArray();
  return;
}


//========================================================================================
//! \fn Reconstruction::ReconstructionFuncX2()
//  \brief Returns L/R interface values in X2-dir constructed using fourth-order PPM and
//         Mignone limiting over [kl,ku][jl,ju][il,iu]

void RadIntegrator::ThirdOrderFluxX2Uniform(Coordinates *pco, const int kl, const int ku, 
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const AthenaArray<Real> &vel, AthenaArray<Real> &flx)
{



  // CS08 constant used in second derivative limiter, >1 , independent of h
  const Real C2 = 1.25;
  Real qa,qb,qc,qd,qe,rho;
  int tot_ang=pmy_rad->n_fre_ang;

  // 1D scratch arrays
  int ncells1 = iu-il + 2*(NGHOST);
  AthenaArray<Real> dph, dph_jp1, qplus, qminus, dqf_plus, dqf_minus, d2qf;
  AthenaArray<Real> d2qc_jm1, d2qc, d2qc_jp1;
  dph.NewAthenaArray(ncells1,tot_ang);
  dph_jp1.NewAthenaArray(ncells1,tot_ang);
  qplus.NewAthenaArray(ncells1,tot_ang);
  qminus.NewAthenaArray(ncells1,tot_ang);
  dqf_plus.NewAthenaArray(ncells1,tot_ang);
  dqf_minus.NewAthenaArray(ncells1,tot_ang);
  d2qf.NewAthenaArray(ncells1,tot_ang);
  d2qc_jm1.NewAthenaArray(ncells1,tot_ang);
  d2qc.NewAthenaArray(ncells1,tot_ang);
  d2qc_jp1.NewAthenaArray(ncells1,tot_ang);

  for (int k=kl; k<=ku; ++k) {

//--- Step 1a. ---------------------------------------------------------------------------
// Reconstruct interface averages <a>_{j-1/2} at j=jl-1 (CW eq 1.6)

    // initialize interface states along 1-D vector at j=jl-1
//#pragma simd
    for (int i=il; i<=iu; ++i) {
      for(int n=0; n<tot_ang;++n){
        dph(i,n) = ( 7.0*(q(k,jl-2,i,n) + q(k,jl-1,i,n)) -
                       (q(k,jl  ,i,n) + q(k,jl-3,i,n)) )/12.0;
        // Approximate second-derivative at cell center using cell averages +/-
        d2qc_jm1(i,n) = q(k,jl-3,i,n) - 2.0*q(k,jl-2,i,n) + q(k,jl-1,i,n);
        d2qc(i,n)     = q(k,jl-2,i,n) - 2.0*q(k,jl-1,i,n) + q(k,jl  ,i,n);
      }
    }

    // Limit interpolated interface states at j=jl-1 as in CD section 4.3.1
    // #pragma simd // poor vectorization efficiency
    for (int i=il; i<=iu; ++i) {
      for(int n=0; n<tot_ang; ++n){
        qa = dph(i,n) - q(k,jl-2,i,n); // (CD eq 84a)
        qb = q(k,jl-1,i,n) - dph(i,n);   // (CD eq 84b)
        if (qa*qb < 0.0) { // Local extrema detected at i-1/2 face
          qa = 3.0*(q(k,jl-2,i,n) - 2.0*dph(i,n) + q(k,jl-1,i,n)); // (CD eq 85b)
          qb = d2qc_jm1(i,n); // (CD eq 85a)(no 1/2)
          qc = d2qc(i,n);     // (CD eq 85c)(no 1/2)
          qd = 0.0;
          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)){
            qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
          }
          dph(i,n) = 0.5*(q(k,jl-2,i,n)+q(k,jl-1,i,n)) - qd/6.0;
        }
      }
    }

// start loop over all j
    for (int j=jl-1; j<=ju; ++j) {

//--- Step 1b. ---------------------------------------------------------------------------
// Reconstruct interface averages <a>_{j-1/2} at j=j+1 (CW eq 1.6)
//#pragma simd
      for (int i=il; i<=iu; ++i) {
        for(int n=0; n<tot_ang; ++n){
          dph_jp1(i,n) = ( 7.0*(q(k,j  ,i,n) + q(k,j+1,i,n)) -
                           (q(k,j+2,i,n) + q(k,j-1,i,n)) )/12.0;
          d2qc_jp1(i,n) = q(k,j,i,n) - 2.0*q(k,j+1,i,n) + q(k,j+2,i,n);
        }
      }

      // Limit interpolated interface states at j=j+1 as in CD section 4.3.1
    // #pragma simd // poor vectorization efficiency
      for (int i=il; i<=iu; ++i) {
        for(int n=0; n<tot_ang; ++n){
          qa = dph_jp1(i,n) - q(k,j,i,n); // (CD eq 84a)
          qb = q(k,j+1,i,n) - dph_jp1(i,n); // (CD eq 84b)
          if (qa*qb < 0.0) { // Local extrema detected at i-1/2 face
            qa = 3.0*(q(k,j,i,n) - 2.0*dph_jp1(i,n) + q(k,j+1,i,n)); // (CD eq 85b)
            qb = d2qc(i,n);      // (CD eq 85a)(no 1/2)
            qc = d2qc_jp1(i,n);  // (CD eq 85c)(no 1/2)
            qd = 0.0;
            if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)){
              qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
            }
            dph_jp1(i,n) = 0.5*(q(k,j,i,n)+q(k,j+1,i,n)) - qd/6.0;
          }
        }
      }

      // Initialize cell-indexed interface states / parabolic coefficients
//#pragma simd
      for (int i=il; i<=iu; ++i) {
        for(int n=0; n<tot_ang; ++n){
          qminus(i,n) = dph(i,n);       // value at j
          qplus(i,n)  = dph_jp1(i,n);   // value at j+1
        }
      }

//--- Step 2. ----------------------------------------------------------------------------
// Compute cell-centered difference stencils (MC section 2.4.1)
//#pragma simd
      for (int i=il; i<=iu; ++i) {
        for(int n=0; n<tot_ang; ++n){
          dqf_minus(i,n) = q(k,j,i,n) - qminus(i,n);
          dqf_plus(i,n) = qplus(i,n) - q(k,j,i,n);
          d2qf(i,n) = 6.0*(dph(i,n) -2.0*q(k,j,i,n) + dph_jp1(i,n)); // a6 coefficient * -2
        }
      }

//--- Step 3. ----------------------------------------------------------------------------
// Apply CS limiters to parabolic interpolant
    // #pragma simd // poor vectorization efficiency
      for (int i=il; i<=iu; ++i) {
        for(int n=0; n<tot_ang; ++n){
          qa = dqf_minus(i,n)*dqf_plus(i,n);
          qb = (q(k,j+1,i,n) - q(k,j,i,n))*(q(k,j,i,n) - q(k,j-1,i,n));

          // check for local extrema
          if (qa <= 0.0 || qb <= 0.0 ) {
            // Check if extrema is smooth
            qa = d2qc_jm1(i,n);
            qb = d2qc(i,n);
            qc = d2qc_jp1(i,n);
            qd = d2qf(i,n);
            if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)) {
              // Extrema is smooth
              qe = SIGN(qd)* std::min(std::min(C2*fabs(qa),C2*fabs(qb)),
                                      std::min(C2*fabs(qc),fabs(qd))); // (CS eq 22)
            } else {
              // Extrema is at interface adjacent to a discontinuity: flatten derivative
              qe = 0.0;
            }

            // Check if 2nd derivative is close to roundoff error
            qa = std::max(fabs(q(k,j-1,i,n)),fabs(q(k,j-2,i,n)));
            qb = std::max(std::max(fabs(q(k,j,i,n)),fabs(q(k,j+1,i,n))),
                          fabs(q(k,j+2,i,n)));

            if (fabs(qd) <= (1.0e-12)*std::max(qa,qb)) {
              // 2nd derivative magnitude is too small: flatten
              rho = 0.0;
            } else {
              // Limiter is not sensitive to roundoff. Use limited ratio (MC eq 27)
              rho = qe/qd;
            }

            // Check if relative change in limited 2nd deriv is > roundoff
            if (rho <= (1.0 - (1.0e-12))) {
              // Limit smooth extrema
              qminus(i,n) = q(k,j,i,n) - rho*dqf_minus(i,n); // (CS eq 23)
              qplus(i,n)  = q(k,j,i,n) + rho*dqf_plus(i,n);
            }

          // No extrema detected
          } else {
            // Overshoot j-1/2,R / j,(-) state
            if (fabs(dqf_minus(i,n)) >= 2.0*fabs(dqf_plus(i,n))) {
              qminus(i,n) = q(k,j,i,n) - 2.0*dqf_plus(i,n);
            }
            // Overshoot j+1/2,L / j,(+) state
            if (fabs(dqf_plus(i,n)) >= 2.0*fabs(dqf_minus(i,n))) {
              qplus(i,n) = q(k,j,i,n) + 2.0*dqf_minus(i,n);
            }
          }
        }
      }

//--- Step 4. ----------------------------------------------------------------------------
// Convert limited cell-centered values to interface-centered L/R Riemann states
// both L/R values defined over [jl,ju]
//#pragma simd
      for (int i=il; i<=iu; ++i) {
        for(int n=0; n<tot_ang;++n){
          if(vel(k,j,i,n) > 0.0){
            flx(k,j,i,n) += qplus(i-1,n) * vel(k,j,i,n);
          }else{
            flx(k,j,i,n) += qminus(i,n) * vel(k,j,i,n);
          }
        }
      }

      // Copy 1D temporary arrays for next value of j unless j-loop finished
      if (j < ju) {
//#pragma simd
        for (int i=il; i<=iu; ++i) {
          for(int n=0; n<tot_ang; ++n){
            dph(i,n) = dph_jp1(i,n);
            d2qc_jm1(i,n) = d2qc    (i,n);
            d2qc    (i,n) = d2qc_jp1(i,n);
          }
        }
      }

    } // end loop over [jl-1,ju]
  }   // end loop over [kl,ku]

  dph.DeleteAthenaArray();
  dph_jp1.DeleteAthenaArray();
  qplus.DeleteAthenaArray();
  qminus.DeleteAthenaArray();
  dqf_plus.DeleteAthenaArray();
  dqf_minus.DeleteAthenaArray();
  d2qf.DeleteAthenaArray();
  d2qc_jm1.DeleteAthenaArray();
  d2qc.DeleteAthenaArray();
  d2qc_jp1.DeleteAthenaArray();
  return;
}










//========================================================================================
//! \fn Reconstruction::ReconstructionFuncX3()
//  \brief Returns L/R interface values in X3-dir constructed using fourth-order PPM and
//         Mignone limiting over [kl,ku][jl,ju][il,iu]

void RadIntegrator::ThirdOrderFluxX3Uniform(Coordinates *pco, const int kl, const int ku, 
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const AthenaArray<Real> &vel, AthenaArray<Real> &flx)
{


  // CS08 constant used in second derivative limiter
  Real C2 = 1.25; // >1 , independent of h
  Real qa,qb,qc,qd,qe,rho;

  int tot_ang=pmy_rad->n_fre_ang;

  // 1D scratch arrays
  int ncells1 = iu-il + 2*(NGHOST);
  AthenaArray<Real> dph, dph_kp1, qplus, qminus, dqf_plus, dqf_minus, d2qf;
  AthenaArray<Real> d2qc_km1, d2qc, d2qc_kp1;
  dph.NewAthenaArray(ncells1,tot_ang);
  dph_kp1.NewAthenaArray(ncells1,tot_ang);
  qplus.NewAthenaArray(ncells1,tot_ang);
  qminus.NewAthenaArray(ncells1,tot_ang);
  dqf_plus.NewAthenaArray(ncells1,tot_ang);
  dqf_minus.NewAthenaArray(ncells1,tot_ang);
  d2qf.NewAthenaArray(ncells1,tot_ang);
  d2qc_km1.NewAthenaArray(ncells1,tot_ang);
  d2qc.NewAthenaArray(ncells1,tot_ang);
  d2qc_kp1.NewAthenaArray(ncells1,tot_ang);

  for (int j=jl; j<=ju; ++j) {

//--- Step 1a. ---------------------------------------------------------------------------
// Reconstruct interface averages <a>_{k-1/2} at k=kl-1 (CW eq 1.6)

    // initialize interface states along 1-D vector at k=kl-1
//#pragma simd
    for (int i=il; i<=iu; ++i) {
      for(int n=0; n<tot_ang; ++n){
        dph(i,n) = ( 7.0*(q(kl-2,j,i,n) + q(kl-1,j,i,n)) -
                       (q(kl  ,j,i,n) + q(kl-3,j,i,n)) )/12.0;
        // Approximate second-derivative at cell center using cell averages +/-
        d2qc_km1(i,n) = q(kl-3,j,i,n) - 2.0*q(kl-2,j,i,n) + q(kl-1,j,i,n);
        d2qc(i,n)     = q(kl-2,j,i,n) - 2.0*q(kl-1,j,i,n) + q(kl  ,j,i,n);
      }
    }

    // Limit interpolated interface states at k=kl-1 as in CD section 4.3.1
    // #pragma simd // poor vectorization efficiency
    for (int i=il; i<=iu; ++i) {
      for(int n=0; n<tot_ang; ++n){
        qa = dph(i,n) - q(kl-2,j,i,n);   // (CD eq 84a)
        qb = q(kl-1,j,i,n) - dph(i,n);   // (CD eq 84b)
        if (qa*qb < 0.0) { // Local extrema detected at i-1/2 face
          qa = 3.0*(q(kl-2,j,i,n) - 2.0*dph(i,n) + q(kl-1,j,i,n)); // (CD eq 85b)
          qb = d2qc_km1(i,n); // (CD eq 85a)(no 1/2)
          qc = d2qc(i,n);     // (CD eq 85c)(no 1/2)
          qd = 0.0;
          if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)){
            qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
          }
          dph(i,n) = 0.5*(q(kl-2,j,i,n)+q(kl-1,j,i,n)) - qd/6.0;
        }
      }
    }

// start loop over all k
    for (int k=kl-1; k<=ku; ++k) {

//--- Step 1b. ---------------------------------------------------------------------------
// Reconstruct interface averages <a>_{k-1/2} at k=k+1 (CW eq 1.6)
//#pragma simd
      for (int i=il; i<=iu; ++i) {
        for(int n=0; n<tot_ang; ++n){
          dph_kp1(i,n) = ( 7.0*(q(k  ,j,i,n) + q(k+1,j,i,n)) -
                             (q(k+2,j,i,n) + q(k-1,j,i,n)) )/12.0;
          d2qc_kp1(i,n) = q(k,j,i,n) - 2.0*q(k+1,j,i,n) + q(k+2,j,i,n);
        }
      }

      // Limit interpolated interface states at k=k+1 as in CD section 4.3.1
    // #pragma simd // poor vectorization efficiency
      for (int i=il; i<=iu; ++i) {
        for(int n=0; n<tot_ang; ++n){
          qa = dph_kp1(i,n) - q(k,j,i,n);   // (CD eq 84a)
          qb = q(k+1,j,i,n) - dph_kp1(i,n); // (CD eq 84b)
          if (qa*qb < 0.0) { // Local extrema detected at i-1/2 face
            qa = 3.0*(q(k,j,i,n) - 2.0*dph_kp1(i,n) + q(k+1,j,i,n)); // (CD eq 85b)
            qb = d2qc(i,n);      // (CD eq 85a)(no 1/2)
            qc = d2qc_kp1(i,n);  // (CD eq 85c)(no 1/2)
            qd = 0.0;
            if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc)){
              qd = SIGN(qa)* std::min(C2*fabs(qb),std::min(C2*fabs(qc),fabs(qa)));
            }
            dph_kp1(i,n) = 0.5*(q(k,j,i,n)+q(k+1,j,i,n)) - qd/6.0;
          }
        }
      }

      // Initialize cell-indexed interface states / parabolic coefficients
      for (int i=il; i<=iu; ++i) {
        for(int n=0; n<tot_ang; ++n){
          qminus(i,n) = dph(i,n);       // value at k
          qplus(i,n)  = dph_kp1(i,n);   // value at k+1
        }
      }

//--- Step 2. ----------------------------------------------------------------------------
// Compute cell-centered difference stencils (MC section 2.4.1)
//#pragma simd
      for (int i=il; i<=iu; ++i) {
        for(int n=0; n<tot_ang; ++n){
          dqf_minus(i,n) = q(k,j,i,n) - qminus(i,n);
          dqf_plus(i,n) = qplus(i,n) - q(k,j,i,n);
          d2qf(i,n) = 6.0*(dph(i,n) -2.0*q(k,j,i,n) + dph_kp1(i,n)); // a6 coefficient * -2
        }
      }

//--- Step 3. ----------------------------------------------------------------------------
// Apply CS limiters to parabolic interpolant
    // #pragma simd // poor vectorization efficiency
      for (int i=il; i<=iu; ++i) {
        for(int n=0; n<tot_ang; ++n){
          qa = dqf_minus(i,n)*dqf_plus(i,n);
          qb = (q(k+1,j,i,n) - q(k,j,i,n))*(q(k,j,i,n) - q(k-1,j,i,n));

          // check for local extrema
          if (qa <= 0.0 || qb <= 0.0 ) {
            // Check if extrema is smooth
            qa = d2qc_km1(i,n);
            qb = d2qc(i,n);
            qc = d2qc_kp1(i,n);
            qd = d2qf(i,n);
            if (SIGN(qa) == SIGN(qb) && SIGN(qa) == SIGN(qc) && SIGN(qa) == SIGN(qd)) {
              // Extrema is smooth
              qe = SIGN(qd)* std::min(std::min(C2*fabs(qa),C2*fabs(qb)),
                                      std::min(C2*fabs(qc),fabs(qd))); // (CS eq 22)
            } else {
              // Extrema is at interface adjacent to a discontinuity: flatten derivative
              qe = 0.0;
            }

            // Check if 2nd derivative is close to roundoff error
            qa = std::max(fabs(q(k-1,j,i,n)),fabs(q(k-2,j,i,n)));
            qb = std::max(std::max(fabs(q(k,j,i,n)),fabs(q(k+1,j,i,n))),
                          fabs(q(k+2,j,i,n)));

            if (fabs(qd) <= (1.0e-12)*std::max(qa,qb)) {
              // 2nd derivative magnitude is too small: flatten
              rho = 0.0;
            } else {
              // Limiter is not sensitive to roundoff. Use limited ratio (MC eq 27)
              rho = qe/qd;
            }

            // Check if relative change in limited 2nd deriv is > roundoff
            if (rho <= (1.0 - (1.0e-12))) {
              // Limit smooth extrema
              qminus(i,n) = q(k,j,i,n) - rho*dqf_minus(i,n); // (CS eq 23)
              qplus(i,n)  = q(k,j,i,n) + rho*dqf_plus(i,n);
            }

          // No extrema detected
          } else {
            // Overshoot k-1/2,R / j,(-) state
            if (fabs(dqf_minus(i,n)) >= 2.0*fabs(dqf_plus(i,n))) {
              qminus(i,n) = q(k,j,i,n) - 2.0*dqf_plus(i,n);
            }
            // Overshoot k+1/2,L / j,(+) state
            if (fabs(dqf_plus(i,n)) >= 2.0*fabs(dqf_minus(i,n))) {
              qplus(i,n) = q(k,j,i,n) + 2.0*dqf_minus(i,n);
            }
          }
        }
      }

//--- Step 4. ----------------------------------------------------------------------------
// Convert limited cell-centered values to interface-centered L/R Riemann states
// both L/R values defined over [jl,ju]
//#pragma simd
      for (int i=il; i<=iu; ++i) {
        for(int n=0; n<tot_ang;++n){
          if(vel(k,j,i,n) > 0.0){
            flx(k,j,i,n) += qplus(i-1,n) * vel(k,j,i,n);
          }else{
            flx(k,j,i,n) += qminus(i,n) * vel(k,j,i,n);
          }
        }
      }

      // Copy 1D temporary arrays for next value of k unless k-loop finished
      if (k < ku) {
//#pragma simd
        for (int i=il; i<=iu; ++i) {
          for(int n=0; n<tot_ang; ++n){
            dph(i,n) = dph_kp1(i,n);
            d2qc_km1(i,n) = d2qc    (i,n);
            d2qc    (i,n) = d2qc_kp1(i,n);
          }
        }
      }

    } // end loop over [kl-1,ku]
  }   // end loop over [jl,ju]

  dph.DeleteAthenaArray();
  dph_kp1.DeleteAthenaArray();
  qplus.DeleteAthenaArray();
  qminus.DeleteAthenaArray();
  dqf_plus.DeleteAthenaArray();
  dqf_minus.DeleteAthenaArray();
  d2qf.DeleteAthenaArray();
  d2qc_km1.DeleteAthenaArray();
  d2qc.DeleteAthenaArray();
  d2qc_kp1.DeleteAthenaArray();
  return;
}


