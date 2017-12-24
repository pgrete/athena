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


void RadIntegrator::ThirdOrderFluxX1(Coordinates *pco, const int kl, const int ku, 
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const AthenaArray<Real> &vel, AthenaArray<Real> &flx)
{
  // 1D scratch arrays
  int ncells1 = iu-il + 2*(NGHOST);
  AthenaArray<Real> dph, qplus, qminus, dqf_plus, dqf_minus;
  AthenaArray<Real> dd, c1, c2, c3, c4, c5, c6;
  int tot_ang=pmy_rad->n_fre_ang;

  dph.NewAthenaArray(ncells1,tot_ang);
  qplus.NewAthenaArray(ncells1,tot_ang);
  qminus.NewAthenaArray(ncells1,tot_ang);
  dqf_plus.NewAthenaArray(ncells1,tot_ang);
  dqf_minus.NewAthenaArray(ncells1,tot_ang);
  dd.NewAthenaArray(ncells1,tot_ang);

  c1.NewAthenaArray(ncells1);
  c2.NewAthenaArray(ncells1);
  c3.NewAthenaArray(ncells1);
  c4.NewAthenaArray(ncells1);
  c5.NewAthenaArray(ncells1);
  c6.NewAthenaArray(ncells1);

  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
  //--- Step 1. ----------------------------------------------------------------------------
  // Reconstruct interface averages <a>_{i-1/2} using PPM for nonuniform mesh (CW eq 1.6)
  // (slow since repeated for all NWAVE interpolations)
  // These coefficients are also independent of x2, x3 position/ j,k indices. Could store
  // for outer loops instead of repeating?

  // Compute average slope in ith zone
    for (int i=il-2; i<=iu+1; ++i){
      // Compute coefficients used in interpolation formulae
      Real& dx_im1 = pco->dx1f(i-1);
      Real& dx_i   = pco->dx1f(i);
      Real& dx_ip1   = pco->dx1f(i+1);
      Real qa = dx_i/(dx_im1 + dx_i + dx_ip1); // Outermost coefficient in CW eq 1.7
      c1(i) = qa*(2.0*dx_im1+dx_i)/(dx_ip1 + dx_i); // First term in CW eq 1.7 bracket
      c2(i) = qa*(2.0*dx_ip1+dx_i)/(dx_im1 + dx_i); // Second term in CW eq 1.7 bracket

      // Above c1, c2 only used in this loop to compute average zone slope, dd(i)
      // So, either turn c1(), c2() to scalars, or compute and store arrays for all NHYDRO
      // and split this loop into two parts: only below this line inside loop over NHYDRO
      for(int n=0; n<tot_ang; ++n){
        Real dplus = q(k,j,i+1,n) - q(k,j,i,n);
        Real dmnus = q(k,j,i,n) - q(k,j,i-1,n);
        dd(i,n) = c1(i)*dplus + c2(i)*dmnus;
      }
    }

    // Construct the interpolation coefficients c*()
    for (int i=il-1; i<=iu+1; ++i){
      Real& dx_im2 = pco->dx1f(i-2);
      Real& dx_im1 = pco->dx1f(i-1);
      Real& dx_i   = pco->dx1f(i);
      Real& dx_ip1   = pco->dx1f(i+1);
      Real qa = dx_im2 + dx_im1 + dx_i + dx_ip1;
      Real qb = dx_im1/(dx_im1 + dx_i);
      Real qc = (dx_im2 + dx_im1)/(2.0*dx_im1 + dx_i);
      Real qd = (dx_ip1 + dx_i)/(2.0*dx_i + dx_im1);
      qb = qb + 2.0*dx_i*qb/qa*(qc-qd);
      c3(i) = 1.0 - qb;
      c4(i) = qb;
      c5(i) = dx_i/qa*qd;
      c6(i) = -dx_im1/qa*qc;
    }

    for (int i=il-1; i<=(iu+1); ++i) {
      // Fourth-order reconstruction function for nonuniform mesh spacing
      for(int n=0; n<tot_ang; ++n){
        dph(i,n)= c3(i)*q(k,j,i-1,n) + c4(i)*q(k,j,i,n) + c5(i)*dd(i-1,n) + c6(i)*dd(i,n);
      }
    }

    for (int i=il-1; i<=iu; ++i) {
      for(int n=0; n<tot_ang;++n){
        qminus(i,n) = dph(i,n);
        qplus(i,n) =  dph(i+1,n);
      // Apply strict monotonicity constraints to interpolated states (Mignone eq 45)
      // (this step could be applied to dph() instead of 2x times to local cell states)
        qminus(i,n) = std::min(qminus(i,n), std::max(q(k,j,i,n),q(k,j,i-1,n)));
        qplus(i,n) = std::min(qplus(i,n), std::max(q(k,j,i,n),q(k,j,i+1,n)));

        qminus(i,n) = std::max(qminus(i,n), std::min(q(k,j,i,n),q(k,j,i-1,n)));
        qplus(i,n) = std::max(qplus(i,n), std::min(q(k,j,i,n),q(k,j,i+1,n)));
      }
    }

    //--- Step 2. ----------------------------------------------------------------------------
    // Compute cell-centered difference stencils (Mignone \delta Q^+_i, -\delta Q^-_i)

    for (int i=il-1; i<=iu; ++i) {
      for(int n=0; n<tot_ang;++n){
        dqf_minus(i,n) = q(k,j,i,n) - qminus(i,n);
        dqf_plus(i,n)  = qplus(i,n) - q(k,j,i,n);
      }
    }

    //--- Step 3. ----------------------------------------------------------------------------
    // Apply Mignone limiters to parabolic interpolant

    for (int i=il-1; i<=iu; ++i) {
      for(int n=0; n<tot_ang; ++n){
        Real qa = dqf_minus(i,n)*dqf_plus(i,n);
        // Mignone limiter does not check for cell-averaged extrema:
        Real qb = (q(k,j,i+1,n) - q(k,j,i,n))*(q(k,j,i,n) - q(k,j,i-1,n));

        // Check for local extrema (Mignone eq 46)
        if (qa <= 0.0) { // Local extrema detected
          qminus(i,n) = q(k,j,i,n);
          qplus(i,n) = q(k,j,i,n);

        } else { // No extrema detected
          // Overshoot i-1/2,R / i,(-) state
          if (fabs(dqf_minus(i,n)) >= 2.0*fabs(dqf_plus(i,n))) {
            qminus(i,n) = q(k,j,i,n) - 2.0*dqf_plus(i,n);
          }
          // Overshoot i+1/2,L / i,(+) state
          if (fabs(dqf_plus(i,n)) >= 2.0*fabs(dqf_minus(i,n))) {
            qplus(i,n) = q(k,j,i,n) + 2.0*dqf_minus(i,n);
          }
        }
      }// end n
    }

    //--- Step 4. ----------------------------------------------------------------------------
    // Convert limited cell-centered values to interface-centered L/R Riemann states
    // both L/R values defined over [il,iu]

    for (int i=il; i<=iu; ++i) {
      for(int n=0; n<tot_ang;++n){
        if(vel(k,j,i,n) > 0.0){
          flx(k,j,i,n) += qplus(i-1,n) * vel(k,j,i,n);
        }else{
          flx(k,j,i,n) += qminus(i,n) * vel(k,j,i,n);
        }
      }
    }

  }}// end k,j


  dph.DeleteAthenaArray();
  qplus.DeleteAthenaArray();
  qminus.DeleteAthenaArray();
  dqf_plus.DeleteAthenaArray();
  dqf_minus.DeleteAthenaArray();
  return;
}


//========================================================================================
//! \fn Reconstruction::ReconstructionFuncX2()
//  \brief Returns L/R interface values in X2-dir constructed using fourth-order PPM and
//         Mignone limiting over [kl,ku][jl,ju][il,iu]

void RadIntegrator::ThirdOrderFluxX2(Coordinates *pco, const int kl, const int ku, 
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const AthenaArray<Real> &vel, AthenaArray<Real> &flx)
{
  Real qa,qb,qc,qd;
  // 1D scratch arrays
  int ncells1 = iu-il + 2*(NGHOST);
  AthenaArray<Real> dph, dph_jp1, qplus, qminus, dqf_plus, dqf_minus;
  AthenaArray<Real> dd, dd_jm1;
  Real c1, c2, c3, c4, c5, c6;
  Real c1_jm1, c2_jm1; // only needed for the initialization of jl-2 ave slopes

  int tot_ang=pmy_rad->n_fre_ang;

  dph.NewAthenaArray(ncells1,tot_ang);
  dph_jp1.NewAthenaArray(ncells1,tot_ang);
  qplus.NewAthenaArray(ncells1,tot_ang);
  qminus.NewAthenaArray(ncells1,tot_ang);
  dqf_plus.NewAthenaArray(ncells1,tot_ang);
  dqf_minus.NewAthenaArray(ncells1,tot_ang);

  // No need to store these x1-sliced arrays, since they are merely intermediate
  // quantities used for the i,j-1 and i,j and i,j+1 entries of dph(), dph_jp1()
  dd.NewAthenaArray(ncells1,tot_ang);
  dd_jm1.NewAthenaArray(ncells1,tot_ang);

  for (int k=kl; k<=ku; ++k) {

//--- Step 1a. ---------------------------------------------------------------------------
// Reconstruct interface averages <a>_{j-1/2} using PPM for nonuniform mesh (CW eq 1.6)
// at j=jl-1 (slow since repeated for all NWAVE interpolations)

    // Compute coefficients used in interpolation formulae
    // These coefficients are independent of x1 position/ i-index
    //    Real& dx_jm3 = pco->dx2f(jl-3);
    Real& dx_jm2 = pco->dx2f(jl-3);
    Real& dx_jm1 = pco->dx2f(jl-2);
    Real& dx_j   = pco->dx2f(jl-1);
    Real& dx_jp1 = pco->dx2f(jl);

    // x1 slice at j=jl-1
    qa = dx_j/(dx_jm1 + dx_j + dx_jp1); // Outermost coefficient in CW eq 1.7
    c1 = qa*(2.0*dx_jm1+dx_j)/(dx_jp1 + dx_j); // First term in CW eq 1.7 bracket
    c2 = qa*(2.0*dx_jp1+dx_j)/(dx_jm1 + dx_j); // Second term in CW eq 1.7 bracket

    // x1 slice at j=jl-2 (only needed for ave slopes in the jl-1 interpolation)
    qa = dx_jm1/(dx_jm2 + dx_jm1 + dx_j);
    c1_jm1 = qa*(2.0*dx_jm2+dx_jm1)/(dx_j + dx_jm1);
    c2_jm1 = qa*(2.0*dx_j+dx_jm1)/(dx_jm2 + dx_jm1);

    // Construct the interpolation coefficients c
    // x1 slice at j=jl-1
    qa = dx_jm2 + dx_jm1 + dx_j + dx_jp1;
    qb = dx_jm1/(dx_jm1 + dx_j);
    qc = (dx_jm2 + dx_jm1)/(2.0*dx_jm1 + dx_j);
    qd = (dx_jp1 + dx_j)/(2.0*dx_j + dx_jm1);
    qb = qb + 2.0*dx_j*qb/qa*(qc-qd);
    c3 = 1.0 - qb;
    c4 = qb;
    c5 = dx_j/qa*qd;
    c6 = -dx_jm1/qa*qc;

    for (int i=il; i<=iu; ++i){
      for(int n=0; n<tot_ang;++n){
        // Above c1, c2 only used in this loop to compute average zone slope, dd(i)
        // compute jl-1 average zone slopes
        Real dplus = q(k,jl,i,n) - q(k,jl-1,i,n);
        Real dmnus = q(k,jl-1,i,n) - q(k,jl-2,i,n);
        dd(i,n) = c1*dplus + c2*dmnus;

        // compute jl-2 average zone slopes
        dplus = dmnus;
        dmnus = q(k,jl-2,i,n) - q(k,jl-3,i,n);
        dd_jm1(i,n) = c1_jm1*dplus + c2_jm1*dmnus;

        // initialize interface states along 1-D vector at j=jl-1
        // Fourth-order reconstruction function for nonuniform mesh spacing
        dph(i,n)= c3*q(k,jl-2,i,n) + c4*q(k,jl-1,i,n) + c5*dd_jm1(i,n) + c6*dd(i,n);
      }
    }

// start loop over all j
    for (int j=jl-1; j<=ju; ++j) {

//--- Step 1b. ---------------------------------------------------------------------------
// Reconstruct interface averages <a>_{j-1/2} at j=j+1 (CW eq 1.6)
      // Compute coefficients used in interpolation formulae
      // for x1 slice at j+1 x2 index
      dx_jm2 = pco->dx2f(j-1);
      dx_jm1 = pco->dx2f(j);
      dx_j   = pco->dx2f(j+1);
      dx_jp1 = pco->dx2f(j+2);

      // reuse coefficients from x1 slice at j
      // Not needed, since we are reusing dd_jm1() array
      c1_jm1 = c1;
      c2_jm1 = c2;

      // compute coefficients for x1 slice at j+1
      qa = dx_j/(dx_jm1 + dx_j + dx_jp1); // Outermost coefficient in CW eq 1.7
      c1 = qa*(2.0*dx_jm1+dx_j)/(dx_jp1 + dx_j); // First term in CW eq 1.7 bracket
      c2 = qa*(2.0*dx_jp1+dx_j)/(dx_jm1 + dx_j); // Second term in CW eq 1.7 bracket

      // Construct the interpolation coefficients c
      // x1 slice at j+1
      qa = dx_jm2 + dx_jm1 + dx_j + dx_jp1;
      qb = dx_jm1/(dx_jm1 + dx_j);
      qc = (dx_jm2 + dx_jm1)/(2.0*dx_jm1 + dx_j);
      qd = (dx_jp1 + dx_j)/(2.0*dx_j + dx_jm1);
      qb = qb + 2.0*dx_j*qb/qa*(qc-qd);
      c3 = 1.0 - qb;
      c4 = qb;
      c5 = dx_j/qa*qd;
      c6 = -dx_jm1/qa*qc;

      for (int i=il; i<=iu; ++i){
        // Above c1, c2 only used in this loop to compute average zone slope, dd(i)
        // reuse j average zone slopes
        for(int n=0; n<tot_ang;++n){
          dd_jm1(i,n) = dd(i,n);
          // compute j+1 average zone slopes
          Real dplus = q(k,j+2,i,n) - q(k,j+1,i,n);
          Real dmnus = q(k,j+1,i,n) - q(k,j,i,n);
          dd(i,n) = c1*dplus + c2*dmnus;

          // initialize interface states along 1-D vector at j
          // Fourth-order reconstruction function for nonuniform mesh spacing
          dph_jp1(i,n)= c3*q(k,j,i,n) + c4*q(k,j+1,i,n) + c5*dd_jm1(i,n) + c6*dd(i,n);
        }
      }
      // Initialize cell-indexed interface states / parabolic coefficients
      for (int i=il; i<=iu; ++i) {
        for(int n=0; n<tot_ang;++n){
          qminus(i,n) = dph(i,n);       // value at j
          qplus(i,n)  = dph_jp1(i,n);   // value at j+1
          // Apply strict monotonicity constraints to interpolated states (Mignone eq 45)
          // (this step could be applied to dph() instead of 2x times to local cell states)
          qminus(i,n) = std::min(qminus(i,n), std::max(q(k,j,i,n),q(k,j-1,i,n)));
          qplus(i,n) = std::min(qplus(i,n), std::max(q(k,j,i,n),q(k,j+1,i,n)));

          qminus(i,n) = std::max(qminus(i,n), std::min(q(k,j,i,n),q(k,j-1,i,n)));
          qplus(i,n) = std::max(qplus(i,n), std::min(q(k,j,i,n),q(k,j+1,i,n)));
        }
      }

//--- Step 2. ----------------------------------------------------------------------------
// Compute cell-centered difference stencils (Mignone \delta Q^+_j, -\delta Q^-_j)
      for (int i=il; i<=iu; ++i) {
        for(int n=0; n<tot_ang;++n){
          dqf_minus(i,n) = q(k,j,i,n) - qminus(i,n);
          dqf_plus(i,n) = qplus(i,n) - q(k,j,i,n);
        }
      }

//--- Step 3. ----------------------------------------------------------------------------
// Apply Mignone limiters to parabolic interpolant

      for (int i=il; i<=iu; ++i) {
        for(int n=0; n<tot_ang;++n){
          qa = dqf_minus(i,n)*dqf_plus(i,n);
          // Mignone limiter does not check for cell-averaged extrema:
          qb = (q(k,j+1,i,n) - q(k,j,i,n))*(q(k,j,i,n) - q(k,j-1,i,n));

          // Check for local extrema (Mignone eq 46)
          if (qa <= 0.0) { // Local extrema detected
            qminus(i,n) = q(k,j,i,n);
            qplus(i,n)  = q(k,j,i,n);
          } else { // No extrema detected
            // Overshoot j-1/2,R / j,(-) state
            if (fabs(dqf_minus(i,n)) >= 2.0*fabs(dqf_plus(i,n))) {
              qminus(i,n) = q(k,j,i,n) - 2.0*dqf_plus(i,n);
            }
            // Overshoot j+1/2,L / j,(+) state
            if (fabs(dqf_plus(i,n)) >= 2.0*fabs(dqf_minus(i,n))) {
              qplus(i,n) = q(k,j,i,n) + 2.0*dqf_minus(i,n);
            }
          }
        }//end n
      }

//--- Step 4. ----------------------------------------------------------------------------
// Convert limited cell-centered values to interface-centered L/R Riemann states
// both L/R values defined over [jl,ju]

      for (int i=il; i<=iu; ++i) {
        for(int n=0; n<tot_ang;++n){
          if(vel(k,j,i,n) > 0.0){
            flx(k,j,i,n) += qplus(i-1,n) * vel(k,j,i,n);
          }else{
            flx(k,j,i,n) += qminus(i,n) * vel(k,j,i,n);
          }
        }
      }// end i

      // Copy 1D temporary arrays for next value of j unless j-loop finished
      if (j < ju) {
        for (int i=il; i<=iu; ++i) {
          for(int n=0; n<tot_ang;++n){
            dph(i,n) = dph_jp1(i,n);
          }
        }
      }

    } // end loop over [jl-1,ju]
  }   // end loop over [kl,ku]
  dd.DeleteAthenaArray();
  dd_jm1.DeleteAthenaArray();
  dph.DeleteAthenaArray();
  dph_jp1.DeleteAthenaArray();
  qplus.DeleteAthenaArray();
  qminus.DeleteAthenaArray();
  dqf_plus.DeleteAthenaArray();
  dqf_minus.DeleteAthenaArray();
  return;
}

//========================================================================================
//! \fn Reconstruction::ReconstructionFuncX3()
//  \brief Returns L/R interface values in X3-dir constructed using fourth-order PPM and
//         Mignone limiting over [kl,ku][jl,ju][il,iu]

void RadIntegrator::ThirdOrderFluxX3(Coordinates *pco, const int kl, const int ku, 
  const int jl, const int ju, const int il, const int iu, const AthenaArray<Real> &q,
  const AthenaArray<Real> &vel, AthenaArray<Real> &flx)
{
  Real qa,qb,qc,qd;
  // 1D scratch arrays
  int ncells1 = iu-il + 2*(NGHOST);
  AthenaArray<Real> dph, dph_kp1, qplus, qminus, dqf_plus, dqf_minus;
  AthenaArray<Real> dd, dd_km1;
  Real c1, c2, c3, c4, c5, c6;
  Real c1_km1, c2_km1; // only needed for the initialization of kl-2 ave slopes

  int tot_ang=pmy_rad->n_fre_ang;

  dph.NewAthenaArray(ncells1,tot_ang);
  dph_kp1.NewAthenaArray(ncells1,tot_ang);
  qplus.NewAthenaArray(ncells1,tot_ang);
  qminus.NewAthenaArray(ncells1,tot_ang);
  dqf_plus.NewAthenaArray(ncells1,tot_ang);
  dqf_minus.NewAthenaArray(ncells1,tot_ang);

  // No need to store these x1-sliced arrays, since they are merely intermediate
  // quantities used for the i,k-1 and i,k and i,k+1 entries of dph(), dph_kp1()
  dd.NewAthenaArray(ncells1,tot_ang);
  dd_km1.NewAthenaArray(ncells1,tot_ang);

  for (int j=jl; j<=ju; ++j) {

//--- Step 1a. ---------------------------------------------------------------------------
// Reconstruct interface averages <a>_{k-1/2} using PPM for nonuniform mesh (CW eq 1.6)
// at k=kl-1 (slow since repeated for all NWAVE interpolations)

    // Compute coefficients used in interpolation formulae
    // These coefficients are independent of x1 position/ i-index
    //    Real& dx_km3 = pco->dx2f(kl-3);
    Real& dx_km2 = pco->dx2f(kl-3);
    Real& dx_km1 = pco->dx2f(kl-2);
    Real& dx_k   = pco->dx2f(kl-1);
    Real& dx_kp1 = pco->dx2f(kl);

    // x1 slice at k=kl-1
    qa = dx_k/(dx_km1 + dx_k + dx_kp1); // Outermost coefficient in CW eq 1.7
    c1 = qa*(2.0*dx_km1+dx_k)/(dx_kp1 + dx_k); // First term in CW eq 1.7 bracket
    c2 = qa*(2.0*dx_kp1+dx_k)/(dx_km1 + dx_k); // Second term in CW eq 1.7 bracket

    // x1 slice at k=kl-2 (only needed for ave slopes in the kl-1 interpolation)
    qa = dx_km1/(dx_km2 + dx_km1 + dx_k);
    c1_km1 = qa*(2.0*dx_km2+dx_km1)/(dx_k + dx_km1);
    c2_km1 = qa*(2.0*dx_k+dx_km1)/(dx_km2 + dx_km1);

    // Construct the interpolation coefficients c
    // x1 slice at k=kl-1
    qa = dx_km2 + dx_km1 + dx_k + dx_kp1;
    qb = dx_km1/(dx_km1 + dx_k);
    qc = (dx_km2 + dx_km1)/(2.0*dx_km1 + dx_k);
    qd = (dx_kp1 + dx_k)/(2.0*dx_k + dx_km1);
    qb = qb + 2.0*dx_k*qb/qa*(qc-qd);
    c3 = 1.0 - qb;
    c4 = qb;
    c5 = dx_k/qa*qd;
    c6 = -dx_km1/qa*qc;

    for (int i=il; i<=iu; ++i){
      for(int n=0; n<tot_ang;++n){
        // Above c1, c2 only used in this loop to compute average zone slope, dd(i)
        // compute kl-1 average zone slopes
        Real dplus = q(kl,j,i,n) - q(kl-1,j,i,n);
        Real dmnus = q(kl-1,j,i,n) - q(kl-2,j,i,n);
        dd(i,n) = c1*dplus + c2*dmnus;

        // compute kl-2 average zone slopes
        dplus = dmnus;
        dmnus = q(kl-2,j,i,n) - q(kl-3,j,i,n);
        dd_km1(i,n) = c1_km1*dplus + c2_km1*dmnus;

        // initialize interface states along 1-D vector at j=jl-1
        // Fourth-order reconstruction function for nonuniform mesh spacing
        dph(i,n)= c3*q(kl-2,j,i,n) + c4*q(kl-1,j,i,n) + c5*dd_km1(i,n) + c6*dd(i,n);
      }
    }

// start loop over all k
    for (int k=kl-1; k<=ku; ++k) {

//--- Step 1b. ---------------------------------------------------------------------------
// Reconstruct interface averages <a>_{k-1/2} at k=k+1 (CW eq 1.6)
      // Compute coefficients used in interpolation formulae
      // for x1 slice at k+1 x2 index
      dx_km2 = pco->dx2f(k-1);
      dx_km1 = pco->dx2f(k);
      dx_k   = pco->dx2f(k+1);
      dx_kp1 = pco->dx2f(k+2);

      // reuse coefficients from x1 slice at k
      // Not needed, since we are reusing dd_km1() array
      c1_km1 = c1;
      c2_km1 = c2;

      // compute coefficients for x1 slice at k+1
      qa = dx_k/(dx_km1 + dx_k + dx_kp1); // Outermost coefficient in CW eq 1.7
      c1 = qa*(2.0*dx_km1+dx_k)/(dx_kp1 + dx_k); // First term in CW eq 1.7 bracket
      c2 = qa*(2.0*dx_kp1+dx_k)/(dx_km1 + dx_k); // Second term in CW eq 1.7 bracket

      // Construct the interpolation coefficients c
      // x1 slice at k+1
      qa = dx_km2 + dx_km1 + dx_k + dx_kp1;
      qb = dx_km1/(dx_km1 + dx_k);
      qc = (dx_km2 + dx_km1)/(2.0*dx_km1 + dx_k);
      qd = (dx_kp1 + dx_k)/(2.0*dx_k + dx_km1);
      qb = qb + 2.0*dx_k*qb/qa*(qc-qd);
      c3 = 1.0 - qb;
      c4 = qb;
      c5 = dx_k/qa*qd;
      c6 = -dx_km1/qa*qc;

      for (int i=il; i<=iu; ++i){
        for(int n=0; n<tot_ang;++n){
          // Above c1, c2 only used in this loop to compute average zone slope, dd(i)
          // reuse k average zone slopes
          dd_km1(i,n) = dd(i,n);
          // compute k+1 average zone slopes
          Real dplus = q(k+2,j,i,n) - q(k+1,j,i,n);
          Real dmnus = q(k+1,j,i,n) - q(k,j,i,n);
          dd(i,n) = c1*dplus + c2*dmnus;

          // initialize interface states along 1-D vector at k
          // Fourth-order reconstruction function for nonuniform mesh spacing
          dph_kp1(i,n)= c3*q(k,j,i,n) + c4*q(k+1,j,i,n) + c5*dd_km1(i,n) + c6*dd(i,n);
        }
      }

      // Initialize cell-indexed interface states / parabolic coefficients
      for (int i=il; i<=iu; ++i) {
        for(int n=0; n<tot_ang;++n){
          qminus(i,n) = dph(i,n);       // value at k
          qplus(i,n)  = dph_kp1(i,n);   // value at k+1
        // Apply strict monotonicity constraints to interpolated states (Mignone eq 45)
        // (this step could be applied to dph() instead of 2x times to local cell states)
          qminus(i,n) = std::min(qminus(i,n), std::max(q(k,j,i,n),q(k-1,j,i,n)));
          qplus(i,n) = std::min(qplus(i,n), std::max(q(k,j,i,n),q(k+1,j,i,n)));

          qminus(i,n) = std::max(qminus(i,n), std::min(q(k,j,i,n),q(k-1,j,i,n)));
          qplus(i,n) = std::max(qplus(i,n), std::min(q(k,j,i,n),q(k+1,j,i,n)));
        }
      }

//--- Step 2. ----------------------------------------------------------------------------
// Compute cell-centered difference stencils (Mignone \delta Q^+_k, -\delta Q^-_k)

      for (int i=il; i<=iu; ++i) {
        for(int n=0; n<tot_ang;++n){
          dqf_minus(i,n) = q(k,j,i,n) - qminus(i,n);
          dqf_plus(i,n) = qplus(i,n) - q(k,j,i,n);
        }
      }

//--- Step 3. ----------------------------------------------------------------------------
// Apply Mignone limiters to parabolic interpolant

      for (int i=il; i<=iu; ++i) {
        for(int n=0; n<tot_ang;++n){
          qa = dqf_minus(i,n)*dqf_plus(i,n);
          // Mignone limiter does not check for cell-averaged extrema:
          qb = (q(k+1,j,i,n) - q(k,j,i,n))*(q(k,j,i,n) - q(k-1,j,i,n));

          // Check for local extrema (Mignone eq 46)
          if (qa <= 0.0) { // Local extrema detected
            qminus(i,n) = q(k,j,i,n);
            qplus(i,n)  = q(k,j,i,n);
          } else { // No extrema detected
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
        for (int i=il; i<=iu; ++i) {
          dph(i) = dph_kp1(i);
        }
      }

    } // end loop over [kl-1,ku]
  }   // end loop over [jl,ju]

  dd.DeleteAthenaArray();
  dd_km1.DeleteAthenaArray();
  dph.DeleteAthenaArray();
  dph_kp1.DeleteAthenaArray();
  qplus.DeleteAthenaArray();
  qminus.DeleteAthenaArray();
  dqf_plus.DeleteAthenaArray();
  dqf_minus.DeleteAthenaArray();

  return;
}


