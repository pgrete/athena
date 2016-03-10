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

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../radiation.hpp"
#include "../../../mesh.hpp"
#include "../../../coordinates/coordinates.hpp"

// this class header
#include "../rad_integrators.hpp"

//--------------------------------------------------------------------------------------
//! \fn RadIntegrator::SecondOrderFluxX1()
//  \brief 

void RadIntegrator::SecondOrderFluxX1(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &q, const AthenaArray<Real> &vel,
  AthenaArray<Real> &flx)
{
  Coordinates *pco = pmy_rad->pmy_block->pcoord;
  Real dql, dqr, dqc, dq2;

  for (int i=il; i<=iu; ++i){
    Real dx_im2i = 1.0/pco->dx1v(i-2);
    Real dx_im1i = 1.0/pco->dx1v(i-1);
    Real dx_i = 1.0/pco->dx1v(i);
    Real dxfr=pco->x1f(i)-pco->x1v(i-1);
    Real dxfl=pco->x1v(i)-pco->x1f(i);
    Real cfl=pco->dx1v(i-1)/dxfr;
    Real cbl=pco->dx1v(i-2)/(pco->x1v(i-1)-pco->x1f(i-1));
    Real cfr=pco->dx1v(i)/(pco->x1f(i+1)-pco->x1v(i));
    Real cbr=pco->dx1v(i-1)/dxfl;
#pragma simd
    for(int n=0; n<pmy_rad->n_fre_ang; ++n){
      
      dqc = (q(k,j,i,  n) - q(k,j,i-1,n))*dx_im1i;
      
      if(vel(i,n) > 0.0){
        dql = (q(k,j,i-1,n) - q(k,j,i-2,n))*dx_im2i;
        dq2 = dql*dqc;
        Real ql = q(k,j,i-1,n);
      // compute ql_(i-1/2) using Mignone 2014's modified van-Leer limiter
        if(dq2>0.0) {
          ql += dxfr*dq2*(cfl*dql+cbl*dqc)/(dql*dql+(cfl+cbl-2.0)*dq2+dqc*dqc);
        }
        flx(i,n) = ql * vel(i,n);
      }
      if(vel(i,n) <=0.0){
        dqr = (q(k,j,i+1,n) - q(k,j,i,  n))*dx_i;
        dq2=dqc*dqr;
        Real qr = q(k,j,i,n);
        if(dq2>0.0) {
          qr -= dxfl*dq2*(cfr*dqc+cbr*dqr)/(dqc*dqc+(cfr+cbr-2.0)*dq2+dqr*dqr);
        }
       flx(i,n) = qr * vel(i,n);
      }
    }
  }


  return;
}

//--------------------------------------------------------------------------------------
//! \fn RadIntegrator::SecondOrderFluxX2()
//  \brief 

void RadIntegrator::SecondOrderFluxX2(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &q, const AthenaArray<Real> &vel,
  AthenaArray<Real> &flx)
{
  Coordinates *pco = pmy_rad->pmy_block->pcoord;
  Real dql, dqr, dqc, dq2;
  
  Real dx_jm2i = 1.0/pco->dx2v(j-2);
  Real dx_jm1i = 1.0/pco->dx2v(j-1);
  Real dx_ji = 1.0/pco->dx2v(j);
  Real dxfr=pco->x2f(j)-pco->x2v(j-1);
  Real dxfl=pco->x2v(j)-pco->x2f(j);
  Real cfl=pco->dx2v(j-1)/dxfr;
  Real cbl=pco->dx2v(j-2)/(pco->x2v(j-1)-pco->x2f(j-1));
  Real cfr=pco->dx2v(j)/(pco->x2f(j+1)-pco->x2v(j));
  Real cbr=pco->dx2v(j-1)/dxfl;

  for (int i=il; i<=iu; ++i){
#pragma simd
    for(int n=0; n<pmy_rad->n_fre_ang; ++n){
      
      dqc = (q(k,j,i,  n) - q(k,j-1,i,n))*dx_jm1i;
      
      if(vel(i,n) > 0.0){
        dql = (q(k,j-1,i,n) - q(k,j-2,i,n))*dx_jm2i;
        dq2 = dql*dqc;
        Real ql = q(k,j-1,i,n);
      // compute ql_(j-1/2) using Mignone 2014's modified van-Leer limiter
        if(dq2>0.0) {
          ql += dxfr*dq2*(cfl*dql+cbl*dqc)/(dql*dql+(cfl+cbl-2.0)*dq2+dqc*dqc);
        }
        flx(i,n) = ql * vel(i,n);
      }
      if(vel(i,n) <=0.0){
        dqr = (q(k,j+1,i,n) - q(k,j,i,  n))*dx_ji;
        dq2=dqc*dqr;
        Real qr = q(k,j,i,n);
        if(dq2>0.0) {
          qr -= dxfl*dq2*(cfr*dqc+cbr*dqr)/(dqc*dqc+(cfr+cbr-2.0)*dq2+dqr*dqr);
        }
       flx(i,n) = qr * vel(i,n);
      }
    }
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn RadIntegrator::SecondOrderFluxX3()
//  \brief 

void RadIntegrator::SecondOrderFluxX3(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &q, const AthenaArray<Real> &vel,
  AthenaArray<Real> &flx)
{
  Coordinates *pco = pmy_rad->pmy_block->pcoord;
  Real dql, dqr, dqc, dq2;
  
  Real dx_km2i = 1.0/pco->dx3v(k-2);
  Real dx_km1i = 1.0/pco->dx3v(k-1);
  Real dx_ki = 1.0/pco->dx3v(k);
  Real dxfr=pco->x3f(k)-pco->x3v(k-1);
  Real dxfl=pco->x3v(k)-pco->x3f(k);
  Real cfl=pco->dx3v(k-1)/dxfr;
  Real cbl=pco->dx3v(k-2)/(pco->x3v(k-1)-pco->x3f(k-1));
  Real cfr=pco->dx3v(k)/(pco->x3f(k+1)-pco->x3v(k));
  Real cbr=pco->dx3v(k-1)/dxfl;

  for (int i=il; i<=iu; ++i){
#pragma simd
    for(int n=0; n<pmy_rad->n_fre_ang; ++n){
      
      dqc = (q(k,j,i,  n) - q(k-1,j,i,n))*dx_km1i;
      
      if(vel(i,n) > 0.0){
        dql = (q(k-1,j,i,n) - q(k-2,j,i,n))*dx_km2i;
        dq2 = dql*dqc;
        Real ql = q(k-1,j,i,n);
      // compute ql_(k-1/2) using Mignone 2014's modified van-Leer limiter
        if(dq2>0.0) {
          ql += dxfr*dq2*(cfl*dql+cbl*dqc)/(dql*dql+(cfl+cbl-2.0)*dq2+dqc*dqc);
        }
        flx(i,n) = ql * vel(i,n);
      }
      if(vel(i,n) <=0.0){
        dqr = (q(k+1,j,i,n) - q(k,j,i,  n))*dx_ki;
        dq2=dqc*dqr;
        Real qr = q(k,j,i,n);
        if(dq2>0.0) {
          qr -= dxfl*dq2*(cfr*dqc+cbr*dqr)/(dqc*dqc+(cfr+cbr-2.0)*dq2+dqr*dqr);
        }
       flx(i,n) = qr * vel(i,n);
      }
    }
  }

  return;
}

