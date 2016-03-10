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
//! \file first_order_flux.cpp
//  \brief  piecewise constant reconstruction
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
//! \fn RadIntegrator::FirstOrderFluxX1()
//  \brief 

void RadIntegrator::FirstOrderFluxX1(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &q, const AthenaArray<Real> &vel,
  AthenaArray<Real> &flx)
{
  for (int i=il; i<=iu; ++i){
#pragma simd
    for(int n=0; n<pmy_rad->n_fre_ang; ++n){
      if(vel(i,n)>0.0){
        flx(i,n) = q(k,j,i-1,n) * vel(i,n);
      }
      if(vel(i,n) <= 0.0){
        flx(i,n) = q(k,j,i,n) * vel(i,n);
      }
    }
  }


  return;
}
//--------------------------------------------------------------------------------------
//! \fn RadIntegrator::FirstOrderFluxX2()
//  \brief 

void RadIntegrator::FirstOrderFluxX2(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &q, const AthenaArray<Real> &vel,
  AthenaArray<Real> &flx)
{
  for (int i=il; i<=iu; ++i){
#pragma simd
    for(int n=0; n<pmy_rad->n_fre_ang; ++n){
      if(vel(i,n)>0.0){
        flx(i,n) = q(k,j-1,i,n) * vel(i,n);
      }
      if(vel(i,n) <= 0.0){
        flx(i,n) = q(k,j,i,n) * vel(i,n);
      }
    }
  }


  return;
}

//--------------------------------------------------------------------------------------
//! \fn RadIntegrator::FirstOrderFluxX3()
//  \brief 

void RadIntegrator::FirstOrderFluxX3(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &q, const AthenaArray<Real> &vel,
  AthenaArray<Real> &flx)
{
  for (int i=il; i<=iu; ++i){
#pragma simd
    for(int n=0; n<pmy_rad->n_fre_ang; ++n){
      if(vel(i,n)>0.0){
        flx(i,n) = q(k-1,j,i,n) * vel(i,n);
      }
      if(vel(i,n) <= 0.0){
        flx(i,n) = q(k,j,i,n) * vel(i,n);
      }
    }
  }


  return;
}
