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
#include "../../../mesh.hpp"
#include "../../../coordinates/coordinates.hpp"

// this class header
#include "../rad_integrators.hpp"

void LRStatePPM(Real imu[5], Real *ileft);

//--------------------------------------------------------------------------------------
//! \fn RadIntegrator::ThirdOrderFluxX1()
//  \brief 

void RadIntegrator::ThirdOrderFluxX1(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &q, const AthenaArray<Real> &vel,
  AthenaArray<Real> &flx)
{
  if(NGHOST < 3){
    std::stringstream msg;
    msg << "### FATAL ERROR in RadIntegrator::ThirdOrderFluxX1" << std::endl
        << "Third order reconstruction requires 3 ghost zones, but NGHOST="
        << NGHOST << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  Real irt[5], ileft;

  for (int i=il; i<=iu; ++i){
    for(int n=0; n<pmy_rad->n_fre_ang; ++n){
      if(vel(i,n) > 0.0){
        irt[0] = q(k,j,i-3,n);
        irt[1] = q(k,j,i-2,n);
        irt[2] = q(k,j,i-1,n);
        irt[3] = q(k,j,i  ,n);
        irt[4] = q(k,j,i+1,n);
      }
      if(vel(i,n) <=0.0){
        irt[0] = q(k,j,i+2,n);
        irt[1] = q(k,j,i+1,n);
        irt[2] = q(k,j,i,  n);
        irt[3] = q(k,j,i-1,n);
        irt[4] = q(k,j,i-2,n);
      }
      LRStatePPM(irt, &ileft);
      flx(i,n) = ileft * vel(i,n);
    }
  }


  return;
}

//--------------------------------------------------------------------------------------
//! \fn RadIntegrator::ThirdOrderFluxX2()
//  \brief 

void RadIntegrator::ThirdOrderFluxX2(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &q, const AthenaArray<Real> &vel,
  AthenaArray<Real> &flx)
{
  if(NGHOST < 3){
    std::stringstream msg;
    msg << "### FATAL ERROR in RadIntegrator::ThirdOrderFluxX2" << std::endl
        << "Third order reconstruction requires 3 ghost zones, but NGHOST="
        << NGHOST << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  Real irt[5], ileft;

  for (int i=il; i<=iu; ++i){
    for(int n=0; n<pmy_rad->n_fre_ang; ++n){
      if(vel(i,n) > 0.0){
        irt[0] = q(k,j-3,i,n);
        irt[1] = q(k,j-2,i,n);
        irt[2] = q(k,j-1,i,n);
        irt[3] = q(k,j,i  ,n);
        irt[4] = q(k,j+1,i,n);
      }
      if(vel(i,n) <=0.0){
        irt[0] = q(k,j+2,i,n);
        irt[1] = q(k,j+1,i,n);
        irt[2] = q(k,j,i,  n);
        irt[3] = q(k,j-1,i,n);
        irt[4] = q(k,j-2,i,n);
      }
      LRStatePPM(irt, &ileft);
      flx(i,n) = ileft * vel(i,n);
    }
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn RadIntegrator::ThirdOrderFluxX3()
//  \brief 

void RadIntegrator::ThirdOrderFluxX3(const int k, const int j,
  const int il, const int iu,
  const AthenaArray<Real> &q, const AthenaArray<Real> &vel,
  AthenaArray<Real> &flx)
{

  if(NGHOST < 3){
    std::stringstream msg;
    msg << "### FATAL ERROR in RadIntegrator::ThirdOrderFluxX2" << std::endl
        << "Third order reconstruction requires 3 ghost zones, but NGHOST="
        << NGHOST << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  Real irt[5], ileft;

  for (int i=il; i<=iu; ++i){
    for(int n=0; n<pmy_rad->n_fre_ang; ++n){
      if(vel(i,n) > 0.0){
        irt[0] = q(k-3,j,i,n);
        irt[1] = q(k-2,j,i,n);
        irt[2] = q(k-1,j,i,n);
        irt[3] = q(k,j,i  ,n);
        irt[4] = q(k+1,j,i,n);
      }
      if(vel(i,n) <=0.0){
        irt[0] = q(k+2,j,i,n);
        irt[1] = q(k+1,j,i,n);
        irt[2] = q(k,j,i,  n);
        irt[3] = q(k-1,j,i,n);
        irt[4] = q(k-2,j,i,n);
      }
      LRStatePPM(irt, &ileft);
      flx(i,n) = ileft * vel(i,n);
    }
  }

  return;
}



// calculate left and right state using third order reconstruction for radiation
// the original function lrstate_PPM
// This assumes uniform spacing, need to extend to arbitrary grids Collallo's paper
// put here for testing purpose
//Only calculate the left state, the upwind state
void LRStatePPM(Real imu[5], Real *ileft)
{
  // imu[0:4] i-3, i-2, i-1, i, i+1
  
  // first, calculate the half state
  //  i-3    | i-2   |    i -1     | i    | i + 1
  //                half[0]    half[1]
  //            dI[0]      dI[1]      dI[2]
  Real &i03 = imu[0];
  Real &i02 = imu[1];
  Real &i01 = imu[2];
  Real &i0  = imu[3];
  Real &i1  = imu[4];

  // calculate Ihalf0
  Real ihalf0 = (7.0 * (i02 + i01) - (i0 + i03)) / 12.0;
  Real dic = 3.0 * (i02 - 2.0 * ihalf0 + i01);
  Real dil = (i03 - 2.0 * i02 + i01);
  Real dir = (i02 - 2.0 * i01 + i0);
  Real dilim = 0.0;

  Real lim_slope = std::min(fabs(dil),fabs(dir));
  
  if((dic > 0.0) && (dil > 0.0) && (dir > 0.0)){
    dilim = SIGN(dic) * std::min(1.25 * lim_slope, fabs(dic));
  }

  if((dic < 0.0) && (dil < 0.0) && (dir < 0.0)){
    dilim = SIGN(dic) * std::min(1.25 * lim_slope, fabs(dic));
  }

  ihalf0 = 0.5 * ((i02 + i01) - dilim/3.0);

  // now calculate the ihalf1
  Real ihalf1 = (7.0 * (i01 + i0) - (i1 + i02)) / 12.0;
  dic = 3.0 * (i01 - 2.0 * ihalf1 + i0);
  dil = (i02 - 2.0 * i01 + i0);
  dir = (i01 - 2.0 * i0 + i1);
  dilim = 0.0;

  lim_slope = std::min(fabs(dil),fabs(dir));

  if((dic > 0.0) && (dil > 0.0) && (dir > 0.0)){
    dilim = SIGN(dic) * std::min(1.25 * lim_slope, fabs(dic));
  }

  if((dic < 0.0) && (dil < 0.0) && (dir < 0.0)){
    dilim = SIGN(dic) * std::min(1.25 * lim_slope, fabs(dic));
  }

  ihalf1 = 0.5 * ((i01 + i0) - dilim/3.0);


  // now calculate the left state
  
  Real il = ihalf0;
  Real ir = ihalf1;
  Real qa = (ir - i01) * (i01 - il);
  Real qb = (i02 - i01) * (i01 - i0);
  Real di, qc;
  if((qa <= 0.0) && (qb <= 0.0)){
    qc = 6.0 * (i01 - 0.5 * (il + ir));
    di = -2.0 * qc;
    dic = i02 - 2.0 * i01 + i0;
    dil = i03 - 2.0 * i02 + i01;
    dir = i01 - 2.0 * i0 + i1;
    dilim = 0.0;
    lim_slope = std::min(fabs(dil), fabs(dir));
    lim_slope = std::min(fabs(dic),lim_slope);

    if((dic > 0.0) && (dil > 0.0) && (dir > 0.0) && (di > 0.0)){
      dilim = SIGN(di) * std::min(1.25 * lim_slope, fabs(di));

    }

    if((dic < 0.0) && (dil < 0.0) && (dir < 0.0) && (di < 0.0)){
      dilim = SIGN(di) * std::min(1.25 * lim_slope, fabs(di));

    }

    if(di == 0.0){
      il = i01;
      ir = i01;
    }else{
      il = i01 + (il - i01) * dilim/di;
      ir = i01 + (ir - i01) * dilim/di;

    }
  }
  
  // constrain the value
  qa = (ir - i01) * (i01 - il);
  qb = ir - ir;
  qc = 6.0 * (i01 - 0.5 * (il + ir));

  if(qa <= 0.0){
    il = i01;
    ir = i01;
  }else if((qb * qc) > (qb * qb)){
    il = 3.0 * i01 - 2.0 * ir;
  } else if ((qb * qc) < -(qb * qb)){
    ir = 3.0 * i01 - 2.0 * il;
  }

  il = std::max(std::min(i01,i02),il);
  il = std::min(std::max(i01,i02),il);

  ir = std::max(std::min(i01,i0),ir);
  ir = std::min(std::max(i01,i0),ir);

  *ileft = ir;
}





