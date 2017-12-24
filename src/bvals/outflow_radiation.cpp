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
//! \file outflow.cpp
//  \brief implementation of outflow BCs in each dimension
//======================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "bvals.hpp"
#include "../radiation/radiation.hpp"

//--------------------------------------------------------------------------------------
//! \fn void RadOutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions for radiation, inner x1 boundary

void RadOutflowInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // copy radiation variables into ghost zones
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=(NGHOST); ++i) {
//#pragma simd
      for(int n=0; n<pmb->prad->n_fre_ang; n++){
        a(k,j,is-i,n) = a(k,j,is,n);
      }
    }}
  }



  return;
}

//--------------------------------------------------------------------------------------
//! \fn void RadOutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions for radiation, outer x1 boundary

void RadOutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // copy radiation variables into ghost zones

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=(NGHOST); ++i) {
//#pragma simd
       for(int n=0; n<pmb->prad->n_fre_ang; n++){
         a(k,j,ie+i,n) = a(k,j,ie,n);
       }
      }
    }
  }


  return;
}

//--------------------------------------------------------------------------------------
//! \fn void RadOutflowInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions for radiation, inner x2 boundary

void RadOutflowInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // copy rad variables into ghost zones
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      for (int i=is; i<=ie; ++i) {
//#pragma simd
        for(int n=0; n<pmb->prad->n_fre_ang; n++){
          a(k,js-j,i,n) = a(k,js,i,n);
        }
      }
    }
  }


  return;
}

//--------------------------------------------------------------------------------------
//! \fn void RadOutflowOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions for radiation, outer x2 boundary

void RadOutflowOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
            Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // copy rad variables into ghost zones
  for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
      for (int i=is; i<=ie; ++i) {
//#pragma simd
        for(int n=0; n<pmb->prad->n_fre_ang; n++){
           a(k,je+j,i,n) = a(k,je,i,n);
        }
      }
    }
  }


  return;
}

//--------------------------------------------------------------------------------------
//! \fn void RadOutflowInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions for radiation, inner x3 boundary

void RadOutflowInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
           Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // copy rad variables into ghost zones
  for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
//#pragma simd
       for(int n=0; n<pmb->prad->n_fre_ang; n++){
         a(ks-k,j,i,n) = a(ks,j,i,n);
       }
      }
    }
  }


  return;
}

//--------------------------------------------------------------------------------------
//! \fn void OutflowOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions for radiation, outer x3 boundary

void RadOutflowOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
           Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // copy hydro variables into ghost zones
  for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
//#pragma simd
       for(int n=0; n<pmb->prad->n_fre_ang; n++){
        a(ke+k,j,i,n) = a(ke,j,i,n);
       }
      }
    }
  }


  return;
}
