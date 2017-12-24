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
//! \file Vacuum.cpp
//  \brief implementation of VACUUM BCs in each dimension
//======================================================================================

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "bvals.hpp"
#include "../radiation/radiation.hpp"

// The angular octant is
//   1  |  0       5  |  4
//   -------      ---------
//   3  |  2       7  |  6
// in radiatin class, n_ang is angles per octant, noct is the number of octant
// Vacuum boundary copy outgoing intensity but set incoming intensity to be zero
//--------------------------------------------------------------------------------------
//! \fn void RadVacuumInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief VACUUM boundary conditions for radiation, inner x1 boundary

void RadVacuumInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
         Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  


  int &nang = pmb->prad->nang;
  int &nfreq = pmb->prad->nfreq; // number of frequency bands

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=1; i<=(NGHOST); ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    for(int n=0; n<nang; ++n){
      int ang=ifr*nang+n;
      Real& miux=pmb->prad->mu(0,k,j,is,n);
      if(miux < 0){
        a(k,j,is-i,ang) = a(k,j,is,ang);
      }else{
        a(k,j,is-i,ang) = 0.0;
      }
    }
  
  }
  }}}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void RadVacuumOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief VACUUM boundary conditions for radiation, outer x1 boundary

void RadVacuumOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
          Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  int &nang = pmb->prad->nang;
  int &nfreq = pmb->prad->nfreq; // number of frequency bands

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=1; i<=(NGHOST); ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    for(int n=0; n<nang; ++n){
      int ang=ifr*nang+n;
      Real& miux=pmb->prad->mu(0,k,j,ie,n);
      if(miux > 0){
        a(k,j,ie+i,ang) = a(k,j,ie,ang);
      }else{
        a(k,j,ie+i,ang) = 0.0;
      }
    }

  
  }
  }}}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void RadReflecInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief VACUUM boundary conditions for x2, inner x2 boundary

void RadVacuumInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
         Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  int &nang = pmb->prad->nang;
  int &nfreq = pmb->prad->nfreq; // number of frequency bands

  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {
  for (int i=is; i<=ie; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
  
    for(int n=0; n<nang; ++n){
      int ang=ifr*nang+n;
      Real& miuy=pmb->prad->mu(1,k,js,i,n);
      if(miuy < 0){
        a(k,js-j,i,ang) = a(k,js,i,ang);
      }else{
        a(k,js-j,i,ang) = 0.0;
      }
    }
  }
  }}}


  return;
}

//--------------------------------------------------------------------------------------
//! \fn void RadVacuumOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief VACUUM boundary radiation conditions, outer x2 boundary

void RadVacuumOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
          Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{

  int &nang = pmb->prad->nang;
  int &nfreq = pmb->prad->nfreq; // number of frequency bands

  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {
  for (int i=is; i<=ie; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    
    for(int n=0; n<nang; ++n){
      int ang=ifr*nang+n;
      Real& miuy=pmb->prad->mu(1,k,je,i,n);
      if(miuy > 0){
        a(k,je+j,i,ang) = a(k,je,i,ang);
      }else{
        a(k,je+j,i,ang) = 0.0;
      }
    }
  }
  }}}


  return;
}

//--------------------------------------------------------------------------------------
//! \fn void RadVacuumInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief VACUUM boundary conditions for radiation, inner x3 boundary

void RadVacuumInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
         Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{

  int &nang = pmb->prad->nang;
  int &nfreq = pmb->prad->nfreq; // number of frequency bands

  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=is; i<=ie; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
  
    for(int n=0; n<nang; ++n){
      int ang=ifr*nang+n;
      Real& miuz=pmb->prad->mu(2,ks,j,i,n);
      if(miuz < 0){
        a(ks-k,j,i,ang) = a(ks,j,i,ang);
      }else{
        a(ks-k,j,i,ang) = 0.0;
      }
    }
    
  
  }
  }}}


  return;
}

//--------------------------------------------------------------------------------------
//! \fn void RadVacuumOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief VACUUM boundary conditions for radiation, outer x3 boundary

void RadVacuumOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
         Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{

  int &nang = pmb->prad->nang;
  int &nfreq = pmb->prad->nfreq; // number of frequency bands
  
  
  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=is; i<=ie; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
  
    for(int n=0; n<nang; ++n){
      int ang=ifr*nang+n;
      Real& miuz=pmb->prad->mu(2,ke,j,i,n);
      if(miuz > 0){
        a(ke+k,j,i,ang) = a(ke,j,i,ang);
      }else{
        a(ke+k,j,i,ang) = 0.0;
      }
    }
    
    
  
  }
  }}}


  return;
}
