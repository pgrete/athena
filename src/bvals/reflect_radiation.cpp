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
//! \file reflect.cpp
//  \brief implementation of reflecting BCs in each dimension
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
// radiation relfection means set specific intensitiy corresponding
// the value in the opposite octant

// Temporary function to copy intensity
void CopyIntensity(Real *iri, Real *iro, int li, int lo, int n_ang)
{
  // here ir is only intensity for each cell and each frequency band
//#pragma simd
  for(int n=0; n<n_ang; ++n){
    int angi = li * n_ang + n;
    int ango = lo * n_ang + n;
    iro[angi] = iri[ango];
    iro[ango] = iri[angi];
  }
}

//--------------------------------------------------------------------------------------
//! \fn void RadReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions for radiation, inner x1 boundary

void RadReflectInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
           Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  
  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx

  int &noct = pmb->prad->noct;
  int n_ang = pmb->prad->nang/noct; // angles per octant
  int &nfreq = pmb->prad->nfreq; // number of frequency bands

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=1; i<=(NGHOST); ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    Real *iri = &(a(k,j,(is+i-1),ifr*pmb->prad->nang));
    Real *iro = &(a(k,j, is-i, ifr*pmb->prad->nang));
    CopyIntensity(iri, iro, 0, 1, n_ang);

    if(noct > 2)
        CopyIntensity(iri, iro, 2, 3, n_ang);
    
    if(noct > 3){
        CopyIntensity(iri, iro, 4, 5, n_ang);
        CopyIntensity(iri, iro, 6, 7, n_ang);
    }
  
  }
  }}}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void RadReflectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions for radiation, outer x1 boundary

void RadReflectOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
          Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // copy radiation variables into ghost zones,
  // reflect rays along angles with opposite nx

  int &noct = pmb->prad->noct;
  int n_ang = pmb->prad->nang/noct; // angles per octant
  int &nfreq = pmb->prad->nfreq; // number of frequency bands

  for (int k=ks; k<=ke; ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=1; i<=(NGHOST); ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    Real *iri = &(a(k,j,(ie-i+1),ifr*pmb->prad->nang));
    Real *iro = &(a(k,j, ie+i, ifr*pmb->prad->nang));
    CopyIntensity(iri, iro, 0, 1, n_ang);

    if(noct > 2)
        CopyIntensity(iri, iro, 2, 3, n_ang);
    
    if(noct > 3){
        CopyIntensity(iri, iro, 4, 5, n_ang);
        CopyIntensity(iri, iro, 6, 7, n_ang);
    }
  
  }
  }}}

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void RadReflecInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions for x2, inner x2 boundary

void RadReflectInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
           Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // copy radiation variables into ghost zones
  // reflect rays along angles with opposite ny


  int &noct = pmb->prad->noct;
  int n_ang = pmb->prad->nang/noct; // angles per octant
  int &nfreq = pmb->prad->nfreq; // number of frequency bands

  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {
  for (int i=is; i<=ie; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    Real *iri = &(a(k,js+j-1,i,ifr*pmb->prad->nang));
    Real *iro = &(a(k,js-j,i,ifr*pmb->prad->nang));
    CopyIntensity(iri, iro, 0, 2, n_ang);
    CopyIntensity(iri, iro, 1, 3, n_ang);
    
    if(noct > 3){
        CopyIntensity(iri, iro, 4, 6, n_ang);
        CopyIntensity(iri, iro, 5, 7, n_ang);
    }
  
  }
  }}}


  return;
}

//--------------------------------------------------------------------------------------
//! \fn void RadReflectOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary radiation conditions, outer x2 boundary

void RadReflectOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
           Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{

  // copy radiation variables into ghost zones
  // reflect rays along angles with opposite ny


  int &noct = pmb->prad->noct;
  int n_ang = pmb->prad->nang/noct; // angles per octant
  int &nfreq = pmb->prad->nfreq; // number of frequency bands

  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {
  for (int i=is; i<=ie; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    Real *iri = &(a(k,je-j+1,i,ifr*pmb->prad->nang));
    Real *iro = &(a(k,je+j,i,ifr*pmb->prad->nang));
    CopyIntensity(iri, iro, 0, 2, n_ang);
    CopyIntensity(iri, iro, 1, 3, n_ang);
    
    if(noct > 3){
        CopyIntensity(iri, iro, 4, 6, n_ang);
        CopyIntensity(iri, iro, 5, 7, n_ang);
    }
  
  }
  }}}


  return;
}

//--------------------------------------------------------------------------------------
//! \fn void RadReflectInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions for radiation, inner x3 boundary

void RadReflectInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
          Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{

  // copy radiation variables into ghost zones
  // reflect rays along angles with opposite nz


  int &noct = pmb->prad->noct;
  int n_ang = pmb->prad->nang/noct; // angles per octant
  int &nfreq = pmb->prad->nfreq; // number of frequency bands

  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=is; i<=ie; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    Real *iri = &(a(ks+k-1,j,i,ifr*pmb->prad->nang));
    Real *iro = &(a(ks-k,j,i,ifr*pmb->prad->nang));
    CopyIntensity(iri, iro, 0, 4, n_ang);
    CopyIntensity(iri, iro, 1, 5, n_ang);
    CopyIntensity(iri, iro, 2, 6, n_ang);
    CopyIntensity(iri, iro, 3, 7, n_ang);
    
  
  }
  }}}


  return;
}

//--------------------------------------------------------------------------------------
//! \fn void RadReflectOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief REFLECTING boundary conditions for radiation, outer x3 boundary

void RadReflectOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
          Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  // copy radiation variables into ghost zones
  // reflect rays along angles with opposite nz


  int &noct = pmb->prad->noct;
  int n_ang = pmb->prad->nang/noct; // angles per octant
  int &nfreq = pmb->prad->nfreq; // number of frequency bands

  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=is; i<=ie; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    Real *iri = &(a(ke-k+1,j,i,ifr*pmb->prad->nang));
    Real *iro = &(a(ke+k,j,i,ifr*pmb->prad->nang));
    CopyIntensity(iri, iro, 0, 4, n_ang);
    CopyIntensity(iri, iro, 1, 5, n_ang);
    CopyIntensity(iri, iro, 2, 6, n_ang);
    CopyIntensity(iri, iro, 3, 7, n_ang);
    
  
  }
  }}}


  return;
}
