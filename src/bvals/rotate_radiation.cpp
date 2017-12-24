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

// The angular octant ( in x-y plane) is
//   1  |  0       5  |  4
//   -------      ---------
//   3  |  2       7  |  6
// in X-Z plane, it is
//   5  |  4       7  |  6
//   -------      ---------
//   1  |  0       3  |  2


// in radiatin class, n_ang is angles per octant, noct is the number of octant


// Temporary function to copy intensity
void CopyIntensity2(Real *ir, int n_ang, int direction)
{
  // here ir is only intensity for each cell and each frequency band
  if(direction == 1){
   for(int n=0; n<n_ang; ++n){
     // from 0 to 4, 4 to 5, 5 to 1, 1 to 0
     int ang1 = 0 * n_ang + n;
     int ang2 = 4 * n_ang + n;
     int ang3 = 5 * n_ang + n;
     int ang4 = 1 * n_ang + n;
     Real temp = ir[ang1];
     ir[ang1] = ir[ang4];
     ir[ang4] = ir[ang3];
     ir[ang3] = ir[ang2];
     ir[ang2] = temp;
   }
   for(int n=0; n<n_ang; ++n){
     // from 0 to 4, 4 to 5, 5 to 1, 1 to 0
     int ang1 = 2 * n_ang + n;
     int ang2 = 6 * n_ang + n;
     int ang3 = 7 * n_ang + n;
     int ang4 = 3 * n_ang + n;
     Real temp = ir[ang1];
     ir[ang1] = ir[ang4];
     ir[ang4] = ir[ang3];
     ir[ang3] = ir[ang2];
     ir[ang2] = temp;
   }
  }
  else if(direction == -1){
   for(int n=0; n<n_ang; ++n){
     // from 0 to 4, 4 to 5, 5 to 1, 1 to 0
     int ang1 = 0 * n_ang + n;
     int ang2 = 4 * n_ang + n;
     int ang3 = 5 * n_ang + n;
     int ang4 = 1 * n_ang + n;
     Real temp = ir[ang1];
     ir[ang1] = ir[ang2];
     ir[ang2] = ir[ang3];
     ir[ang3] = ir[ang4];
     ir[ang4] = temp;
   }
   for(int n=0; n<n_ang; ++n){
     // from 0 to 4, 4 to 5, 5 to 1, 1 to 0
     int ang1 = 2 * n_ang + n;
     int ang2 = 6 * n_ang + n;
     int ang3 = 7 * n_ang + n;
     int ang4 = 3 * n_ang + n;
     Real temp = ir[ang1];
     ir[ang1] = ir[ang2];
     ir[ang2] = ir[ang3];
     ir[ang3] = ir[ang4];
     ir[ang4] = temp;
   }
  }
}

void CopyIntensity3(Real *ir, int n_ang, int direction)
{
  // here ir is only intensity for each cell and each frequency band
   for(int n=0; n<n_ang; ++n){
     // switch 0 and 4, switch 5 and 1
     int ang1 = 0 * n_ang + n;
     int ang2 = 4 * n_ang + n;
     int ang3 = 5 * n_ang + n;
     int ang4 = 1 * n_ang + n;
     Real temp = ir[ang1];
     ir[ang1] = ir[ang2];
     ir[ang2] = temp;
     
     temp = ir[ang4];
     ir[ang4] = ir[ang3];
     ir[ang3] = temp;
   }
   for(int n=0; n<n_ang; ++n){
     // switch 2 and 6, switch 3 and 7
     int ang1 = 2 * n_ang + n;
     int ang2 = 6 * n_ang + n;
     int ang3 = 7 * n_ang + n;
     int ang4 = 3 * n_ang + n;
     Real temp = ir[ang1];
     ir[ang1] = ir[ang2];
     ir[ang2] = temp;
     
     temp = ir[ang4];
     ir[ang4] = ir[ang3];
     ir[ang3] = temp;
   }
}

// The angular octant ( in x-y plane) is
//   1  |  0       5  |  4
//   -------      ---------
//   3  |  2       7  |  6

void CopyIntensity4(Real *ir, int n_ang, int direction)
{
  // here ir is only intensity for each cell and each frequency band
  if(direction == 1){
   for(int n=0; n<n_ang; ++n){
     // from 0 to 1, 1 to 3, 3 to 2, 2 to 0
     int ang1 = 0 * n_ang + n;
     int ang2 = 1 * n_ang + n;
     int ang3 = 3 * n_ang + n;
     int ang4 = 2 * n_ang + n;
     Real temp = ir[ang1];
     ir[ang1] = ir[ang4];
     ir[ang4] = ir[ang3];
     ir[ang3] = ir[ang2];
     ir[ang2] = temp;
     
   }
   for(int n=0; n<n_ang; ++n){
     // from 4 to 5, 5 to 7, 7 to 6, 6 to 1
     int ang1 = 4 * n_ang + n;
     int ang2 = 5 * n_ang + n;
     int ang3 = 7 * n_ang + n;
     int ang4 = 6 * n_ang + n;
     
     Real temp = ir[ang1];
     ir[ang1] = ir[ang4];
     ir[ang4] = ir[ang3];
     ir[ang3] = ir[ang2];
     ir[ang2] = temp;
     
   }
 }else if(direction == -1){
   for(int n=0; n<n_ang; ++n){
     // from 0 to 2, 2 to 3, 3 to 1, 1 to 0
     int ang1 = 0 * n_ang + n;
     int ang2 = 2 * n_ang + n;
     int ang3 = 3 * n_ang + n;
     int ang4 = 1 * n_ang + n;
     Real temp = ir[ang1];
     ir[ang1] = ir[ang4];
     ir[ang4] = ir[ang3];
     ir[ang3] = ir[ang2];
     ir[ang2] = temp;
     
   }
   for(int n=0; n<n_ang; ++n){
     // from 4 to 6, 6 to 7, 7 to 5, 5 to 4
     int ang1 = 4 * n_ang + n;
     int ang2 = 6 * n_ang + n;
     int ang3 = 7 * n_ang + n;
     int ang4 = 5 * n_ang + n;
     
     Real temp = ir[ang1];
     ir[ang1] = ir[ang4];
     ir[ang4] = ir[ang3];
     ir[ang3] = ir[ang2];
     ir[ang2] = temp;
     
   }

 }
}

void CopyIntensity5(Real *ir, int n_ang, int direction)
{
  // here ir is only intensity for each cell and each frequency band
   for(int n=0; n<n_ang; ++n){
     // swap 0 and 2, swap 1 and 3
     int ang1 = 0 * n_ang + n;
     int ang2 = 2 * n_ang + n;
     int ang3 = 1 * n_ang + n;
     int ang4 = 3 * n_ang + n;
     Real temp = ir[ang1];
     ir[ang1] = ir[ang2];
     ir[ang2] = temp;
     
     
     temp = ir[ang3];
     ir[ang3] = ir[ang4];
     ir[ang4] = temp;
     
     
   }
   for(int n=0; n<n_ang; ++n){
     // swap 4 and 6, 5 and 7
     int ang1 = 4 * n_ang + n;
     int ang2 = 6 * n_ang + n;
     int ang3 = 5 * n_ang + n;
     int ang4 = 7 * n_ang + n;
     
     Real temp = ir[ang1];
     ir[ang1] = ir[ang2];
     ir[ang2] = temp;
     
     temp = ir[ang3];
     ir[ang3] = ir[ang4];
     ir[ang4] = temp;
     
   }
}

//--------------------------------------------------------------------------------------
//! \fn void RotateHPi_InnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief ROTATE boundary conditions for x2, inner x2 boundary by Pi/2

// anti-clockwise rotation
// This function is used after periodic copy

void RotateHPi_InnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                   int is, int ie, int js, int je, int ks, int ke)
{
  // copy radiation variables into ghost zones


  int &noct = pmb->prad->noct;
  int n_ang = pmb->prad->nang/noct; // angles per octant
  int &nfreq = pmb->prad->nfreq; // number of frequency bands

  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {
  for (int i=is; i<=ie; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    Real *ir = &(a(k,js-j,i,ifr*pmb->prad->nang));
    CopyIntensity3(ir, n_ang, 1);
  }
  }}}


  return;
}

//--------------------------------------------------------------------------------------
//! \fn void RotateHPi_OuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief  ROTATE boundary radiation conditions, outer x2 boundary

void RotateHPi_OuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                    int is, int ie, int js, int je, int ks, int ke)
{



  int &noct = pmb->prad->noct;
  int n_ang = pmb->prad->nang/noct; // angles per octant
  int &nfreq = pmb->prad->nfreq; // number of frequency bands

  for (int k=ks; k<=ke; ++k) {
  for (int j=1; j<=(NGHOST); ++j) {
  for (int i=is; i<=ie; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    Real *ir = &(a(k,je+j,i,ifr*pmb->prad->nang));
    CopyIntensity3(ir, n_ang, -1);
  }
  }}}

  return;
}


//--------------------------------------------------------------------------------------
//! \fn void RotateHPi_InnerX3((MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                     int is, int ie, int js, int je, int ks, int ke)
//  \brief ROTATE boundary conditions for x3, inner x3 boundary by Pi/2

// anti-clockwise rotation
// This function is used after periodic copy

void RotateHPi_InnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                     int is, int ie, int js, int je, int ks, int ke)
{
  // copy radiation variables into ghost zones


  int &noct = pmb->prad->noct;
  int n_ang = pmb->prad->nang/noct; // angles per octant
  int &nfreq = pmb->prad->nfreq; // number of frequency bands

  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=is; i<=ie; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    Real *ir = &(a(ks-k,j,i,ifr*pmb->prad->nang));
    CopyIntensity4(ir, n_ang, -1);
  }
  }}}


  return;
}


//--------------------------------------------------------------------------------------
//! \fn void RotateHPi_OuterX3((MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                     int is, int ie, int js, int je, int ks, int ke)
//  \brief ROTATE boundary conditions for x3, outer x3 boundary by Pi/2

// anti-clockwise rotation
// This function is used after periodic copy

void RotateHPi_OuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                     int is, int ie, int js, int je, int ks, int ke)
{
  // copy radiation variables into ghost zones


  int &noct = pmb->prad->noct;
  int n_ang = pmb->prad->nang/noct; // angles per octant
  int &nfreq = pmb->prad->nfreq; // number of frequency bands

  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=is; i<=ie; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    Real *ir = &(a(ke+k,j,i,ifr*pmb->prad->nang));
    CopyIntensity4(ir, n_ang, 1);
  }
  }}}


  return;
}


//--------------------------------------------------------------------------------------
//! \fn void RotatePi_InnerX3((MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                     int is, int ie, int js, int je, int ks, int ke)
//  \brief ROTATE boundary conditions for x3, inner x3 boundary by Pi

// anti-clockwise rotation
// This function is used after periodic copy

void RotatePi_InnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                     int is, int ie, int js, int je, int ks, int ke)
{
  // copy radiation variables into ghost zones


  int &noct = pmb->prad->noct;
  int n_ang = pmb->prad->nang/noct; // angles per octant
  int &nfreq = pmb->prad->nfreq; // number of frequency bands

  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=is; i<=ie; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    Real *ir = &(a(ks-k,j,i,ifr*pmb->prad->nang));
    CopyIntensity5(ir, n_ang, 1);
  }
  }}}


  return;
}


//--------------------------------------------------------------------------------------
//! \fn void RotateHPi_OuterX3((MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                     int is, int ie, int js, int je, int ks, int ke)
//  \brief ROTATE boundary conditions for x3, outer x3 boundary by Pi

// anti-clockwise rotation
// This function is used after periodic copy

void RotatePi_OuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                     int is, int ie, int js, int je, int ks, int ke)
{
  // copy radiation variables into ghost zones


  int &noct = pmb->prad->noct;
  int n_ang = pmb->prad->nang/noct; // angles per octant
  int &nfreq = pmb->prad->nfreq; // number of frequency bands

  for (int k=1; k<=(NGHOST); ++k) {
  for (int j=js; j<=je; ++j) {
  for (int i=is; i<=ie; ++i) {
  for (int ifr=0; ifr<nfreq; ++ifr){
    Real *ir = &(a(ke+k,j,i,ifr*pmb->prad->nang));
    CopyIntensity5(ir, n_ang, 1);
  }
  }}}


  return;
}



