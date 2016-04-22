//======================================================================================
/* Athena++ astrophysical MHD code
 * Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
 *
 * This program is free software: you can redistribute and/or modify it under the terms
 * of the GNU General Public License (GPL) as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of GNU GPL in the file LICENSE included in the code
 * distribution.  If not see <http://www.gnu.org/licenses/>.
 *====================================================================================*/

// C++ headers
#include <iostream>   // endl
#include <fstream>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cmath>      // sqrt
#include <algorithm>  // min

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh.hpp"
#include "../parameter_input.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/eos/eos.hpp"
#include "../bvals/bvals.hpp"
#include "../hydro/srcterms/srcterms.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"
#include "../radiation/radiation.hpp"
#include "../radiation/integrators/rad_integrators.hpp"

// The global parameters
static Real kappaes = 28.702;
static Real kappaffp = 2.7297e-2;
static Real kappaffr = 7.37756e-4;
static Real tinject = 0.1;
static Real rinject = 1.0;
static Real mdotin = 1.86797e5;
static Real vinject = -132.189;
static Real inb0 = 50.0;
static Real tfloor;
static Real dfloor;
//======================================================================================
/*! \file radshock.cpp
 *  \brief radiation shock test for the radiative transfer module
 *
 *====================================================================================*/

void Inject_zo(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
               int is, int ie, int js, int je, int ks, int ke);


void Inject_zi(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
               int is, int ie, int js, int je, int ks, int ke);

void Inject_xo(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
               int is, int ie, int js, int je, int ks, int ke);


void Inject_xi(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
               int is, int ie, int js, int je, int ks, int ke);

void Inject_yo(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
               int is, int ie, int js, int je, int ks, int ke);


void Inject_yi(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
               int is, int ie, int js, int je, int ks, int ke);

void DiskOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);


void Inject_rad_zo(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                     int is, int ie, int js, int je, int ks, int ke);

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  
    // Enroll boundary functions

  EnrollUserBoundaryFunction(INNER_X3, Inject_zi);
  EnrollUserBoundaryFunction(OUTER_X3, Inject_zo);
  EnrollUserBoundaryFunction(INNER_X1, Inject_xi);
  EnrollUserBoundaryFunction(OUTER_X1, Inject_xo);
  EnrollUserBoundaryFunction(INNER_X2, Inject_yi);
  EnrollUserBoundaryFunction(OUTER_X2, Inject_yo);
  EnrollUserRadBoundaryFunction(OUTER_X3, Inject_rad_zo);
  
  tfloor = pin->GetOrAddReal("radiation", "tfloor", 0.01);
  dfloor = pin->GetOrAddReal("hydro", "dfloor", 1.e-4);
  
  

  return;
}


void MeshBlock::InitUserMeshBlockProperties(ParameterInput *pin)
{
  
  
  if(RADIATION_ENABLED){
      prad->EnrollOpacityFunction(DiskOpacity);
  }
  

  return;
}


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  
  Real gamma = phydro->peos->GetGamma();

  AthenaArray<Real> ir_cm;
  Real *ir_lab;
  
  if(RADIATION_ENABLED)
    ir_cm.NewAthenaArray(prad->n_fre_ang);

  Real rhoinject = fabs(mdotin/(PI*rinject*rinject*vinject));
  
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      Real &x2 = pcoord->x2v(j);
      for (int i=is; i<=ie; ++i) {
        Real &x1 = pcoord->x1v(i);
        Real dis=sqrt(x1*x1+x2*x2);
        Real rho=exp(-pow(dis/rinject,2.0)) * rhoinject;
        Real vel=exp(-pow(dis/rinject,2.0)) * vinject;
        rho=std::max(rho,dfloor);
        Real pg = rho * tinject;
        

        
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = rho * vel;
        if (NON_BAROTROPIC_EOS){

          phydro->u(IEN,k,j,i) = pg/(gamma-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }
        
        if(RADIATION_ENABLED){
          for(int n=0; n<prad->n_fre_ang; ++n)
             ir_cm(n) = tinject * tinject * tinject * tinject;
          
          Real *mux = &(prad->mu(0,k,j,i,0));
          Real *muy = &(prad->mu(1,k,j,i,0));
          Real *muz = &(prad->mu(2,k,j,i,0));

          
          ir_lab = &(prad->ir(k,j,i,0));
          
          prad->pradintegrator->ComToLab(0,0,vinject,mux,muy,muz,ir_cm,ir_lab);

        }
      }
    }
  }

  // Opacity will be set during initialization

  if(RADIATION_ENABLED)
    ir_cm.DeleteAthenaArray();
  
// initialize interface B

  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie+1; ++i) {
      pfield->b.x1f(k,j,i) = 0.0;
    }}}
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
    for (int i=is; i<=ie; ++i) {
      pfield->b.x2f(k,j,i) = 0.0;
    }}}
    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
    for (int i=is; i<=ie; ++i) {
      pfield->b.x3f(k,j,i) = inb0;
    }}}
    if (NON_BAROTROPIC_EOS) {
      for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        phydro->u(IEN,k,j,i) += 0.5*(SQR(inb0));
      }}}
    }
  }

  
  return;
}



void DiskOpacity(MeshBlock *pmb, AthenaArray<Real> &prim)
{
  Radiation *prad = pmb->prad;
  int il = pmb->is; int jl = pmb->js; int kl = pmb->ks;
  int iu = pmb->ie; int ju = pmb->je; int ku = pmb->ke;
  il -= NGHOST;
  iu += NGHOST;
  if(ju > jl){
    jl -= NGHOST;
    ju += NGHOST;
  }
  if(ku > kl){
    kl -= NGHOST;
    ku += NGHOST;
  }
  
  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
  for (int i=il; i<=iu; ++i) {
  for (int ifr=0; ifr<prad->nfreq; ++ifr){
    Real rho  = prim(IDN,k,j,i);
    Real gast = std::max(prim(IEN,k,j,i)/rho,tfloor);
    Real tpower= 1.0/(gast*gast*gast*sqrt(gast));

    prad->sigma_s(k,j,i,ifr) = kappaes * rho;
    prad->sigma_a(k,j,i,ifr) = kappaffr * rho * rho * tpower;
    prad->sigma_ae(k,j,i,ifr) = prad->sigma_a(k,j,i,ifr);
  }
  }}}

}

// This function sets boundary condition for primitive variables


void Inject_zi(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
               int is, int ie, int js, int je, int ks, int ke)
{

  
  for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
          a(IDN,ks-k,j,i) = a(IDN,ks+k-1,j,i);
          a(IVX,ks-k,j,i) = a(IVX,ks+k-1,j,i);
          a(IVY,ks-k,j,i) = a(IVY,ks+k-1,j,i);
          if(a(IDN,ks,j,i) > 1.01 * dfloor)
             a(IVZ,ks-k,j,i) = -a(IVZ,ks+k-1,j,i);
          else
             a(IVZ,ks-k,j,i) = 0.0;
          a(IEN,ks-k,j,i) = a(IEN,ks+k-1,j,i);
        
      }
    }
  }
   // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int k=1; k<=NGHOST; ++k){
    for(int j=js; j<=je; ++j){
#pragma simd
      for(int i=is; i<=ie+1; ++i){
        b.x1f(ks-k,j,i) = b.x1f(ks-k+1,j,i);
      }
    }}

    for(int k=1; k<=NGHOST; ++k){
    for(int j=js; j<=je+1; ++j){
#pragma simd
      for(int i=is; i<=ie; ++i){
        b.x2f(ks-k,j,i) = b.x2f(ks-k+1,j,i);
      }
    }}

    for(int k=1; k<=NGHOST; ++k){
    for(int j=js; j<=je; ++j){
#pragma simd
      for(int i=is; i<=ie; ++i){
        b.x3f(ks-k,j,i) = inb0;
      }
    }}
    
  }


  return;
}

void Inject_zo(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
               int is, int ie, int js, int je, int ks, int ke)
{

  Real rhoinject = fabs(mdotin/(PI * rinject * rinject * vinject));

  for (int k=1; k<=(NGHOST); ++k) {
    for (int j=js; j<=je; ++j) {
          Real &x2 = pco->x2v(j);
      for (int i=is; i<=ie; ++i) {
          Real &x1 = pco->x1v(i);
          Real dis=sqrt(x1*x1+x2*x2);
          Real rho=exp(-pow(dis/rinject,2.0)) * rhoinject;
          Real vel=exp(-pow(dis/rinject,2.0)) * vinject;
          rho=std::max(rho,dfloor);
        
        
          a(IDN,ke+k,j,i) = rho;
          a(IVX,ke+k,j,i) = 0.0 * a(IVX,ke,j,i);
          a(IVY,ke+k,j,i) = 0.0 * a(IVY,ke,j,i);
          a(IVZ,ke+k,j,i) = vel;
          a(IEN,ke+k,j,i) = rho * tinject;
        
      }
    }
  }
   // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int k=1; k<=NGHOST; ++k){
    for(int j=js; j<=je; ++j){
#pragma simd
      for(int i=is; i<=ie+1; ++i){
        b.x1f(ke+k,j,i) = 0.0 * b.x1f(ke+k-1,j,i);
      }
    }}

    for(int k=1; k<=NGHOST; ++k){
    for(int j=js; j<=je+1; ++j){
#pragma simd
      for(int i=is; i<=ie; ++i){
        b.x2f(ke+k,j,i) = 0.0 * b.x2f(ke+k-1,j,i);
      }
    }}

    for(int k=1; k<=NGHOST; ++k){
    for(int j=js; j<=je; ++j){
#pragma simd
      for(int i=is; i<=ie; ++i){
        b.x3f(ke+k+1,j,i) = inb0;
      }
    }}
    
  }


  return;
}


void Inject_rad_zo(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                     int is, int ie, int js, int je, int ks, int ke)
{
  Radiation *prad=pmb->prad;
  

  AthenaArray<Real> ir_cm;
  Real *ir_lab;
  
  ir_cm.NewAthenaArray(prad->n_fre_ang);
  
  // set intensity in the co-moving frame to T^4
  // And Lorentz tranform to the lab frame

  
  for (int k=1; k<=NGHOST; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        for(int n=0; n<prad->n_fre_ang; ++n){
          ir_cm(n) = tinject * tinject * tinject * tinject;
        }
        
          Real *mux = &(prad->mu(0,ke+k,j,i,0));
          Real *muy = &(prad->mu(1,ke+k,j,i,0));
          Real *muz = &(prad->mu(2,ke+k,j,i,0));
        
          ir_lab = &(a(ke+k,j,i,0));
        
          Real vel = pmb->phydro->w(IVZ,ke+k,j,i);
        
          prad->pradintegrator->ComToLab(0,0,vel,mux,muy,muz,ir_cm,ir_lab);
      }
    }
  }
  
  ir_cm.DeleteAthenaArray();

  return;
}

void Inject_xi(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                    int is, int ie, int js, int je, int ks, int ke)
{
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        a(n,k,j,is-i) = a(n,k,j,is);
        if(n==IVX) a(n,k,j,is-i) = std::min(a(n,k,j,is-i),0.0);
      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x1f(k,j,(is-i)) = b.x1f(k,j,is);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x2f(k,j,(is-i)) = b.x2f(k,j,is);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x3f(k,j,(is-i)) = b.x3f(k,j,is);
      }
    }}
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void OutflowOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions, outer x1 boundary

void Inject_xo(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                    int is, int ie, int js, int je, int ks, int ke)
{
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        a(n,k,j,ie+i) = a(n,k,j,ie);
        if(n==IVX) a(n,k,j,ie+i) = std::max(a(n,k,j,ie+i),0.0);

      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x1f(k,j,(ie+i+1)) = b.x1f(k,j,(ie+1));
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je+1; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x2f(k,j,(ie+i)) = b.x2f(k,j,ie);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
#pragma simd
      for (int i=1; i<=(NGHOST); ++i) {
        b.x3f(k,j,(ie+i)) = b.x3f(k,j,ie);
      }
    }}
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void OutflowInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions, inner x2 boundary

void Inject_yi(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                    int is, int ie, int js, int je, int ks, int ke)
{
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        a(n,k,js-j,i) = a(n,k,js,i);
        if(n==IVY) a(n,k,js-j,i) = std::min(a(n,k,js-j,i),0.0);

      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f(k,(js-j),i) = b.x1f(k,js,i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,(js-j),i) = b.x2f(k,js,i);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x3f(k,(js-j),i) = b.x3f(k,js,i);
      }
    }}
  }

  return;
}

//--------------------------------------------------------------------------------------
//! \fn void OutflowOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
//                          int is, int ie, int js, int je, int ks, int ke)
//  \brief OUTFLOW boundary conditions, outer x2 boundary

void Inject_yo(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
                    int is, int ie, int js, int je, int ks, int ke)
{
  // copy hydro variables into ghost zones
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        a(n,k,je+j,i) = a(n,k,je,i);
        if(n==IVY) a(n,k,je+j,i) = std::max(a(n,k,je+j,i),0.0);

      }
    }}
  }

  // copy face-centered magnetic fields into ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie+1; ++i) {
        b.x1f(k,(je+j  ),i) = b.x1f(k,(je  ),i);
      }
    }}

    for (int k=ks; k<=ke; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x2f(k,(je+j+1),i) = b.x2f(k,(je+1),i);
      }
    }}

    for (int k=ks; k<=ke+1; ++k) {
    for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
      for (int i=is; i<=ie; ++i) {
        b.x3f(k,(je+j  ),i) = b.x3f(k,(je  ),i);
      }
    }}
  }

  return;
}




