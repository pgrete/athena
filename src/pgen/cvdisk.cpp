///======================================================================================
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
#include <cstdlib>    // srand

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
#include "../utils/utils.hpp"

// The global parameters
static Real kappaes = 9.38e3;
static Real kappaffp = 8.51705e9;
static Real kappaffr = 2.30191e8;
static Real rho0 = 0.298994; // Density normalization of torus
static Real inib0 = 0.0; // initial magnetic field strength
static Real tfloor; // temperature floor
static Real rhofloor; // density floor
static int bconf = 0; // bconf=1: pure B_phi
                      // bconf=0: vector potential proportional to density
                      // bconf=2: two loops
static int iniflag=8;


static Real gm1 = 1.02737e5;
static Real gm2 = 1.02737e4;
static Real qm = gm2/gm1;
static Real omega0 = 1.79925;
static Real rm2 = 32.6823;

//the initial thickness of the injected stream
static Real wheight=0.248556;
//the initial injected radial velocity
static Real vr_l1=-0.447214;
static Real vphi_l1=-vr_l1*0.40035;
static Real t_l1=0.2;
static Real beta=0.0;
static Real amp=1.e-3;


static Real lunit;
static Real rhounit;
static Real tunit;


// the opacity table

static AthenaArray<Real> opacitytable;
static AthenaArray<Real> logttable;
static AthenaArray<Real> logrhottable;


//======================================================================================
/*! \file globaldisk.cpp
 *  \brief global accretion disk problem with radiation
 *
 *====================================================================================*/

void Inflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
               int is, int ie, int js, int je, int ks, int ke);


void Outflow_X2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
               int is, int ie, int js, int je, int ks, int ke);

void Outflow_rad_X2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                     int is, int ie, int js, int je, int ks, int ke);

void DiskOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);


void LoadRadVariable(MeshBlock *pmb);

void TidalPotential(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, 
  AthenaArray<Real> &bcc, AthenaArray<Real> &cons);




void Mesh::InitUserMeshData(ParameterInput *pin)
{
  
    // Enroll boundary functions

  EnrollUserBoundaryFunction(OUTER_X1, Outflow_X2);
  
  tfloor = pin->GetOrAddReal("radiation", "tfloor", 0.01);
  rhofloor = pin->GetOrAddReal("hydro", "dfloor", 1.e-7);
  rhounit = pin->GetOrAddReal("radiation", "rhounit", 1.e-4);
  tunit = pin->GetOrAddReal("radiation", "Tunit", 5.e4);
  lunit = pin->GetOrAddReal("radiation", "lunit", 4.69e8);
  
  EnrollUserSourceTermFunction(TidalPotential);
  

  if(RADIATION_ENABLED){
  
    EnrollUserRadBoundaryFunction(OUTER_X1, Outflow_rad_X2);
  
    // the opacity table
    opacitytable.NewAthenaArray(212,46);
    logttable.NewAthenaArray(212);
    logrhottable.NewAthenaArray(46);
    
    // read in the opacity table
    FILE *fkappa, *flogt, *flogrhot;
    if ( (fkappa=fopen("./aveopacity.txt","r"))==NULL )
    {
      printf("Open input file error");
      return;
    }
  
    if ( (flogt=fopen("./logT.txt","r"))==NULL )
    {
      printf("Open input file error");
      return;
    }
  
    if ( (flogrhot=fopen("./logRhoT.txt","r"))==NULL )
    {
      printf("Open input file error");
      return;
    }

    int i, j;
    for(j=0; j<212; j++){
      for(i=0; i<46; i++){
          fscanf(fkappa,"%lf",&(opacitytable(j,i)));
      }
    }
  
    for(i=0; i<46; i++){
      fscanf(flogrhot,"%lf",&(logrhottable(i)));
    }
  
    for(i=0; i<212; i++){
      fscanf(flogt,"%lf",&(logttable(i)));
    }
  
    fclose(fkappa);
    fclose(flogt);
    fclose(flogrhot);
  
  }
  

  return;
}



//======================================================================================
//! \fn void Mesh::TerminateUserMeshProperties(void)
//  \brief Clean up the Mesh properties
//======================================================================================
void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  // free memory
  if(RADIATION_ENABLED){
  
    opacitytable.DeleteAthenaArray();
    logttable.DeleteAthenaArray();
    logrhottable.DeleteAthenaArray();
  
  }
  
  return;
}

void MeshBlock::InitUserMeshBlockProperties(ParameterInput *pin)
{
  
  
  if(RADIATION_ENABLED){
    
      prad->EnrollOpacityFunction(DiskOpacity);
    
    

      if(NRADFOV > 0)
        prad->EnrollInternalVariableFunction(LoadRadVariable);

      prad->set_source_flag = 0;


  }else{

  }
  

  return;
}


void MeshBlock::UserWorkInLoop(void)
{
  if(RADIATION_ENABLED){
  
    if(prad->set_source_flag > 0)
       prad->set_source_flag--;
    
    int il=is, iu=ie, jl=js, ju=je, kl=ks, ku=ke;
    il -= NGHOST;
    iu += NGHOST;
    if(ju>jl){
       jl -= NGHOST;
       ju += NGHOST;
    }
    if(ku>kl){
      kl -= NGHOST;
      ku += NGHOST;
    }
    Real gammma1 = phydro->peos->GetGamma() - 1.0;
    
     for (int k=kl; k<=ku; ++k){
      for (int j=jl; j<=ju; ++j){
       for (int i=il; i<=iu; ++i){
         
          Real& vx=phydro->w(IVX,k,j,i);
          Real& vy=phydro->w(IVY,k,j,i);
          Real& vz=phydro->w(IVZ,k,j,i);
         
          Real& rho=phydro->w(IDN,k,j,i);

/*         if(iniflag){
              if(pcoord->x1v(i) > 50.0 && rho < 1.e-6){
                  phydro->w(IDN,k,j,i) = 1.e-7;
                  phydro->u(IDN,k,j,i) = 1.e-7;

              }

         }
*/
         
          Real vel = sqrt(vx*vx+vy*vy+vz*vz);

          if(vel > prad->vmax * prad->crat){
            Real ratio = prad->vmax * prad->crat / vel;
            vx *= ratio;
            vy *= ratio;
            vz *= ratio;
            
            phydro->u(IM1,k,j,i) = rho*vx;
            phydro->u(IM2,k,j,i) = rho*vy;
            phydro->u(IM3,k,j,i) = rho*vz;

            Real ke = 0.5 * rho * (vx*vx+vy*vy+vz*vz);
            
            Real pb=0.0;
            if(MAGNETIC_FIELDS_ENABLED){
               pb = 0.5*(SQR(pfield->bcc(IB1,k,j,i))+SQR(pfield->bcc(IB2,k,j,i))
                     +SQR(pfield->bcc(IB3,k,j,i)));
            }
            
            Real  eint = phydro->w(IEN,k,j,i)/gammma1;
            
            phydro->u(IEN,k,j,i) = eint + ke + pb;

          }
  
      }}}
    }else{
      int il=is, iu=ie, jl=js, ju=je, kl=ks, ku=ke;
      il -= NGHOST;
      iu += NGHOST;
      if(ju>jl){
        jl -= NGHOST;
        ju += NGHOST;
      }
      if(ku>kl){
        kl -= NGHOST;
        ku += NGHOST;
      }
      Real gamma1 = phydro->peos->GetGamma() - 1.0;
    
      for (int k=kl; k<=ku; ++k){
       for (int j=jl; j<=ju; ++j){
        for (int i=il; i<=iu; ++i){
         
          Real& vx=phydro->w(IVX,k,j,i);
          Real& vy=phydro->w(IVY,k,j,i);
          Real& vz=phydro->w(IVZ,k,j,i);
         
          Real& rho=phydro->w(IDN,k,j,i);
          
          Real tgas=phydro->w(IEN,k,j,i)/rho;
/*
          if(tgas > 10.0 && rho < 1.e-5){
            phydro->w(IEN,k,j,i) = rho * tfloor;
            Real eint = rho * tfloor/gamma1;
            Real ke = 0.5 * rho * (vx*vx+vy*vy+vz*vz);
            Real pb=0.0;
            if(MAGNETIC_FIELDS_ENABLED){
               pb = 0.5*(SQR(pfield->bcc(IB1,k,j,i))+SQR(pfield->bcc(IB2,k,j,i))
                     +SQR(pfield->bcc(IB3,k,j,i)));
            }
            
            phydro->u(IEN,k,j,i) = eint + ke + pb;
          
          }
*/     
        }
       }
      }
    
    }
//    if(iniflag > 0) iniflag--;
  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  
  Real gamma = phydro->peos->GetGamma();
  Real gm = pin->GetReal("problem", "GM");

  //initialize random number
  

  AthenaArray<Real> ir_cm;
  Real *ir_lab;
  
  Real crat, prat;
  if(RADIATION_ENABLED){
    ir_cm.NewAthenaArray(prad->n_fre_ang);
    crat = prad->crat;
    prat = prad->prat;
  }else{
    crat = 17178.9;
    prat = 0.0;
  }
  


  
  int kl=ks, ku=ke;
  if(ku > kl){
    ku += NGHOST;
    kl -= NGHOST;
  }
  int jl=js, ju=je;
  if(ju > jl){
    ju += NGHOST;
    jl -= NGHOST;
  }
  int il = is-NGHOST, iu=ie+NGHOST;

  
  // Initialize hydro variable
  for(int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=il; i<=iu; ++i) {

        // Initialize the hydro quantity
        phydro->u(IDN,k,j,i) = rhofloor;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = 0.1*rhofloor/(gamma-1.0);
        phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        
        // initialize radiation quantity
        if(RADIATION_ENABLED){
          for(int n=0; n<prad->n_fre_ang; ++n)
             prad->ir(k,j,i,n) = 1.e-4;
        
        }// End Rad
        
      }// i
    }// j
  }// k

  // Opacity will be set during initialization

  if(RADIATION_ENABLED)
    ir_cm.DeleteAthenaArray();
  
// initialize interface B

  if (MAGNETIC_FIELDS_ENABLED) {

    for (int k=ks; k<=ke; ++k) {
      // reset loop limits for polar boundary
      jl=js; ju=je+1;
      if (block_bcs[INNER_X2] == POLAR_BNDRY) jl=js+1;
      if (block_bcs[OUTER_X2] == POLAR_BNDRY) ju=je;
      for (int j=jl; j<=ju; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfield->b.x2f(k,j,i) = 0.0;
        }
      }
    }

    for (int k=ks; k<=ke+1; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfield->b.x3f(k,j,i) = 0.0;
        }
      }
    }

    if (block_size.nx2 > 1) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie+1; ++i) {
            pfield->b.x1f(k,j,i) = 0.0;
          }
        }
      }
  
  
  }


      
     // Update total energy with mangefew
    if(RADIATION_ENABLED){
     // Get cell-centered magnetic field
     pfield->CalculateCellCenteredField(pfield->b,pfield->bcc,pcoord,is,ie,js,je,ks,ke);
    
    
      for(int k=ks; k<=ke; ++k){
      for(int j=js; j<=je; ++j){
      for(int i=is; i<=ie; ++i){
          phydro->u(IEN,k,j,i) +=
            0.5*(SQR((pfield->bcc(IB1,k,j,i)))
               + SQR((pfield->bcc(IB2,k,j,i)))
               + SQR((pfield->bcc(IB3,k,j,i))));
        
         }
      }}
      
    }

 
  }// End MHD

  
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
   
    Real logt = log10(gast * tunit);
    Real logrhot = log10(rho* rhounit) - 3.0* logt + 18.0;
    
    
    int nrhot1 = 0;
    int nrhot2 = 0;
    
    while((logrhot > logrhottable(nrhot2)) && (nrhot2 < 45)){
      nrhot1 = nrhot2;
      nrhot2++;
    }
  
  /* The data point should between NrhoT1 and NrhoT2 */
    int nt1 = 0;
    int nt2 = 0;
    while((logt > logttable(nt2)) && (nt2 < 211)){
      nt1 = nt2;
      nt2++;
    }
    
    Real kappa;
    Real kappas = 0.2 * (1.0 + 0.6);
    Real sqrtgasT = sqrt(gast);
    Real kappff = (4908.12) * rho /(sqrtgasT * gast * gast* gast);
    

    Real kappa_t1_rho1=opacitytable(nt1,nrhot1);
    Real kappa_t1_rho2=opacitytable(nt1,nrhot2);
    Real kappa_t2_rho1=opacitytable(nt2,nrhot1);
    Real kappa_t2_rho2=opacitytable(nt2,nrhot2);

    Real rho_1 = logrhottable(nrhot1);
    Real rho_2 = logrhottable(nrhot2);
    Real t_1 = logttable(nt1);
    Real t_2 = logttable(nt2);

    
    if(nrhot1 == nrhot2){
      if(nt1 == nt2){
        kappa = kappa_t1_rho1;
      }else{
        kappa = kappa_t1_rho1 + (kappa_t2_rho1 - kappa_t1_rho1) *
                                (logt - t_1)/(t_2 - t_1);
      }/* end same T*/
    }else{
      if(nt1 == nt2){
        kappa = kappa_t1_rho1 + (kappa_t1_rho2 - kappa_t1_rho1) *
                                (logrhot - rho_1)/(rho_2 - rho_1);
      }else{
        kappa = kappa_t1_rho1 * (t_2 - logt) * (rho_2 - logrhot)/
                                ((t_2 - t_1) * (rho_2 - rho_1))
              + kappa_t2_rho1 * (logt - t_1) * (rho_2 - logrhot)/
                                ((t_2 - t_1) * (rho_2 - rho_1))
              + kappa_t1_rho2 * (t_2 - logt) * (logrhot - rho_1)/
                                ((t_2 - t_1) * (rho_2 - rho_1))
              + kappa_t2_rho2 * (logt - t_1) * (logrhot - rho_1)/
                                ((t_2 - t_1) * (rho_2 - rho_1));
      }
    }/* end same rhoT */

 
    if(kappa < kappas){
     if(gast > 1.0 && kappff < kappa){
       kappas = kappa - kappff;
       kappa = kappff;
     }else{
       kappas = 0.0;
     }
    }else{
      kappa -= kappas;
    }
    
    
    prad->sigma_s(k,j,i,ifr) = kappas * rho * rhounit * lunit;
    prad->sigma_a(k,j,i,ifr) = kappa * rho * rhounit * lunit;
    prad->sigma_ae(k,j,i,ifr) = prad->sigma_a(k,j,i,ifr);
  }
  }}}

}

// This function sets boundary condition for primitive variables


void Inflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
               int is, int ie, int js, int je, int ks, int ke)
{

  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=NGHOST; ++i) {
          a(IDN,k,j,is-i) = a(IDN,k,j,is);
          a(IVX,k,j,is-i) = std::min(a(IVX,k,j,is),0.0);
          a(IVY,k,j,is-i) = a(IVY,k,j,is);
          a(IVZ,k,j,is-i) = a(IVZ,k,j,is);
          a(IEN,k,j,is-i) = a(IEN,k,j,is);
        
      }
    }
  }
   // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){
#pragma simd
      for(int i=1; i<=NGHOST; ++i){
        b.x1f(k,j,is-i) = b.x1f(k,j,is);
      }
    }}

    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je+1; ++j){
#pragma simd
      for(int i=1; i<=NGHOST; ++i){
        b.x2f(k,j,is-i) = b.x2f(k,j,is);
      }
    }}

    for(int k=ks; k<=ke+1; ++k){
    for(int j=js; j<=je; ++j){
#pragma simd
      for(int i=1; i<=NGHOST; ++i){
        b.x3f(k,j,is-i) = b.x3f(k,j,is);
      }
    }}
    
  }


  return;
}

void Outflow_X2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
               int is, int ie, int js, int je, int ks, int ke)
{
  //initialize random number
  std::srand(pmb->gid);

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=NGHOST; ++i) {
          Real radius=pmb->pcoord->x1v(ie+i);
          Real theta=pmb->pcoord->x2v(j);
          Real phi=pmb->pcoord->x3v(k);
          Real xpos=radius*sin(theta)*cos(phi);
          Real ypos=radius*sin(theta)*sin(phi);
          Real zpos=radius*cos(theta);
          Real dis=xpos*xpos+(ypos-radius)*(ypos-radius)+zpos*zpos;
          Real rhostream=rho0*exp(-dis/(wheight*wheight));

          if(rhostream < rhofloor){
            a(IDN,k,j,ie+i) = a(IDN,k,j,ie);
            a(IVX,k,j,ie+i) = std::max(a(IVX,k,j,ie),0.0);
            a(IVY,k,j,ie+i) = a(IVY,k,j,ie);
            a(IVZ,k,j,ie+i) = a(IVZ,k,j,ie);
            a(IEN,k,j,ie+i) = a(IEN,k,j,ie);
          }else{
            a(IDN,k,j,ie+i) = rhostream;
            a(IVX,k,j,ie+i) = vr_l1;
            a(IVY,k,j,ie+i) = 0.0;
            a(IVZ,k,j,ie+i) = vphi_l1;
            a(IEN,k,j,ie+i) = t_l1 * rhostream;
                
             // add random perturbation
            if(pmb->pmy_mesh->time < 1.e-12){
                a(IDN,k,j,ie+i) *= (1.0 + amp *
                    ((double)rand()/(double)RAND_MAX-0.5));
                
            }
                
          }
      }
    }
  }
   // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){
#pragma simd
      for(int i=1; i<=NGHOST; ++i){
        b.x1f(k,j,ie+i+1) = b.x1f(k,j,ie+1);
      }
    }}

    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je+1; ++j){
#pragma simd
      for(int i=1; i<=NGHOST; ++i){
        b.x2f(k,j,ie+i) = b.x2f(k,j,ie);
      }
    }}

    for(int k=ks; k<=ke+1; ++k){
    for(int j=js; j<=je; ++j){
#pragma simd
      for(int i=1; i<=NGHOST; ++i){
          Real radius=pmb->pcoord->x1v(ie+i);
          Real theta=pmb->pcoord->x2v(j);
          Real phi=pmb->pcoord->x3v(k);
          Real xpos=radius*sin(theta)*cos(phi);
          Real ypos=radius*sin(theta)*sin(phi);
          Real zpos=radius*cos(theta);
          Real dis=xpos*xpos+(ypos-radius)*(ypos-radius)+zpos*zpos;
          Real rhostream=rho0*exp(-dis/(wheight*wheight));
          if(rhostream < rhofloor){
            b.x3f(k,j,ie+i) = b.x3f(k,j,ie);
          }else{
             b.x3f(k,j,ie+i)=0.0*sqrt(2.0*rho0*t_l1/beta);
          }
      }
    }}
    
  }

  return;
}


void Outflow_rad_X2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
                     int is, int ie, int js, int je, int ks, int ke)
{
  Radiation *prad=pmb->prad;

  AthenaArray<Real> ir_cm;
  Real *ir_lab;
  
  ir_cm.NewAthenaArray(prad->n_fre_ang);
  
  Hydro *phydro = pmb->phydro;
  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=NGHOST; ++i) {
         Real radius=pmb->pcoord->x1v(ie+i);
         Real theta=pmb->pcoord->x2v(j);
         Real phi=pmb->pcoord->x3v(k);
         Real xpos=radius*sin(theta)*cos(phi);
         Real ypos=radius*sin(theta)*sin(phi);
         Real zpos=radius*cos(theta);
         Real dis=xpos*xpos+(ypos-radius)*(ypos-radius)+zpos*zpos;
         Real rhostream=rho0*exp(-dis/(wheight*wheight));

         if(rhostream > rhofloor){
           
           for(int n=0; n<prad->n_fre_ang; ++n)
             ir_cm(n) = t_l1*t_l1*t_l1*t_l1;
           
           Real vr=phydro->w(IVX,k,j,ie+i);
           Real vtheta=phydro->w(IVY,k,j,ie+i);
           Real vphi=phydro->w(IVZ,k,j,ie+i);
           
           Real *mux = &(prad->mu(0,k,j,ie+i,0));
           Real *muy = &(prad->mu(1,k,j,ie+i,0));
           Real *muz = &(prad->mu(2,k,j,ie+i,0));

           ir_lab = &(prad->ir(k,j,ie+i,0));
          
           prad->pradintegrator->ComToLab(vr,vtheta,vphi,mux,muy,muz,ir_cm,ir_lab);
              
         }else{
           for(int ifr=0; ifr<prad->nfreq; ++ifr){
              for(int n=0; n<prad->nang; ++n){
                Real miuz = prad->mu(0,k,j,ie+i,ifr*prad->nang+n);
                if(miuz > 0.0){
                  prad->ir(k,j,ie+i,ifr*prad->nang+n)
                                = prad->ir(k,j,ie+i-1,ifr*prad->nang+n);
                }else{
                  prad->ir(k,j,ie+i,ifr*prad->nang+n) = 0.0;
                }
         
            }
         
           }
         }
        
      }//i
    }//j
  }//k
  

  return;
}




void LoadRadVariable(MeshBlock *pmb)
{

  int il=pmb->is, iu=pmb->ie;
  int jl=pmb->js, ju=pmb->je;
  int kl=pmb->ks, ku=pmb->ke;
  
  for(int k=kl; k<=ku; ++k)
    for(int j=jl; j<=ju; ++j)
      for(int i=il; i<=iu; ++i){
          pmb->prad->rad_ifov(0,k,j,i) =  pmb->prad->mu(0,k,j,i,0);
          pmb->prad->rad_ifov(1,k,j,i) =  pmb->prad->mu(1,k,j,i,0);
          pmb->prad->rad_ifov(2,k,j,i) =  pmb->prad->mu(2,k,j,i,0);
        
          pmb->prad->rad_ifov(3,k,j,i) =  pmb->prad->mu(0,k,j,i,3);
          pmb->prad->rad_ifov(4,k,j,i) =  pmb->prad->mu(1,k,j,i,3);
          pmb->prad->rad_ifov(5,k,j,i) =  pmb->prad->mu(2,k,j,i,3);
        
      }
  
  
  
  return;
}


Real grav_pot(const Real radius, const Real theta, const Real phi)
{
  // the companion is located at \theta=90, phi=0, r=rm2
  //x=rm2, y=0, z=0
  // current point r\sin\theta \cosphi, r\sin\theta\sin phi, r\sin\theta
  Real dist_r2=sqrt(radius*radius+rm2*rm2-2.0*radius*rm2*sin(theta)*cos(phi));
  
  Real potphi=-gm1/radius-gm2/dist_r2-0.5*omega0*omega0*radius*radius
        *sin(theta)*sin(theta)+gm2*radius*sin(theta)*cos(phi)/(rm2*rm2);
  return potphi;
}


void TidalPotential(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, 
  AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{

// add the effective tidal potential in the co-rotating frame
// -GM1/r-GM2/(r-R2)+GM2(r\cos\theta)/R2^2
  AthenaArray<Real> &x1flux=pmb->phydro->flux[x1face];
  AthenaArray<Real> &x2flux=pmb->phydro->flux[x2face];
  AthenaArray<Real> &x3flux=pmb->phydro->flux[x3face];
  
  for(int k=pmb->ks; k<=pmb->ke; ++k){
    for(int j=pmb->js; j<=pmb->je; ++j){
      for(int i=pmb->is; i<=pmb->ie; ++i){
        Real rho = prim(IDN,k,j,i);
        Real rcen = pmb->pcoord->x1v(i);
        Real rleft = pmb->pcoord->x1f(i);
        Real rright = pmb->pcoord->x1f(i+1);
        
        Real thetacen = pmb->pcoord->x2v(j);
        Real thetaleft = pmb->pcoord->x2f(j);
        Real thetaright = pmb->pcoord->x2f(j+1);
        
        Real phicen = pmb->pcoord->x3v(k);
        Real phileft = pmb->pcoord->x3f(k);
        Real phiright = pmb->pcoord->x3f(k+1);
        
        Real vol=pmb->pcoord->GetCellVolume(k,j,i);
        Real phic = grav_pot(rcen,thetacen,phicen);
        
        // radial direction
        
        Real phil = grav_pot(rleft,thetacen,phicen);
        Real phir = grav_pot(rright,thetacen,phicen);
        
        Real areal=pmb->pcoord->GetFace1Area(k,j,i);
        Real arear=pmb->pcoord->GetFace1Area(k,j,i+1);

        Real src = - dt * rho * (phir - phil)/pmb->pcoord->dx1f(i);
        cons(IM1,k,j,i) += src;
        Real phidivrhov = (arear*x1flux(IDN,k,j,i+1) -
                           areal*x1flux(IDN,k,j,i))*phic/vol;
        Real divrhovphi = (arear*x1flux(IDN,k,j,i+1)*phir -
                           areal*x1flux(IDN,k,j,i)*phil)/vol;
        cons(IEN,k,j,i) += (dt*(phidivrhov - divrhovphi));
        
        //theta direction

        phil = grav_pot(rcen,thetaleft,phicen);
        phir = grav_pot(rcen,thetaright,phicen);
        
        areal=0.5*(rright*rright-rleft*rleft)*fabs(sin(thetaleft))*
                   pmb->pcoord->dx3f(k);
        arear=0.5*(rright*rright-rleft*rleft)*fabs(sin(thetaright))*
                   pmb->pcoord->dx3f(k);
        
        src = - dt * rho * (phir - phil)/(rcen*pmb->pcoord->dx2f(j));
        cons(IM2,k,j,i) += src;
        phidivrhov = (arear*x2flux(IDN,k,j+1,i) -
                           areal*x2flux(IDN,k,j,i))*phic/vol;
        divrhovphi = (arear*x2flux(IDN,k,j+1,i)*phir -
                           areal*x2flux(IDN,k,j,i)*phil)/vol;
        cons(IEN,k,j,i) += (dt*(phidivrhov - divrhovphi));
        
        //phi direction
        
        phil = grav_pot(rcen,thetacen,phileft);
        phir = grav_pot(rcen,thetacen,phiright);
        
        areal=0.5*(rright*rright-rleft*rleft)*pmb->pcoord->dx2f(j);
        arear=areal;
        
        src = - dt * rho * (phir - phil)/(rcen*fabs(sin(thetacen))*
                                pmb->pcoord->dx3f(k));
        cons(IM3,k,j,i) += src;
        phidivrhov = (arear*x3flux(IDN,k+1,j,i) -
                           areal*x3flux(IDN,k,j,i))*phic/vol;
        divrhovphi = (arear*x3flux(IDN,k+1,j,i)*phir -
                           areal*x3flux(IDN,k,j,i)*phil)/vol;
        cons(IEN,k,j,i) += (dt*(phidivrhov - divrhovphi));
        
        // Add the coriolis force
       //dM/dt=-2\rho \Omega_0\times V
       // Omega_0=(\Omega_0\cos\theta,-\Omega_0\sin\theta,0)
       // because we use semi-implicit method, we need velocity
       // from conservative quantities
          rho = cons(IDN,k,j,i);
        Real vr=cons(IVX,k,j,i)/rho;
        Real vtheta=cons(IVY,k,j,i)/rho;
        Real vphi=cons(IVZ,k,j,i)/rho;
        Real dtomega = dt*omega0;
        Real sintheta=sin(thetacen);
        Real costheta=cos(thetacen);
        
        
        Real vphinew = -2.0 * sintheta*vr - 2.0*costheta*vtheta-(dtomega-1.0/dtomega)*vphi;
        vphinew /= (dtomega+1.0/dtomega);
        
        Real vrnew = dtomega * sintheta*vphinew + vr + dtomega*sintheta*vphi;
        Real vthetanew = dtomega * costheta*vphinew + vtheta + dtomega*costheta*vphi;
        
        cons(IM1,k,j,i) = vrnew * rho;
        
        cons(IM2,k,j,i) = vthetanew * rho;
        
        cons(IM3,k,j,i) = vphinew * rho;
        
      }
    }
  }

}




