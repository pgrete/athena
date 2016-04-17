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
static Real kappaes = 5.04655e5;
static Real kappaffp = 1.18573e4;
static Real kappaffr = 320.469;
static Real rho0 = 10.0; // Density normalization of torus
static Real inib0 = 1.0; // initial magnetic field strength
static Real r0 = 25.0; // Center of initial torus
static Real tfloor; // temperature floor
static Real rhofloor; // density floor
static int bconf = 0; // bconf=1: pure B_phi
                      // bconf=0: vector potential proportional to density
                      // bconf=2: two loops
static Real lprofile= 0.4;
static Real vs0 = 100.0;

static Real gm;

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


void DiskOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);


void LoadRadVariable(MeshBlock *pmb);

void PseudoNewtonian( MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim,
  AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

void Mesh::InitUserMeshProperties(ParameterInput *pin)
{
  
    // Enroll boundary functions

  EnrollUserBoundaryFunction(INNER_X1, Inflow_X1);
  EnrollUserBoundaryFunction(OUTER_X1, Outflow_X2);
  
  tfloor = pin->GetOrAddReal("radiation", "tfloor", 0.01);
  rhofloor = pin->GetOrAddReal("hydro", "dfloor", 1.e-5);
  rhounit = pin->GetOrAddReal("radiation", "rhounit", 1.e-8);
  tunit = pin->GetOrAddReal("radiation", "Tunit", 2.e5);
  lunit = pin->GetOrAddReal("radiation", "lunit", 1.48428e14);
  
  EnrollUserSourceTermFunction(PseudoNewtonian);
  
  if(RADIATION_ENABLED){
  
    // the opacity table
    opacitytable.NewAthenaArray(138,37);
    logttable.NewAthenaArray(138);
    logrhottable.NewAthenaArray(37);
    
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
    for(j=0; j<138; j++){
      for(i=0; i<37; i++){
          fscanf(fkappa,"%lf",&(opacitytable(j,i)));
      }
    }
  
    for(i=0; i<37; i++){
      fscanf(flogrhot,"%lf",&(logrhottable(i)));
    }
  
    for(i=0; i<138; i++){
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
void Mesh::TerminateUserMeshProperties(ParameterInput *pin)
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
    
      gm = 0.5 * prad->crat * prad->crat;
    
    

      if(NRADFOV > 0)
        prad->EnrollInternalVariableFunction(LoadRadVariable);

  }else{
      gm = 0.5 * 5694.76 * 5694.76;
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

  //initialize random number
  std::srand(gid);
  

  AthenaArray<Real> ir_cm;
  Real *ir_lab;
  
  Real crat, prat;
  if(RADIATION_ENABLED){
    ir_cm.NewAthenaArray(prad->n_fre_ang);
    crat = prad->crat;
    prat = prad->prat;
  }else{
    crat = 5694.76;
    prat = 0.0;
  }
  

  Real l0 = sqrt(0.5 * r0) * r0 * crat/(r0 - 1.0);
  Real nindex = 1.0/(gamma-1.0);
  Real amp;
  
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
      Real &x2 = pcoord->x2v(j);
      for (int i=il; i<=iu; ++i) {
        Real &x1 = pcoord->x1v(i);
        Real angradius = x1 *sin(x2);
        Real langular = l0 * pow(angradius/r0,lprofile);
        Real vphi = langular/angradius;
        Real effphi = -crat*crat/(2.0*(x1-1.0))
                  +pow((langular/angradius),2.0)/(2.0*(1.0-lprofile));
        Real effphi0 = -crat*crat/(2.0*(r0-1.0))
                  +pow((l0/r0),2.0)/(2.0*(1.0-lprofile));
        
        Real tempphi = ((effphi - effphi0)/nindex)/(vs0 * vs0);
        
        Real rho, press;
        
        if((fabs(tempphi)<1.0) && (x1 > 8.0)){
          rho = rho0*pow(fabs(1.0-tempphi),nindex);
          rho = std::max(rho,rhofloor);
          press = rho0*vs0*vs0*pow(rho/rho0,gamma)/gamma;
          amp = 0.01;

        }else{
          rho = rhofloor;
          press = rhofloor * 1.0;
          amp = 0.0;
          vphi = 0.0;
        }
        
        
        Real temp0 = press/rho;
        Real coef1 = prat/3.0;
        Real coef2 = rho;
        Real coef3 = -press;
        
        Real gast;
        if(RADIATION_ENABLED){
          gast = Rtsafe(Tequilibrium, 0.0, temp0, 1.e-12, coef1, coef2, coef3, 0.0);
          if(gast < 1.0) gast = 1.0;
        }else{
          gast = press/rho;
        }
        
        press = gast * rho;
        
        // Add perturbation
        rho *= (1.0 + amp * ((double)rand()/(double)RAND_MAX-0.5));
        rho = std::max(rho, rhofloor);
        
        // Initialize the hydro quantity
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = vphi * rho;
        phydro->u(IEN,k,j,i) = press/(gamma-1.0);
        phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        
        // initialize radiation quantity
        if(RADIATION_ENABLED){
          for(int n=0; n<prad->n_fre_ang; ++n)
             ir_cm(n) = gast * gast * gast * gast;

          Real *mux = &(prad->mu(0,k,j,i,0));
          Real *muy = &(prad->mu(1,k,j,i,0));
          Real *muz = &(prad->mu(2,k,j,i,0));

          ir_lab = &(prad->ir(k,j,i,0));
          
          prad->pradintegrator->ComToLab(0,0,vphi,mux,muy,muz,ir_cm,ir_lab);
        
        }// End Rad
        
      }// i
    }// j
  }// k

  // Opacity will be set during initialization

  if(RADIATION_ENABLED)
    ir_cm.DeleteAthenaArray();
  
// initialize interface B

  if (MAGNETIC_FIELDS_ENABLED) {
  
      int nx1 = ie-is+1+2*NGHOST;
      int nx2 = 1;
      if(je > js) nx2 = je-js+1+2*NGHOST;
      int nx3 = 1;
      if(ke > ks) nx3 = ke-ks+1+2*NGHOST;
    
      AthenaArray<Real> baphi;

      baphi.NewAthenaArray(nx3,nx2,nx1);

      AthenaArray<Real> area, len, len_p1;
      
      area.NewAthenaArray(nx1);
      len.NewAthenaArray(nx1);
      len_p1.NewAthenaArray(nx1);
    



    // need vector potential
    if(bconf == 0 || bconf == 2){
      Real aphi;

      for(int k=kl; k<=ku; ++k){
      for(int j=jl; j<=ju; ++j){
      for(int i=il; i<=iu; ++i){
       Real &x1 = pcoord->x1v(i);
       if(x1 > 4.0){
         if(bconf == 0){
           if(phydro->u(IDN,k,j,i) > 1.01 *rhofloor){
              aphi = inib0 * phydro->u(IDN,k,j,i);
           }else{
              aphi = 0.0;
           }
         }else if(bconf == 2){
           Real rholimit = 1.0;
           if(phydro->u(IDN,k,j,i) > rholimit){
             aphi = inib0 * pow(((phydro->u(IDN,k,j,i)-rholimit)/rho0)*pow(x1,0.75)
                    ,2.0)*sin(log(x1/(0.4*r0))/0.16);
           }else{
             aphi = 0.0;
           }
         }
       }// end x1 > 4
       else{
         aphi = 0.0;
       }
       baphi(k,j,i) = aphi;
      }}}
    
    
    }// end bconf=0 and bconf=2

       
    //B=div X Phi
    // vector potential only has non-zero phi component
    // in spherical polar coordinate system
    // B_r= 1/rsintheta\partial (sintheta A_phi)/partial theta
    //      - 1/rsintheta \partial A_theta/\partial phi
    
    // B_theta=1/rsintheta\partial A_r/\partial phi - 1/r\partial (rA_phi)/\partial r
    // B_phi=1/r\partial rA_theta/\partial r - 1/r \partial A_r/\partial \theta
    
    // For non-zero A_phi component, B_r, B_theta poloidal component are non-zero
    

    for (int k=ks; k<=ke; ++k) {
      // reset loop limits for polar boundary
      jl=js; ju=je+1;
      if (block_bcs[INNER_X2] == POLAR_BNDRY) jl=js+1;
      if (block_bcs[OUTER_X2] == POLAR_BNDRY) ju=je;
      for (int j=jl; j<=ju; ++j) {
        pcoord->Face2Area(k,j,is,ie,area);
        pcoord->Edge3Length(k,j,is,ie+1,len);
        for (int i=is; i<=ie; ++i) {
          pfield->b.x2f(k,j,i) = -1.0*(len(i+1)*baphi(k,j,i+1)
                               -len(i)*baphi(k,j,i))/area(i);
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
          pcoord->Face1Area(k,j,is,ie+1,area);
          pcoord->Edge3Length(k,j  ,is,ie+1,len);
          pcoord->Edge3Length(k,j+1,is,ie+1,len_p1);
          for (int i=is; i<=ie+1; ++i) {
            pfield->b.x1f(k,j,i) = (len_p1(i)*baphi(k,j+1,i) -
                                    len(i)*baphi(k,j,i))/area(i);
          }
        }
      }

    }// end nx2 > 1


      
     // Update total energy with mangefew
    if(RADIATION_ENABLED){
      for(int k=ks; k<=ke; ++k){
      for(int j=js; j<=je; ++j){
      for(int i=is; i<=ie; ++i){
          phydro->u(IEN,k,j,i) +=
            0.5*(SQR(0.5*(pfield->b.x1f(k,j,i+1) + pfield->b.x1f(k,j,i)))
               + SQR(0.5*(pfield->b.x2f(k,j+1,i) + pfield->b.x2f(k,j,i)))
               + SQR(0.5*(pfield->b.x3f(k+1,j,i) + pfield->b.x3f(k,j,i))));
        
         }
      }}
      
    }

    baphi.DeleteAthenaArray();
    area.DeleteAthenaArray();
    len.DeleteAthenaArray();
    len_p1.DeleteAthenaArray();
 
  }// End MHD

  
  return;
}

//======================================================================================
//! \fn void MeshBlock::UserWorkInLoop(void)
//  \brief User-defined work function for every time step
//======================================================================================

void MeshBlock::UserWorkInLoop(void)
{
  // nothing to do
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
    
    while((logrhot > logrhottable(nrhot2)) && (nrhot2 < 36)){
      nrhot1 = nrhot2;
      nrhot2++;
    }
  
  /* The data point should between NrhoT1 and NrhoT2 */
    int nt1 = 0;
    int nt2 = 0;
    while((logt > logttable(nt2)) && (nt2 < 137)){
      nt1 = nt2;
      nt2++;
    }
    
    Real kappa;
    Real kappas = 0.2 * (1.0 + 0.6);
    
    if(nrhot1 == nrhot2){
      if(nt1 == nt2){
        kappa = opacitytable(nt1,nrhot1);
      }else{
        kappa = opacitytable(nt1,nrhot1) + (opacitytable(nt2,nrhot1)
              - opacitytable(nt1,nrhot1)) * (logt - logttable(nt1))/(logttable(nt2)
              - logttable(nt1));
      }/* end same T*/
    }else{
      if(nt1 == nt2){
        kappa = opacitytable(nt1,nrhot1) + (opacitytable(nt1,nrhot2)
              - opacitytable(nt1,nrhot1)) * (logrhot
              - logrhottable(nrhot1))/(logrhottable(nrhot2) - logrhottable(nrhot1));
      }else{
        kappa = opacitytable(nt1,nrhot1) * (logttable(nt2) - logt) *
                  (logrhottable(nrhot2)
                - logrhot)/((logttable(nt2) - logttable(nt1)) * (logrhottable(nrhot2)
                - logrhottable(nrhot1))) + opacitytable(nt2,nrhot1) * (logt
                - logttable(nt1)) * (logrhottable(nrhot2) - logrhot)/((logttable(nt2)
                - logttable(nt1)) * (logrhottable(nrhot2) - logrhottable(nrhot1)))
                + opacitytable(nt1,nrhot2) * (logttable(nt2) - logt) * (logrhot
                - logrhottable(nrhot1))/((logttable(nt2) - logttable(nt1))
                * (logrhottable(nrhot2) - logrhottable(nrhot1)))
                + opacitytable(nt2,nrhot2) * (logt - logttable(nt1)) * (logrhot
                - logrhottable(nrhot1))/((logttable(nt2)
                - logttable(nt1)) * (logrhottable(nrhot2) - logrhottable(nrhot1)));
      }
    }/* end same rhoT */
  
    if(kappa < kappas){
      kappas = kappa;
      kappa = 0.0;
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

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=NGHOST; ++i) {
          a(IDN,k,j,ie+i) = a(IDN,k,j,ie);
          a(IVX,k,j,ie+i) = std::max(a(IVX,k,j,ie),0.0);
          a(IVY,k,j,ie+i) = a(IVY,k,j,ie);
          a(IVZ,k,j,ie+i) = a(IVZ,k,j,ie);
          a(IEN,k,j,ie+i) = a(IEN,k,j,ie);
        
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
        b.x3f(k,j,ie+i) = b.x3f(k,j,ie);
      }
    }}
    
  }

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



void PseudoNewtonian(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, 
  AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{

  AthenaArray<Real> &x1flux=pmb->phydro->flux[x1face];

  for(int k=pmb->ks; k<=pmb->ke; ++k){
    for(int j=pmb->js; j<=pmb->je; ++j){
      for(int i=pmb->is; i<=pmb->ie; ++i){
        Real rho = prim(IDN,k,j,i);
        Real phic = -gm/(pmb->pcoord->x1v(i)-1.0);
        Real phil = -gm/(pmb->pcoord->x1f(i)-1.0);
        Real phir = -gm/(pmb->pcoord->x1f(i+1)-1.0);
        Real rr = pmb->pcoord->x1f(i+1);
        Real rl = pmb->pcoord->x1f(i);
        
        Real areal = rl * rl;
        Real arear = rr * rr;
        Real vol = (rr*rr*rr-rl*rl*rl)/3.0;
        Real src = - dt * rho * (phir - phil)/pmb->pcoord->dx1f(i);
        cons(IM1,k,j,i) += src;
        Real phidivrhov = (arear*x1flux(IDN,k,j,i+1) -
                           areal*x1flux(IDN,k,j,i))*phic/vol;
        Real divrhovphi = (arear*x1flux(IDN,k,j,i+1)*phir -
                           areal*x1flux(IDN,k,j,i)*phil)/vol;
        cons(IEN,k,j,i) += (dt*(phidivrhov - divrhovphi));
        
      }
    }
  }

}
