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
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../bvals/bvals.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"
#include "../radiation/radiation.hpp"
#include "../radiation/integrators/rad_integrators.hpp"
#include "../utils/utils.hpp"

// The global parameters
static Real kappaes = 6681.63;
static Real kappaffp = 177.58;
static Real kappaffr = 4.79946;
static Real rho0 = 1.0; // Density normalization of torus
static Real inib0 = 0.1; // initial magnetic field strength
static Real r0 = 40.0; // Center of initial torus
static Real tfloor; // temperature floor
static Real rhofloor; // density floor
static int bconf = 2; // bconf=1: pure B_phi
                      // bconf=0: vector potential proportional to density
                      // bconf=2: two loops
static int nloop = 1;
static Real lprofile= 0.4;
static Real vs0 = 7.07;
static int iniflag=0;

static Real gm;


//======================================================================================
/*! \file globaldisk.cpp
 *  \brief global accretion disk problem with radiation
 *
 *====================================================================================*/

void Inflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);


void Outflow_X2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);


void DiskOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);


void PseudoNewtonian( MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim,
  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons);


void Mesh::InitUserMeshData(ParameterInput *pin)
{
  
    // Enroll boundary functions

  EnrollUserBoundaryFunction(INNER_X1, Inflow_X1);
  EnrollUserBoundaryFunction(OUTER_X1, Outflow_X2);
  
  tfloor = pin->GetOrAddReal("radiation", "tfloor", 0.005);
  rhofloor = pin->GetOrAddReal("hydro", "dfloor", 1.e-9);
  EnrollUserExplicitSourceFunction(PseudoNewtonian);
  

  return;
}



void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  
  
  if(RADIATION_ENABLED){
      prad->EnrollOpacityFunction(DiskOpacity);
    
      gm = 0.5 * prad->crat * prad->crat;

  
      prad->set_source_flag = 0;
  
  }else{
      gm = 0.5 * 805.338 * 805.338;
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
    Real gm1 = peos->GetGamma() - 1.0;

    Real wi[(NWAVE)];
    Real cfmax;
    
    
    for (int k=kl; k<=ku; ++k){
      for (int j=jl; j<=ju; ++j){
       for (int i=il; i<=iu; ++i){

          Real& vx=phydro->w(IVX,k,j,i);
          Real& vy=phydro->w(IVY,k,j,i);
          Real& vz=phydro->w(IVZ,k,j,i);

          Real& rho=phydro->w(IDN,k,j,i);

//          if(iniflag){
//             if(rho < 1.01e-6){
//                 phydro->w(IDN,k,j,i) = 5.e-7;
//                 phydro->u(IDN,k,j,i) = 5.e-7;
//             }

//          }



          Real vel = sqrt(vx*vx+vy*vy+vz*vz);
          Real ke = 0.5 * rho * vel * vel;
          int flag=0;

          Real pb=0.0;

          // case 1, check superlum velocity
          if(vel > prad->vmax * prad->crat){
            Real ratio = prad->vmax * prad->crat / vel;
            vx *= ratio;
            vy *= ratio;
            vz *= ratio;

            phydro->u(IM1,k,j,i) = rho*vx;
            phydro->u(IM2,k,j,i) = rho*vy;
            phydro->u(IM3,k,j,i) = rho*vz;

            ke = 0.5 * rho * (vx*vx+vy*vy+vz*vz);
           
            flag=1;
          }
         
          // case 2, check fast speed too large
          if(MAGNETIC_FIELDS_ENABLED){
            wi[IDN]=rho;
            wi[IVX]=vx;
            wi[IVY]=vy;
            wi[IVZ]=vz;

            Real bx = pfield->bcc(IB1,k,j,i) +
                      fabs(pfield->b.x1f(k,j,i)-pfield->bcc(IB1,k,j,i));
            
            wi[IBY] = pfield->bcc(IB2,k,j,i);
            wi[IBZ] = pfield->bcc(IB3,k,j,i);
            Real cf = peos->FastMagnetosonicSpeed(wi,bx);
            cfmax = (fabs(wi[IVX]) + cf);

            bx = pfield->bcc(IB2,k,j,i) +
                      fabs(pfield->b.x2f(k,j,i)-pfield->bcc(IB2,k,j,i));

            wi[IBY] = pfield->bcc(IB3,k,j,i);
            wi[IBZ] = pfield->bcc(IB1,k,j,i);
            cf = peos->FastMagnetosonicSpeed(wi,bx);
            cf = (fabs(wi[IVY]) + cf);

            if(cf> cfmax) cfmax = cf;

            bx = pfield->bcc(IB3,k,j,i) +
                fabs(pfield->b.x3f(k,j,i)-pfield->bcc(IB3,k,j,i));

            wi[IBY] = pfield->bcc(IB1,k,j,i);
            wi[IBZ] = pfield->bcc(IB2,k,j,i);
            cf = peos->FastMagnetosonicSpeed(wi,bx);
            cf = (fabs(wi[IVZ]) + cf);
               
            if(cf> cfmax) cfmax = cf;
            
            Real ratio = cfmax/(prad->vmax * prad->crat);
/*            if(ratio > 1.0){
               ratio = ratio * ratio;
               phydro->w(IDN,k,j,i) *= ratio;
               phydro->u(IDN,k,j,i) *= ratio;

               phydro->u(IM1,k,j,i) *= ratio;
               phydro->u(IM2,k,j,i) *= ratio;
               phydro->u(IM3,k,j,i) *= ratio;

               ke *= ratio;
              
               flag=1;

            }
*/            
            pb = 0.5*(SQR(pfield->bcc(IB1,k,j,i))+SQR(pfield->bcc(IB2,k,j,i))
                 +SQR(pfield->bcc(IB3,k,j,i)));
          }

          if(flag > 0){
             // Only do this for bad cells
             Real  eint = phydro->w(IEN,k,j,i)/gm1;
             phydro->u(IEN,k,j,i) = eint + ke + pb;

          }

      }}}
    }

  if(iniflag > 0) iniflag=iniflag-1;
  return;

}


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  
  Real gamma = peos->GetGamma();

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
    crat = 805.338;
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

  if(RADIATION_ENABLED){
    
    ir_cm.DeleteAthenaArray();

  }
  
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
    


   if(nloop == 1){
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
           Real rholimit = 0.001;
           Real press=rho0*vs0*vs0*pow(phydro->u(IDN,k,j,i)/rho0,gamma)/gamma;
//           if(phydro->u(IDN,k,j,i) > rholimit){
//             aphi = inib0 * pow(((phydro->u(IDN,k,j,i)-rholimit)/rho0)*pow(x1,0.75)
//                    ,2.0)*sin(log(x1/(0.5*r0))/0.01);
           if(press > rholimit){
             aphi = inib0 * pow(((press-rholimit)/(rho0*vs0*vs0))*pow(x1,0.75)
                    ,2.0)*sin(log(x1/(0.5*r0))/0.01);
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
   }else if(nloop==3){
     Real lambdaB = 2.0;
     Real ri = 20.0;
     Real ro = 60.0;
     Real aphi;

     for(int k=ks; k<=ku; ++k){
     for(int j=jl; j<=ju; ++j){
     for(int i=il; i<=iu; ++i){
      Real rho = phydro->u(IDN,k,j,i);
    
      // reconstruct the density at (r, pi/2)
      Real &x1 = pcoord->x1v(i);
      // the cylindrical radius
      Real angradius = x1;
      Real langular = l0 * pow(angradius/r0,lprofile);
      Real vphi = langular/angradius;
      Real effphi = -crat*crat/(2.0*(x1-1.0))
                  +pow((langular/angradius),2.0)/(2.0*(1.0-lprofile));
      Real effphi0 = -crat*crat/(2.0*(r0-1.0))
                  +pow((l0/r0),2.0)/(2.0*(1.0-lprofile));
    
      Real tempphi = ((effphi-effphi0)/nindex)/(vs0*vs0);
      Real rhomid;
      if((fabs(tempphi)<1.0) && (x1 > 8.0)){
        rhomid = rho0*pow(fabs(1.0-tempphi),nindex);
        rhomid = std::max(rhomid,rhofloor);
      }else{
        rhomid = rhofloor;
      }
      
      // reconstruct the density at (ro, pi/2)
      angradius = ro;
      langular = l0 * pow(angradius/r0,lprofile);
      vphi = langular/angradius;
      effphi = -crat*crat/(2.0*(ro-1.0))
                  +pow((langular/angradius),2.0)/(2.0*(1.0-lprofile));
      effphi0 = -crat*crat/(2.0*(r0-1.0))
                  +pow((l0/r0),2.0)/(2.0*(1.0-lprofile));
    
      tempphi = ((effphi-effphi0)/nindex)/(vs0*vs0);
      Real rhoend;
      if((fabs(tempphi)<1.0) && (ro > 8.0)){
        rhoend = rho0*pow(fabs(1.0-tempphi),nindex);
        rhoend = std::max(rhoend,rhofloor);
      }else{
        rhoend = rhofloor;
      }
      
      Real qfactor;
      Real rhoratio = (rho-rhoend)/(rhomid-rhoend);
      if(x1 > ri && x1 < ro && rho > 1000*rhofloor){
        Real sintheta = sin(pcoord->x2v(j));
        qfactor = sintheta*sintheta*sintheta*(rhoratio - 0.2)/0.8;
      }else{
        qfactor = 0.0;
      }
      Real ffactor = (pow(x1,2.0/3.0)+15.0*pow(x1,-0.4)/8.0)/lambdaB;
      Real fifactor = (pow(ri,2.0/3.0)+15.0*pow(ri,-0.4)/8.0)/lambdaB;
      
      aphi = inib0 * qfactor * sin(ffactor - fifactor)*sin(0.5*PI-pcoord->x2v(j));
      
      baphi(k,j,i) = aphi;

  
  
     }}}// end k, j, i
   
   
   
   }
       
    //B=div X Phi
    // vector potential only has non-zero phi component
    // in spherical polar coordinate system
    // B_r= 1/rsintheta\partial (sintheta A_phi)/partial theta
    //      - 1/rsintheta \partial A_theta/\partial phi
    
    // B_theta=1/rsintheta\partial A_r/\partial phi - 1/r\partial (rA_phi)/\partial r
    // B_phi=1/r\partial rA_theta/\partial r - 1/r \partial A_r/\partial \theta
    
    // For non-zero A_phi component, B_r, B_theta poloidal component are non-zero
    
  if(bconf ==0 || bconf == 2){
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

  }else if(bconf == 1){

    for (int k=ks; k<=ke; ++k) {
      // reset loop limits for polar boundary
      jl=js; ju=je+1;
      if (block_bcs[INNER_X2] == POLAR_BNDRY) jl=js+1;
      if (block_bcs[OUTER_X2] == POLAR_BNDRY) ju=je;
      for (int j=jl; j<=ju; ++j) {
        pcoord->Face2Area(k,j,is,ie,area);
        pcoord->Edge3Length(k,j,is,ie+1,len);
        for (int i=is; i<=ie; ++i) {
          pfield->b.x2f(k,j,i) = inib0;
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
            pfield->b.x1f(k,j,i) = 0.0;
          }
        }
      }

    }// end nx2 > 1


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

    baphi.DeleteAthenaArray();
    area.DeleteAthenaArray();
    len.DeleteAthenaArray();
    len_p1.DeleteAthenaArray();
 
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
    Real tpower= 1.0/(gast*gast*gast*sqrt(gast));

    prad->sigma_s(k,j,i,ifr) = kappaes * rho;
    prad->sigma_a(k,j,i,ifr) = kappaffr * rho * rho * tpower;
    prad->sigma_ae(k,j,i,ifr) = prad->sigma_a(k,j,i,ifr);
    prad->sigma_planck(k,j,i,ifr) = (kappaffp-kappaffr)*rho*rho*tpower;
  }
  }}}

}

// This function sets boundary condition for primitive variables


void Inflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
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
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
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



void PseudoNewtonian(MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim, 
  const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{

  AthenaArray<Real> &x1flux=pmb->phydro->flux[X1DIR];

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
