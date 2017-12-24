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


//======================================================================================
/*! \file beam.cpp
 *  \brief Beam test for the radiative transfer module
 *
 *====================================================================================*/

// The space for opacity table

static AthenaArray<Real> opacitytable;
static AthenaArray<Real> planckopacity;
static AthenaArray<Real> logttable;
static AthenaArray<Real> logrhottable;

static AthenaArray<Real> ir_cm;


// The global variable

static Real gm;
static Real consFr = 6.34448e-4;
static Real kappaes = 111.28;

static Real rhounit = 5.0e-9;
static Real tunit;
static Real lunit = 6.955e10;
static Real tfloor;

static Real lbottom=25.0;

static int ninputline = 167938;

static int addperturbation=0;

void StarOpacity(MeshBlock *pmb, AthenaArray<Real> &prim);

//provide density and temperature, returns the opacity
void rossopacity(const Real rho, const Real tgas, Real &kappa, Real &kappa_planck);


void Inflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);


void Outflow_X2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);


void Inflow_rad_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);


void Outflow_rad_X2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);


void Mesh::InitUserMeshData(ParameterInput *pin)
{

   tfloor = pin->GetOrAddReal("radiation", "tfloor", 0.005);
   tunit = pin->GetOrAddReal("radiation","Tunit",1.6668e5);
   gm = pin->GetOrAddReal("problem","GM",0.0);
  
   EnrollUserBoundaryFunction(INNER_X1, Inflow_X1);
   EnrollUserBoundaryFunction(OUTER_X1, Outflow_X2);

  
   if(RADIATION_ENABLED){
   
     EnrollUserRadBoundaryFunction(INNER_X1, Inflow_rad_X1);
//     EnrollUserRadBoundaryFunction(OUTER_X1, Outflow_rad_X2);
   
     
      // create the memory and read in the opacity table
     opacitytable.NewAthenaArray(138,37);
     planckopacity.NewAthenaArray(138,37);
     
     logttable.NewAthenaArray(138);
     logrhottable.NewAthenaArray(37);
     
      FILE *fkappa, *flogT, *flogrhoT, *fplanck;
     
      if ( (fkappa=fopen("./aveopacity.txt","r"))==NULL )
      {
         printf("Open input file error");
         return;
      }

      if ( (fplanck=fopen("./PlanckOpacity.txt","r"))==NULL )
      {
         printf("Open input file error");
         return;
      }
  
      if ( (flogT=fopen("./logT.txt","r"))==NULL )
      {
         printf("Open input file error");
         return;
      }
  
      if ( (flogrhoT=fopen("./logRhoT.txt","r"))==NULL )
      {
         printf("Open input file error");
         return;
      }
  
     
      for(int j=0; j<138; j++){
        for(int i=0; i<37; i++){
          fscanf(fkappa,"%lf",&(opacitytable(j,i)));
        }
      }

      for(int j=0; j<138; j++){
        for(int i=0; i<37; i++){
          fscanf(fplanck,"%lf",&(planckopacity(j,i)));
        }
      }
   
      for(int i=0; i<37; i++){
        fscanf(flogrhoT,"%lf",&(logrhottable(i)));
      }
  
      for(int i=0; i<138; i++){
        fscanf(flogT,"%lf",&(logttable(i)));
      }
   
  
  
     fclose(fkappa);
     fclose(fplanck);
     fclose(flogT);
     fclose(flogrhoT);
   
   }// End Radiation

  return;
}

void MeshBlock::UserWorkInLoop()
{


   if(prad->set_source_flag > 0)
     prad->set_source_flag--;

   if(addperturbation==1){
       //initialize random number
     std::srand(gid);
   
     for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for(int i=is; i<=ie; ++i) {
          Real &x1 = pcoord->x1v(i);
          Real amp=0.0;
          if(x1 > 25.0){
            amp = 5.e-2;
          }else{
            amp = 0.0;
          }
          phydro->u(IDN,k,j,i) *= (1.0 + amp * ((double)rand()/(double)RAND_MAX-0.5));
        }
      }
    }
    addperturbation = 0;
  }
    
  return;
}



//======================================================================================
//! \fn void Mesh::TerminateUserMeshProperties(void)
//  \brief Clean up the Mesh properties
//======================================================================================
void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  
     if(RADIATION_ENABLED){
     
       opacitytable.DeleteAthenaArray();
       planckopacity.DeleteAthenaArray();
       logttable.DeleteAthenaArray();
       logrhottable.DeleteAthenaArray();
       
       ir_cm.DeleteAthenaArray();
       
     }
  
  
  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  

  
  if(RADIATION_ENABLED){
     prad->EnrollOpacityFunction(StarOpacity);

    
      // the bottom value of z coordinate in the block
     AthenaArray<Real> height, tgas, density;
    
     prad->set_source_flag = 0;
     // first, get the coordinate at the bottom of box
    
     ir_cm.NewAthenaArray(prad->n_fre_ang);
    

    
     // Initial profile from input data
     // Only the block with physical boundary condition needs this
     Real &x1 = pcoord->x1v(is-1);
     if(x1 < pmy_mesh->mesh_size.x1min){
     
        height.NewAthenaArray(ninputline);
        tgas.NewAthenaArray(ninputline);
        density.NewAthenaArray(ninputline);
    
        FILE *finput;
        if ( (finput=fopen("./Input.txt","r"))==NULL )
	      {
		      printf("Open input file error");
		      return;
	      }
        for(int i=0; i<ninputline; ++i){
          fscanf(finput,"%lf",&(height(i)));
          fscanf(finput,"%lf",&(tgas(i)));
	 	      fscanf(finput,"%lf",&(density(i)));
        }
  
        fclose(finput);
       
        // Get the temperature at the ghozt zones
        for(int i=1; i<=NGHOST; ++i){
           Real &x1 = pcoord->x1v(is-i);
         // get the position
           int lleft=0;

          int lright=1;
          while((x1 > height(lright)) && (lright < ninputline-1)){
              lright = lright+1;
          }
          if(lright - lleft > 1) lleft = lright -1;
    
          prad->tbot(i-1) = tgas(lleft) + (x1 - height(lleft)) *
                                (tgas(lright) - tgas(lleft))
                               /(height(lright) - height(lleft));
          
          prad->rhobot(i-1) = density(lleft) + (x1 - height(lleft)) *
                               (density(lright) - density(lleft))
                               /(height(lright) - height(lleft));

        }
    
    
        height.DeleteAthenaArray();
        tgas.DeleteAthenaArray();
        density.DeleteAthenaArray();
       
       
     
     
      }
    

   }else{

   }
  

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
  
  // the bottom value of z coordinate in the block
  AthenaArray<Real> height, tgas, density;
  
  Real amp = 0.0;
  
  
  height.NewAthenaArray(ninputline);
  tgas.NewAthenaArray(ninputline);
  density.NewAthenaArray(ninputline);
  
  
  // Initial profile from input data
  
  FILE *finput;
  if ( (finput=fopen("./Input.txt","r"))==NULL )
	{   
		printf("Open input file error");
		return;
	}
  for(int i=0; i<ninputline; ++i){
    fscanf(finput,"%lf",&(height(i)));
    fscanf(finput,"%lf",&(tgas(i)));
	 	fscanf(finput,"%lf",&(density(i)));
  }
  
  fclose(finput);



  //////////////////////////////////////////////
  
  // Initialize hydro variable
  for(int i=is; i<=ie; ++i) {
    Real &x1 = pcoord->x1v(i);
   
    
    // get the position
    int lleft=0;

    int lright=1;
    while((x1 > height(lright)) && (lright < ninputline-1)){
       lright = lright+1;
    }
    if(lright - lleft > 1) lleft = lright -1;
    
    Real rho = density(lleft) + (x1 - height(lleft)) *
                                (density(lright) - density(lleft))
                               /(height(lright) - height(lleft));
    Real tem = tgas(lleft) + (x1 - height(lleft)) *
                                (tgas(lright) - tgas(lleft))
                               /(height(lright) - height(lleft));
    
    if(rho > 5.e-2 && x1 > 22.0 && x1 < 35.0) amp = 5.e-2;
    else amp = 0.0;
    
 //   rho *= (1.0 + amp * ((double)rand()/(double)RAND_MAX-0.5));
    
    
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = rho * (amp * ((double)rand()/(double)RAND_MAX-0.5));
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS){

          phydro->u(IEN,k,j,i) = tem * rho/(gamma-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }
        
        if(RADIATION_ENABLED){
          Real er = tem * tem * tem * tem;
          // geometric dialution
          Real radratio = (lbottom/x1) * (lbottom/x1);
          
          for(int ifr=0; ifr<prad->nfreq; ++ifr){
            Real coefa = 0.0, coefb = 0.0;
            for(int n=0; n<prad->nang; ++n){
              // spherical polar coordinate
              Real &miuz = prad->mu(0,k,j,i,n);
              Real &weight = prad->wmu(n);
              if(miuz > 0.0){
                coefa += weight;
                coefb += (miuz * weight);
              }
            }
            
            for(int n=0; n<prad->nang; ++n){
              Real &miuz = prad->mu(0,k,j,i,n);
            
              if(miuz > 0.0){
                prad->ir(k,j,i,ifr*prad->nang+n) = 0.5 *
                                       (er/coefa + consFr * radratio/coefb);
              }else{
                prad->ir(k,j,i,ifr*prad->nang+n) = 0.5 *
                                       (er/coefa - consFr * radratio/coefb);
              
              }
            
            }
            
          }
        }// End Rad
        
      }// end j
    }// end k
  }// end i

  
  
  height.DeleteAthenaArray();
  tgas.DeleteAthenaArray();
  density.DeleteAthenaArray();
  
  return;
}



void StarOpacity(MeshBlock *pmb, AthenaArray<Real> &prim)
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

    // electron scattering opacity
  Real kappas = 0.2 * (1.0 + 0.6);
  Real kappaa = 0.0;
  
  for (int k=kl; k<=ku; ++k) {
  for (int j=jl; j<=ju; ++j) {
  for (int i=il; i<=iu; ++i) {
  for (int ifr=0; ifr<prad->nfreq; ++ifr){
    Real rho  = prim(IDN,k,j,i);
    Real gast = std::max(prim(IEN,k,j,i)/rho,tfloor);
    Real kappa, kappa_planck;
    rossopacity(rho, gast, kappa, kappa_planck);
    
    if(kappa < kappas){
      if(gast < 0.1){
        kappaa = kappa;
        kappa = 0.0;
      }else{
        kappaa = 0.0;
      }
    }else{
      kappaa = kappa - kappas;
      kappa = kappas;
    }

    prad->sigma_s(k,j,i,ifr) = kappa * rho * rhounit * lunit;
    prad->sigma_a(k,j,i,ifr) = kappaa * rho * rhounit * lunit;
    prad->sigma_ae(k,j,i,ifr) = prad->sigma_a(k,j,i,ifr);
    if(kappaa < kappa_planck)
      prad->sigma_planck(k,j,i,ifr) = (kappa_planck-kappaa)*rho*rhounit*lunit;
    else
      prad->sigma_planck(k,j,i,ifr) = 0.0;
  }
  }}}

}



void rossopacity(const Real rho, const Real tgas, Real &kappa, Real &kappa_planck)
{
  
    
    Real logt = log10(tgas * tunit);
    Real logrhot = log10(rho* rhounit) - 3.0* logt + 18.0;
    int nrhot1 = 0;
    int nrhot2 = 0;
    
    while((logrhot > logrhottable(nrhot2)) && (nrhot2 < 36)){
      nrhot1 = nrhot2;
      nrhot2++;
    }
    if(nrhot2==36 && (logrhot > logrhottable(nrhot2)))
      nrhot1=nrhot2;
  
  /* The data point should between NrhoT1 and NrhoT2 */
    int nt1 = 0;
    int nt2 = 0;
    while((logt > logttable(nt2)) && (nt2 < 137)){
      nt1 = nt2;
      nt2++;
    }
    if(nt2==137 && (logt > logttable(nt2)))
      nt1=nt2;
  

    Real kappa_t1_rho1=opacitytable(nt1,nrhot1);
    Real kappa_t1_rho2=opacitytable(nt1,nrhot2);
    Real kappa_t2_rho1=opacitytable(nt2,nrhot1);
    Real kappa_t2_rho2=opacitytable(nt2,nrhot2);

    Real planck_t1_rho1=planckopacity(nt1,nrhot1);
    Real planck_t1_rho2=planckopacity(nt1,nrhot2);
    Real planck_t2_rho1=planckopacity(nt2,nrhot1);
    Real planck_t2_rho2=planckopacity(nt2,nrhot2);


    Real rho_1 = logrhottable(nrhot1);
    Real rho_2 = logrhottable(nrhot2);
    Real t_1 = logttable(nt1);
    Real t_2 = logttable(nt2);

    
    if(nrhot1 == nrhot2){
      if(nt1 == nt2){
        kappa = kappa_t1_rho1;
        kappa_planck = planck_t1_rho1;
      }else{
        kappa = kappa_t1_rho1 + (kappa_t2_rho1 - kappa_t1_rho1) *
                                (logt - t_1)/(t_2 - t_1);
        kappa_planck = planck_t1_rho1 + (planck_t2_rho1 - planck_t1_rho1) *
                                (logt - t_1)/(t_2 - t_1);
      }/* end same T*/
    }else{
      if(nt1 == nt2){
        kappa = kappa_t1_rho1 + (kappa_t1_rho2 - kappa_t1_rho1) *
                                (logrhot - rho_1)/(rho_2 - rho_1);
        kappa_planck = planck_t1_rho1 + (planck_t1_rho2 - planck_t1_rho1) *
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
        
        kappa_planck = planck_t1_rho1 * (t_2 - logt) * (rho_2 - logrhot)/
                                ((t_2 - t_1) * (rho_2 - rho_1))
                     + planck_t2_rho1 * (logt - t_1) * (rho_2 - logrhot)/
                                ((t_2 - t_1) * (rho_2 - rho_1))
                     + planck_t1_rho2 * (t_2 - logt) * (logrhot - rho_1)/
                                ((t_2 - t_1) * (rho_2 - rho_1))
                     + planck_t2_rho2 * (logt - t_1) * (logrhot - rho_1)/
                                ((t_2 - t_1) * (rho_2 - rho_1));
      }
    }/* end same rhoT */
  
    return;

}



// This function sets boundary condition for primitive variables


void Inflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a, FaceField &b,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
   // set density according to force balance
  Radiation *prad=pmb->prad;
  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      Real jtop = 0.0;
      for(int n=0; n<prad->nang; ++n){
          Real& weight = prad->wmu(n);
          jtop += weight * prad->ir(k,j,is,0*prad->nang+n);
      }
      Real rho = pmb->phydro->w(IDN,k,j,is);
      Real tgas = pmb->phydro->w(IEN,k,j,is)/rho;
      Real kappa, kappa_planck;
      rossopacity(rho, tgas, kappa, kappa_planck);
      kappa *= (rho * rhounit * lunit);
      Real drsum = 0.0;
      
      for (int i=1; i<=NGHOST; ++i) {
        
//        kappa1 *= (rho1 * rhounit * lunit);

          Real dr = -pco->x1v(is-i) + pco->x1v(is-i+1);
          Real x1 = pmb->pcoord->x1f(is-i+1);
          Real radratio = (lbottom/x1) * (lbottom/x1);
          drsum += dr*radratio * kappa * consFr;

        
          Real jlocal = jtop + 3.0 * drsum;
          tgas=sqrt(sqrt(jlocal));
        

      
//          a(IDN,k,j,is-i) = a(IDN,k,j,is-i+1);
//         a(IDN,k,j,is-i) = a(IDN,k,j,is+i-1);
           a(IDN,k,j,is-i) = prad->rhobot(i-1);
/*          if(a(IVX,k,j,is) < 0.0)
            a(IVX,k,j,is-i) = -a(IVX,k,j,is);
          else
            a(IVX,k,j,is-i) = a(IVX,k,j,is-i+1)*rho*x1*x1/(x1g*x1g*a(IDN,k,j,is-i));
*/
          a(IVX,k,j,is-i) = -a(IVX,k,j,is+i-1);
          a(IVY,k,j,is-i) = a(IVY,k,j,is+i-1);
          a(IVZ,k,j,is-i) = a(IVZ,k,j,is+i-1);
        
          a(IEN,k,j,is-i) = prad->tbot(i-1)*a(IDN,k,j,is-i);
        
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
          Real &x1g = pco->x1v(ie+i);
          Real &x1 = pco->x1v(ie+i-1);
          if(a(IVX,k,j,ie) < 0.0){
            a(IDN,k,j,ie+i) = a(IDN,k,j,ie);
            a(IVX,k,j,ie+i) = 0.0;
          }else{
            a(IDN,k,j,ie+i) = a(IDN,k,j,ie+i-1) * x1*x1/(x1g*x1g);
            a(IVX,k,j,ie+i) = a(IVX,k,j,ie);
          }
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

void Inflow_rad_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  Radiation *prad=pmb->prad;
  
  Hydro *phydro = pmb->phydro;
  
  // assume T^4=Er at the bottom
  // And diffusiion equation d(Er)/dz=-3\sigma consFr
  
  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {

      Real jtop = 0.0;
        for(int n=0; n<prad->nang; ++n){
          Real& weight = prad->wmu(n);
          jtop += weight * prad->ir(k,j,is,0*prad->nang+n);
      }
      Real rho = phydro->w(IDN,k,j,is);
      Real tgas = phydro->w(IEN,k,j,is)/rho;
      Real kappa, kappa_planck;
      rossopacity(rho, tgas, kappa, kappa_planck);
      kappa *= (rho * rhounit * lunit);
      Real drsum = 0.0;
    
    
      for (int i=1; i<=NGHOST; ++i) {
      
//        Real rho = phydro->w(IDN,k,j,is-i+1);
//        Real tgas = phydro->w(IEN,k,j,is-i+1)/rho;
        
//        Real kappa = rossopacity(rho, tgas);
        
//        Real rho1 = phydro->w(IDN,k,j,is-i);
//        Real tgas1 = phydro->w(IEN,k,j,is-i)/rho1;
        
//        Real kappa1 = rossopacity(rho1, tgas1);

//        Real rho=phydro->w(IDN,k,j,is);
//        Real rho=rhobot[i-1];
//        Real kappa = rossopacity(rhobot[i-1], tbot[i-1]);
//        kappa *= (rho * rhounit * lunit);
        
        
//        kappa1 *= (rho1 * rhounit * lunit);

        Real dr = -pco->x1v(is-i) + pco->x1v(is-i+1);
        Real x1v = pmb->pcoord->x1v(is-i);
        Real radratio = (lbottom/x1v) * (lbottom/x1v);
        drsum += dr*radratio * kappa * consFr;


        Real r1 = pco->x1v(is-i+1);
        Real r2 = pco->x1v(is-i);
        
        
        for(int ifr=0; ifr<prad->nfreq; ++ifr){
          
//           Real jlocal = jtop + 3.0 * drsum;
          
           Real jlocal=prad->tbot(i-1) * prad->tbot(i-1) * prad->tbot(i-1) * prad->tbot(i-1);
/*
           for(int n=0; n<prad->nang; ++n){
             Real irtop = prad->ir(k,j,is-i+1,ifr*prad->nang+n);
             Real miur1 = prad->mu(0,k,j,is-i+1,n);
             Real miur2 = prad->mu(0,k,j,is-i,n);
             Real temptop = miur1 * irtop * r1 * r1 + 0.5*dr * kappa *
                            (jlocal + jtop - irtop)*x1*x1;
             Real tempbot = miur2 * r2 * r2 + 0.5 * dr * kappa*x1*x1;
             
             prad->ir(k,j,is-i,ifr*prad->nang+n) = temptop/tempbot;

           }

*/
         // construct Ir to satisfy jlocal and Fr
//         x1 = pmb->pcoord->x1v(is-i);
//         radratio = (lbottom/x1) * (lbottom/x1);
         Real coefa = 0.0, coefb = 0.0;
         Real coefa1 = 0.0, coefb1 = 0.0;
         for(int n=0; n<prad->nang; ++n){
           Real &miuz = prad->mu(0,k,j,is-i,n);
           Real &weight = prad->wmu(n);
           if(miuz > 0.0){
              coefa += weight;
              coefb += (miuz * weight);
           }else{
              coefa1 += weight;
              coefb1 += (miuz * weight);
           }

         }
        
         for(int n=0; n<prad->nang; ++n){
            Real &miuz = prad->mu(0,k,j,i,n);

            if(miuz > 0.0){
              prad->ir(k,j,is-i,ifr*prad->nang+n) = 0.5 *
                              (jlocal/coefa + consFr * radratio/coefb);
            }else{
              prad->ir(k,j,is-i,ifr*prad->nang+n) = 0.5 *
                              (jlocal/coefa1 + consFr * radratio/coefb1);
            }
          }

        
       }// end ifr
     
        
     }// end i
    }// end j
  }// end k
  

  return;
}

void Outflow_rad_X2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
    Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  Radiation *prad=pmb->prad;

  
  Hydro *phydro = pmb->phydro;
  
  Real kappas = 0.2 * (1 + 0.6);

  
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=NGHOST; ++i) {
      
         Real rho = phydro->w(IDN,k,j,ie+i);
         Real tgas = phydro->w(IEN,k,j,ie+i)/rho;
        
         Real kappa, kappa_planck;
         rossopacity(rho, tgas, kappa, kappa_planck);
//         if(kappa < kappas) kappa = kappas;
        
         kappa *= (rho * rhounit * lunit);
       
         Real dr = pmb->pcoord->dx1v(ie+i-1);
        
         Real tausq = dr * kappa * dr * kappa;
        
        
         if(tausq > 1.e-3) tausq = 1.0 - exp(-tausq);
        
        
         for(int ifr=0; ifr<prad->nfreq; ++ifr){
            for(int n=0; n<prad->nang; ++n){
                Real miuz = prad->mu(0,k,j,ie+i,ifr*prad->nang+n);
                if(miuz > 0.0){
                  prad->ir(k,j,ie+i,ifr*prad->nang+n)
                                = prad->ir(k,j,ie+i-1,ifr*prad->nang+n);
                }else{
                  prad->ir(k,j,ie+i,ifr*prad->nang+n) = tausq *
                                prad->ir(k,j,ie+i-1,ifr*prad->nang+n);
                }
         
            }
         }
        
        
      }
    }
  }
  

  return;
}



