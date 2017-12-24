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
#include "../utils/utils.hpp"

static Real inib0 = 0.5;
static Real rhoinf;
static Real rang;
static Real csinf;
static Real coolcoef;
static Real dfloor;
static Real tfloor;

static AthenaArray<Real> bondisol;

void BondiRho(double rho, double coef1, double coef2, double coef3,
    double gamma, double * fval, double *dfval);


void CoolingFunc( MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim,
  AthenaArray<Real> &bcc, AthenaArray<Real> &cons);

void BondiInflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
               FaceField &b, int is, int ie, int js, int je, int ks, int ke);

void BondiOutflow_X2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
               FaceField &b, int is, int ie, int js, int je, int ks, int ke);

//======================================================================================
/*! \file beam.cpp
 *  \brief Beam test for the radiative transfer module
 *
 *====================================================================================*/


void Mesh::InitUserMeshData(ParameterInput *pin)
{
   rhoinf=pin->GetReal("problem","rhoinf");
   rang=pin->GetReal("problem","Jinf");
   coolcoef=pin->GetReal("problem","Cooling");
   dfloor=pin->GetReal("hydro","dfloor");
   tfloor=pin->GetReal("hydro","tfloor");
   csinf = 1.0;
  
   // Enrol boundary condition
   EnrollUserBoundaryFunction(INNER_X1, BondiInflow_X1);
   EnrollUserBoundaryFunction(OUTER_X1, BondiOutflow_X2);
  
   EnrollUserSourceTermFunction(CoolingFunc);


  return;
}

// construct the Bondi solution for boundary condition
// particularly needed after restart

void MeshBlock::InitUserMeshBlockProperties(ParameterInput *pin)
{

   bondisol.NewAthenaArray(NHYDRO,je-js+1+2*NGHOST,NGHOST);
   // This is for boundary condition
   // We only need the primitive variables
  
   Real gamma = phydro->peos->GetGamma();
   Real fgam = pow(2.0/(5.0 - 3.0 * gamma), (5.0 - 3.0 * gamma)/(2*(gamma - 1.0)));
   Real mdot = PI * fgam * 0.25 * rhoinf  * csinf;
   Real rsonic = (5.0 - 3.0 * gamma) * 0.125;
   Real coef3 = fgam * fgam /(32.0 * 16.0);
   
   Real ang = sqrt(rang * 0.5); /* Then J = ang * RB *Cs */
  
  // Only for the ghost zones
   for (int k=ks; k<=ke; ++k) {
    for (int j=0; j<=je+NGHOST; ++j) {
    for (int i=ie+1; i<=ie+NGHOST; ++i) {
      Real radius = pcoord->x1v(i);
      Real theta = pcoord->x2v(j);
      Real coef1 =  radius * radius * radius * radius/(gamma - 1.0);
	    Real coef2 = -(0.5/radius + 1.0/(gamma - 1.0)) *
                 radius * radius * radius * radius;

  	  /* The solution we need to solve to find density is *
	   * coef1 * rho^(gamma+1) + coef2 * rho^2 + coef3 = 0
	   */
	    Real rho1 = pow(0.5 * (gamma - 1.0)/radius + 1.0, 1.0/(gamma - 1.0));
	    Real rho2 = pow((gamma - 1.0)/((gamma + 1.0) * radius)
            + 2.0/(gamma + 1.0), 1.0/(gamma - 1.0));
	    Real rho3 = (fgam /(SQRT2 * 16.0 * radius * radius)) *
                sqrt(2.0 * radius * (gamma - 1.0)/(gamma - 1.0 + 2.0 * radius));
	    Real rho;
      if(radius < rsonic)
	      rho = Rtsafe(BondiRho, rho3, rho2, 1.e-12, coef1, coef2, coef3, gamma);
	    else
	      rho = Rtsafe(BondiRho, rho2, rho1, 1.e-12, coef1, coef2, coef3, gamma);

	    rho *= rhoinf;

      
      Real rcyl = radius * sin(theta);
      Real rscale = pow(rcyl/(0.5*rang), 4.0);


      bondisol(IDN,j,i-ie-1) = rho;
      // velocity

      bondisol(IVX,j,i-ie-1) = -mdot/(4.0 * PI * radius * radius * rho);
      bondisol(IVY,j,i-ie-1) = 0.0;

      if(radius > rang)
        bondisol(IVZ,j,i-ie-1) = ang * sqrt(rscale/(1.0 + rscale))/rcyl;
      else
        bondisol(IVZ,j,i-ie-1) = 0.0;

      // pressure
      bondisol(IEN,j,i-ie-1) = rhoinf * (pow(rho/rhoinf, gamma))
                                /(gamma * (gamma - 1.0));
    
    }}}

  return;
}

void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  bondisol.DeleteAthenaArray();
  
  return;
}


//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  
  // The problem generator needs conservative variables
  
  Real gamma = phydro->peos->GetGamma();
  Real fgam = pow(2.0/(5.0 - 3.0 * gamma), (5.0 - 3.0 * gamma)/(2*(gamma - 1.0)));
  Real mdot = PI * fgam * 0.25 * rhoinf  * csinf;
  Real rsonic = (5.0 - 3.0 * gamma) * 0.125;
  Real coef3 = fgam * fgam /(32.0 * 16.0);
  
  Real ang = sqrt(rang * 0.5); /* Then J = ang * RB *Cs */
  
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
      
        Real radius = pcoord->x1v(i);
        Real theta = pcoord->x2v(j);
        Real coef1 =  radius * radius * radius * radius/(gamma - 1.0);
	      Real coef2 = -(0.5/radius + 1.0/(gamma - 1.0)) *
                 radius * radius * radius * radius;

  	  /* The solution we need to solve to find density is *
	   * coef1 * rho^(gamma+1) + coef2 * rho^2 + coef3 = 0
	   */
	      Real rho1 = pow(0.5 * (gamma - 1.0)/radius + 1.0, 1.0/(gamma - 1.0));
	      Real rho2 = pow((gamma - 1.0)/((gamma + 1.0) * radius)
              + 2.0/(gamma + 1.0), 1.0/(gamma - 1.0));
	      Real rho3 = (fgam /(SQRT2 * 16.0 * radius * radius)) *
                sqrt(2.0 * radius * (gamma - 1.0)/(gamma - 1.0 + 2.0 * radius));
	      Real rho;
        if(radius < rsonic)
	        rho = Rtsafe(BondiRho, rho3, rho2, 1.e-12, coef1, coef2, coef3, gamma);
	      else
	        rho = Rtsafe(BondiRho, rho2, rho1, 1.e-12, coef1, coef2, coef3, gamma);

	      rho *= rhoinf;
        
      
        Real rcyl = radius * sin(theta);
        Real rscale = pow(rcyl/(0.5*rang), 4.0);

        phydro->u(IDN,k,j,i) = rho;
        // velocity

        phydro->u(IM1,k,j,i) = -mdot/(4.0 * PI * radius * radius);
        phydro->u(IM2,k,j,i) = 0.0;

        if(radius > 0.0 * rang)
          phydro->u(IM3,k,j,i) = rho * ang * sqrt(rscale/(1.0 + rscale))/rcyl;
        else
          phydro->u(IM3,k,j,i) = 0.0;

        // pressure
        phydro->u(IEN,k,j,i) = rhoinf * (pow(rho/rhoinf, gamma))
                                /(gamma * (gamma - 1.0));

        phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
        phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);

      
      }
    }
  }
  
  // Magnetic field is carried inward from outer boundary condition
  // No need to initialize magnetic field inside the domain

  if (MAGNETIC_FIELDS_ENABLED) {
  
  
    for (int k=ks; k<=ke; ++k) {
      // reset loop limits for polar boundary
      int jl=js; int ju=je+1;

      for (int j=jl; j<=ju; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfield->b.x2f(k,j,i) = 0.0;
        }
      }
    }

    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          pfield->b.x3f(k,j,i) = 0.0;
        }
      }
    }


      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie+1; ++i) {
            pfield->b.x1f(k,j,i) = 0.0;
          }
        }
      }




  }

 
  
  
  
  
  return;
}


// out boundary at out boundary
// Fix the outer boundary to be Bondi solution with uniform B
void BondiOutflow_X2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
               FaceField &b, int is, int ie, int js, int je, int ks, int ke)
{

  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=1; i<=NGHOST; ++i) {
         if(a(IVX,k,j,ie) > 0.0){
            a(IDN,k,j,ie+i) = a(IDN,k,j,ie);
            a(IVX,k,j,ie+i) = a(IVX,k,j,ie);
            a(IVY,k,j,ie+i) = a(IVY,k,j,ie);
            a(IVZ,k,j,ie+i) = a(IVZ,k,j,ie);
            a(IEN,k,j,ie+i) = a(IEN,k,j,ie);
         }else{
            a(IDN,k,j,ie+i)=bondisol(IDN,j,i-1);
            a(IVX,k,j,ie+i)=bondisol(IVX,j,i-1);
            a(IVY,k,j,ie+i)=bondisol(IVY,j,i-1);
            a(IVZ,k,j,ie+i)=bondisol(IVZ,j,i-1);
            a(IEN,k,j,ie+i)=bondisol(IEN,j,i-1);
         }
      }
    }
  }
   // set magnetic field in inlet ghost zones
  if (MAGNETIC_FIELDS_ENABLED) {
    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){
      for(int i=1; i<=NGHOST; ++i){
        if(a(IVX,k,j,ie) > 0.0)
          b.x1f(k,j,ie+i+1) = b.x1f(k,j,ie+1);
        else
          b.x1f(k,j,ie+i+1) = 0.0;
      }
    }}

    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je+1; ++j){
      for(int i=1; i<=NGHOST; ++i){
        if(a(IVX,k,j,ie) > 0.0)
          b.x2f(k,j,ie+i) = b.x2f(k,j,ie);
        else
          b.x2f(k,j,ie+i) = 0.0;
      }
    }}

    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je; ++j){
      for(int i=1; i<=NGHOST; ++i){
        if(a(IVX,k,j,ie) > 0.0)
           b.x3f(k,j,ie+i) = b.x3f(k,j,ie);
        else
           b.x3f(k,j,ie+i) = inib0;
      }
    }}
    
  }

  return;
}


void BondiInflow_X1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &a,
               FaceField &b, int is, int ie, int js, int je, int ks, int ke)
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
      for(int i=1; i<=NGHOST; ++i){
        b.x1f(k,j,is-i) = b.x1f(k,j,is);
      }
    }}

    for(int k=ks; k<=ke; ++k){
    for(int j=js; j<=je+1; ++j){
      for(int i=1; i<=NGHOST; ++i){
        b.x2f(k,j,is-i) = b.x2f(k,j,is);
      }
    }}
    
    int kl=ks, ku=ke;
    if(ku>kl) ku++;

    for(int k=kl; k<=ku; ++k){
    for(int j=js; j<=je; ++j){
      for(int i=1; i<=NGHOST; ++i){
        b.x3f(k,j,is-i) = b.x3f(k,j,is);
      }
    }}
    
  }


  return;
}




/* Function to find the density for Bondi solution */
void BondiRho(double rho, double coef1, double coef2, double coef3,
    double gamma, double * fval, double *dfval)
{

/* function is
 *  coef1 * rho^(gamma+1) + coef2 * rho^2 + coef3 == 0 *
 */

  *fval = coef1 * pow(rho, gamma + 1.0) + coef2 * rho * rho + coef3;
  *dfval = (gamma + 1.0) * coef1 * pow(rho, gamma) + 2.0 * coef2 * rho;

  return;
}


void CoolingFunc( MeshBlock *pmb, const Real time, const Real dt,
  const AthenaArray<Real> &prim,
  AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{

  Real gamma = pmb->phydro->peos->GetGamma();

  if (coolcoef < TINY_NUMBER) return;


  for (int k=pmb->ks; k<=pmb->ke; ++k) {
  for (int j=pmb->js; j<=pmb->je; ++j) {
    for (int i=pmb->is; i<=pmb->ie; ++i) {
	// first, set density and temperature floor
      if(cons(IDN,k,j,i) < dfloor) cons(IDN,k,j,i) = dfloor;

      Real rho = cons(IDN,k,j,i);

      Real solcoef = dt * coolcoef * rho * sqrt(rho);
      Real eint = prim(IEN,k,j,i); // pressure
      if(eint/rho < tfloor) eint = rho * tfloor; // set temperature floor
      eint /= (gamma - 1.0);
      // internal energy = pressure/(gamma - 1)

      Real enew = 0.5 * solcoef * (sqrt(1.0 + 4.0 * eint/(solcoef * solcoef)) - 1.0);
      enew = (enew * enew);


      Real den = (enew - eint);
      if(fabs(den/eint) > 0.05) den *= fabs((0.05 * eint)/(enew - eint));
      
      if(rho > 1.e4 || den > 0.) den = 0.0;

      cons(IEN,k,j,i) += den;

    }
  }}
  
  
}

