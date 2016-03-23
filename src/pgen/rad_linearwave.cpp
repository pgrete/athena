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

static Real amp = 1.e-6;
static Real sigma0 = 1.0;

//======================================================================================
/*! \file beam.cpp
 *  \brief Beam test for the radiative transfer module
 *
 *====================================================================================*/


void Mesh::InitUserMeshProperties(ParameterInput *pin)
{

  return;
}

//======================================================================================
//! \fn void Mesh::TerminateUserMeshProperties(void)
//  \brief Clean up the Mesh properties
//======================================================================================
void Mesh::TerminateUserMeshProperties(void)
{
  // nothing to do
  return;
}



//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief beam test
//======================================================================================
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{

  Real gamma = phydro->peos->GetGamma();
  Real knum = 2.0 * PI;
  Real omegareal = 6.398479314825398;
  Real omegaimg =  0.7044806086435424;
  Real rho0 = 1.0, p0 = 1.0;
  Real e0 = p0/(gamma-1.0);


  
  // Initialize hydro variable
  for(int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        Real &x1 = pcoord->x1v(i);
        Real &x2 = pcoord->x2v(j);
        Real &x3 = pcoord->x3v(k);
        Real theta = knum * x1;
        Real delv = amp* (1.0183496112257058 * cos(theta)
                    + 0.1121215711780068 * sin(theta));
        Real delp = amp * (1.0220380692314723 * cos(theta)
                    + 0.18993018794365163 * sin(theta));
      
      
        phydro->u(IDN,k,j,i) = rho0 + amp * cos(theta);
        phydro->u(IM1,k,j,i) = delv;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        if (NON_BAROTROPIC_EOS){

          phydro->u(IEN,k,j,i) = (p0 + delp)/(gamma-1.0);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM1,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM2,k,j,i))/phydro->u(IDN,k,j,i);
          phydro->u(IEN,k,j,i) += 0.5*SQR(phydro->u(IM3,k,j,i))/phydro->u(IDN,k,j,i);
        }
        
        if(RADIATION_ENABLED){
          Real der = amp * (-0.026018127896336885 * cos(theta)
                  + 0.12095401964915764 * sin(theta));
          Real dfr = amp * (-0.10566859341556321 * cos(theta)
                  + 0.030196412832965945 * sin(theta));
        
          Real jr = (1.0+der);
          Real hr = dfr;
          for(int ifr=0; ifr<prad->nfreq; ++ifr){
            for(int n=0; n<prad->nang; ++n){
               Real& weight = prad->wmu(n);
               Real& miux = prad->mu(0,k,j,i,n);
               prad->ir(k,j,i,ifr*prad->nang+n)=(jr/(4.0*weight) + hr/(4.0*weight*miux));
            }
            
            prad->sigma_s(k,j,i,ifr) = 0.0;
            prad->sigma_a(k,j,i,ifr) = sigma0;
            prad->sigma_ae(k,j,i,ifr) = sigma0;
          }
        
        }// End rad
      }// end i
    }
  }
  
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





