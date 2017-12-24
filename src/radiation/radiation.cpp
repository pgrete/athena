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
//! \file radiation.cpp
//  \brief implementation of functions in class Radiation
//======================================================================================


#include <sstream>  // msg
#include <iostream>  // cout
#include <stdexcept> // runtime erro
#include <stdio.h>  // fopen and fwrite


// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp" 
#include "radiation.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../globals.hpp"
#include "integrators/rad_integrators.hpp"

// constructor, initializes data structures and parameters

// The default opacity function.
// Do nothing. Keep the opacity as the initial value
inline void DefaultOpacity(MeshBlock *pmb, AthenaArray<Real> &prim)
{
  
}


Radiation::Radiation(MeshBlock *pmb, ParameterInput *pin)
{
  // read in the parameters
  int nmu = pin->GetInteger("radiation","nmu");
  int angle_flag = pin->GetOrAddInteger("radiation","angle_flag",0);
  prat = pin->GetReal("radiation","Prat");
  crat = pin->GetReal("radiation","Crat");
  rotate_theta=pin->GetOrAddInteger("radiation","rotate_theta",0);
  rotate_phi=pin->GetOrAddInteger("radiation","rotate_phi",0);
  reduced_c  = crat * pin->GetOrAddReal("radiation","reduced_factor",1.0);
  nfreq = pin->GetOrAddInteger("radiation","n_frequency",1);
  vmax = pin->GetOrAddReal("radiation","vmax",0.9);
  tunit = pin->GetOrAddReal("radiation","Tunit",1.e7);
  t_floor_ = pin->GetOrAddReal("radiation", "tfloor", TINY_NUMBER);

  ir_output=pin->GetOrAddInteger("radiation","ir_output",0);
  
  set_source_flag = 0;

  // equivalent temperature for electron
  telectron = 5.94065e9;
  telectron /= tunit;

  
  pmy_block = pmb;
  
  // calculate noct based on dimension
  int ndim = 1;
  if(pmb->block_size.nx2 > 1) ndim = 2;
  if(pmb->block_size.nx3 > 1) ndim = 3;
  
  int n1z = pmb->block_size.nx1 + 2*(NGHOST);
  int n2z = 1;
  int n3z = 1;
  
  
  
  int n_ang; // number of angles per octant and number of octant
  // total calculate total number of angles based on dimensions
  if(ndim == 1){
    n_ang = nmu;
    noct = 2;
  }else if(ndim == 2){
    noct = 4;
    n2z = pmb->block_size.nx2 + 2*(NGHOST);
    if(angle_flag == 0){
      n_ang = nmu * (nmu + 1)/2;
    }else if(angle_flag == 10){
      n_ang = nmu;
    }
  }else if(ndim == 3){
    noct = 8;
    n2z = pmb->block_size.nx2 + 2*(NGHOST);
    n3z = pmb->block_size.nx3 + 2*(NGHOST);
    if(angle_flag == 0){
      n_ang = nmu * (nmu + 1)/2;
    }else if(angle_flag == 10){
      n_ang = nmu * nmu/2;
    }
  }// end 3D
  
  nang = n_ang * noct;
  
  n_fre_ang = nang * nfreq;

  if(ir_output > n_fre_ang){    

    std::stringstream msg;
    msg << "### FATAL ERROR in Radiation Class" << std::endl
        << "number of output specific intensity is too large!";
    throw std::runtime_error(msg.str().c_str());
  }
  
  if(ir_output > 0){
    ir_index.NewAthenaArray(ir_output);
    dump_ir.NewAthenaArray(ir_output,n3z,n2z,n1z);
  }
  
  // allocate arrays
  // store frequency and angles as [nfre][ang]
  ir.NewAthenaArray(n3z,n2z,n1z,n_fre_ang);
  ir1.NewAthenaArray(n3z,n2z,n1z,n_fre_ang);
  
  rad_mom.NewAthenaArray(13,n3z,n2z,n1z);
  rad_mom_cm.NewAthenaArray(4,n3z,n2z,n1z);
  sigma_s.NewAthenaArray(n3z,n2z,n1z,nfreq);
  sigma_a.NewAthenaArray(n3z,n2z,n1z,nfreq);
  sigma_ae.NewAthenaArray(n3z,n2z,n1z,nfreq);
  sigma_planck.NewAthenaArray(n3z,n2z,n1z,nfreq);
  
  grey_sigma.NewAthenaArray(3,n3z,n2z,n1z);

  
  mu.NewAthenaArray(3,n3z,n2z,n1z,nang);
  wmu.NewAthenaArray(nang);

  wfreq.NewAthenaArray(nfreq);
  
  AngularGrid(angle_flag, nmu);
  
  //allocate memory to store the flux
  flux[X1DIR].NewAthenaArray(n3z,n2z,n1z,n_fre_ang);
  if(n2z > 1) flux[X2DIR].NewAthenaArray(n3z,n2z,n1z,n_fre_ang);
  if(n3z > 1) flux[X3DIR].NewAthenaArray(n3z,n2z,n1z,n_fre_ang);
  
  // Initialize the frequency weight
  FrequencyGrid();
  
  // set a default opacity function
  UpdateOpacity = DefaultOpacity;
  
  
  pradintegrator = new RadIntegrator(this, pin);
  
  // dump the angular grid and radiation parameters in a file
  if(Globals::my_rank ==0){
    FILE *pfile;
    std::stringstream msg;
    if((pfile = fopen("Rad_angles.txt","w")) == NULL){
        msg << "### FATAL ERROR in Radiation Class" << std::endl
            << "Output file Rad_angles.txt could not be opened";
        throw std::runtime_error(msg.str().c_str());
    }
      // damp the angular grid in one cell
    // The angles are still in cartesian
    fprintf(pfile,"Prat          %4.2e \n",prat);
    fprintf(pfile,"Crat          %4.2e \n",crat);
    fprintf(pfile,"reduced_c     %4.2e \n",reduced_c);
    fprintf(pfile,"Vmax          %4.2e \n",vmax);
    fprintf(pfile,"Tunit         %4.2e \n",tunit);
    fprintf(pfile,"Compt         %d  \n",pradintegrator->compton_flag_);
    fprintf(pfile,"Planck        %d  \n",pradintegrator->planck_flag_);
    fprintf(pfile,"Tfloor        %4.2e \n",t_floor_);
    fprintf(pfile,"rotate_theta  %d  \n",rotate_theta);
    fprintf(pfile,"rotate_phi    %d  \n",rotate_phi);
    fprintf(pfile,"adv_flag:     %d  \n",pradintegrator->adv_flag_);
    
    for(int n=0; n<nang; ++n){
      fprintf(pfile,"%2d   %e   %e   %e    %e\n",n,mu(0,0,0,0,n),mu(1,0,0,0,n),
             mu(2,0,0,0,n), wmu(n));
    }

    
    fclose(pfile);
  
  }
  
  

}

// destructor

Radiation::~Radiation()
{
  ir.DeleteAthenaArray();
  ir1.DeleteAthenaArray();
  rad_mom.DeleteAthenaArray();
  rad_mom_cm.DeleteAthenaArray();
  sigma_s.DeleteAthenaArray();
  sigma_a.DeleteAthenaArray();
  sigma_ae.DeleteAthenaArray();
  sigma_planck.DeleteAthenaArray();
  grey_sigma.DeleteAthenaArray();
  
  if(ir_output > 0){
    ir_index.DeleteAthenaArray();
    dump_ir.DeleteAthenaArray();
  }
  
  mu.DeleteAthenaArray();
  wmu.DeleteAthenaArray();
  wfreq.DeleteAthenaArray();
  
  flux[X1DIR].DeleteAthenaArray();
  if(pmy_block->block_size.nx2 > 1) flux[X2DIR].DeleteAthenaArray();
  if(pmy_block->block_size.nx3 > 1) flux[X3DIR].DeleteAthenaArray();
  
  delete pradintegrator;
  
}


//Enrol the function to update opacity

void Radiation::EnrollOpacityFunction(Opacity_t MyOpacityFunction)
{
  UpdateOpacity = MyOpacityFunction;
  
}


