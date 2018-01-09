//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particle-gas.cpp
//  \brief tests the implementation of particles.

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../particles/particles.hpp"

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Sets the initial conditions.
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  // Set a uniform, steady gas.
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        phydro->u(IDN,k,j,i) = 1.0;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = 1.0 / (peos->GetGamma() - 1.0);
      }
    }
  }

  // Set the particles.
  if (block_size.x1min <= 1.0 && 1.0 < block_size.x1max) {
    ppar->npar = 1;
    ppar->pid(0) = 0;
    ppar->xp(0) = 1.0;
    ppar->yp(0) = 0.0;
    ppar->zp(0) = 0.0;
    ppar->vpx(0) = 0.0;
    ppar->vpy(0) = 0.0;
    ppar->vpz(0) = 0.0;
  } else
    ppar->npar = 0;
}
