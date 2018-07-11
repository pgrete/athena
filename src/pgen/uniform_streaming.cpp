//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file one_particle.cpp
//  \brief tests one particle.

// C++ standard libraries
#include <sstream>

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
  // Get the (uniform) velocity of the gas.
  Real ux0, uy0, uz0;
  ux0 = pin->GetOrAddReal("problem", "ux0", 0.0);
  uy0 = pin->GetOrAddReal("problem", "uy0", 0.0);
  uz0 = pin->GetOrAddReal("problem", "uz0", 0.0);

  // Set a uniform, steady gas.
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        phydro->u(IDN,k,j,i) = 1.0;
        phydro->u(IM1,k,j,i) = ux0;
        phydro->u(IM2,k,j,i) = uy0;
        phydro->u(IM3,k,j,i) = uz0;
        phydro->u(IEN,k,j,i) = 1.0 / (peos->GetGamma() - 1.0);
      }
    }
  }

  // Get the dust-to-gas ratio and the velocity of the particles.
  Real dtog, vpx0, vpy0, vpz0;
  dtog = pin->GetOrAddReal("problem", "dtog", 1.0);
  vpx0 = pin->GetOrAddReal("problem", "vpx0", 0.0);
  vpy0 = pin->GetOrAddReal("problem", "vpy0", 0.0);
  vpz0 = pin->GetOrAddReal("problem", "vpz0", 0.0);

  // Find the mass of each particle.
  long npar = 1;
  Real dvol = 1.0;
  if (block_size.nx1 > 1) {
    npar *= block_size.nx1;
    dvol *= block_size.x1max - block_size.x1min;
  }
  if (block_size.nx2 > 1) {
    npar *= block_size.nx2;
    dvol *= block_size.x2max - block_size.x2min;
  }
  if (block_size.nx3 > 1) {
    npar *= block_size.nx3;
    dvol *= block_size.x3max - block_size.x3min;
  }
  ppar->mass = dtog * dvol / npar;

  // Check the memory allocation.
  if (npar > ppar->nparmax) {
    std::ostringstream msg;
    msg << "### FATAL ERROR in function [MeshBlock::ProblemGenerator]" << std::endl
        << "ncells = " << npar << " > nparmax = " << ppar->nparmax << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  // Assign the particles.
  long ipar = 0, ipbase = gid * npar;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        ppar->xp(ipar) = pcoord->x1v(i);
        ppar->yp(ipar) = pcoord->x2v(j);
        ppar->zp(ipar) = pcoord->x3v(k);
        ppar->vpx(ipar) = vpx0;
        ppar->vpy(ipar) = vpy0;
        ppar->vpz(ipar) = vpz0;
        ppar->pid(ipar) = ipbase + ipar;
        ++ipar;
      }

  if (ipar != npar) {
    std::ostringstream msg;
    msg << "### FATAL ERROR in function [MeshBlock::ProblemGenerator]" << std::endl
        << "index error" << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  ppar->npar = npar;
}
