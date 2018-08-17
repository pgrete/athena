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

  // Find the total number of particles in each direction.
  RegionSize& mesh_size = pmy_mesh->mesh_size;
  int npx1 = (block_size.nx1 > 1) ?
                  pin->GetOrAddInteger("problem", "npx1", mesh_size.nx1) : 1,
      npx2 = (block_size.nx2 > 1) ?
                  pin->GetOrAddInteger("problem", "npx2", mesh_size.nx2) : 1,
      npx3 = (block_size.nx3 > 1) ?
                  pin->GetOrAddInteger("problem", "npx3", mesh_size.nx3) : 1;

  // Find the mass of each particle and the distance between adjacent particles.
  Real vol = 1.0;
  Real dx1 = 0.0, dx2 = 0.0, dx3 = 0.0, length;

  length = mesh_size.x1max - mesh_size.x1min;
  dx1 = length / npx1;
  if (block_size.nx1 > 1) vol *= length;

  length = mesh_size.x2max - mesh_size.x2min;
  dx2 = length / npx2;
  if (block_size.nx2 > 1) vol *= length;

  length = mesh_size.x3max - mesh_size.x3min;
  dx3 = length / npx3;
  if (block_size.nx3 > 1) vol *= length;

  ppar->mass = dtog * vol / (npx1 * npx2 * npx3);

  // Find the local number of particles and their beginning index.
  int ix1 = 0, ix2 = 0, ix3 = 0;
  int npx1_loc = 1, npx2_loc = 1, npx3_loc = 1;
  if (block_size.nx1 > 1) {
    ix1 = lround((block_size.x1min - mesh_size.x1min) / dx1);
    npx1_loc = lround((block_size.x1max - block_size.x1min) / dx1);
  }
  if (block_size.nx2 > 1) {
    ix2 = lround((block_size.x2min - mesh_size.x2min) / dx2);
    npx2_loc = lround((block_size.x2max - block_size.x2min) / dx2);
  }
  if (block_size.nx3 > 1) {
    ix3 = lround((block_size.x3min - mesh_size.x3min) / dx3);
    npx3_loc = lround((block_size.x3max - block_size.x3min) / dx3);
  }

  // Check the memory allocation.
  if (npx1_loc * npx2_loc * npx3_loc > ppar->nparmax) {
    std::ostringstream msg;
    msg << "### FATAL ERROR in function [MeshBlock::ProblemGenerator]" << std::endl
        << "nparmax is too small" << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  // Assign the particles.
  int ipar = 0, id = ix1 + npx1 * (ix2 + npx2 * ix3);
  for (int k = 0; k < npx3_loc; ++k) {
    Real zp1 = block_size.x3min + (k + 0.5) * dx3;
    for (int j = 0; j < npx2_loc; ++j) {
      Real yp1 = block_size.x2min + (j + 0.5) * dx2;
      for (int i = 0; i < npx1_loc; ++i) {
        Real xp1 = block_size.x1min + (i + 0.5) * dx1;
        ppar->xp(ipar) = xp1;
        ppar->yp(ipar) = yp1;
        ppar->zp(ipar) = zp1;
        ppar->vpx(ipar) = vpx0;
        ppar->vpy(ipar) = vpy0;
        ppar->vpz(ipar) = vpz0;
        ppar->pid(ipar) = id++;
        ++ipar;
      }
      id += npx1 - npx1_loc;
    }
    id += npx1 * (npx2 - npx2_loc);
  }
  ppar->npar = ipar;
}
