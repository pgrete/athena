//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file streaming_instability.cpp
//  \brief sets up a linear mode for the streaming instability between gas and particles.

// NOTE: In this setup, Y <-> Z.

// C++ standard libraries
#include <cmath>
#include <sstream>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../particles/particles.hpp"

// Global constants
static Real OMEGA = 1.0;

// Global parameters
static Real two_omega = 2.0 * OMEGA;
static Real omega_half = 0.5 * OMEGA;
static Real gas_accel_x = 0.0;

//======================================================================================
//! \fn void SourceTermsForGas(MeshBlock *pmb, const Real time, const Real dt,
//               const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
//               AthenaArray<Real> &cons) {
//  \brief Adds source terms to the gas.
//======================================================================================
void SourceTermsForGas(MeshBlock *pmb, const Real time, const Real dt,
         const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc,
         AthenaArray<Real> &cons) {
  // Apply the Coriolis and centrifugal forces, and linear gravity from the star, and
  // the background radial pressure gradient.
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real rho_dt = prim(IDN,k,j,i) * dt;
        cons(IM1,k,j,i) += rho_dt * (two_omega * prim(IVZ,k,j,i) + gas_accel_x);
        cons(IM3,k,j,i) -= rho_dt * omega_half * prim(IVX,k,j,i);
      }
    }
  }
}

//======================================================================================
//! \fn void DustParticles::AddSourceTerms(Real t, Real dt,
//                              const AthenaArray<Real>& meshsrc)
//  \brief Adds source terms to the particles.
//======================================================================================
void DustParticles::AddSourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc) {
  // Preprocess the constants.
  Real dt_two_omega = dt * two_omega,
       dt_omega_half = dt * omega_half;

  // Apply the Coriolis and centrifugal forces, and linear gravity from the star.
  for (int k = 0; k < npar; ++k) {
    vpx(k) += dt_two_omega * vpz0(k);
    vpz(k) -= dt_omega_half * vpx0(k);
  }
}

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Initializes problem-specific data in Mesh class.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Preprocess constants.
  Real cs0 = pin->GetReal("hydro", "iso_sound_speed");
  Real omega = pin->GetOrAddReal("problem", "omega", OMEGA);
  Real dux0 = pin->GetReal("problem", "dux0");
  two_omega = 2.0 * omega;
  omega_half = 0.5 * omega;
  gas_accel_x = 2.0 * dux0 * cs0 * omega;

  // Enroll source terms.
  EnrollUserExplicitSourceFunction(SourceTermsForGas);
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief Sets the initial conditions.
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Get the dust-to-gas density ratio.
  Real dtog = pin->GetReal("problem", "dtog");

  // Get the stopping time.
  Real taus = pin->GetReal("particles", "taus");
  Real omega = pin->GetOrAddReal("problem", "omega", OMEGA);
  taus *= omega;

  // Find the Nakagawa-Sekiya-Hayashi (1986) equilibrium solution.
  Real dux0 = pin->GetReal("problem", "dux0");
  Real cs0 = pin->GetReal("hydro", "iso_sound_speed");
  dux0 *= cs0 / (std::pow(1.0 + dtog, 2) + std::pow(taus, 2));
  Real ux0 = 2.0 * dtog * taus * dux0,
       uy0 = -((1.0 + dtog) + std::pow(taus, 2)) * dux0,
       vpx0 = -2.0 * taus * dux0,
       vpy0 = -(1.0 + dtog) * dux0;

  // Set a uniform, steady gas.
  for (int k = ks; k <= ke; ++k) {
    for (int j = js; j <= je; ++j) {
      for (int i = is; i <= ie; ++i) {
        phydro->u(IDN,k,j,i) = 1.0;
        phydro->u(IM1,k,j,i) = ux0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = uy0;
      }
    }
  }

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
    ix1 = static_cast<int>(std::round((block_size.x1min - mesh_size.x1min) / dx1));
    npx1_loc = static_cast<int>(std::round((block_size.x1max - block_size.x1min) / dx1));
  }
  if (block_size.nx2 > 1) {
    ix2 = static_cast<int>(std::round((block_size.x2min - mesh_size.x2min) / dx2));
    npx2_loc = static_cast<int>(std::round((block_size.x2max - block_size.x2min) / dx2));
  }
  if (block_size.nx3 > 1) {
    ix3 = static_cast<int>(std::round((block_size.x3min - mesh_size.x3min) / dx3));
    npx3_loc = static_cast<int>(std::round((block_size.x3max - block_size.x3min) / dx3));
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
        ppar->vpy(ipar) = 0.0;
        ppar->vpz(ipar) = vpy0;
        ppar->pid(ipar) = id++;
        ++ipar;
      }
      id += npx1 - npx1_loc;
    }
    id += npx1 * (npx2 - npx2_loc);
  }
  ppar->npar = ipar;
}
