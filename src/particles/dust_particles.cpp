//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file dust_particles.cpp
//  \brief implements functions in the DustParticles class

#include "particles.hpp"
#include "../hydro/hydro.hpp"

// Class variable initialization
bool DustParticles::initialized = false;
int DustParticles::iux = -1, DustParticles::iuy = -1, DustParticles::iuz = -1;
AthenaArray<int> DustParticles::pm_meshindices, DustParticles::pm_auxindices;

//--------------------------------------------------------------------------------------
//! \fn DustParticles::Initialize()
//  \brief initializes the class.

void DustParticles::Initialize()
{
  // Initialize first the parent class.
  if (!Particles::initialized) Particles::Initialize();

  // Add gas velocity at each particle.
  iux = AddAuxProperty();
  iuy = AddAuxProperty();
  iuz = AddAuxProperty();

  // Assign indice mapping for particle-mesh.
  pm_meshindices.NewAthenaArray(3);
  pm_auxindices.NewAthenaArray(3);
  pm_meshindices(0) = IVX;
  pm_meshindices(1) = IVY;
  pm_meshindices(2) = IVZ;
  pm_auxindices(0) = iux;
  pm_auxindices(1) = iuy;
  pm_auxindices(2) = iuz;

  initialized = true;
}

//--------------------------------------------------------------------------------------
//! \fn DustParticles::DustParticles(MeshBlock *pmb, ParameterInput *pin)
//  \brief constructs a DustParticles instance.

DustParticles::DustParticles(MeshBlock *pmb, ParameterInput *pin)
  : Particles(pmb, pin)
{
  // Initialize the class when first called.
  if (!initialized) Initialize();

  // Define mass.
  mass = pin->GetOrAddInteger("particles", "mass", 1);

  // Define stopping time.
  taus = pin->GetOrAddInteger("particles", "taus", 1);
}

//--------------------------------------------------------------------------------------
//! \fn DustParticles::~DustParticles()
//  \brief destroys a DustParticles instance.

DustParticles::~DustParticles()
{
}

//--------------------------------------------------------------------------------------
//! \fn void DustParticles::AssignShorthands()
//  \brief assigns shorthands by shallow coping slices of the data.

void DustParticles::AssignShorthands()
{
  Particles::AssignShorthands();
  ux.InitWithShallowSlice(auxprop, 2, iux, 1);
  uy.InitWithShallowSlice(auxprop, 2, iuy, 1);
  uz.InitWithShallowSlice(auxprop, 2, iuz, 1);
}

//--------------------------------------------------------------------------------------
//! \fn void DustParticles::AddAcceleration()
//  \brief adds acceleration to particles.

void DustParticles::AddAcceleration(Real t, Real dt)
{
  // Interpolate gas velocity onto particles.
  InterpolateMeshToParticles(pmy_block->phydro->w, pm_meshindices, pm_auxindices);

  // Add drag force to particles.
  if (taus == 0.0) {
    for (long k = 0; k < npar; ++k) {
      vpx(k) = ux(k);
      vpy(k) = uy(k);
      vpz(k) = uz(k);
    }
  } else {
    Real taus1 = 1.0 / taus;
    for (long k = 0; k < npar; ++k) {
      apx(k) += taus1 * (ux(k) - vpx(k));
      apy(k) += taus1 * (uy(k) - vpy(k));
      apz(k) += taus1 * (uz(k) - vpz(k));
    }
  }
}
