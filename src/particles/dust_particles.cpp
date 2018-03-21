//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file dust_particles.cpp
//  \brief implements functions in the DustParticles class

#include "particles.hpp"

// Class variable initialization
bool DustParticles::initialized = false;
int DustParticles::iux = -1, DustParticles::iuy = -1, DustParticles::iuz = -1;

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
}
