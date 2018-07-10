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
int DustParticles::iwx = -1, DustParticles::iwy = -1, DustParticles::iwz = -1;
AthenaArray<int> DustParticles::imeshsrc, DustParticles::iwork;

bool DustParticles::backreaction = false;
Real DustParticles::mass = 1.0, DustParticles::taus = 0.0;

//--------------------------------------------------------------------------------------
//! \fn DustParticles::Initialize()
//  \brief initializes the class.

void DustParticles::Initialize(ParameterInput *pin)
{
  // Initialize first the parent class.
  if (!Particles::initialized) Particles::Initialize(pin);

  if (!initialized) {
    // Add gas velocity at each particle.
    iwx = AddWorkingArray();
    iwy = AddWorkingArray();
    iwz = AddWorkingArray();

    // Assign indice mapping for particle-mesh.
    imeshsrc.NewAthenaArray(3);
    iwork.NewAthenaArray(3);
    imeshsrc(0) = IVX;
    imeshsrc(1) = IVY;
    imeshsrc(2) = IVZ;
    iwork(0) = iwx;
    iwork(1) = iwy;
    iwork(2) = iwz;

    // Define mass.
    mass = pin->GetOrAddReal("particles", "mass", 1.0);

    // Define stopping time.
    taus = pin->GetOrAddReal("particles", "taus", 0.0);

    // Turn on/off back reaction.
    backreaction = pin->GetOrAddBoolean("particles", "backreaction", false);
    if (taus == 0.0) backreaction = false;

    initialized = true;
  }
}

//--------------------------------------------------------------------------------------
//! \fn DustParticles::DustParticles(MeshBlock *pmb, ParameterInput *pin)
//  \brief constructs a DustParticles instance.

DustParticles::DustParticles(MeshBlock *pmb, ParameterInput *pin)
  : Particles(pmb, pin)
{
  // Assign shorthands (need to do this for every constructor of a derived class)
  AssignShorthands();
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
  wx.InitWithShallowSlice(work, 2, iwx, 1);
  wy.InitWithShallowSlice(work, 2, iwy, 1);
  wz.InitWithShallowSlice(work, 2, iwz, 1);
}

//--------------------------------------------------------------------------------------
//! \fn void DustParticles::AddAcceleration()
//  \brief adds acceleration to particles.

void DustParticles::AddAcceleration(Real t, Real dt, const AthenaArray<Real>& meshsrc)
{
  // Interpolate gas velocity onto particles.
  ppm->InterpolateMeshToParticles(meshsrc, imeshsrc, work, iwork);

  // Add drag force to particles.
  if (taus > 0.0) {
    Real taus1 = 1.0 / taus;
    for (long k = 0; k < npar; ++k) {
      wx(k) = taus1 * (vpx(k) - wx(k));
      wy(k) = taus1 * (vpy(k) - wy(k));
      wz(k) = taus1 * (vpz(k) - wz(k));
      apx(k) -= wx(k);
      apy(k) -= wy(k);
      apz(k) -= wz(k);
    }
  } else if (taus == 0.0) {
    for (long k = 0; k < npar; ++k) {
      vpx(k) = wx(k);
      vpy(k) = wy(k);
      vpz(k) = wz(k);
    }
  }
}
