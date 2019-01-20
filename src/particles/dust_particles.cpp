//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file dust_particles.cpp
//  \brief implements functions in the DustParticles class

// Athena++ headers
#include "../athena.hpp"
#include "particles.hpp"

// Class variable initialization
bool DustParticles::initialized = false;
int DustParticles::iwx = -1, DustParticles::iwy = -1, DustParticles::iwz = -1;
int DustParticles::idpx = -1, DustParticles::idpy = -1, DustParticles::idpz = -1;

bool DustParticles::backreaction = false;
Real DustParticles::mass = 1.0, DustParticles::taus = 0.0;

//--------------------------------------------------------------------------------------
//! \fn DustParticles::Initialize()
//  \brief initializes the class.

void DustParticles::Initialize(ParameterInput *pin) {
  // Initialize first the parent class.
  if (!Particles::initialized) Particles::Initialize(pin);

  if (!initialized) {
    // Add gas velocity at each particle.
    iwx = AddWorkingArray();
    iwy = AddWorkingArray();
    iwz = AddWorkingArray();

    // Define mass.
    mass = pin->GetOrAddReal("particles", "mass", 1.0);

    // Define stopping time.
    taus = pin->GetOrAddReal("particles", "taus", 0.0);

    // Turn on/off back reaction.
    backreaction = pin->GetOrAddBoolean("particles", "backreaction", false);
    if (taus == 0.0) backreaction = false;

    if (backreaction) {
      idpx = imvpx;
      idpy = imvpy;
      idpz = imvpz;
    }

    initialized = true;
  }
}

//--------------------------------------------------------------------------------------
//! \fn DustParticles::DustParticles(MeshBlock *pmb, ParameterInput *pin)
//  \brief constructs a DustParticles instance.

DustParticles::DustParticles(MeshBlock *pmb, ParameterInput *pin)
  : Particles(pmb, pin) {
  // Assign shorthands (need to do this for every constructor of a derived class)
  AssignShorthands();

  if (backreaction) {
    dpx.InitWithShallowSlice(ppm->meshaux, 4, idpx, 1);
    dpy.InitWithShallowSlice(ppm->meshaux, 4, idpy, 1);
    dpz.InitWithShallowSlice(ppm->meshaux, 4, idpz, 1);
  }
}

//--------------------------------------------------------------------------------------
//! \fn DustParticles::~DustParticles()
//  \brief destroys a DustParticles instance.

DustParticles::~DustParticles() {
  wx.DeleteAthenaArray();
  wy.DeleteAthenaArray();
  wz.DeleteAthenaArray();

  if (backreaction) {
    dpx.DeleteAthenaArray();
    dpy.DeleteAthenaArray();
    dpz.DeleteAthenaArray();
  }
}

//--------------------------------------------------------------------------------------
//! \fn void DustParticles::AssignShorthands()
//  \brief assigns shorthands by shallow coping slices of the data.

void DustParticles::AssignShorthands() {
  Particles::AssignShorthands();
  wx.InitWithShallowSlice(work, 2, iwx, 1);
  wy.InitWithShallowSlice(work, 2, iwy, 1);
  wz.InitWithShallowSlice(work, 2, iwz, 1);
}

//--------------------------------------------------------------------------------------
//! \fn void DustParticles::AddAcceleration()
//  \brief adds acceleration to particles.

void DustParticles::AddAcceleration(Real t, Real dt, const AthenaArray<Real>& meshsrc) {
  // Interpolate gas velocity onto particles.
  if (backreaction)
    ppm->InterpolateMeshAndAssignParticles(meshsrc, IVX, work, iwx, 3,
                                           realprop, ivpx, idpx, 3);
  else
    ppm->InterpolateMeshToParticles(meshsrc, IVX, work, iwx, 3);

  // Add drag force to particles.
  if (taus > 0.0) {
    Real taus1 = 1.0 / taus;
    for (int k = 0; k < npar; ++k) {
      apx(k) -= taus1 * (vpx(k) - wx(k));
      apy(k) -= taus1 * (vpy(k) - wy(k));
      apz(k) -= taus1 * (vpz(k) - wz(k));
    }
  } else if (taus == 0.0) {
    for (int k = 0; k < npar; ++k) {
      vpx(k) = wx(k);
      vpy(k) = wy(k);
      vpz(k) = wz(k);
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn void DustParticles::AddSourceTerms(Real t, Real dt,
//                                         const AthenaArray<Real>& meshsrc)
//  \brief adds additional source terms to particles, overloaded by the user.

void __attribute__((weak)) DustParticles::AddSourceTerms(
    Real t, Real dt, const AthenaArray<Real>& meshsrc) {
}

//--------------------------------------------------------------------------------------
//! \fn void DustParticles::DepositToMesh(Real t, Real dt,
//               const AthenaArray<Real>& meshsrc, AthenaArray<Real>& meshdst);
//  \brief Deposits meshaux to Mesh.

void DustParticles::DepositToMesh(
         Real t, Real dt, const AthenaArray<Real>& meshsrc, AthenaArray<Real>& meshdst) {
  if (!backreaction) return;

  const int ias = ppm->is, jas = ppm->js, kas = ppm->ks;
  const int ibs = pmy_block->is, jbs = pmy_block->js, kbs = pmy_block->ks;
  const int nx1 = pmy_block->block_size.nx1,
            nx2 = pmy_block->block_size.nx2,
            nx3 = pmy_block->block_size.nx3;

  // Compute the momentum change.
  Real c = dt * mass / taus;
  #pragma ivdep
  for (int k = 0; k < nx3; ++k)
    for (int j = 0; j < nx2; ++j)
      for (int i = 0; i < nx1; ++i) {
        int ia = ias + i, ja = jas + j, ka = kas + k;
        int ib = ibs + i, jb = jbs + j, kb = kbs + k;
        Real w = ppm->weight(ka,ja,ia);
        dpx(ka,ja,ia) = c * (dpx(ka,ja,ia) - w * meshsrc(IVX,kb,jb,ib));
        dpy(ka,ja,ia) = c * (dpy(ka,ja,ia) - w * meshsrc(IVY,kb,jb,ib));
        dpz(ka,ja,ia) = c * (dpz(ka,ja,ia) - w * meshsrc(IVZ,kb,jb,ib));
      }

  // Deposit it to the mesh.
  ppm->DepositMeshAux(meshdst, idpx, IM1, 3);
}
