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
AthenaArray<int> DustParticles::imeshsrc, DustParticles::iwork,
                 DustParticles::irealprop, DustParticles::imeshaux,
                 DustParticles::imeshdst;

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

    if (backreaction) {
      idpx = ParticleMesh::AddMeshAux();
      idpy = ParticleMesh::AddMeshAux();
      idpz = ParticleMesh::AddMeshAux();

      irealprop.NewAthenaArray(3);
      irealprop(0) = ivpx;
      irealprop(1) = ivpy;
      irealprop(2) = ivpz;

      imeshaux.NewAthenaArray(3);
      imeshaux(0) = idpx;
      imeshaux(1) = idpy;
      imeshaux(2) = idpz;

      imeshdst.NewAthenaArray(3);
      imeshdst(0) = IM1;
      imeshdst(1) = IM2;
      imeshdst(2) = IM3;
    }

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

  if (backreaction) {
    dpx.InitWithShallowSlice(ppm->meshaux, 4, idpx, 1);
    dpy.InitWithShallowSlice(ppm->meshaux, 4, idpy, 1);
    dpz.InitWithShallowSlice(ppm->meshaux, 4, idpz, 1);
  }
}

//--------------------------------------------------------------------------------------
//! \fn DustParticles::~DustParticles()
//  \brief destroys a DustParticles instance.

DustParticles::~DustParticles()
{
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
  if (backreaction)
    ppm->InterpolateMeshAndAssignParticles(meshsrc, IVX, IVZ, work, iwx,
                                           realprop, ivpx, ivpz, idpx);
  else
    ppm->InterpolateMeshToParticles(meshsrc, IVX, IVZ, work, iwx);

  // Add drag force to particles.
  if (taus > 0.0) {
    Real taus1 = 1.0 / taus;
    for (long k = 0; k < npar; ++k) {
      apx(k) -= taus1 * (vpx(k) - wx(k));
      apy(k) -= taus1 * (vpy(k) - wy(k));
      apz(k) -= taus1 * (vpz(k) - wz(k));
    }
  } else if (taus == 0.0) {
    for (long k = 0; k < npar; ++k) {
      vpx(k) = wx(k);
      vpy(k) = wy(k);
      vpz(k) = wz(k);
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn void DustParticles::DepositToMesh(Real t, Real dt,
//               const AthenaArray<Real>& meshsrc, AthenaArray<Real>& meshdst);
//  \brief Deposits meshaux to Mesh.

void DustParticles::DepositToMesh(
         Real t, Real dt, const AthenaArray<Real>& meshsrc, AthenaArray<Real>& meshdst)
{
  if (!backreaction) return;

  // Compute the momentum change.
  Real c = dt * mass / taus;
  for (int ka = ppm->ks, kb = pmy_block->ks; ka <= ppm->ke; ++ka, ++kb)
    for (int ja = ppm->js, jb = pmy_block->js; ja <= ppm->je; ++ja, ++jb)
      for (int ia = ppm->is, ib = pmy_block->is; ia <= ppm->ie; ++ia, ++ib) {
        Real w = ppm->weight(ka,ja,ia);
        dpx(ka,ja,ia) = c * (dpx(ka,ja,ia) - w * meshsrc(IVX,kb,jb,ib));
        dpy(ka,ja,ia) = c * (dpy(ka,ja,ia) - w * meshsrc(IVY,kb,jb,ib));
        dpz(ka,ja,ia) = c * (dpz(ka,ja,ia) - w * meshsrc(IVZ,kb,jb,ib));
      }

  // Deposit it to the mesh.
  ppm->DepositMeshAux(meshdst, idpx, idpz, IM1);
}
