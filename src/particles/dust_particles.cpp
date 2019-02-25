//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file dust_particles.cpp
//  \brief implements functions in the DustParticles class

// C++ headers
#include <algorithm>  // min()

// Athena++ headers
#include "../athena.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "particles.hpp"

// Class variable initialization
bool DustParticles::initialized = false;
int DustParticles::iwx = -1, DustParticles::iwy = -1, DustParticles::iwz = -1;
int DustParticles::idpx1 = -1, DustParticles::idpx2 = -1, DustParticles::idpx3 = -1;

bool DustParticles::backreaction = false;
Real DustParticles::mass = 1.0, DustParticles::taus = 0.0;

//--------------------------------------------------------------------------------------
//! \fn void DustParticles::Initialize(Mesh *pm, ParameterInput *pin)
//  \brief initializes the class.

void DustParticles::Initialize(Mesh *pm, ParameterInput *pin) {
  // Initialize first the parent class.
  Particles::Initialize(pm, pin);

  if (!initialized) {
    // Add working array at particles for gas velocity/particle momentum change.
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
      idpx1 = imvpx;
      idpx2 = imvpy;
      idpx3 = imvpz;
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
    dpx1.InitWithShallowSlice(ppm->meshaux, 4, idpx1, 1);
    dpx2.InitWithShallowSlice(ppm->meshaux, 4, idpx2, 1);
    dpx3.InitWithShallowSlice(ppm->meshaux, 4, idpx3, 1);
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
    dpx1.DeleteAthenaArray();
    dpx2.DeleteAthenaArray();
    dpx3.DeleteAthenaArray();
  }
}

//--------------------------------------------------------------------------------------
//! \fn Real DustParticles::NewBlockTimeStep();
//  \brief returns the time step required by particles in the block.

Real DustParticles::NewBlockTimeStep() {
  // Run first the parent class.
  Real dt = Particles::NewBlockTimeStep();

  // Nothing to do for tracer particles.
  if (taus <= 0.0) return dt;

  Real epsmax = 0;
  if (backreaction) {
    // Find the maximum local solid-to-gas density ratio.
    Coordinates *pc = pmy_block->pcoord;
    Hydro *phydro = pmy_block->phydro;
    const int is = ppm->is, js = ppm->js, ks = ppm->ks;
    const int ie = ppm->ie, je = ppm->je, ke = ppm->ke;
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          Real epsilon = ppm->weight(k,j,i) / (
                         pc->GetCellVolume(k,j,i) * phydro->u(IDN,k,j,i));
          epsmax = std::max(epsmax, epsilon);
        }
    epsmax *= mass;
  }

  // Return the drag timescale.
  return std::min(dt, static_cast<Real>(cfl_par * taus / (1.0 + epsmax)));
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
//! \fn void DustParticles::SourceTerms()
//  \brief adds acceleration to particles.

void DustParticles::SourceTerms(Real t, Real dt, const AthenaArray<Real>& meshsrc) {
  // Transform the particle velocity in mesh coordinates.
  const Coordinates *pc = pmy_block->pcoord;
  for (int k = 0; k < npar; ++k)
    pc->CartesianToMeshCoordsVector(xp(k), yp(k), zp(k), vpx(k), vpy(k), vpz(k),
                                    vpx(k), vpy(k), vpz(k));

  // Communicate gas and particle velocities.
  if (backreaction)
    ppm->InterpolateMeshAndAssignParticles(meshsrc, IVX, work, iwx, 3,
                                           realprop, ivpx, idpx1, 3);
  else
    ppm->InterpolateMeshToParticles(meshsrc, IVX, work, iwx, 3);

  if (taus > 0.0) {
    // Compute the drag force.
    Real c = dt / taus;
    for (int k = 0; k < npar; ++k) {
      wx(k) = c * (wx(k) - vpx(k));
      wy(k) = c * (wy(k) - vpy(k));
      wz(k) = c * (wz(k) - vpz(k));
    }
  }

  // Transform the drag force or the gas velocity back in Cartesian.
  for (int k = 0; k < npar; ++k) {
    Real x1, x2, x3;
    pc->CartesianToMeshCoords(xp(k), yp(k), zp(k), x1, x2, x3);
    pc->MeshCoordsToCartesianVector(x1, x2, x3, wx(k), wy(k), wz(k), wx(k), wy(k), wz(k));
  }

  if (taus > 0.0) {
    // Add the drag force to particles.
    for (int k = 0; k < npar; ++k) {
      vpx(k) = vpx0(k) + wx(k);
      vpy(k) = vpy0(k) + wy(k);
      vpz(k) = vpz0(k) + wz(k);
    }
  } else if (taus == 0.0) {
    // Assign gas velocity to particle velocity for tracer particles.
    for (int k = 0; k < npar; ++k) {
      vpx(k) = wx(k);
      vpy(k) = wy(k);
      vpz(k) = wz(k);
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn void DustParticles::UserSourceTerms(Real t, Real dt,
//                                          const AthenaArray<Real>& meshsrc)
//  \brief adds additional source terms to particles, overloaded by the user.

void __attribute__((weak)) DustParticles::UserSourceTerms(
    Real t, Real dt, const AthenaArray<Real>& meshsrc) {
}

//--------------------------------------------------------------------------------------
//! \fn void DustParticles::DepositToMesh(Real t, Real dt,
//               const AthenaArray<Real>& meshsrc, AthenaArray<Real>& meshdst);
//  \brief Deposits meshaux to Mesh.

void DustParticles::DepositToMesh(
         Real t, Real dt, const AthenaArray<Real>& meshsrc, AthenaArray<Real>& meshdst) {
  if (!backreaction) return;

  // Compute the momentum change.
  const Real c = dt * mass / taus;
  const int is = ppm->is, js = ppm->js, ks = ppm->ks;
  const int ie = ppm->ie, je = ppm->je, ke = ppm->ke;
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        Real w = ppm->weight(k,j,i);
        dpx1(k,j,i) = c * (dpx1(k,j,i) - w * meshsrc(IVX,k,j,i));
        dpx2(k,j,i) = c * (dpx2(k,j,i) - w * meshsrc(IVY,k,j,i));
        dpx3(k,j,i) = c * (dpx3(k,j,i) - w * meshsrc(IVZ,k,j,i));
      }

  // Deposit particle momentum changes to the gas.
  ppm->DepositMeshAux(meshdst, idpx1, IM1, 3);
}
