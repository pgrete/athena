//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particles.cpp
//  \brief implements functions in particle classes

// C++ Standard Libraries
#include <cmath>
#include <sstream>
#include <string>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../globals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "particles.hpp"

// Class variable initialization
bool Particles::initialized = false;
int Particles::nint = 0;
int Particles::nreal = 0;
int Particles::naux = 0;
int Particles::nwork = 0;
int Particles::ipid = -1;
int Particles::ixp = -1, Particles::iyp = -1, Particles::izp = -1;
int Particles::ivpx = -1, Particles::ivpy = -1, Particles::ivpz = -1;
int Particles::ixp0 = -1, Particles::iyp0 = -1, Particles::izp0 = -1;
int Particles::ivpx0 = -1, Particles::ivpy0 = -1, Particles::ivpz0 = -1;
int Particles::ixi1 = -1, Particles::ixi2 = -1, Particles::ixi3 = -1;
int Particles::iapx = -1, Particles::iapy = -1, Particles::iapz = -1;

// Local function prototypes
static void _CartesianToMeshCoords(Real x, Real y, Real z, Real& x1, Real& x2, Real& x3);
static void _MeshCoordsToCartesian(Real x1, Real x2, Real x3, Real& x, Real& y, Real& z);
static void _MeshCoordsToIndices(MeshBlock *pmb, Real x1, Real x2, Real x3,
                                 Real& xi1, Real& xi2, Real& xi3);
static void _IndicesToMeshCoords(MeshBlock *pmb, Real xi1, Real xi2, Real xi3,
                                 Real& x1, Real& x2, Real& x3);
static int CheckSide(int xi, int xi1, int xi2);

//--------------------------------------------------------------------------------------
//! \fn Particles::Initialize(ParameterInput *pin)
//  \brief initializes the class.

void Particles::Initialize(ParameterInput *pin)
{
  if (initialized) return;

  // Add particle ID.
  ipid = AddIntProperty();

  // Add particle position.
  ixp = AddRealProperty();
  iyp = AddRealProperty();
  izp = AddRealProperty();

  // Add particle velocity.
  ivpx = AddRealProperty();
  ivpy = AddRealProperty();
  ivpz = AddRealProperty();

  // Add old particle position.
  ixp0 = AddAuxProperty();
  iyp0 = AddAuxProperty();
  izp0 = AddAuxProperty();

  // Add old particle velocity.
  ivpx0 = AddAuxProperty();
  ivpy0 = AddAuxProperty();
  ivpz0 = AddAuxProperty();

  // Add particle position indices.
  ixi1 = AddWorkingArray();
  ixi2 = AddWorkingArray();
  ixi3 = AddWorkingArray();

  // Add acceleration components.
  iapx = AddWorkingArray();
  iapy = AddWorkingArray();
  iapz = AddWorkingArray();

  // Initiate ParticleMesh class.
  ParticleMesh::Initialize(pin);

  initialized = true;
}

//--------------------------------------------------------------------------------------
//! \fn Particles::Particles(MeshBlock *pmb, ParameterInput *pin)
//  \brief constructs a Particles instance.

Particles::Particles(MeshBlock *pmb, ParameterInput *pin)
{
  // Point to the calling MeshBlock.
  pmy_block = pmb;
  pmy_mesh = pmb->pmy_mesh;
  pbval_ = pmb->pbval;
  nparmax = pin->GetOrAddInteger("particles", "nparmax", 1);
  npar = 0;

  // Check active dimensions.
  active1_ = pmy_mesh->mesh_size.nx1 > 1;
  active2_ = pmy_mesh->mesh_size.nx2 > 1;
  active3_ = pmy_mesh->mesh_size.nx3 > 1;

  // Allocate integer properties.
  intprop.NewAthenaArray(nint,nparmax);

  // Allocate integer properties.
  realprop.NewAthenaArray(nreal,nparmax);

  // Allocate auxiliary properties.
  if (naux > 0) auxprop.NewAthenaArray(naux,nparmax);

  // Allocate working arrays.
  if (nwork > 0) work.NewAthenaArray(nwork,nparmax);

  // Allocate mesh auxiliaries.
  ppm = new ParticleMesh(this);

  // Shallow copy to shorthands.
  AssignShorthands();

  // Initiate ParticleBuffer class.
  ParticleBuffer::SetNumberOfProperties(nint, nreal + naux);
}

//--------------------------------------------------------------------------------------
//! \fn Particles::~Particles()
//  \brief destroys a Particles instance.

Particles::~Particles()
{
  // Delete integer properties.
  intprop.DeleteAthenaArray();

  // Delete real properties.
  realprop.DeleteAthenaArray();

  // Delete auxiliary properties.
  if (naux > 0) auxprop.DeleteAthenaArray();

  // Delete working arrays.
  if (nwork > 0) work.DeleteAthenaArray();

  // Delete mesh auxiliaries.
  delete ppm;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::ApplyBoundaryConditions(long k, Real &x1, Real &x2, Real &x3)
//  \brief applies boundary conditions to particle k and returns its updated mesh
//         coordinates (x1,x2,x3).

void Particles::ApplyBoundaryConditions(long k, Real &x1, Real &x2, Real &x3)
{
  bool flag = false;
  RegionSize& mesh_size = pmy_mesh->mesh_size;

  // Find the mesh coordinates.
  Real x10, x20, x30;
  _IndicesToMeshCoords(pmy_block, xi1(k), xi2(k), xi3(k), x1, x2, x3);
  _CartesianToMeshCoords(xp0(k), yp0(k), zp0(k), x10, x20, x30);

  if (active1_) {
    if (x1 < mesh_size.x1min) {
      // Inner x1
      if (pmy_mesh->mesh_bcs[INNER_X1] == PERIODIC_BNDRY) {
        x1 += mesh_size.x1len;
        x10 += mesh_size.x1len;
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::ApplyBoundaryConditions]"
            << std::endl
            << "Non-periodic boundary for inner x1 not supported. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      flag = true;
    } else if (x1 >= mesh_size.x1max) {
      // Outer x1
      if (pmy_mesh->mesh_bcs[OUTER_X1] == PERIODIC_BNDRY) {
        x1 -= mesh_size.x1len;
        x10 -= mesh_size.x1len;
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::ApplyBoundaryConditions]"
            << std::endl
            << "Non-periodic boundary for outer x1 not supported. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      flag = true;
    }
  }

  if (active2_) {
    if (x2 < mesh_size.x2min) {
      // Inner x2
      if (pmy_mesh->mesh_bcs[INNER_X2] == PERIODIC_BNDRY) {
        x2 += mesh_size.x2len;
        x20 += mesh_size.x2len;
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::ApplyBoundaryConditions]"
            << std::endl
            << "Non-periodic boundary for inner x2 not supported. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      flag = true;
    } else if (x2 >= mesh_size.x2max) {
      // Outer x2
      if (pmy_mesh->mesh_bcs[OUTER_X2] == PERIODIC_BNDRY) {
        x2 -= mesh_size.x2len;
        x20 -= mesh_size.x2len;
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::ApplyBoundaryConditions]"
            << std::endl
            << "Non-periodic boundary for outer x2 not supported. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      flag = true;
    }
  }

  if (active3_) {
    if (x3 < mesh_size.x3min) {
      // Inner x3
      if (pmy_mesh->mesh_bcs[INNER_X3] == PERIODIC_BNDRY) {
        x3 += mesh_size.x3len;
        x30 += mesh_size.x3len;
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::ApplyBoundaryConditions]"
            << std::endl
            << "Non-periodic boundary for inner x3 not supported. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      flag = true;
    } else if (x3 >= mesh_size.x3max) {
      // Outer x3
      if (pmy_mesh->mesh_bcs[OUTER_X3] == PERIODIC_BNDRY) {
        x3 -= mesh_size.x3len;
        x30 -= mesh_size.x3len;
      } else {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::ApplyBoundaryConditions]"
            << std::endl
            << "Non-periodic boundary for outer x3 not supported. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      flag = true;
    }
  }

  if (flag) {
    _MeshCoordsToCartesian(x1, x2, x3, xp(k), yp(k), zp(k));
    _MeshCoordsToCartesian(x10, x20, x30, xp0(k), yp0(k), zp0(k));
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::GetPositionIndices(MeshBlock *pmb, long npar,
//                                         const AthenaArray<Real>& xp,
//                                         const AthenaArray<Real>& yp,
//                                         const AthenaArray<Real>& zp,
//                                         AthenaArray<Real>& xi1,
//                                         AthenaArray<Real>& xi2,
//                                         AthenaArray<Real>& xi3)
//  \brief finds the position indices of each particle with respect to the local grid.

void Particles::GetPositionIndices(MeshBlock *pmb, long npar,
                                   const AthenaArray<Real>& xp,
                                   const AthenaArray<Real>& yp,
                                   const AthenaArray<Real>& zp,
                                   AthenaArray<Real>& xi1,
                                   AthenaArray<Real>& xi2,
                                   AthenaArray<Real>& xi3)
{
  for (long k = 0; k < npar; ++k) {
    // Convert to the Mesh coordinates.
    Real x1, x2, x3;
    _CartesianToMeshCoords(xp(k), yp(k), zp(k), x1, x2, x3);

    // Convert to the index space.
    _MeshCoordsToIndices(pmb, x1, x2, x3, xi1(k), xi2(k), xi3(k));
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::SendParticlesAndMesh(int step)
//  \brief send particles and meshaux near boundaries to neighbors.

void Particles::SendParticlesAndMesh(int step)
{
  // Send particles.
  if (npar > 0) {
    SetPositionIndices();
    SendToNeighbors();
  }

  // Send MeshAux boundary.
  if (step > 0 && ppm->nmeshaux > 0)
    ppm->SendBoundary();
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::ReceiveParticlesAndMesh(int step)
//  \brief receives particles and meshaux near boundaries from neighbors.

void Particles::ReceiveParticlesAndMesh(int step)
{
  // Receive particles from neighbor blocks.
  ReceiveFromNeighbors();

  // Flush ParticleMesh receive buffers and deposit MeshAux to MeshBlock.
  if (ppm->nmeshaux > 0) {
    ppm->ReceiveBoundary();

    Hydro *phydro = pmy_block->phydro;
    Real t, dt;

    switch (step) {

    case 1:
      t = pmy_mesh->time;
      dt = 0.5 * pmy_mesh->dt;
      DepositToMesh(t, dt, phydro->u, phydro->u1);
      break;

    case 2:
      t = pmy_mesh->time + 0.5 * pmy_mesh->dt;
      dt = pmy_mesh->dt;
      DepositToMesh(t, dt, phydro->u1, phydro->u);
      break;
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::EulerStep(Real t, Real dt, const AthenaArray<Real>& meshsrc)
//  \brief evolves the particle positions and velocities by one Euler step.

void Particles::EulerStep(Real t, Real dt, const AthenaArray<Real>& meshsrc)
{
  // Get the accelerations.
  for (long k = 0; k < npar; ++k) {
    apx(k) = 0.0;
    apy(k) = 0.0;
    apz(k) = 0.0;
  }
  AddAcceleration(t, dt, meshsrc);

  // Update the positions and velocities **from the beginning of the time step**.
  for (long k = 0; k < npar; ++k) {
    xp(k) = xp0(k) + dt * vpx(k);
    yp(k) = yp0(k) + dt * vpy(k);
    zp(k) = zp0(k) + dt * vpz(k);
    vpx(k) = vpx0(k) + dt * apx(k);
    vpy(k) = vpy0(k) + dt * apy(k);
    vpz(k) = vpz0(k) + dt * apz(k);
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::SaveStatus()
//  \brief saves the current positions and velocities for later use.

void Particles::SaveStatus()
{
  for (long k = 0; k < npar; ++k) {
    // Save current positions.
    xp0(k) = xp(k);
    yp0(k) = yp(k);
    zp0(k) = zp(k);

    // Save current velocities.
    vpx0(k) = vpx(k);
    vpy0(k) = vpy(k);
    vpz0(k) = vpz(k);
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::SendToNeighbors()
//  \brief sends particles outside boundary to the buffers of neighboring meshblocks.

void Particles::SendToNeighbors()
{
  const int IS = pmy_block->is;
  const int IE = pmy_block->ie;
  const int JS = pmy_block->js;
  const int JE = pmy_block->je;
  const int KS = pmy_block->ks;
  const int KE = pmy_block->ke;

  // TODO: Currently only works for Cartesian.
  if (COORDINATE_SYSTEM != "cartesian") {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [Particles::SendToNeighbors]" << std::endl
        << "Non-Cartesian coordinates not yet implemented. " << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  for (long k = 0; k < npar; ) {
    // Check if a particle is outside the boundary.
    int xi1i = int(xi1(k)), xi2i = int(xi2(k)), xi3i = int(xi3(k));
    int ox1 = active1_ ? CheckSide(xi1i, IS, IE) : 0,
        ox2 = active2_ ? CheckSide(xi2i, JS, JE) : 0,
        ox3 = active3_ ? CheckSide(xi3i, KS, KE) : 0;
    if (ox1 == 0 && ox2 == 0 && ox3 == 0) {
      ++k;
      continue;
    }

    // Apply boundary conditions and find the mesh coordinates.
    Real x1, x2, x3;
    ApplyBoundaryConditions(k, x1, x2, x3);

    // Find the neighbor MeshBlock to send it to.
    Neighbor *pn = FindTargetNeighbor(ox1, ox2, ox3, xi1i, xi2i, xi3i);
    if (pn == NULL) {
      std::stringstream msg;
      msg << "### FATAL ERROR in function [Particles::SendToNeighbors]" << std::endl
          << "cannot find the neighbor block to send the particle to. " << std::endl;
      throw std::runtime_error(msg.str().data());
      continue;
    }
    MeshBlock *pnmb = pmy_mesh->FindMeshBlock(pn->pnb->gid);
    Particles *pnp = pnmb->ppar;

    // No need to send if back to the same block.
    if (pnmb == pmy_block) {
      _MeshCoordsToIndices(pmy_block, x1, x2, x3, xi1(k), xi2(k), xi3(k));
      ++k;
      continue;
    }

    // Check the buffer size of the target MeshBlock.
    ParticleBuffer& nrecv = pnp->recv_[pn->pnb->targetid];
    if (nrecv.npar >= nrecv.nparmax)
      nrecv.Reallocate((nrecv.nparmax > 0) ? 2 * nrecv.nparmax : 1);

    // Copy the properties of the particle to the neighbor.
    long *pi = nrecv.ibuf + ParticleBuffer::nint * nrecv.npar;
    for (int j = 0; j < nint; ++j)
      *pi++ = intprop(j,k);
    Real *pr = nrecv.rbuf + ParticleBuffer::nreal * nrecv.npar;
    for (int j = 0; j < nreal; ++j)
      *pr++ = realprop(j,k);
    for (int j = 0; j < naux; ++j)
      *pr++ = auxprop(j,k);
    ++nrecv.npar;

    // Pop the particle from the current MeshBlock.
    if (--npar != k) {
      xi1(k) = xi1(npar);
      xi2(k) = xi2(npar);
      xi3(k) = xi3(npar);
      for (int j = 0; j < nint; ++j)
        intprop(j,k) = intprop(j,npar);
      for (int j = 0; j < nreal; ++j)
        realprop(j,k) = realprop(j,npar);
      for (int j = 0; j < naux; ++j)
        auxprop(j,k) = auxprop(j,npar);
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::ReceiveFromNeighbors()
//  \brief receives particles from neighboring meshblocks.

void Particles::ReceiveFromNeighbors()
{
  for (int i = 0; i < pbval_->nneighbor; ++i) {
    NeighborBlock& nb = pbval_->neighbor[i];
    ParticleBuffer& recv = recv_[nb.bufid];
    if (recv.npar > 0) FlushReceiveBuffer(recv);
  }
}

//--------------------------------------------------------------------------------------
//! \fn MeshBlock* Particles::FindTargetNeighbor(
//          int ox1, int ox2, int ox3, int xi1, int xi2, int xi3)
//  \brief finds the neighbor to send a particle to.

struct Neighbor* Particles::FindTargetNeighbor(
    int ox1, int ox2, int ox3, int xi1, int xi2, int xi3)
{
  // Find the head of the linked list.
  Neighbor *pn = &neighbor_[ox1+1][ox2+1][ox3+1];

  // Search down the list if the neighbor is at a finer level.
  if (pn->pnb->level > pmy_block->loc.level) {
    RegionSize& bs = pmy_block->block_size;
    int fi[2] = {0, 0}, i = 0;
    if (active1_ && ox1 == 0) fi[i++] = 2 * (xi1 - pmy_block->is) / bs.nx1;
    if (active2_ && ox2 == 0) fi[i++] = 2 * (xi2 - pmy_block->js) / bs.nx2;
    if (active3_ && ox3 == 0) fi[i++] = 2 * (xi3 - pmy_block->ks) / bs.nx3;
    while (pn != NULL) {
      NeighborBlock *pnb = pn->pnb;
      if (pnb->fi1 == fi[0] && pnb->fi2 == fi[1]) break;
      pn = pn->next;
    }
  }

  // Return the target neighbor.
  return pn;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::FlushReceiveBuffer(ParticleBuffer& recv)
//  \brief adds particles from the receive buffer.

void Particles::FlushReceiveBuffer(ParticleBuffer& recv)
{
  // Check the memory size.
  int nprecv = recv.npar;
  if (npar + nprecv > nparmax) {
    // Increase maximum number of particles allowed.
    nparmax += 2 * (npar + nprecv - nparmax);

    // Increase size of property arrays
    intprop.ResizeLastDimension(nparmax);
    realprop.ResizeLastDimension(nparmax);
    if (naux > 0) auxprop.ResizeLastDimension(nparmax);
    if (nwork > 0) work.ResizeLastDimension(nparmax);

    // Reassign the shorthands.
    AssignShorthands();
  }

  // Flush the receive buffers.
  long *pi = recv.ibuf;
  Real *pr = recv.rbuf;
  for (long k = npar; k < npar + nprecv; ++k) {
    for (int j = 0; j < nint; ++j)
      intprop(j,k) = *pi++;
    for (int j = 0; j < nreal; ++j)
      realprop(j,k) = *pr++;
    for (int j = 0; j < naux; ++j)
      auxprop(j,k) = *pr++;
  }

  // Find their position indices.
  AthenaArray<Real> xps, yps, zps, xi1s, xi2s, xi3s;
  xps.InitWithShallowSlice(xp, 1, npar, nprecv);
  yps.InitWithShallowSlice(yp, 1, npar, nprecv);
  zps.InitWithShallowSlice(zp, 1, npar, nprecv);
  xi1s.InitWithShallowSlice(xi1, 1, npar, nprecv);
  xi2s.InitWithShallowSlice(xi2, 1, npar, nprecv);
  xi3s.InitWithShallowSlice(xi3, 1, npar, nprecv);
  GetPositionIndices(pmy_block, nprecv, xps, yps, zps, xi1s, xi2s, xi3s);

  // Clear the receive buffers.
  npar += nprecv;
  recv.npar = 0;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::LinkNeighbors()
//  \brief fetches neighbor information for later communication.

void Particles::LinkNeighbors()
{
  neighbor_[1][1][1].pmb = pmy_block;

  for (int i = 0; i < pbval_->nneighbor; ++i) {
    NeighborBlock& nb = pbval_->neighbor[i];
    Neighbor *pn = &neighbor_[nb.ox1+1][nb.ox2+1][nb.ox3+1];
    while (pn->next != NULL)
      pn = pn->next;
    if (pn->pnb != NULL) {
      pn->next = new Neighbor;
      pn = pn->next;
    }
    pn->pnb = &nb;
    if (nb.rank == Globals::my_rank)
      pn->pmb = pmy_mesh->FindMeshBlock(nb.gid);
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::SetPositionIndices()
//  \brief updates position indices of particles.

void Particles::SetPositionIndices()
{
  GetPositionIndices(pmy_block, npar, xp, yp, zp, xi1, xi2, xi3);
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::Integrate(int step)
//  \brief updates all particle positions and velocities from t to t + dt.

void Particles::Integrate(int step)
{
  Real t, dt;

  switch (step) {

  case 1:
    t = pmy_mesh->time;
    dt = 0.5 * pmy_mesh->dt;
    SaveStatus();
    EulerStep(t, dt, pmy_block->phydro->w);
    ReactToMeshAux(t, dt, pmy_block->phydro->w);
    break;

  case 2:
    t = pmy_mesh->time + 0.5 * pmy_mesh->dt;
    dt = pmy_mesh->dt;
    EulerStep(t, dt, pmy_block->phydro->w1);
    ReactToMeshAux(t, dt, pmy_block->phydro->w1);
    break;
  }
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::AddIntProperty()
//  \brief adds one integer property to the particles and returns the index.

int Particles::AddIntProperty()
{
  return nint++;
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::AddRealProperty()
//  \brief adds one real property to the particles and returns the index.

int Particles::AddRealProperty()
{
  return nreal++;
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::AddAuxProperty()
//  \brief adds one auxiliary property to the particles and returns the index.

int Particles::AddAuxProperty()
{
  return naux++;
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::AddWorkingArray()
//  \brief adds one working array to the particles and returns the index.

int Particles::AddWorkingArray()
{
  return nwork++;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::AssignShorthands()
//  \brief assigns shorthands by shallow coping slices of the data.

void Particles::AssignShorthands()
{
  pid.InitWithShallowSlice(intprop, 2, ipid, 1);

  xp.InitWithShallowSlice(realprop, 2, ixp, 1);
  yp.InitWithShallowSlice(realprop, 2, iyp, 1);
  zp.InitWithShallowSlice(realprop, 2, izp, 1);
  vpx.InitWithShallowSlice(realprop, 2, ivpx, 1);
  vpy.InitWithShallowSlice(realprop, 2, ivpy, 1);
  vpz.InitWithShallowSlice(realprop, 2, ivpz, 1);

  xp0.InitWithShallowSlice(auxprop, 2, ixp0, 1);
  yp0.InitWithShallowSlice(auxprop, 2, iyp0, 1);
  zp0.InitWithShallowSlice(auxprop, 2, izp0, 1);
  vpx0.InitWithShallowSlice(auxprop, 2, ivpx0, 1);
  vpy0.InitWithShallowSlice(auxprop, 2, ivpy0, 1);
  vpz0.InitWithShallowSlice(auxprop, 2, ivpz0, 1);

  xi1.InitWithShallowSlice(work, 2, ixi1, 1);
  xi2.InitWithShallowSlice(work, 2, ixi2, 1);
  xi3.InitWithShallowSlice(work, 2, ixi3, 1);
  apx.InitWithShallowSlice(work, 2, iapx, 1);
  apy.InitWithShallowSlice(work, 2, iapy, 1);
  apz.InitWithShallowSlice(work, 2, iapz, 1);
}

//--------------------------------------------------------------------------------------
//! \fn Particles::GetSizeInBytes()
//  \brief returns the data size in bytes in the meshblock.

size_t Particles::GetSizeInBytes()
{
  size_t size = sizeof(npar);
  if (npar > 0) size += npar * (nint * sizeof(long) + nreal * sizeof(Real));
  return size;
}

//--------------------------------------------------------------------------------------
//! \fn Particles::ReadRestart()
//  \brief reads the particle data from the restart file.

#include <cstring>

void Particles::ReadRestart(char *mbdata, int &os)
{
  // Read number of particles.
  std::memcpy(&npar, &(mbdata[os]), sizeof(npar));
  os += sizeof(npar);
  if (npar > nparmax)
  {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [Particles::ReadRestart]" << std::endl
        << "npar = " << npar << " > nparmax = " << nparmax << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  if (npar > 0) {
    // Read integer properties.
    size_t size = npar * sizeof(long);
    for (int k = 0; k < nint; ++k) {
      std::memcpy(&(intprop(k,0)), &(mbdata[os]), size);
      os += size;
    }

    // Read real properties.
    size = npar * sizeof(Real);
    for (int k = 0; k < nreal; ++k) {
      std::memcpy(&(realprop(k,0)), &(mbdata[os]), size);
      os += size;
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn Particles::WriteRestart()
//  \brief writes the particle data to the restart file.

void Particles::WriteRestart(char *&pdata)
{
  // Write number of particles.
  memcpy(pdata, &npar, sizeof(npar));
  pdata += sizeof(npar);

  if (npar > 0) {
    // Write integer properties.
    size_t size = npar * sizeof(long);
    for (int k = 0; k < nint; ++k) {
      std::memcpy(pdata, &(intprop(k,0)), size);
      pdata += size;
    }

    // Write real properties.
    size = npar * sizeof(Real);
    for (int k = 0; k < nreal; ++k) {
      std::memcpy(pdata, &(realprop(k,0)), size);
      pdata += size;
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn Particles::FormattedTableOutput()
//  \brief outputs the particle data in tabulated format.

#include <fstream>
#include <iomanip>

void Particles::FormattedTableOutput(Mesh *pm, OutputParameters op)
{
  MeshBlock *pmb = pm->pblock;
  Particles *ppar;
  std::stringstream fname, msg;
  std::ofstream os;

  // Loop over MeshBlocks
  while (pmb != NULL) {
    ppar = pmb->ppar;

    // Create the filename.
    fname << op.file_basename
          << ".block" << pmb->gid << '.' << op.file_id
          << '.' << std::setw(5) << std::right << std::setfill('0') << op.file_number
          << '.' << "par.tab";

    // Open the file for write.
    os.open(fname.str().data());
    if (!os.is_open()) {
      msg << "### FATAL ERROR in function [Particles::FormattedTableOutput]"
          << std::endl << "Output file '" << fname.str() << "' could not be opened"
          << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }

    // Write the time.
    os << std::setprecision(18);
    os << "# Athena++ particle data at time = " << pm->time << std::endl;

    // Write the particle data in the meshblock.
    for (long k = 0; k < ppar->npar; ++k)
      os << ppar->pid(k) << "  "
         << ppar->xp(k) << "  " << ppar->yp(k) << "  " << ppar->zp(k) << "  "
         << ppar->vpx(k) << "  " << ppar->vpy(k) << "  " << ppar->vpz(k) << std::endl;

    // Close the file and get the next meshblock.
    os.close();
    fname.str("");
    pmb = pmb->next;
  }
}

//--------------------------------------------------------------------------------------
//! \fn void _CartesianToMeshCoords(x, y, z, x1, x2, x3)
//  \brief returns in (x1, x2, x3) the coordinates used by the mesh from Cartesian 
//         coordinates (x, y, z).
// TODO: Currently only supports Cartesian to Cartensian.
// TODO: Generalize and move this to the Coordinates class.

inline void _CartesianToMeshCoords(Real x, Real y, Real z, Real& x1, Real& x2, Real& x3)
{
  x1 = x;
  x2 = y;
  x3 = z;
}

//--------------------------------------------------------------------------------------
//! \fn void _MeshCoordsToCartesian(x1, x2, x3, x, y, z)
//  \brief returns in Cartesian coordinates (x, y, z) from (x1, x2, x3) the coordinates
//         used by the mesh.
// TODO: Currently only supports Cartesian to Cartensian.
// TODO: Generalize and move this to the Coordinates class.

inline void _MeshCoordsToCartesian(Real x1, Real x2, Real x3, Real& x, Real& y, Real& z)
{
  x = x1;
  y = x2;
  z = x3;
}

//--------------------------------------------------------------------------------------
//! \fn void _MeshCoordsToIndices(MeshBlock *pmb, Real x1, Real x2, Real x3,
//                                Real& xi1, Real& xi2, Real& xi3)
//  \brief returns in index coordinates (xi1, xi2, xi3) with respect to the local
//         grid of MeshBlock pmb from the physical coordinates (x1, x2, x3).
// TODO: Currently only supports uniform mesh.
// TODO: Generalize and move this to the Coordinates class.

void _MeshCoordsToIndices(MeshBlock *pmb, Real x1, Real x2, Real x3,
                          Real& xi1, Real& xi2, Real& xi3)
{
  // Get the meshblock info.
  const int IS = pmb->is;
  const int JS = pmb->js;
  const int KS = pmb->ks;
  const RegionSize& block_size = pmb->block_size;
  const Coordinates *pcoord = pmb->pcoord;

  // Make the conversion.
  xi1 = (block_size.nx1 > 1) ? IS + (x1 - block_size.x1min) / pcoord->dx1f(IS) : IS;
  xi2 = (block_size.nx2 > 1) ? JS + (x2 - block_size.x2min) / pcoord->dx2f(JS) : JS;
  xi3 = (block_size.nx3 > 1) ? KS + (x3 - block_size.x3min) / pcoord->dx3f(KS) : KS;
}

//--------------------------------------------------------------------------------------
//! \fn void _IndicesToMeshCoords(MeshBlock *pmb, Real xi1, Real xi2, Real xi3,
//                                Real& x1, Real& x2, Real& x3)
//  \brief returns in mesh coordinates (x1, x2, x3) from index coordinates
//         (xi1, xi2, xi3) with respect to the local grid of MeshBlock pmb.
// TODO: Currently only supports uniform mesh.
// TODO: Generalize and move this to the Coordinates class.

void _IndicesToMeshCoords(MeshBlock *pmb, Real xi1, Real xi2, Real xi3,
                          Real& x1, Real& x2, Real& x3)
{
  // Get the meshblock info.
  const int IS = pmb->is;
  const int JS = pmb->js;
  const int KS = pmb->ks;
  const RegionSize& block_size = pmb->block_size;
  const Coordinates *pcoord = pmb->pcoord;

  // Make the conversion.
  x1 = (block_size.nx1 > 1) ?
           block_size.x1min + (xi1 - IS) * pcoord->dx1f(IS) : pcoord->x1v(IS);
  x2 = (block_size.nx2 > 1) ?
           block_size.x2min + (xi2 - JS) * pcoord->dx2f(JS) : pcoord->x2v(JS);
  x3 = (block_size.nx3 > 1) ?
           block_size.x3min + (xi3 - KS) * pcoord->dx3f(KS) : pcoord->x3v(KS);
}

//--------------------------------------------------------------------------------------
//! \fn int CheckSide(int xi, nx, int xi1, int xi2)
//  \brief returns -1 if xi < xi1, +1 if xi > xi2, or 0 otherwise.

inline int CheckSide(int xi, int xi1, int xi2)
{
   if (xi < xi1) return -1;
   if (xi > xi2) return +1;
   return 0;
}
