//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particles.cpp
//  \brief implements functions in particle classes

// C++ Standard Libraries
#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
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
#ifdef MPI_PARALLEL
MPI_Comm Particles::my_comm = MPI_COMM_NULL;
#endif

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

void Particles::Initialize(ParameterInput *pin) {
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

#ifdef MPI_PARALLEL
  // Get my MPI communicator.
  MPI_Comm_dup(MPI_COMM_WORLD, &my_comm);
#endif

  initialized = true;
}

//--------------------------------------------------------------------------------------
//! \fn Particles::Particles(MeshBlock *pmb, ParameterInput *pin)
//  \brief constructs a Particles instance.

Particles::Particles(MeshBlock *pmb, ParameterInput *pin) {
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
  ClearBoundary();
}

//--------------------------------------------------------------------------------------
//! \fn Particles::~Particles()
//  \brief destroys a Particles instance.

Particles::~Particles() {
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
//! \fn void Particles::ClearBoundary()
//  \brief resets boundary for particle transportation.

void Particles::ClearBoundary() {
  for (int i = 0; i < pbval_->nneighbor; ++i) {
    NeighborBlock& nb = pbval_->neighbor[i];
    bstatus_[nb.bufid] = BNDRY_WAITING;
#ifdef MPI_PARALLEL
    if (nb.rank != Globals::my_rank) {
      ParticleBuffer& recv = recv_[nb.bufid];
      recv.flagn = recv.flagi = recv.flagr = 0;
      send_[nb.bufid].npar = 0;
    }
#endif
  }

  ppm->ClearBoundary();
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::Integrate(int step)
//  \brief updates all particle positions and velocities from t to t + dt.

void Particles::Integrate(int stage) {
  Real t, dt;

  switch (stage) {

  case 1:
    t = pmy_mesh->time;
    dt = 0.5 * pmy_mesh->dt;
    SaveStatus();
    break;

  case 2:
    t = pmy_mesh->time + 0.5 * pmy_mesh->dt;
    dt = pmy_mesh->dt;
    break;
  }

  EulerStep(t, dt, pmy_block->phydro->w);
  ReactToMeshAux(t, dt, pmy_block->phydro->w);
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::LinkNeighbors()
//  \brief fetches neighbor information for later communication.

void Particles::LinkNeighbors() {
  // Construct links to neighbors.
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
    if (nb.rank == Globals::my_rank) {
      pn->pmb = pmy_mesh->FindMeshBlock(nb.gid);
    } else {
#ifdef MPI_PARALLEL
      send_[nb.bufid].tag = (nb.gid<<8) | (nb.targetid<<2),
      recv_[nb.bufid].tag = (pmy_block->gid<<8) | (nb.bufid<<2);
#endif
    }
  }

  // Initiate ParticleMesh boundary data.
  ppm->SetBoundaryAttributes();
  ppm->InitiateBoundaryData();
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::SendParticleMesh()
//  \brief send ParticleMesh meshaux near boundaries to neighbors.

void Particles::SendParticleMesh() {
  if (ppm->nmeshaux > 0)
    ppm->SendBoundary();
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::SendToNeighbors()
//  \brief sends particles outside boundary to the buffers of neighboring meshblocks.

void Particles::SendToNeighbors() {
  const int IS = pmy_block->is;
  const int IE = pmy_block->ie;
  const int JS = pmy_block->js;
  const int JE = pmy_block->je;
  const int KS = pmy_block->ks;
  const int KE = pmy_block->ke;

  // TODO(ccyang): Currently only works for Cartesian.
  if (COORDINATE_SYSTEM != "cartesian") {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [Particles::SendToNeighbors]" << std::endl
        << "Non-Cartesian coordinates not yet implemented. " << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  for (int k = 0; k < npar; ) {
    // Check if a particle is outside the boundary.
    int xi1i = static_cast<int>(xi1(k)),
        xi2i = static_cast<int>(xi2(k)),
        xi3i = static_cast<int>(xi3(k));
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

    // Find the neighbor block to send it to.
    Neighbor *pn = FindTargetNeighbor(ox1, ox2, ox3, xi1i, xi2i, xi3i);
    if (pn == NULL) {
      std::stringstream msg;
      msg << "### FATAL ERROR in function [Particles::SendToNeighbors]" << std::endl
          << "cannot find the neighbor block to send the particle to. " << std::endl;
      throw std::runtime_error(msg.str().data());
      continue;
    }
    NeighborBlock *pnb = pn->pnb;

    // Determine which particle buffer to use.
    ParticleBuffer *ppb = NULL;
    if (pnb->rank == Globals::my_rank) {
      // No need to send if back to the same block.
      if (pnb->gid == pmy_block->gid) {
        _MeshCoordsToIndices(pmy_block, x1, x2, x3, xi1(k), xi2(k), xi3(k));
        ++k;
        continue;
      }
      // Use the target receive buffer.
      ppb = &pn->pmb->ppar->recv_[pnb->targetid];

    } else {
#ifdef MPI_PARALLEL
      // Use the send buffer.
      ppb = &send_[pnb->bufid];
#endif
    }

    // Check the buffer size.
    if (ppb->npar >= ppb->nparmax)
      ppb->Reallocate((ppb->nparmax > 0) ? 2 * ppb->nparmax : 1);

    // Copy the properties of the particle to the buffer.
    int *pi = ppb->ibuf + ParticleBuffer::nint * ppb->npar;
    for (int j = 0; j < nint; ++j)
      *pi++ = intprop(j,k);
    Real *pr = ppb->rbuf + ParticleBuffer::nreal * ppb->npar;
    for (int j = 0; j < nreal; ++j)
      *pr++ = realprop(j,k);
    for (int j = 0; j < naux; ++j)
      *pr++ = auxprop(j,k);
    ++ppb->npar;

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

  // Send to neighbor processes and update boundary status.
  for (int i = 0; i < pbval_->nneighbor; ++i) {
    NeighborBlock& nb = pbval_->neighbor[i];
    int dst = nb.rank;
    if (dst == Globals::my_rank) {
      Particles *ppar = pmy_mesh->FindMeshBlock(nb.gid)->ppar;
      ppar->bstatus_[nb.targetid] =
          (ppar->recv_[nb.targetid].npar > 0) ? BNDRY_ARRIVED : BNDRY_COMPLETED;
    } else {
#ifdef MPI_PARALLEL
      ParticleBuffer& send = send_[nb.bufid];
      int npsend = send.npar;
      MPI_Send(&npsend, 1, MPI_INT, nb.rank, send.tag, my_comm);
      if (npsend > 0) {
        MPI_Request req = MPI_REQUEST_NULL;
        MPI_Isend(send.ibuf, npsend * ParticleBuffer::nint, MPI_LONG,
                  dst, send.tag + 1, my_comm, &req);
        MPI_Request_free(&req);
        MPI_Isend(send.rbuf, npsend * ParticleBuffer::nreal, MPI_ATHENA_REAL,
                  dst, send.tag + 2, my_comm, &req);
        MPI_Request_free(&req);
      }
#endif
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::SetPositionIndices()
//  \brief updates position indices of particles.

void Particles::SetPositionIndices() {
  GetPositionIndices(pmy_block, npar, xp, yp, zp, xi1, xi2, xi3);
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::StartReceiving()
//  \brief starts receiving ParticleMesh meshaux near boundary from neighbor processes.

void Particles::StartReceiving() {
  ppm->StartReceiving();
}

//--------------------------------------------------------------------------------------
//! \fn bool Particles::ReceiveFromNeighbors()
//  \brief receives particles from neighboring meshblocks and returns a flag indicating
//         if all receives are completed.

bool Particles::ReceiveFromNeighbors() {
  bool flag = true;

  for (int i = 0; i < pbval_->nneighbor; ++i) {
    NeighborBlock& nb = pbval_->neighbor[i];
    enum BoundaryStatus& bstatus = bstatus_[nb.bufid];

#ifdef MPI_PARALLEL
    // Communicate with neighbor processes.
    if (nb.rank != Globals::my_rank && bstatus == BNDRY_WAITING) {
      ParticleBuffer& recv = recv_[nb.bufid];
      if (!recv.flagn) {
        // Get the number of incoming particles.
        if (recv.reqi == MPI_REQUEST_NULL)
          MPI_Irecv(&recv.npar, 1, MPI_INT, nb.rank, recv.tag, my_comm, &recv.reqi);
        else
          MPI_Test(&recv.reqi, &recv.flagn, MPI_STATUS_IGNORE);
        if (recv.flagn) {
          if (recv.npar > 0) {
            // Check the buffer size.
            int nprecv = recv.npar;
            if (nprecv > recv.nparmax) {
              recv.npar = 0;
              recv.Reallocate(2 * nprecv - recv.nparmax);
              recv.npar = nprecv;
            }
          } else {
            // No incoming particles.
            bstatus = BNDRY_COMPLETED;
          }
        }
      }
      if (recv.flagn && recv.npar > 0) {
        // Receive data from the neighbor.
        if (!recv.flagi) {
          if (recv.reqi == MPI_REQUEST_NULL)
            MPI_Irecv(recv.ibuf, recv.npar * ParticleBuffer::nint, MPI_LONG,
                      nb.rank, recv.tag + 1, my_comm, &recv.reqi);
          else
            MPI_Test(&recv.reqi, &recv.flagi, MPI_STATUS_IGNORE);
        }
        if (!recv.flagr) {
          if (recv.reqr == MPI_REQUEST_NULL)
            MPI_Irecv(recv.rbuf, recv.npar * ParticleBuffer::nreal, MPI_ATHENA_REAL,
                      nb.rank, recv.tag + 2, my_comm, &recv.reqr);
          else
            MPI_Test(&recv.reqr, &recv.flagr, MPI_STATUS_IGNORE);
        }
        if (recv.flagi && recv.flagr)
          bstatus = BNDRY_ARRIVED;
      }
    }
#endif

    switch (bstatus) {

      case BNDRY_COMPLETED:
        break;

      case BNDRY_WAITING:
        flag = false;
        break;

      case BNDRY_ARRIVED:
        ParticleBuffer& recv = recv_[nb.bufid];
        FlushReceiveBuffer(recv);
        bstatus = BNDRY_COMPLETED;
        break;
    }
  }

  return flag;
}

//--------------------------------------------------------------------------------------
//! \fn bool Particles::ReceiveParticleMesh(int step)
//  \brief receives ParticleMesh meshaux near boundaries from neighbors and returns a
//         flag indicating if all receives are completed.

bool Particles::ReceiveParticleMesh(int stage) {
  if (ppm->nmeshaux <= 0) return true;

  // Flush ParticleMesh receive buffers.
  bool flag = ppm->ReceiveBoundary();

  if (flag) {
    // Deposit ParticleMesh meshaux to MeshBlock.
    Hydro *phydro = pmy_block->phydro;
    Real t, dt;

    switch (stage) {

    case 1:
      t = pmy_mesh->time;
      dt = 0.5 * pmy_mesh->dt;
      break;

    case 2:
      t = pmy_mesh->time + 0.5 * pmy_mesh->dt;
      dt = pmy_mesh->dt;
      break;
    }

    DepositToMesh(t, dt, phydro->w, phydro->u);
  }

  return flag;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::ApplyBoundaryConditions(int k, Real &x1, Real &x2, Real &x3)
//  \brief applies boundary conditions to particle k and returns its updated mesh
//         coordinates (x1,x2,x3).

void Particles::ApplyBoundaryConditions(int k, Real &x1, Real &x2, Real &x3) {
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
//! \fn void Particles::GetPositionIndices(MeshBlock *pmb, int npar,
//                                         const AthenaArray<Real>& xp,
//                                         const AthenaArray<Real>& yp,
//                                         const AthenaArray<Real>& zp,
//                                         AthenaArray<Real>& xi1,
//                                         AthenaArray<Real>& xi2,
//                                         AthenaArray<Real>& xi3)
//  \brief finds the position indices of each particle with respect to the local grid.

void Particles::GetPositionIndices(MeshBlock *pmb, int npar,
                                   const AthenaArray<Real>& xp,
                                   const AthenaArray<Real>& yp,
                                   const AthenaArray<Real>& zp,
                                   AthenaArray<Real>& xi1,
                                   AthenaArray<Real>& xi2,
                                   AthenaArray<Real>& xi3) {
  for (int k = 0; k < npar; ++k) {
    // Convert to the Mesh coordinates.
    Real x1, x2, x3;
    _CartesianToMeshCoords(xp(k), yp(k), zp(k), x1, x2, x3);

    // Convert to the index space.
    _MeshCoordsToIndices(pmb, x1, x2, x3, xi1(k), xi2(k), xi3(k));
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::EulerStep(Real t, Real dt, const AthenaArray<Real>& meshsrc)
//  \brief evolves the particle positions and velocities by one Euler step.

void Particles::EulerStep(Real t, Real dt, const AthenaArray<Real>& meshsrc) {
  // Get the accelerations.
  ZeroAcceleration();
  AddAcceleration(t, dt, meshsrc);

  // Update the positions and velocities **from the beginning of the time step**.
  for (int k = 0; k < npar; ++k) {
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

void Particles::SaveStatus() {
  for (int k = 0; k < npar; ++k) {
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
//! \fn void Particles::ZeroAcceleration()
//  \brief initializes acceleration with zeros.

void Particles::ZeroAcceleration() {
  for (int k = 0; k < npar; ++k) {
    apx(k) = 0.0;
    apy(k) = 0.0;
    apz(k) = 0.0;
  }
}

//--------------------------------------------------------------------------------------
//! \fn MeshBlock* Particles::FindTargetNeighbor(
//          int ox1, int ox2, int ox3, int xi1, int xi2, int xi3)
//  \brief finds the neighbor to send a particle to.

struct Neighbor* Particles::FindTargetNeighbor(
    int ox1, int ox2, int ox3, int xi1, int xi2, int xi3) {
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

void Particles::FlushReceiveBuffer(ParticleBuffer& recv) {
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
  int *pi = recv.ibuf;
  Real *pr = recv.rbuf;
  for (int k = npar; k < npar + nprecv; ++k) {
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
//! \fn int Particles::AddIntProperty()
//  \brief adds one integer property to the particles and returns the index.

int Particles::AddIntProperty() {
  return nint++;
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::AddRealProperty()
//  \brief adds one real property to the particles and returns the index.

int Particles::AddRealProperty() {
  return nreal++;
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::AddAuxProperty()
//  \brief adds one auxiliary property to the particles and returns the index.

int Particles::AddAuxProperty() {
  return naux++;
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::AddWorkingArray()
//  \brief adds one working array to the particles and returns the index.

int Particles::AddWorkingArray() {
  return nwork++;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::AssignShorthands()
//  \brief assigns shorthands by shallow coping slices of the data.

void Particles::AssignShorthands() {
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
//! \fn void Particles::NewBlockTimeStep();
//  \brief constrains the time step by particles in the block.

void Particles::NewBlockTimeStep() {
  Coordinates *pc = pmy_block->pcoord;

  // Find the allowed time step for each particle.
  Real dt_inv2_max = 0.0, dt_inv2;
  for (int k = 0; k < npar; ++k) {
    dt_inv2 = 0.0;
    if (active1_) dt_inv2 += std::pow(vpx(k) / pc->dx1f(static_cast<int>(xi1(k))), 2);
    if (active2_) dt_inv2 += std::pow(vpy(k) / pc->dx2f(static_cast<int>(xi2(k))), 2);
    if (active3_) dt_inv2 += std::pow(vpz(k) / pc->dx3f(static_cast<int>(xi3(k))), 2);
    dt_inv2_max = std::max(dt_inv2_max, dt_inv2);
  }

  // Constrain the time step by the minimum.
  if (dt_inv2_max > 0.0) {
    Real dt_min = 1.0 / std::sqrt(dt_inv2_max);
    if (dt_min < pmy_block->new_block_dt) pmy_block->new_block_dt = dt_min;
  }
}

//--------------------------------------------------------------------------------------
//! \fn Particles::GetSizeInBytes()
//  \brief returns the data size in bytes in the meshblock.

size_t Particles::GetSizeInBytes() {
  size_t size = sizeof(npar);
  if (npar > 0) size += npar * (nint * sizeof(int) + nreal * sizeof(Real));
  return size;
}

//--------------------------------------------------------------------------------------
//! \fn Particles::ReadRestart()
//  \brief reads the particle data from the restart file.

#include <cstring>

void Particles::ReadRestart(char *mbdata, int &os) {
  // Read number of particles.
  std::memcpy(&npar, &(mbdata[os]), sizeof(npar));
  os += sizeof(npar);
  if (npar > nparmax) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [Particles::ReadRestart]" << std::endl
        << "npar = " << npar << " > nparmax = " << nparmax << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  if (npar > 0) {
    // Read integer properties.
    size_t size = npar * sizeof(int);
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

void Particles::WriteRestart(char *&pdata) {
  // Write number of particles.
  memcpy(pdata, &npar, sizeof(npar));
  pdata += sizeof(npar);

  if (npar > 0) {
    // Write integer properties.
    size_t size = npar * sizeof(int);
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

void Particles::FormattedTableOutput(Mesh *pm, OutputParameters op) {
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
    for (int k = 0; k < ppar->npar; ++k)
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
// TODO(ccyang): Currently only supports Cartesian to Cartensian.
// TODO(ccyang): Generalize and move this to the Coordinates class.

inline void _CartesianToMeshCoords(Real x, Real y, Real z, Real& x1, Real& x2, Real& x3) {
  x1 = x;
  x2 = y;
  x3 = z;
}

//--------------------------------------------------------------------------------------
//! \fn void _MeshCoordsToCartesian(x1, x2, x3, x, y, z)
//  \brief returns in Cartesian coordinates (x, y, z) from (x1, x2, x3) the coordinates
//         used by the mesh.
// TODO(ccyang): Currently only supports Cartesian to Cartensian.
// TODO(ccyang): Generalize and move this to the Coordinates class.

inline void _MeshCoordsToCartesian(Real x1, Real x2, Real x3, Real& x, Real& y, Real& z) {
  x = x1;
  y = x2;
  z = x3;
}

//--------------------------------------------------------------------------------------
//! \fn void _MeshCoordsToIndices(MeshBlock *pmb, Real x1, Real x2, Real x3,
//                                Real& xi1, Real& xi2, Real& xi3)
//  \brief returns in index coordinates (xi1, xi2, xi3) with respect to the local
//         grid of MeshBlock pmb from the physical coordinates (x1, x2, x3).
// TODO(ccyang): Currently only supports uniform mesh.
// TODO(ccyang): Generalize and move this to the Coordinates class.

void _MeshCoordsToIndices(MeshBlock *pmb, Real x1, Real x2, Real x3,
                          Real& xi1, Real& xi2, Real& xi3) {
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
// TODO(ccyang): Currently only supports uniform mesh.
// TODO(ccyang): Generalize and move this to the Coordinates class.

void _IndicesToMeshCoords(MeshBlock *pmb, Real xi1, Real xi2, Real xi3,
                          Real& x1, Real& x2, Real& x3) {
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

inline int CheckSide(int xi, int xi1, int xi2) {
   if (xi < xi1) return -1;
   if (xi > xi2) return +1;
   return 0;
}
