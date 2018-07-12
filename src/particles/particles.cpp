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
#include "../mesh/meshblock_tree.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "particles.hpp"

// Class variable initialization
bool Particles::initialized = false;
int Particles::nint = 0;
int Particles::nreal = 0;
int Particles::naux = 0;
int Particles::nwork = 0;
int Particles::nmeshaux = 0;
int Particles::ipid = -1;
int Particles::ixp = -1, Particles::iyp = -1, Particles::izp = -1;
int Particles::ivpx = -1, Particles::ivpy = -1, Particles::ivpz = -1;
int Particles::ixp0 = -1, Particles::iyp0 = -1, Particles::izp0 = -1;
int Particles::ivpx0 = -1, Particles::ivpy0 = -1, Particles::ivpz0 = -1;
int Particles::ixi1 = -1, Particles::ixi2 = -1, Particles::ixi3 = -1;
int Particles::iapx = -1, Particles::iapy = -1, Particles::iapz = -1;

// Local function prototypes
void _CartesianToMeshCoords(Real x, Real y, Real z, Real& x1, Real& x2, Real& x3);
void _MeshCoordsToCartesian(Real x1, Real x2, Real x3, Real& x, Real& y, Real& z);
void _MeshCoordsToIndices(MeshBlock *pmb, Real x1, Real x2, Real x3,
                          Real& xi1, Real& xi2, Real& xi3);
void _IndicesToMeshCoords(MeshBlock *pmb, Real xi1, Real xi2, Real xi3,
                          Real& x1, Real& x2, Real& x3);
int _CheckSide(Real xi, int nx, int xi1, int xi2);

//--------------------------------------------------------------------------------------
//! \fn Particles::Initialize(ParameterInput *pin)
//  \brief initializes the class.

void Particles::Initialize(ParameterInput *pin)
{
  if (!initialized) {
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

    initialized = true;
  }
}

//--------------------------------------------------------------------------------------
//! \fn Particles::Particles(MeshBlock *pmb, ParameterInput *pin)
//  \brief constructs a Particles instance.

Particles::Particles(MeshBlock *pmb, ParameterInput *pin)
{
  // Point to the calling MeshBlock.
  pmy_block = pmb;
  pmy_mesh = pmb->pmy_mesh;
  nparmax = pin->GetOrAddInteger("particles", "nparmax", 1);
  npar = 0;

  // Allocate integer properties.
  intprop.NewAthenaArray(nint,nparmax);

  // Allocate integer properties.
  realprop.NewAthenaArray(nreal,nparmax);

  // Allocate auxiliary properties.
  if (naux > 0) auxprop.NewAthenaArray(naux,nparmax);

  // Allocate working arrays.
  if (nwork > 0) work.NewAthenaArray(nwork,nparmax);

  // Allocate mesh auxiliaries.
  ppm = new ParticleMesh(this, nmeshaux);

  // Shallow copy to shorthands.
  AssignShorthands();

  // Allocate buffers.
  nprecvmax = 1;
  nprecv = 0;
  irecv.NewAthenaArray(nint,nprecvmax);
  rrecv.NewAthenaArray(nreal+naux,nprecvmax);
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

  // Delete buffers.
  irecv.DeleteAthenaArray();
  rrecv.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
//! \fn bool Particles::ApplyBoundaryConditions(Mesh *pm, Real &x1, Real &x2, Real &x3)
//  \brief applies boundary conditions to (x1,x2,x3) and return true if (x1,x2,x3) is
//         outside the mesh.  Otherwise, (x1,x2,x3) is unchanged and false is returned.

bool Particles::ApplyBoundaryConditions(Mesh *pm, Real &x1, Real &x2, Real &x3)
{
  bool flag = false;

  if (pm->mesh_size.nx1 > 1) {
    if (x1 <= pm->mesh_size.x1min) {
      // Inner x1
      if (pm->mesh_bcs[INNER_X1] == PERIODIC_BNDRY)
        x1 += pm->mesh_size.x1max - pm->mesh_size.x1min;
      else {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::ApplyBoundaryConditions]"
            << std::endl
            << "Non-periodic boundary for inner x1 not supported. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      flag = true;
    } else if (x1 >= pm->mesh_size.x1max) {
      // Outer x1
      if (pm->mesh_bcs[OUTER_X1] == PERIODIC_BNDRY)
        x1 -= pm->mesh_size.x1max - pm->mesh_size.x1min;
      else {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::ApplyBoundaryConditions]"
            << std::endl
            << "Non-periodic boundary for outer x1 not supported. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      flag = true;
    }
  }

  if (pm->mesh_size.nx2 > 1) {
    if (x2 <= pm->mesh_size.x2min) {
      // Inner x2
      if (pm->mesh_bcs[INNER_X2] == PERIODIC_BNDRY)
        x2 += pm->mesh_size.x2max - pm->mesh_size.x2min;
      else {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::ApplyBoundaryConditions]"
            << std::endl
            << "Non-periodic boundary for inner x2 not supported. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      flag = true;
    } else if (x2 >= pm->mesh_size.x2max) {
      // Outer x2
      if (pm->mesh_bcs[OUTER_X2] == PERIODIC_BNDRY)
        x2 -= pm->mesh_size.x2max - pm->mesh_size.x2min;
      else {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::ApplyBoundaryConditions]"
            << std::endl
            << "Non-periodic boundary for outer x2 not supported. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      flag = true;
    }
  }

  if (pm->mesh_size.nx3 > 1) {
    if (x3 <= pm->mesh_size.x3min) {
      // Inner x3
      if (pm->mesh_bcs[INNER_X3] == PERIODIC_BNDRY)
        x3 += pm->mesh_size.x3max - pm->mesh_size.x3min;
      else {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::ApplyBoundaryConditions]"
            << std::endl
            << "Non-periodic boundary for inner x3 not supported. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      flag = true;
    } else if (x3 >= pm->mesh_size.x3max) {
      // Outer x3
      if (pm->mesh_bcs[OUTER_X3] == PERIODIC_BNDRY)
        x3 -= pm->mesh_size.x3max - pm->mesh_size.x3min;
      else {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::ApplyBoundaryConditions]"
            << std::endl
            << "Non-periodic boundary for outer x3 not supported. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
      flag = true;
    }
  }

  return flag;
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
    GetPositionIndices(pmy_block, npar, xp, yp, zp, xi1, xi2, xi3);
    SendToNeighbors();
  }

  // Send MeshAux boundary.
  if (step > 0 && nmeshaux > 0)
    ppm->SendBoundary();
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::ReceiveParticlesAndMesh(int step)
//  \brief receives particles and meshaux near boundaries from neighbors.

void Particles::ReceiveParticlesAndMesh(int step)
{
  // Flush Particles receive buffer.
  if (nprecv > 0) FlushReceiveBuffer();

  // Flush ParticleMesh receive buffers and deposit MeshAux to MeshBlock.
  if (nmeshaux > 0) {
    ppm->ReceiveBoundary();
    switch (step) {
    case 1:
      DepositToMesh(pmy_block->phydro->u1);
      break;
    case 2:
      DepositToMesh(pmy_block->phydro->u);
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
  const int NX1 = pmy_block->block_size.nx1;
  const int NX2 = pmy_block->block_size.nx2;
  const int NX3 = pmy_block->block_size.nx3;
  const int IS = pmy_block->is;
  const int IE = pmy_block->ie;
  const int JS = pmy_block->js;
  const int JE = pmy_block->je;
  const int KS = pmy_block->ks;
  const int KE = pmy_block->ke;

  Mesh *pm = pmy_block->pmy_mesh;

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
    int ox1 = _CheckSide(xi1(k), NX1, IS, IE),
        ox2 = _CheckSide(xi2(k), NX2, JS, JE),
        ox3 = _CheckSide(xi3(k), NX3, KS, KE);
    if (ox1 == 0 && ox2 == 0 && ox3 == 0) {
      ++k;
      continue;
    }

    // Convert to the MeshBlock coordinates.
    Real x1, x2, x3;
    _IndicesToMeshCoords(pmy_block, xi1(k), xi2(k), xi3(k), x1, x2, x3);

    // Find the neighbor MeshBlock to send it to.
    MeshBlockTree *pnmbt =
        pm->tree.FindNeighbor(pmy_block->loc, ox1, ox2, ox3, pm->mesh_bcs,
                              pm->nrbx1, pm->nrbx2, pm->nrbx3, pm->root_level);
    if (pnmbt == NULL) {
      std::stringstream msg;
      msg << "### FATAL ERROR in function [Particles::SendToNeighbors]" << std::endl
          << "cannot find the neighboring MeshBlock. " << std::endl;
      throw std::runtime_error(msg.str().c_str());
      continue;
    }

    MeshBlock *pnmb;
    if (pnmbt->flag)  // Neighbor is on the same or a courser level:
      pnmb = pm->FindMeshBlock(pnmbt->gid);
    else {          // Neighbor is on a finer level:
      bool flag = true;
      for (int i = 0; flag && i < 2; ++i)
        for (int j = 0; flag && j < 2; ++j)
          for (int k = 0; flag && k < 2; ++k) {
            int gid = pnmbt->pleaf[i][j][k]->gid;
            if (gid < 0) continue;
            pnmb = pm->FindMeshBlock(gid);
            flag = false;
            if (pm->mesh_size.nx1 > 1)
              flag = flag || x1 < pnmb->block_size.x1min || x1 > pnmb->block_size.x1max;
            if (pm->mesh_size.nx2 > 1)
              flag = flag || x2 < pnmb->block_size.x2min || x2 > pnmb->block_size.x2max;
            if (pm->mesh_size.nx3 > 1)
              flag = flag || x3 < pnmb->block_size.x3min || x3 > pnmb->block_size.x3max;
          }
      if (flag) {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::SendToNeighbors]" << std::endl
            << "cannot find the neighboring MeshBlock. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
        continue;
      }
    }
    Particles *pnp = pnmb->ppar;

    // Apply boundary conditions and convert back to Cartesian coordinates.
    if (ApplyBoundaryConditions(pm, x1, x2, x3))
      _MeshCoordsToCartesian(x1, x2, x3, xp(k), yp(k), zp(k));

    // Check the buffer size of the target MeshBlock.
    if (pnp->nprecv >= pnp->nprecvmax) {
      pnp->nprecvmax *= 2;
      pnp->irecv.ResizeLastDimension(pnp->nprecvmax);
      pnp->rrecv.ResizeLastDimension(pnp->nprecvmax);
    }

    // Copy the properties of the particle to the neighbor.
    for (int j = 0; j < nint; ++j)
      pnp->irecv(j,pnp->nprecv) = intprop(j,k);
    for (int j = 0; j < nreal; ++j)
      pnp->rrecv(j,pnp->nprecv) = realprop(j,k);
    for (int i = nreal, j = 0; j < naux; ++i, ++j)
      pnp->rrecv(i,pnp->nprecv) = auxprop(j,k);
    ++pnp->nprecv;

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
//! \fn void Particles::FlushReceiveBuffer()
//  \brief adds particles from the receive buffer.

void Particles::FlushReceiveBuffer()
{
  // Check the memory size.
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
  for (int j = 0; j < nint; ++j)
    for (long ip = npar, k = 0; k < nprecv; ++ip, ++k)
      intprop(j,ip) = irecv(j,k);
  for (int j = 0; j < nreal; ++j)
    for (long ip = npar, k = 0; k < nprecv; ++ip, ++k)
      realprop(j,ip) = rrecv(j,k);
  for (int i = nreal, j = 0; j < naux; ++i, ++j)
    for (long ip = npar, k = 0; k < nprecv; ++ip, ++k)
      auxprop(j,ip) = rrecv(i,k);

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
  nprecv = 0;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::Integrate(int step)
//  \brief updates all particle positions and velocities from t to t + dt.
//======================================================================================

void Particles::Integrate(int step)
{
  Real t, dt;

  switch (step) {

  case 1:
    t = pmy_mesh->time;
    dt = 0.5 * pmy_mesh->dt;
    SaveStatus();
    EulerStep(t, dt, pmy_block->phydro->w);
    break;

  case 2:
    t = pmy_mesh->time + 0.5 * pmy_mesh->dt;
    dt = pmy_mesh->dt;
    EulerStep(t, dt, pmy_block->phydro->w1);
    break;
  }

  ReactToMeshAux(t, dt);
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
//! \fn int Particles::AddMeshAux()
//  \brief adds one auxiliary to the mesh and returns the index.

int Particles::AddMeshAux()
{
  return nmeshaux++;
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
  xi1 = (pmb->block_size.nx1 > 1) ?
            IS + (x1 - block_size.x1min) / pcoord->dx1f(IS) : IS;
  xi2 = (pmb->block_size.nx2 > 1) ?
            JS + (x2 - block_size.x2min) / pcoord->dx2f(JS) : JS;
  xi3 = (pmb->block_size.nx3 > 1) ?
            KS + (x3 - block_size.x3min) / pcoord->dx3f(KS) : KS;
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
  x1 = (pmb->block_size.nx1 > 1) ?
           block_size.x1min + (xi1 - IS) * pcoord->dx1f(IS) : pcoord->x1v(IS);
  x2 = (pmb->block_size.nx2 > 1) ?
           block_size.x2min + (xi2 - JS) * pcoord->dx2f(JS) : pcoord->x2v(JS);
  x3 = (pmb->block_size.nx3 > 1) ?
           block_size.x3min + (xi3 - KS) * pcoord->dx3f(KS) : pcoord->x3v(KS);
}

//--------------------------------------------------------------------------------------
//! \fn int _CheckSide(Real xi, int nx, int xi1, int xi2)
//  \brief returns -1 if int(xi) < xi1 and nx > 1, +1 if int(xi) > xi2 and nx > 1,
//         and 0 otherwise.

inline int _CheckSide(Real xi, int nx, int xi1, int xi2)
{
   if (nx > 1) {
     if (int(xi) < xi1) return -1;
     if (int(xi) > xi2) return +1;
   }
   return 0;
}
