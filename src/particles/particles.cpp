//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particles.cpp
//  \brief implements functions in particle classes

#include <sstream>
#include <string>
#include "../athena_arrays.hpp"
#include "../mesh/meshblock_tree.hpp"
#include "../coordinates/coordinates.hpp"
#include "particles.hpp"

bool Particles::initialized = false;
int Particles::nint = 0;
int Particles::nreal = 0;
int Particles::ipid = -1;
int Particles::ixp = -1, Particles::iyp = -1, Particles::izp = -1;
int Particles::ivpx = -1, Particles::ivpy = -1, Particles::ivpz = -1;

void _ErrorIfInitialized(const std::string& calling_function, bool initialized);
void _CartesianToMeshCoords(Real x, Real y, Real z, Real& x1, Real& x2, Real& x3);
void _MeshCoordsToCartesian(Real x1, Real x2, Real x3, Real& x, Real& y, Real& z);
void _MeshCoordsToIndices(MeshBlock *pmb, Real x1, Real x2, Real x3,
                          Real& xi1, Real& xi2, Real& xi3);

//--------------------------------------------------------------------------------------
//! \fn Particles::Initialize()
//  \brief initializes the class.

void Particles::Initialize()
{
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

  initialized = true;
}

//--------------------------------------------------------------------------------------
//! \fn Particles::Particles(MeshBlock *pmb, ParameterInput *pin)
//  \brief constructs a Particles instance.

Particles::Particles(MeshBlock *pmb, ParameterInput *pin)
{
  // Initialize the class if first called.
  if (!initialized) Particles::Initialize();

  // Point to the calling MeshBlock.
  pmy_block = pmb;
  nparmax = pin->GetOrAddInteger("particles", "nparmax", 1);
  npar = 0;

  // Allocate integer properties.
  intprop.NewAthenaArray(nint,nparmax);

  // Allocate integer properties.
  realprop.NewAthenaArray(nreal,nparmax);

  // Shallow copy to shorthands.
  AssignShorthands();

  // Allocate position indices.
  xi1.NewAthenaArray(nparmax);
  xi2.NewAthenaArray(nparmax);
  xi3.NewAthenaArray(nparmax);

  // Allocate buffers.
  nrecvmax = 1;
  nrecv = 0;
  irecv.NewAthenaArray(nint,nrecvmax);
  rrecv.NewAthenaArray(nreal,nrecvmax);
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

  // Delete position indices.
  xi1.DeleteAthenaArray();
  xi2.DeleteAthenaArray();
  xi3.DeleteAthenaArray();

  // Delete buffers.
  irecv.DeleteAthenaArray();
  rrecv.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::ApplyBoundaryConditions(Mesh *pm, Real &x1, Real &x2, Real &x3)
//  \brief applies boundary conditions if (x1,x2,x3) is outside the mesh.

void Particles::ApplyBoundaryConditions(Mesh *pm, Real &x1, Real &x2, Real &x3)
{
  if (pm->mesh_size.nx1 > 1) {
    // Inner x1
    if (x1 <= pm->mesh_size.x1min) {
      if (pm->mesh_bcs[INNER_X1] == PERIODIC_BNDRY)
        x1 += pm->mesh_size.x1max - pm->mesh_size.x1min;
      else {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::ApplyBoundaryConditions]"
            << std::endl
            << "Non-periodic boundary for inner x1 not supported. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
    }

    // Outer x1
    if (x1 >= pm->mesh_size.x1max) {
      if (pm->mesh_bcs[OUTER_X1] == PERIODIC_BNDRY)
        x1 -= pm->mesh_size.x1max - pm->mesh_size.x1min;
      else {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::ApplyBoundaryConditions]"
            << std::endl
            << "Non-periodic boundary for outer x1 not supported. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
    }
  }

  if (pm->mesh_size.nx2 > 1) {
    // Inner x2
    if (x2 <= pm->mesh_size.x2min) {
      if (pm->mesh_bcs[INNER_X2] == PERIODIC_BNDRY)
        x2 += pm->mesh_size.x2max - pm->mesh_size.x2min;
      else {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::ApplyBoundaryConditions]"
            << std::endl
            << "Non-periodic boundary for inner x2 not supported. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
    }

    // Outer x2
    if (x2 >= pm->mesh_size.x2max) {
      if (pm->mesh_bcs[OUTER_X2] == PERIODIC_BNDRY)
        x2 -= pm->mesh_size.x2max - pm->mesh_size.x2min;
      else {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::ApplyBoundaryConditions]"
            << std::endl
            << "Non-periodic boundary for outer x2 not supported. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
    }
  }

  if (pm->mesh_size.nx3 > 1) {
    // Inner x3
    if (x3 <= pm->mesh_size.x3min) {
      if (pm->mesh_bcs[INNER_X3] == PERIODIC_BNDRY)
        x3 += pm->mesh_size.x3max - pm->mesh_size.x3min;
      else {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::ApplyBoundaryConditions]"
            << std::endl
            << "Non-periodic boundary for inner x3 not supported. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
    }

    // Outer x3
    if (x3 >= pm->mesh_size.x3max) {
      if (pm->mesh_bcs[OUTER_X3] == PERIODIC_BNDRY)
        x3 -= pm->mesh_size.x3max - pm->mesh_size.x3min;
      else {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::ApplyBoundaryConditions]"
            << std::endl
            << "Non-periodic boundary for outer x3 not supported. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
    }
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
//! \fn void Particles::Drift(Real t, Real dt)
//  \brief updates the particle positions from t to t + dt given velocities at t.

void Particles::Drift(Real t, Real dt)
{
  for (long k = 0; k < npar; ++k) {
    xp(k) += dt * vpx(k);
    yp(k) += dt * vpy(k);
    zp(k) += dt * vpz(k);
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::Kick(Real t, Real dt)
//  \brief updates the particle velocities from t to t + dt given accelerations at t.

#include <cmath>

void Particles::Kick(Real t, Real dt)
{
  Real a1 = 0.0, a2 = 0.0, a3 = 0.0;  // TODO: might need to be put inside the loop for
                                      // vectorization.

  for (long k = 0; k < npar; ++k) {
    vpx(k) += dt * a1;
    vpy(k) += dt * a2;
    vpz(k) += dt * a3;
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::Migrate(Mesh *pm)
//  \brief migrates particles that are outside of the MeshBlock boundary.
int _CheckSide(int nx, Real x, Real xmin, Real xmax);

void Particles::Migrate(Mesh *pm)
{
  // Send particles.
  MeshBlock *pmb = pm->pblock;
  while (pmb != NULL) {
    pmb->ppar->SendToNeighbors();
    pmb = pmb->next;
  }

  // Flush the receive buffers.
  pmb = pm->pblock;
  Particles *ppar;
  while (pmb != NULL) {
    ppar = pmb->ppar;
    if (ppar->nrecv > 0) ppar->FlushReceiveBuffer();
    pmb = pmb->next;
  }

  // Update the position indices.
  pmb = pm->pblock;
  while (pmb != NULL) {
    ppar = pmb->ppar;
    long npar = ppar->npar;
    if (npar > 0) GetPositionIndices(pmb, npar, ppar->xp, ppar->yp, ppar->zp,
                                                ppar->xi1, ppar->xi2, ppar->xi3);
    pmb = pmb->next;
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
  const Real X1MIN = pmy_block->block_size.x1min;
  const Real X1MAX = pmy_block->block_size.x1max;
  const Real X2MIN = pmy_block->block_size.x2min;
  const Real X2MAX = pmy_block->block_size.x2max;
  const Real X3MIN = pmy_block->block_size.x3min;
  const Real X3MAX = pmy_block->block_size.x3max;

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
    // Convert to the MeshBlock coordinates.
    Real x1, x2, x3;
    _CartesianToMeshCoords(xp(k), yp(k), zp(k), x1, x2, x3);

    // Check if a particle is outside the boundary.
    int ox1 = _CheckSide(NX1, x1, X1MIN, X1MAX),
        ox2 = _CheckSide(NX2, x2, X2MIN, X2MAX),
        ox3 = _CheckSide(NX3, x3, X3MIN, X3MAX);
    if (ox1 == 0 && ox2 == 0 && ox3 == 0) {
      ++k;
      continue;
    }

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
            pnmb = pm->FindMeshBlock(pnmbt->pleaf[i][j][k]->gid);
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

    // Apply boundary conditions.
    ApplyBoundaryConditions(pm, x1, x2, x3);

    // Convert back to Cartesian coordinates.
    _MeshCoordsToCartesian(x1, x2, x3, xp(k), yp(k), zp(k));

    // Check the buffer size of the target MeshBlock.
    if (pnp->nrecv >= pnp->nrecvmax) {
      pnp->nrecvmax *= 2;
      pnp->irecv.ResizeLastDimension(pnp->nrecvmax);
      pnp->rrecv.ResizeLastDimension(pnp->nrecvmax);
    }

    // Copy the properties of the particle to the neighbor.
    for (int j = 0; j < nint; ++j)
      pnp->irecv(j,pnp->nrecv) = intprop(j,k);
    for (int j = 0; j < nreal; ++j)
      pnp->rrecv(j,pnp->nrecv) = realprop(j,k);
    ++pnp->nrecv;

    // Pop the particle from the current MeshBlock.
    if (--npar != k) {
      for (int j = 0; j < nint; ++j)
        intprop(j,k) = intprop(j,npar);
      for (int j = 0; j < nreal; ++j)
        realprop(j,k) = realprop(j,npar);
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::FlushReceiveBuffer()
//  \brief adds particles from the receive buffer.

void Particles::FlushReceiveBuffer()
{
  // Check the memory size.
  if (npar + nrecv > nparmax) {
    nparmax += 2 * (npar + nrecv - nparmax);  // Increase maximum number of particles
    intprop.ResizeLastDimension(nparmax);     // Increase size of property arrays
    realprop.ResizeLastDimension(nparmax);
    AssignShorthands();
    xi1.ResizeLastDimension(nparmax);         // Increase size of index arrays
    xi2.ResizeLastDimension(nparmax);
    xi3.ResizeLastDimension(nparmax);
  }

  // Flush the receive buffers.
  for (int j = 0; j < nint; ++j) {
    long ip = npar, k = 0;
    while (k < nrecv)
      intprop(j,ip++) = irecv(j,k++);
  }
  for (int j = 0; j < nreal; ++j) {
    long ip = npar, k = 0;
    while (k < nrecv)
      realprop(j,ip++) = rrecv(j,k++);
  }
  npar += nrecv;

  // Clear the receive buffers.
  nrecv = 0;
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::Update(Mesh *pm)
//  \brief updates all particle positions and velocities from t to t + dt.
//======================================================================================

void Particles::Update(Mesh *pm)
{
  MeshBlock *pmb = pm->pblock;
  Real dth = 0.5 * pm->dt;

  // Drift particles.
  while (pmb != NULL) {
    pmb->ppar->Drift(pm->time, dth);
    pmb = pmb->next;
  }
  Migrate(pm);

  // Kick particles.
  pmb = pm->pblock;
  while (pmb != NULL) {
    pmb->ppar->Kick(pm->time, pm->dt);
    pmb = pmb->next;
  }

  // Drift particles.
  pmb = pm->pblock;
  while (pmb != NULL) {
    pmb->ppar->Drift(pm->time + dth, dth);
    pmb = pmb->next;
  }
  Migrate(pm);
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::AddIntProperty()
//  \brief adds one integer property to the particles and returns the index.

int Particles::AddIntProperty()
{
  _ErrorIfInitialized("Particles::AddIntProperty", initialized);
  return nint++;
}

//--------------------------------------------------------------------------------------
//! \fn int Particles::AddRealProperty()
//  \brief adds one real property to the particles and returns the index.

int Particles::AddRealProperty()
{
  _ErrorIfInitialized("Particles::AddRealProperty", initialized);
  return nreal++;
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
    os.open(fname.str());
    if (!os.is_open()) {
      msg << "### FATAL ERROR in function [Particles::FormattedTableOutput]"
          << std::endl << "Output file '" << fname.str() << "' could not be opened"
          << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }

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
//! \fn void _ErrorIfInitialized(const string& calling_function)
//  \brief throws an error when the class has already been initialized.

void _ErrorIfInitialized(const std::string& calling_function, bool initialized)
{
  if (initialized)
  {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [" << calling_function << "]" << std::endl << "The Particles class has already been initialized. " << std::endl;
    throw std::runtime_error(msg.str().c_str());
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
  xi1 = IS + (x1 - block_size.x1min) / pcoord->dx1f(IS);
  xi2 = JS + (x2 - block_size.x2min) / pcoord->dx2f(JS);
  xi3 = KS + (x3 - block_size.x3min) / pcoord->dx3f(KS);
}

//--------------------------------------------------------------------------------------
//! \fn int _CheckSide(int nx, Real x, Real xmin, Real xmax)
//  \brief returns -1 if x < xmin and nx > 1, +1 if x > xmax and nx > 1,
//         and 0 otherwise.

inline int _CheckSide(int nx, Real x, Real xmin, Real xmax)
{
   if (nx > 1) {
     if (x < xmin) return -1;
     if (x > xmax) return +1;
   }
   return 0;
}
