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
#include "particles.hpp"

bool Particles::initialized = false;
int Particles::nint = 0;
int Particles::nreal = 0;
int Particles::ipid = -1;
int Particles::ixp = -1, Particles::iyp = -1, Particles::izp = -1;
int Particles::ivpx = -1, Particles::ivpy = -1, Particles::ivpz = -1;

void _ErrorIfInitialized(const std::string& calling_function, bool initialized);

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
  npar = 1;
  nparmax = 1;  // TODO: dynamically adjust npar and nparmax.

  // Allocate integer properties.
  intprop.NewAthenaArray(nint,nparmax);
  pid.InitWithShallowSlice(intprop, 1, ipid, 1);

  // Allocate integer properties.
  realprop.NewAthenaArray(nreal,nparmax);
  xp.InitWithShallowSlice(realprop, 1, ixp, 1);
  yp.InitWithShallowSlice(realprop, 1, iyp, 1);
  zp.InitWithShallowSlice(realprop, 1, izp, 1);
  vpx.InitWithShallowSlice(realprop, 1, ivpx, 1);
  vpy.InitWithShallowSlice(realprop, 1, ivpy, 1);
  vpz.InitWithShallowSlice(realprop, 1, ivpz, 1);

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

  // Delete buffers.
  irecv.DeleteAthenaArray();
  rrecv.DeleteAthenaArray();
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
  while (pmb != NULL) {
    if (pmb->ppar->nrecv > 0)
      pmb->ppar->FlushReceiveBuffer();
    pmb = pmb->next;
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::SendToNeighbors()
//  \brief sends particles outside boundary to the buffers of neighboring meshblocks.

void Particles::SendToNeighbors()
{
  Mesh *pm = pmy_block->pmy_mesh;
  MeshBlock *pnmb;
  MeshBlockTree *pnmbt;
  Particles *pnp;
  int ox1, ox2, ox3;

  // TODO: Currently only works for Cartesian.
  if (COORDINATE_SYSTEM != "cartesian") {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [Particles::SendToNeighbors]" << std::endl
        << "Non-Cartesian coordinates not yet implemented. " << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  for (long k = 0; k < npar; ++k) {
    // Check if a particle is outside the boundary.
    ox1 = _CheckSide(pmy_block->block_size.nx1, xp(k),
                     pmy_block->block_size.x1min, pmy_block->block_size.x1max);
    ox2 = _CheckSide(pmy_block->block_size.nx2, yp(k),
                     pmy_block->block_size.x2min, pmy_block->block_size.x2max);
    ox3 = _CheckSide(pmy_block->block_size.nx3, zp(k),
                     pmy_block->block_size.x3min, pmy_block->block_size.x3max);
    if (ox1 != 0 || ox2 != 0 || ox3 != 0) {

      // Find the neighbor MeshBlock to send it to.
      pnmbt = pm->tree.FindNeighbor(pmy_block->loc, ox1, ox2, ox3, pm->mesh_bcs,
                                    pm->nrbx1, pm->nrbx2, pm->nrbx3, pm->root_level);
      if (pnmbt == NULL) {
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::SendToNeighbors]" << std::endl
            << "cannot find the neighboring MeshBlock. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
        continue;
      }

      if (pnmbt->flag)  // Neighbor is on the same or a courser level:
        pnmb = pm->FindMeshBlock(pnmbt->gid);
      else {          // Neighbor is on a finer level:
        std::stringstream msg;
        msg << "### FATAL ERROR in function [Particles::SendToNeighbors]" << std::endl
            << "coarse-to-fine migration not yet implemented. " << std::endl;
        throw std::runtime_error(msg.str().c_str());
        continue;
      }
      pnp = pnmb->ppar;

      // Check the buffer size of the target MeshBlock.
      if (pnp->nrecv >= nrecvmax) {
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
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::FlushReceiveBuffer()
//  \brief adds particles from the receive buffer.

void Particles::FlushReceiveBuffer()
{
  // Check the memory size.
  if (npar + nrecv > nparmax) {
    nparmax += 2 * (npar + nrecv - nparmax);
    intprop.ResizeLastDimension(nparmax);
    realprop.ResizeLastDimension(nparmax);
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
    msg << "### FATAL ERROR in function [" << calling_function << "]" << std::endl
        << "The Particles class has already been initialized. " << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
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
