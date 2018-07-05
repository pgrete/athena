//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file pm_bvals.cpp
//  \brief implements ParticleMeshBoundaryValues class used for communication between
//         meshblocks needed by particle-mesh methods.

// Standard library
#include <cstring>

// Athena++ classes headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../utils/buffer_utils.hpp"
#include "particles.hpp"

//--------------------------------------------------------------------------------------
//! \fn ParticleMesh::ParticleMesh(int nmeshaux, MeshBlock *pmb)
//  \brief constructs a new ParticleMesh instance.

ParticleMesh::ParticleMesh(int nmeshaux, MeshBlock *pmb)
{
  // Save the inputs.
  pmb_ = pmb;
  pbval_ = pmb->pbval;
  nmeshaux_ = nmeshaux;

  // Determine the dimensions of the block needed.
  RegionSize block_size = pmb->block_size;
  int dim = 0, ncells1 = 1, ncells2 = 1, ncells3 = 1;

  if (block_size.nx1 > 1) {
    ++dim;
    is_ = NGPM;
    ie_ = NGPM + block_size.nx1 - 1;
    ncells1 = block_size.nx1 + 2 * NGPM;
  } else
    is_ = ie_ = 0;

  if (block_size.nx2 > 1) {
    ++dim;
    js_ = NGPM;
    je_ = NGPM + block_size.nx2 - 1;
    ncells2 = block_size.nx2 + 2 * NGPM;
  } else
    js_ = je_ = 0;

  if (block_size.nx3 > 1) {
    ++dim;
    ks_ = NGPM;
    ke_ = NGPM + block_size.nx3 - 1;
    ncells3 = block_size.nx3 + 2 * NGPM;
  } else
    ks_ = ke_ = 0;

  // Allocate the block for particle-mesh.
  meshaux_.NewAthenaArray(nmeshaux, ncells3, ncells2, ncells1);

  // Find the number of neighbors.
  bd_.nbmax = BoundaryBase::BufferID(dim, pmb->pmy_mesh->multilevel);

  // Initialize boundary data.
  for (int n = 0; n < bd_.nbmax; n++) {
    bd_.flag[n] = BNDRY_WAITING;
    bd_.send[n] = NULL;
    bd_.recv[n] = NULL;

    int size = ((pbval_->ni[n].ox1 == 0) ? pmb->block_size.nx1 : NGPM) *
               ((pbval_->ni[n].ox2 == 0) ? pmb->block_size.nx2 : NGPM) *
               ((pbval_->ni[n].ox3 == 0) ? pmb->block_size.nx3 : NGPM) * nmeshaux;
    bd_.send[n] = new Real [size];
    bd_.recv[n] = new Real [size];
  }
}

//--------------------------------------------------------------------------------------
//! \fn ParticleMesh::~ParticleMesh()
//  \brief destructs a ParticleMesh instance.

ParticleMesh::~ParticleMesh()
{
  // Destroy the particle meshblock.
  meshaux_.DeleteAthenaArray();

  // Destroy boundary data.
  for (int n = 0; n < bd_.nbmax; n++) {
    delete [] bd_.send[n];
    delete [] bd_.recv[n];
  }
}

//--------------------------------------------------------------------------------------
//! \fn int ParticleMesh::LoadBoundaryBufferSameLevel(
//              Real *buf, const NeighborBlock& nb)
//  \brief Fill boundary buffers for sending to a block on the same level

int ParticleMesh::LoadBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb)
{
  // Determine the chunk of the block to be communicated.
  int si = is_, sj = js_, sk = ks_, ei = ie_, ej = je_, ek = ke_;

  if (nb.ox1 > 0) {
    si = ie_ + 1;
    ei += NGPM;
  } else if (nb.ox1 < 0) {
    si -= NGPM;
    ei = is_ - 1;
  }

  if (nb.ox2 > 0) {
    sj = je_ + 1;
    ej += NGPM;
  } else if (nb.ox2 < 0) {
    sj -= NGPM;
    ej = js_ - 1;
  }

  if (nb.ox3 > 0) {
    sk = ke_ + 1;
    ek += NGPM;
  } else if (nb.ox3 < 0) {
    sk -= NGPM;
    ek = ks_ - 1;
  }

  // Load the data to the buffer.
  int p = 0;
  BufferUtility::Pack4DData(meshaux_, buf, 0, nmeshaux_ - 1, si, ei, sj, ej, sk, ek, p);
  return p;
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::SendBoundary()
//  \brief Send boundary values to neighboring blocks.

void ParticleMesh::SendBoundary()
{
  int mylevel = pmb_->loc.level;

  for (int n = 0; n < pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];

    // Load boundary values to send buffer.
    int ssize;
    if (nb.level == mylevel)
      ssize = LoadBoundaryBufferSameLevel(bd_.send[nb.bufid], nb);

    // Receive boundary values from neighboring blocks.
    if (nb.rank == Globals::my_rank) {
      MeshBlock *pnb = pmb_->pmy_mesh->FindMeshBlock(nb.gid);
      BoundaryData *ptarget = &(pnb->ppar->ppm->bd_);
      std::memcpy(ptarget->recv[nb.targetid], bd_.send[nb.bufid], ssize*sizeof(Real));
      ptarget->flag[nb.targetid] = BNDRY_ARRIVED;
    }
  }
}
