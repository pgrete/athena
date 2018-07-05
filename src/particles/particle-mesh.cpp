//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particle-mesh.cpp
//  \brief implements ParticleMesh class used for operations involved in particle-mesh 
//         methods.

// Standard library
#include <cstring>
#include <sstream>

// Athena++ classes headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../utils/buffer_utils.hpp"
#include "particles.hpp"

// Local function prototypes.
Real _ParticleMeshWeightFunction(Real dxi);

//--------------------------------------------------------------------------------------
//! \fn ParticleMesh::ParticleMesh(Particles *ppar, int nmeshaux)
//  \brief constructs a new ParticleMesh instance.

ParticleMesh::ParticleMesh(Particles *ppar, int nmeshaux)
{
  // Save some inputs.
  ppar_ = ppar;
  pmb_ = ppar->pmy_block;
  pbval_ = pmb_->pbval;
  nmeshaux_ = nmeshaux;

  // Determine active dimensions.
  RegionSize block_size = pmb_->block_size;
  active1_ = block_size.nx1 > 1;
  active2_ = block_size.nx2 > 1;
  active3_ = block_size.nx3 > 1;

  // Determine the range of a particle cloud.
  dxi1_ = active1_ ? RINF : 0;
  dxi2_ = active2_ ? RINF : 0;
  dxi3_ = active3_ ? RINF : 0;

  // Determine the dimensions of the block for boundary communication.
  int dim = 0, ncells1 = 1, ncells2 = 1, ncells3 = 1;

  if (active1_) {
    ++dim;
    is_ = NGPM;
    ie_ = NGPM + block_size.nx1 - 1;
    ncells1 = block_size.nx1 + 2 * NGPM;
  } else
    is_ = ie_ = 0;

  if (active2_) {
    ++dim;
    js_ = NGPM;
    je_ = NGPM + block_size.nx2 - 1;
    ncells2 = block_size.nx2 + 2 * NGPM;
  } else
    js_ = je_ = 0;

  if (active3_) {
    ++dim;
    ks_ = NGPM;
    ke_ = NGPM + block_size.nx3 - 1;
    ncells3 = block_size.nx3 + 2 * NGPM;
  } else
    ks_ = ke_ = 0;

  // Allocate the block for particle-mesh.
  meshaux_.NewAthenaArray(nmeshaux, ncells3, ncells2, ncells1);

  // Find the number of neighbors.
  bd_.nbmax = BoundaryBase::BufferID(dim, pmb_->pmy_mesh->multilevel);

  // Initialize boundary data.
  for (int n = 0; n < bd_.nbmax; n++) {
    bd_.flag[n] = BNDRY_WAITING;
    bd_.send[n] = NULL;
    bd_.recv[n] = NULL;

    int size = ((pbval_->ni[n].ox1 == 0) ? block_size.nx1 : NGPM) *
               ((pbval_->ni[n].ox2 == 0) ? block_size.nx2 : NGPM) *
               ((pbval_->ni[n].ox3 == 0) ? block_size.nx3 : NGPM) * nmeshaux;
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
//! \fn void ParticleMesh::InterpolateMeshToParticles(
//               const AthenaArray<int>& auxindices,
//               const AthenaArray<Real>& meshprop,
//               const AthenaArray<int>& meshindices)
//  \brief interpolates meshprop at specified property indices meshindices onto
//         auxprop of particles at the corresponding property indices auxindices.

void ParticleMesh::InterpolateMeshToParticles(
         const AthenaArray<int>& auxindices,
         const AthenaArray<Real>& meshprop,
         const AthenaArray<int>& meshindices)
{
  // Check the index mapping.
  int nprop = meshindices.GetSize();
  if (nprop <= 0 || auxindices.GetSize() != nprop) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [Particles::InterpolateMeshToParticles]"
        << std::endl
        << "index arrays meshindices and auxindices does not match." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  // Loop over each particle.
  for (long k = 0; k < ppar_->npar; ++k) {
    for (int i = 0; i < nprop; ++i)
      ppar_->auxprop(auxindices(i),k) = 0;

    // Find the domain the particle influences.
    Real xi1 = ppar_->xi1(k), xi2 = ppar_->xi2(k), xi3 = ppar_->xi3(k);
    int ix1s = int(xi1 - dxi1_), ix1e = int(xi1 + dxi1_);
    int ix2s = int(xi2 - dxi2_), ix2e = int(xi2 + dxi2_);
    int ix3s = int(xi3 - dxi3_), ix3e = int(xi3 + dxi3_);

    // Weight each cell and accumulate the mesh properties onto the particles.
    for (int ix3 = ix3s; ix3 <= ix3e; ++ix3) {
      Real w3 = active3_ ? _ParticleMeshWeightFunction(ix3 + 0.5 - xi3) : 1.0;

      for (int ix2 = ix2s; ix2 <= ix2e; ++ix2) {
        Real w23 = w3 * (active2_ ? _ParticleMeshWeightFunction(ix2 + 0.5 - xi2) : 1.0);

        for (int ix1 = ix1s; ix1 <= ix1e; ++ix1) {
          Real weight = w23 * (active1_ ?
                                    _ParticleMeshWeightFunction(ix1 + 0.5 - xi1) : 1.0);

          for (int i = 0; i < nprop; ++i)
            ppar_->auxprop(auxindices(i),k) +=
                weight * meshprop(meshindices(i),ix3,ix2,ix1);
        }
      }
    }
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
//! \fn void ParticleMesh::AddBoundaryBufferSameLevel(
//              Real *buf, const NeighborBlock& nb)
//  \brief Add boundary buffers from a neighbor block to mesh.

void ParticleMesh::AddBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb)
{
  // Determine the chunk of the block to be added to.
  int si = (nb.ox1 > 0) ? (ie_ - NGPM + 1) : is_,
      ei = (nb.ox1 < 0) ? (is_ + NGPM - 1) : ie_,
      sj = (nb.ox2 > 0) ? (je_ - NGPM + 1) : js_,
      ej = (nb.ox2 < 0) ? (js_ + NGPM - 1) : je_,
      sk = (nb.ox3 > 0) ? (ke_ - NGPM + 1) : ks_,
      ek = (nb.ox3 < 0) ? (ks_ + NGPM - 1) : ke_;

  // Add the data to the mesh.
  for (int n = 0; n < nmeshaux_; ++n)
    for (int k = sk; k <= ek; ++k)
      for (int j = sj; j <= ej; ++j)
        for (int i = si; i <= ei; ++i)
          meshaux_(n,k,j,i) += *buf++;
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

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::ReceiveBoundary()
//  \brief Receive boundary values from neighboring blocks and add to my block.

void ParticleMesh::ReceiveBoundary()
{
  int mylevel = pmb_->loc.level;

  for (int n = 0; n < pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (bd_.flag[nb.bufid] == BNDRY_COMPLETED) continue;

    if (nb.level == mylevel)
      AddBoundaryBufferSameLevel(bd_.recv[nb.bufid], nb);

    bd_.flag[nb.bufid] = BNDRY_COMPLETED;
  }
}

//--------------------------------------------------------------------------------------
//! \fn Real _ParticleMeshWeightFunction(Real dxi)
//  \brief evaluates the weight function given index distance.

Real _ParticleMeshWeightFunction(Real dxi)
{
  if (dxi < 0) dxi = -dxi;

  if (dxi < 0.5)
    return 0.75 - dxi * dxi;

  if (dxi < 1.5) {
    dxi = 1.5 - dxi;
    return 0.5 * (dxi * dxi);
  }

  return 0;
}
