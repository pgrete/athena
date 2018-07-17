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
#include "../coordinates/coordinates.hpp"
#include "../utils/buffer_utils.hpp"
#include "particles.hpp"

// Local function prototypes.
Real _ParticleMeshWeightFunction(Real dxi);

// Local constants.
const int OFFSET = NGHOST - NGPM;  // offset between meshblock and meshaux

//--------------------------------------------------------------------------------------
//! \fn ParticleMesh::ParticleMesh(Particles *ppar, int nmeshaux)
//  \brief constructs a new ParticleMesh instance.

ParticleMesh::ParticleMesh(Particles *ppar, int nmeshaux_)
{
  // Save some inputs.
  ppar_ = ppar;
  pmb_ = ppar->pmy_block;
  pbval_ = pmb_->pbval;
  nmeshaux = nmeshaux_;

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
  int dim = 0, nx1 = 1, nx2 = 1, nx3 = 1;

  if (active1_) {
    ++dim;
    is = NGPM;
    ie = NGPM + block_size.nx1 - 1;
    nx1 = block_size.nx1 + 2 * NGPM;
  } else
    is = ie = 0;

  if (active2_) {
    ++dim;
    js = NGPM;
    je = NGPM + block_size.nx2 - 1;
    nx2 = block_size.nx2 + 2 * NGPM;
  } else
    js = je = 0;

  if (active3_) {
    ++dim;
    ks = NGPM;
    ke = NGPM + block_size.nx3 - 1;
    nx3 = block_size.nx3 + 2 * NGPM;
  } else
    ks = ke = 0;

  // Allocate the block for particle-mesh.
  meshaux.NewAthenaArray(nmeshaux, nx3, nx2, nx1);
  ncells = nx1 * nx2 * nx3;

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
  meshaux.DeleteAthenaArray();

  // Destroy boundary data.
  for (int n = 0; n < bd_.nbmax; n++) {
    delete [] bd_.send[n];
    delete [] bd_.recv[n];
  }
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::InterpolateMeshToParticles(
//               const AthenaArray<Real>& meshsrc, const AthenaArray<int>& imeshsrc,
//               AthenaArray<Real>& par, const AthenaArray<int>& ipar)
//  \brief interpolates meshsrc at specified indices imeshsrc onto particle array par
//         (realprop, auxprop, or work in Particles class) at the corresponding indices
//         ipar.

void ParticleMesh::InterpolateMeshToParticles(
         const AthenaArray<Real>& meshsrc, const AthenaArray<int>& imeshsrc,
         AthenaArray<Real>& par, const AthenaArray<int>& ipar)
{
  // Check the index mapping.
  int nprop = imeshsrc.GetSize();
  if (nprop <= 0 || ipar.GetSize() != nprop) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [Particles::InterpolateMeshToParticles]"
        << std::endl
        << "index arrays imeshsrc and ipar do not match." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  // Zero out the particle arrays.
  for (int n = 0; n < nprop; ++n) {
    Real* pdata = &par(ipar(n),0);
    for (long k = 0; k < ppar_->npar; ++k)
      *pdata++ = 0.0;
  }

  // Loop over each particle.
  for (long k = 0; k < ppar_->npar; ++k) {
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

          for (int n = 0; n < nprop; ++n)
            par(ipar(n),k) += weight * meshsrc(imeshsrc(n),ix3,ix2,ix1);
        }
      }
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::AssignParticlesToMeshAux(
//               const AthenaArray<Real>& par, const AthenaArray<int>& ipar,
//               const AthenaArray<int>& imeshaux)
//  \brief assigns par (realprop, auxprop, or work in Particles class) at specified
//         indices ipar onto meshaux at the corresponding indices imeshaux.

void ParticleMesh::AssignParticlesToMeshAux(
         const AthenaArray<Real>& par, const AthenaArray<int>& ipar,
         const AthenaArray<int>& imeshaux)
{
  // Check the index mapping.
  int nprop = ipar.GetSize();
  if (nprop <= 0 || imeshaux.GetSize() != nprop) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [Particles::AssignParticlesToMeshAux]"
        << std::endl
        << "index arrays ipar and imeshaux does not match." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  // Zero out meshaux.
  for (int n = 0; n < nprop; ++n) {
    Real* pdata = &meshaux(imeshaux(n),0,0,0);
    for (int i = 0; i < ncells; ++i)
      *pdata++ = 0.0;
  }

  // Allocate working array.
  AthenaArray<Real> p;
  p.NewAthenaArray(nprop);

  for (long k = 0; k < ppar_->npar; ++k) {
    // Find the domain the particle influences.
    Real xi1 = ppar_->xi1(k) - (active1_ ? OFFSET : 0),
         xi2 = ppar_->xi2(k) - (active2_ ? OFFSET : 0),
         xi3 = ppar_->xi3(k) - (active3_ ? OFFSET : 0);
    int ix1s = int(xi1 - dxi1_), ix1e = int(xi1 + dxi1_);
    int ix2s = int(xi2 - dxi2_), ix2e = int(xi2 + dxi2_);
    int ix3s = int(xi3 - dxi3_), ix3e = int(xi3 + dxi3_);

    // Copy the particle properties.
    for (int n = 0; n < nprop; ++n)
      p(n) = par(ipar(n),k);

    // Weight each cell and accumulate particle property onto meshaux.
    for (int ix3 = ix3s; ix3 <= ix3e; ++ix3) {
      Real w3 = active3_ ? _ParticleMeshWeightFunction(ix3 + 0.5 - xi3) : 1.0;

      for (int ix2 = ix2s; ix2 <= ix2e; ++ix2) {
        Real w23 = w3 * (active2_ ? _ParticleMeshWeightFunction(ix2 + 0.5 - xi2) : 1.0);

        for (int ix1 = ix1s; ix1 <= ix1e; ++ix1) {
          Real weight = w23 * (active1_ ?
                                    _ParticleMeshWeightFunction(ix1 + 0.5 - xi1) : 1.0);

          for (int n = 0; n < nprop; ++n)
            meshaux(imeshaux(n),ix3,ix2,ix1) += weight * p(n);
        }
      }
    }
  }

  // Release working array.
  p.DeleteAthenaArray();
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::DepositMeshAux(AthenaArray<Real>& u,
//               const AthenaArray<int>& imeshaux, const AthenaArray<int>& imeshblock)
//  \brief deposits data in meshaux at specified indices imeshaux to meshblock u at the
//         corresponding indices imeshblock, divided by cell volume.

void ParticleMesh::DepositMeshAux(AthenaArray<Real>& u,
         const AthenaArray<int>& imeshaux, const AthenaArray<int>& imeshblock)
{
  // Check the index mapping.
  const int nprop = imeshaux.GetSize();
  if (nprop <= 0 || imeshblock.GetSize() != nprop) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [Particles::AssignParticlesToMeshAux]"
        << std::endl
        << "index arrays ipar and imeshaux does not match." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  // Add meshaux to meshblock.
  Coordinates *pc = pmb_->pcoord;
  for (int n = 0; n < nprop; ++n) {
    int ima = imeshaux(n), imb = imeshblock(n);

    for (int ka = ks, kb = pmb_->ks; ka <= ke; ++ka, ++kb)
      for (int ja = js, jb = pmb_->js; ja <= je; ++ja, ++jb)
        for (int ia = is, ib = pmb_->is; ia <= ie; ++ia, ++ib)
          u(imb,kb,jb,ib) += meshaux(ima,ka,ja,ia) / pc->GetCellVolume(kb,jb,ib);
  }
}

//--------------------------------------------------------------------------------------
//! \fn int ParticleMesh::LoadBoundaryBufferSameLevel(
//              Real *buf, const NeighborBlock& nb)
//  \brief Fill boundary buffers for sending to a block on the same level

int ParticleMesh::LoadBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb)
{
  // Determine the chunk of the block to be communicated.
  int si = is, sj = js, sk = ks, ei = ie, ej = je, ek = ke;

  if (nb.ox1 > 0) {
    si = ie + 1;
    ei += NGPM;
  } else if (nb.ox1 < 0) {
    si -= NGPM;
    ei = is - 1;
  }

  if (nb.ox2 > 0) {
    sj = je + 1;
    ej += NGPM;
  } else if (nb.ox2 < 0) {
    sj -= NGPM;
    ej = js - 1;
  }

  if (nb.ox3 > 0) {
    sk = ke + 1;
    ek += NGPM;
  } else if (nb.ox3 < 0) {
    sk -= NGPM;
    ek = ks - 1;
  }

  // Load the data to the buffer.
  int p = 0;
  BufferUtility::Pack4DData(meshaux, buf, 0, nmeshaux - 1, si, ei, sj, ej, sk, ek, p);
  return p;
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::AddBoundaryBufferSameLevel(
//              Real *buf, const NeighborBlock& nb)
//  \brief Add boundary buffers from a neighbor block to mesh.

void ParticleMesh::AddBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb)
{
  // Determine the chunk of the block to be added to.
  int si = (nb.ox1 > 0) ? (ie - NGPM + 1) : is,
      ei = (nb.ox1 < 0) ? (is + NGPM - 1) : ie,
      sj = (nb.ox2 > 0) ? (je - NGPM + 1) : js,
      ej = (nb.ox2 < 0) ? (js + NGPM - 1) : je,
      sk = (nb.ox3 > 0) ? (ke - NGPM + 1) : ks,
      ek = (nb.ox3 < 0) ? (ks + NGPM - 1) : ke;

  // Add the data to the mesh.
  for (int n = 0; n < nmeshaux; ++n)
    for (int k = sk; k <= ek; ++k)
      for (int j = sj; j <= ej; ++j)
        for (int i = si; i <= ei; ++i)
          meshaux(n,k,j,i) += *buf++;
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
