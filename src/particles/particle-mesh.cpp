//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file pm_bvals.cpp
//  \brief implements ParticleMeshBoundaryValues class used for communication between
//         meshblocks needed by particle-mesh methods.

// Athena++ classes headers
#include "../athena.hpp"
#include "particle-mesh.hpp"

//--------------------------------------------------------------------------------------
//! \fn ParticleMeshBoundaryValues::ParticleMeshBoundaryValues(
//          MeshBlock *pmb, enum BoundaryFlag *input_bcs)
//  \brief constructs a new ParticleMeshBoundaryValues instance.

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
  meshaux.NewAthenaArray(nmeshaux, ncells3, ncells2, ncells1);

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
//! \fn ParticleMeshBoundaryValues::~ParticleMeshBoundaryValues()
//  \brief destructs a ParticleMeshBoundaryValues instance.

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
