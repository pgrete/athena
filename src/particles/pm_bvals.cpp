//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file pm_bvals.cpp
//  \brief implements ParticleMeshBoundaryValues class used for communication between
//         meshblocks needed by particle-mesh methods.

// Athena++ classes headers
#include "pm_bvals.hpp"

//--------------------------------------------------------------------------------------
//! \fn ParticleMeshBoundaryValues::ParticleMeshBoundaryValues(
//          MeshBlock *pmb, enum BoundaryFlag *input_bcs)
//  \brief constructs a new ParticleMeshBoundaryValues instance.

ParticleMeshBoundaryValues::ParticleMeshBoundaryValues(int nval,
    MeshBlock *pmb, enum BoundaryFlag *input_bcs)
 : BoundaryBase(pmb->pmy_mesh, pmb->loc, pmb->block_size, input_bcs)
{
  // Initialize boundary data.
  bd_.nbmax = maxneighbor_;

  for (int n = 0; n < bd_.nbmax; n++) {
    bd_.flag[n] = BNDRY_WAITING;
    bd_.send[n] = NULL;
    bd_.recv[n] = NULL;

    int size = ((ni[n].ox1 == 0) ? pmb->block_size.nx1 : NGPM) *
               ((ni[n].ox2 == 0) ? pmb->block_size.nx2 : NGPM) *
               ((ni[n].ox3 == 0) ? pmb->block_size.nx3 : NGPM) * nval;
    bd_.send[n] = new Real [size];
    bd_.recv[n] = new Real [size];
  }
}

//--------------------------------------------------------------------------------------
//! \fn ParticleMeshBoundaryValues::~ParticleMeshBoundaryValues()
//  \brief destructs a ParticleMeshBoundaryValues instance.

ParticleMeshBoundaryValues::~ParticleMeshBoundaryValues()
{
  // Destroy boundary data.
  for (int n = 0; n < bd_.nbmax; n++) {
    delete [] bd_.send[n];
    delete [] bd_.recv[n];
  }
}
