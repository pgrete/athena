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

ParticleMeshBoundaryValues::ParticleMeshBoundaryValues(
    MeshBlock *pmb, enum BoundaryFlag *input_bcs)
 : BoundaryBase(pmb->pmy_mesh, pmb->loc, pmb->block_size, input_bcs)
{
}

//--------------------------------------------------------------------------------------
//! \fn ParticleMeshBoundaryValues::~ParticleMeshBoundaryValues()
//  \brief destructs a ParticleMeshBoundaryValues instance.

ParticleMeshBoundaryValues::~ParticleMeshBoundaryValues()
{
}
