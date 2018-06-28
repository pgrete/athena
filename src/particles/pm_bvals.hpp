//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file pm_bvals.hpp
//  \brief defines ParticleMeshBoundaryValues class used for communication between
//         meshblocks needed by particle-mesh methods.

// Athena++ classes headers
#include "../bvals/bvals.hpp"
#include "../mesh/mesh.hpp"

//--------------------------------------------------------------------------------------
//! \class ParticleMeshBoundaryValues
//  \brief BoundaryValues class for particle-mesh methods

class ParticleMeshBoundaryValues : public BoundaryBase {

public:
  // Constructor and destructor
  ParticleMeshBoundaryValues(MeshBlock *pmb, enum BoundaryFlag *input_bcs);
  ~ParticleMeshBoundaryValues();

};
