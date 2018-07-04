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

// Particle-mesh constants.
const Real RINF = 1;  // radius of influence
const int NGPM = 1;   // number of ghost cells needed.

//--------------------------------------------------------------------------------------
//! \class ParticleMeshBoundaryValues
//  \brief defines BoundaryValues class for particle-mesh methods

class ParticleMeshBoundaryValues : public BoundaryBase {

public:
  // Constructor and destructor
  ParticleMeshBoundaryValues(int nval, MeshBlock *pmb, enum BoundaryFlag *input_bcs);
  ~ParticleMeshBoundaryValues();

private:
  // Instance Variables
  MeshBlock *pmy_block_;  // ptr to my meshblock
  BoundaryData bd_;       // boundary data

};
