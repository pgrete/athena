//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particle-mesh.hpp
//  \brief defines ParticleMesh class used for communication between meshblocks needed 
//         by particle-mesh methods.

// C++ standard library
#include <cmath>

// Athena++ classes headers
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../mesh/mesh.hpp"

// Particle-mesh constants.
const Real RINF = 1;  // radius of influence
const int NGPM = int(std::ceil(RINF));   // number of ghost cells needed.

// Forward declaration
class Particles;

//--------------------------------------------------------------------------------------
//! \class ParticleMesh
//  \brief defines the class for particle-mesh methods

class ParticleMesh {

friend class Particles;
friend class DustParticles;

public:
  // Constructor and destructor
  ParticleMesh(Particles *ppar, int nmeshaux);
  ~ParticleMesh();

protected:
  // Instance Variables
  AthenaArray<Real> meshaux;   // auxiliaries to the meshblock
  int nmeshaux;                // number of auxiliaries to the meshblock
  int is, ie, js, je, ks, ke;  // beginning and ending indices
  int ncells;                  // number of cells in meshaux

  // Instance methods
  void InterpolateMeshToParticles(
           const AthenaArray<Real>& meshsrc, const AthenaArray<int>& imeshsrc,
           AthenaArray<Real>& par, const AthenaArray<int>& ipar);
  void AssignParticlesToMeshAux(
           const AthenaArray<Real>& par, const AthenaArray<int>& ipar,
           const AthenaArray<int>& imeshaux);
  void InterpolateMeshAndAssignParticles(
           const AthenaArray<Real>& meshsrc, const AthenaArray<int>& imeshsrc,
           AthenaArray<Real>& pardst, const AthenaArray<int>& ipardst,
           const AthenaArray<Real>& parsrc, const AthenaArray<int>& iparsrc,
           const AthenaArray<int>& imeshaux);
  void DepositMeshAux(AthenaArray<Real>& u,
           const AthenaArray<int>& imeshaux, const AthenaArray<int>& imeshblock);
  void SendBoundary();
  void ReceiveBoundary();

private:
  // Instance Variables
  bool active1_, active2_, active3_;  // active dimensions
  Real dxi1_, dxi2_, dxi3_;           // range of influence from a particle cloud

  Particles *ppar_;        // ptr to my Particles instance
  MeshBlock *pmb_;         // ptr to my MeshBlock
  Mesh *pmesh_;            // ptr to my Mesh
  BoundaryValues *pbval_;  // ptr to my BoundaryValues
  BoundaryData bd_;        // boundary data

  // Instance methods
  void AssignParticlesToDifferentLevels(
           const AthenaArray<Real>& par, const AthenaArray<int>& ipar,
           const AthenaArray<int>& imeshaux);
  int LoadBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb);
  void AddBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb);

};
