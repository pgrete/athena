#ifndef PARTICLE_HPP
#define PARTICLE_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
//======================================================================================
//! \file particles.hpp
//  \brief defines classes for particle dynamics.
//======================================================================================

// Athena headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../outputs/outputs.hpp"

class ParameterInput;

//--------------------------------------------------------------------------------------
//! \class Particles
//  \brief defines the bass class for all implementations of particles.

class Particles {
public:
  // Class methods
  static void FormattedTableOutput(Mesh *pm, OutputParameters op);
  static void Update(Mesh *pm);  // master integrator

  // Constructor
  Particles(MeshBlock *pmb, ParameterInput *pin);

  // Destructor
  ~Particles();

  // Instance methods
  void Drift(Real t, Real dt);
  void Kick(Real t, Real dt);

  // Particle properties
  AthenaArray<long> id;            // particle id
  AthenaArray<Real> x1, x2, x3;    // particle position
  AthenaArray<Real> v1, v2, v3;    // particle velocity

  // Bookkeeping indices
  MeshBlock* pmy_block;    // MeshBlock pointer
  long npar;               // number of particles
  long nparmax;            // maximum number of particles per meshblock
};

#endif
