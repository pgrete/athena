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
  static void Initialize();
  static void FormattedTableOutput(Mesh *pm, OutputParameters op);
  static void Update(Mesh *pm);  // master integrator

  // Constructor
  Particles(MeshBlock *pmb, ParameterInput *pin);

  // Destructor
  ~Particles();

  // Instance methods
  void Drift(Real t, Real dt);
  void Kick(Real t, Real dt);

  // Instance variables
  AthenaArray<long> pid;            // shorthand for particle ID
  AthenaArray<Real> xp1, xp2, xp3;  // shorthand for position components
  AthenaArray<Real> vp1, vp2, vp3;  // shorthand for velocity components

  long npar;  // number of particles

private:
  // Class variables
  static int ipid;              // index for the particle ID
  static int ixp1, ixp2, ixp3;  // indices for the position components
  static int ivp1, ivp2, ivp3;  // indices for the velocity components

  // Class methods
  static int AddIntProperty();
  static int AddRealProperty();

  // Class variables
  static bool initialized;      // whether or not the class is initialized
  static int nint, nreal;       // numbers of integer and real properties

  // Instance variables
  AthenaArray<long> intprop;   // integer particle properties
  AthenaArray<Real> realprop;  // real particle properties

  long nparmax;            // maximum number of particles per meshblock
  MeshBlock* pmy_block;    // MeshBlock pointer
};

#endif
