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

friend class MeshBlock;  // Make writing initial conditions possible.

public:
  // Class methods
  static void Initialize();
  static void Integrate(Mesh *pm, int step);  // master integrator
  static void FormattedTableOutput(Mesh *pm, OutputParameters op); 
  static void Migrate(Mesh *pm);

  // Constructor
  Particles(MeshBlock *pmb, ParameterInput *pin);

  // Destructor
  ~Particles();

  // Instance methods
  size_t GetSizeInBytes();
  void ReadRestart(char *mbdata, int &os);
  void WriteRestart(char *&pdata);

protected:
  // Class methods
  static int AddIntProperty();
  static int AddRealProperty();
  static int AddAuxProperty();

  // Class variables
  static bool initialized;  // whether or not the class is initialized
  static int nint, nreal;   // numbers of integer and real properties
  static int naux;          // number of auxiliary properties

  static int ipid;                 // index for the particle ID
  static int ixp, iyp, izp;        // indices for the position components
  static int ivpx, ivpy, ivpz;     // indices for the velocity components

  static int iapx, iapy, iapz;     // indices for the acceleration compoenents
  static int ixp0, iyp0, izp0;     // indices for beginning position components
  static int ivpx0, ivpy0, ivpz0;  // indices for beginning velocity components

  // Instance methods
  virtual void AssignShorthands();  // Needs to be called everytime
                                    // intprop, realprop, & auxprop are resized
                                    // Be sure to call back when derived.
  virtual void AddAcceleration(Real t, Real dt) = 0;
  void InterpolateMeshToParticles(
           const AthenaArray<Real>& meshprop, const AthenaArray<int>& meshindices,
           const AthenaArray<int>& auxindices);

  // Instance variables
  AthenaArray<long> intprop;   // integer particle properties
  AthenaArray<Real> realprop;  // real particle properties
  AthenaArray<Real> auxprop;   // auxiliary particle properties (for intermediate
                               //     computations)

  AthenaArray<Real> xi1, xi2, xi3;  // position indices of each particle
                                    // in local meshblock

  AthenaArray<long> pid;            // shorthands for particle ID
  AthenaArray<Real> xp, yp, zp;     // shorthands for position components
  AthenaArray<Real> vpx, vpy, vpz;  // shorthands for velocity components

  AthenaArray<Real> apx, apy, apz;     // shorthands for acceleration components
  AthenaArray<Real> xp0, yp0, zp0;     // shorthands for beginning position components
  AthenaArray<Real> vpx0, vpy0, vpz0;  // shorthands for beginning velocity components

  long npar;             // number of particles
  long nparmax;          // maximum number of particles per meshblock
  MeshBlock* pmy_block;  // MeshBlock pointer

private:
  // Class methods
  static bool ApplyBoundaryConditions(Mesh *pm, Real &x1, Real &x2, Real &x3);
  static void GetPositionIndices(MeshBlock *pmb, long npar,
                                 const AthenaArray<Real>& xp,
                                 const AthenaArray<Real>& yp,
                                 const AthenaArray<Real>& zp,
                                 AthenaArray<Real>& xi1,
                                 AthenaArray<Real>& xi2,
                                 AthenaArray<Real>& xi3);

  // Instance methods
  void EulerStep(Real t, Real dt);
  void FlushReceiveBuffer();
  void SaveStatus();
  void SendToNeighbors();

  // MeshBlock-to-MeshBlock communication:
  AthenaArray<long> irecv;  //   integer receive buffers
  AthenaArray<Real> rrecv;  //   real receive buffers
  int nrecvmax;             //   maximum number of particles per receive buffer
  int nrecv;                //   actual number of particles per receive buffer

  // Particle-mesh related properties:
  Real pm_dxi1, pm_dxi2, pm_dxi3;  // maximum index distance from a particle in each dir
};

//--------------------------------------------------------------------------------------
//! \class DustParticles
//  \brief defines the class for dust particles that interact with the gas via drag
//         force.

class DustParticles : public Particles {

public:
  // Class methods
  static void Initialize();

  // Constructor
  DustParticles(MeshBlock *pmb, ParameterInput *pin);

  // Destructor
  ~DustParticles();

private:
  // Class variables
  static bool initialized;   // whether or not the class is initialized
  static int iux, iuy, iuz;  // indices for the gas velocity
  static AthenaArray<int> pm_meshindices, pm_auxindices;
                             // Array of indices for particle-mesh mapping

  // Instance methods.
  void AssignShorthands();
  void AddAcceleration(Real t, Real dt);

  // Instance variables
  Real mass;  // mass of each particle
  Real taus;  // stopping time (in code units)

  AthenaArray<Real> ux, uy, uz;  // shorthand for gas velocity components
};

#endif
