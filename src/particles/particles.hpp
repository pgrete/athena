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
  static void Update(Mesh *pm);  // master integrator
  static void FormattedTableOutput(Mesh *pm, OutputParameters op); 

  // Constructor
  Particles(MeshBlock *pmb, ParameterInput *pin);

  // Destructor
  ~Particles();

  // Instance methods
  void Drift(Real t, Real dt);
  void Kick(Real t, Real dt);

  size_t GetSizeInBytes();
  void ReadRestart(char *mbdata, int &os);
  void WriteRestart(char *&pdata);

protected:
  // Class methods
  static int AddIntProperty();
  static int AddRealProperty();
  static int AddAuxProperty();
  static void Migrate(Mesh *pm);

  // Class variables
  static bool initialized;  // whether or not the class is initialized
  static int nint, nreal;   // numbers of integer and real properties
  static int naux;          // number of auxiliary properties

  static int ipid;              // index for the particle ID
  static int ixp, iyp, izp;     // indices for the position components
  static int ivpx, ivpy, ivpz;  // indices for the velocity components

  // Instance methods
  void AssignShorthands();  // Needs to be called everytime
                            // intprop & realprop are resized

  // Instance variables
  AthenaArray<long> intprop;   // integer particle properties
  AthenaArray<Real> realprop;  // real particle properties
  AthenaArray<Real> auxprop;   // auxiliary particle properties (for intermediate
                               //     computations)

  AthenaArray<Real> xi1, xi2, xi3;  // position indices of each particle
                                    // in local meshblock

  AthenaArray<long> pid;            // shorthand for particle ID
  AthenaArray<Real> xp, yp, zp;     // shorthand for position components
  AthenaArray<Real> vpx, vpy, vpz;  // shorthand for velocity components

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
  void SendToNeighbors();
  void FlushReceiveBuffer();

                            // MeshBlock-to-MeshBlock communication:
  AthenaArray<long> irecv;  //   integer receive buffers
  AthenaArray<Real> rrecv;  //   real receive buffers
  int nrecvmax;             //   maximum number of particles per receive buffer
  int nrecv;                //   actual number of particles per receive buffer
};

#endif
