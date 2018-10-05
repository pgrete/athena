#ifndef PARTICLES_PARTICLES_HPP_
#define PARTICLES_PARTICLES_HPP_
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
#include "particle_buffer.hpp"
#include "particle-mesh.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

// Forward declarations
class ParameterInput;

//--------------------------------------------------------------------------------------
//! \struct Neighbor
//  \brief defines a structure for links to neighbors

struct Neighbor {
  NeighborBlock *pnb;
  MeshBlock *pmb;
  Neighbor *next;

  Neighbor() : pnb(NULL), pmb(NULL), next(NULL) {}
};

//--------------------------------------------------------------------------------------
//! \class Particles
//  \brief defines the bass class for all implementations of particles.

class Particles {

friend class MeshBlock;  // Make writing initial conditions possible.
friend class ParticleMesh;

public:
  // Class methods
  static void Initialize(ParameterInput *pin);
  static void FormattedTableOutput(Mesh *pm, OutputParameters op); 

  // Constructor
  Particles(MeshBlock *pmb, ParameterInput *pin);

  // Destructor
  ~Particles();

  // Instance methods
  void ClearBoundary();
  void Integrate(int step);
  void LinkNeighbors();
  void NewBlockTimeStep();
  void SendParticleMesh();
  void SendToNeighbors();
  void SetPositionIndices();
  void StartReceiving();
  bool ReceiveFromNeighbors();
  bool ReceiveParticleMesh(int step);

  size_t GetSizeInBytes();
  void ReadRestart(char *mbdata, int &os);
  void WriteRestart(char *&pdata);

protected:
  // Class methods
  static int AddIntProperty();
  static int AddRealProperty();
  static int AddAuxProperty();
  static int AddWorkingArray();

  // Class variables
  static bool initialized;  // whether or not the class is initialized
  static int nint, nreal;   // numbers of integer and real particle properties
  static int naux;          // number of auxiliary particle properties
  static int nwork;         // number of working arrays for particles

  static int ipid;                 // index for the particle ID
  static int ixp, iyp, izp;        // indices for the position components
  static int ivpx, ivpy, ivpz;     // indices for the velocity components

  static int ixp0, iyp0, izp0;     // indices for beginning position components
  static int ivpx0, ivpy0, ivpz0;  // indices for beginning velocity components

  static int ixi1, ixi2, ixi3;     // indices for position indices
  static int iapx, iapy, iapz;     // indices for acceleration components

  // Instance methods
  virtual void AssignShorthands();  // Needs to be called everytime
                                    // intprop, realprop, & auxprop are resized
                                    // Be sure to call back when derived.
  virtual void AddAcceleration(Real t, Real dt, const AthenaArray<Real>& meshsrc) {}
  virtual void ReactToMeshAux(Real t, Real dt, const AthenaArray<Real>& meshsrc) {}
  virtual void DepositToMesh(Real t, Real dt, const AthenaArray<Real>& meshsrc,
                             AthenaArray<Real>& meshdst) {}

  // Instance variables
  int npar;     // number of particles
  int nparmax;  // maximum number of particles per meshblock

                               // Data attached to the particles:
  AthenaArray<int> intprop;   //   integer properties
  AthenaArray<Real> realprop;  //   real properties
  AthenaArray<Real> auxprop;   //   auxiliary properties (communicated when
                               //     particles moving to another meshblock)
  AthenaArray<Real> work;      //   working arrays (not communicated)

  ParticleMesh *ppm;  // ptr to particle-mesh

                                       // Shorthands:
  AthenaArray<int> pid;               //   particle ID
  AthenaArray<Real> xp, yp, zp;        //   position
  AthenaArray<Real> vpx, vpy, vpz;     //   velocity
  AthenaArray<Real> xi1, xi2, xi3;     //   position indices in local meshblock
  AthenaArray<Real> apx, apy, apz;     //   acceleration
  AthenaArray<Real> xp0, yp0, zp0;     //   beginning position
  AthenaArray<Real> vpx0, vpy0, vpz0;  //   beginning velocity

  MeshBlock* pmy_block;  // MeshBlock pointer
  Mesh* pmy_mesh;        // Mesh pointer

private:
  // Class methods
  static void GetPositionIndices(MeshBlock *pmb, int npar,
                                 const AthenaArray<Real>& xp,
                                 const AthenaArray<Real>& yp,
                                 const AthenaArray<Real>& zp,
                                 AthenaArray<Real>& xi1,
                                 AthenaArray<Real>& xi2,
                                 AthenaArray<Real>& xi3);

  // Instance methods
  void ApplyBoundaryConditions(int k, Real &x1, Real &x2, Real &x3);
  void EulerStep(Real t, Real dt, const AthenaArray<Real>& meshsrc);
  void FlushReceiveBuffer(ParticleBuffer& recv);
  void SaveStatus();
  void ZeroAcceleration();
  struct Neighbor* FindTargetNeighbor(
      int ox1, int ox2, int ox3, int xi1, int xi2, int xi3);

  // Instance variables
  bool active1_, active2_, active3_;  // active dimensions

  // MeshBlock-to-MeshBlock communication:
  BoundaryValues *pbval_;            // ptr to my BoundaryValues
  Neighbor neighbor_[3][3][3];       // links to neighbors
  ParticleBuffer recv_[56];          // particle receive buffers
  enum BoundaryStatus bstatus_[56];  // boundary status
#ifdef MPI_PARALLEL
  static MPI_Comm my_comm;   // my MPI communicator
  ParticleBuffer send_[56];  // particle send buffers
#endif
};

//--------------------------------------------------------------------------------------
//! \class DustParticles
//  \brief defines the class for dust particles that interact with the gas via drag
//         force.

class DustParticles : public Particles {

friend class MeshBlock;

public:
  // Class methods
  static void Initialize(ParameterInput *pin);

  // Constructor
  DustParticles(MeshBlock *pmb, ParameterInput *pin);

  // Destructor
  ~DustParticles();

protected:
  static bool backreaction;  // on/off of back reaction
  static Real mass;          // mass of each particle
  static Real taus;          // stopping time (in code units)

private:
  // Class variables
  static bool initialized;      // whether or not the class is initialized
  static int iwx, iwy, iwz;     // indices for working arrays
  static int idpx, idpy, idpz;  // indices for momentum change

  // Instance methods.
  void AssignShorthands();
  void AddAcceleration(Real t, Real dt, const AthenaArray<Real>& meshsrc);
  void DepositToMesh(Real t, Real dt, const AthenaArray<Real>& meshsrc,
                     AthenaArray<Real>& meshdst);

  // Instance variables
  AthenaArray<Real> wx, wy, wz;     // shorthand for working arrays
  AthenaArray<Real> dpx, dpy, dpz;  // shorthand for momentum change
};

#endif  // PARTICLES_PARTICLES_HPP_
