//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
//======================================================================================
//! \file particle_buffer.hpp
//  \brief defines ParticleBuffer class for communication of particles.
//======================================================================================

// Athena++ headers
#include "../athena.hpp"

class ParticleBuffer {

friend class Particles;

public:
  // Class methods
  static void SetNumberOfProperties(int nint0, int nreal0);

  // Constructors and destructor
  ParticleBuffer();
  ParticleBuffer(int nparmax0);
  ~ParticleBuffer();

  // Instance method
  void Reallocate(int new_nparmax);

protected:
  // Class variables
  static int nint;   // number of integer properties per particle
  static int nreal;  // number of real properties per particle

  // Instance variables
  long* ibuf;   // ptr to integer buffer
  Real* rbuf;   // ptr to real buffer
  int nparmax;  // maximum number of particles
  int npar;     // actual number of particles in the buffer
};
