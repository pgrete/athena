//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
//======================================================================================
//! \file particle_buffer.cpp
//  \brief implements ParticleBuffer class for communication of particles.
//======================================================================================

// C++ standard libraries
#include <cstring>
#include <sstream>

// Athena++ headers
#include "particle_buffer.hpp"

// Class variable initialization
int ParticleBuffer::nint = 0, ParticleBuffer::nreal = 0;

//--------------------------------------------------------------------------------------
//! \fn void ParticleBuffer::SetNumberOfProperties(int nint0, int nreal0)
//  \brief sets the number of integer and real propertes.

void ParticleBuffer::SetNumberOfProperties(int nint0, int nreal0)
{
  // Sanity check
  if (nint0 < 0 || nreal0 < 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [ParticleBuffer::SetNumberOfProperties]"
        << std::endl
        << "Invalid nint0 = " << nint0 << " or nreal0 = " << nreal0 << std::endl;
    throw std::runtime_error(msg.str().data());
    return;
  }

  // Set the numbers.
  nint = nint0;
  nreal = nreal0;
}

//--------------------------------------------------------------------------------------
//! \fn ParticleBuffer::ParticleBuffer()
//  \brief initiates a default instance of ParticleBuffer.

ParticleBuffer::ParticleBuffer()
{
  ibuf = NULL;
  rbuf = NULL;
  nparmax = npar = 0;
}

//--------------------------------------------------------------------------------------
//! \fn ParticleBuffer::ParticleBuffer(int nparmax0)
//  \brief initiates a new instance of ParticleBuffer with nparmax = nparmax0.

ParticleBuffer::ParticleBuffer(int nparmax0)
{
  // Sanity check
  if (nparmax0 <= 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [ParticleBuffer::ParticleBuffer]" << std::endl
        << "Invalid nparmax0 = " << nparmax0 << std::endl;
    throw std::runtime_error(msg.str().data());

    ibuf = NULL;
    rbuf = NULL;
    nparmax = npar = 0;
    return;
  }

  // Initialize the instance variables.
  nparmax = nparmax0;
  ibuf = new long [nint * nparmax];
  rbuf = new Real [nreal * nparmax];
  npar = 0;
}

//--------------------------------------------------------------------------------------
//! \fn ParticleBuffer::ParticleBuffer(int nparmax0)
//  \brief destroys an instance of ParticleBuffer.

ParticleBuffer::~ParticleBuffer()
{
  if (ibuf != NULL) delete [] ibuf;
  if (rbuf != NULL) delete [] rbuf;
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleBuffer::Reallocate(int new_nparmax)
//  \brief reallocates the buffers; the old content is preserved.

void ParticleBuffer::Reallocate(int new_nparmax)
{
  // Sanity check
  if (new_nparmax <= 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [ParticleBuffer::Reallocate]" << std::endl
        << "Invalid new_nparmax = " << new_nparmax << std::endl;
    throw std::runtime_error(msg.str().data());
    return;
  }
  if (new_nparmax < npar) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [ParticleBuffer::Reallocate]" << std::endl
        << "new_nparmax = " << new_nparmax << " < npar = " << npar << std::endl;
    throw std::runtime_error(msg.str().data());
    return;
  }

  // Allocate new space.
  nparmax = new_nparmax;
  long *ibuf_new = new long [nint * nparmax];
  Real *rbuf_new = new Real [nreal * nparmax];

  // Move existing data.
  if (npar > 0) {
    std::memcpy(ibuf_new, ibuf, nint * npar * sizeof(long));
    std::memcpy(rbuf_new, rbuf, nreal * npar * sizeof(Real));
  }

  // Delete old space.
  if (ibuf != NULL) delete [] ibuf;
  if (rbuf != NULL) delete [] rbuf;
  ibuf = ibuf_new;
  rbuf = rbuf_new;
}