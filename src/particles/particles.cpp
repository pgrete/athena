//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particles.cpp
//  \brief implements functions in particle classes

#include "particles.hpp"

//--------------------------------------------------------------------------------------
//! \fn Particles::Particles(MeshBlock *pmb, ParameterInput *pin)
//  \brief constructs a Particles instance.

Particles::Particles(MeshBlock *pmb, ParameterInput *pin)
{
  // Point to the calling MeshBlock.
  pmy_block = pmb;
  npar = 1;
  nparmax = 1;

  // Allocate IDs.
  id.NewAthenaArray(nparmax);

  // Allocate positions.
  x1.NewAthenaArray(nparmax);
  x2.NewAthenaArray(nparmax);
  x3.NewAthenaArray(nparmax);

  // Allocate velocities.
  v1.NewAthenaArray(nparmax);
  v2.NewAthenaArray(nparmax);
  v3.NewAthenaArray(nparmax);
}

//--------------------------------------------------------------------------------------
//! \fn Particles::~Particles()
//  \brief destroys a Particles instance.

Particles::~Particles()
{
  // Delete IDs.
  id.DeleteAthenaArray();

  // Delete positions.
  x1.DeleteAthenaArray();
  x2.DeleteAthenaArray();
  x3.DeleteAthenaArray();

  // Delete velocities.
  v1.DeleteAthenaArray();
  v2.DeleteAthenaArray();
  v3.DeleteAthenaArray();
}
