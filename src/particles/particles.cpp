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

//--------------------------------------------------------------------------------------
//! \fn void Particles::Drift(Real t, Real dt)
//  \brief updates the particle positions from t to t + dt given velocities at t.

void Particles::Drift(Real t, Real dt)
{
  for (long k = 0; k < npar; ++k) {
    x1(k) += dt * v1(k);
    x2(k) += dt * v2(k);
    x3(k) += dt * v3(k);
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::Kick(Real t, Real dt)
//  \brief updates the particle velocities from t to t + dt given accelerations at t.

#include <cmath>

void Particles::Kick(Real t, Real dt)
{
  Real a1 = 0.0, a2 = 0.0, a3 = 0.0;  // TODO: might need to be put inside the loop for
                                      // vectorization.

  for (long k = 0; k < npar; ++k) {
    v1(k) += dt * a1;
    v2(k) += dt * a2;
    v3(k) += dt * a3;
  }
}

//--------------------------------------------------------------------------------------
//! \fn void Particles::Update(Mesh *pm)
//  \brief updates all particle positions and velocities from t to t + dt.
//======================================================================================
#include <iostream>

void Particles::Update(Mesh *pm)
{
  MeshBlock *pmb = pm->pblock;
  Real dth = 0.5 * pm->dt;

  // Loop over MeshBlocks
  while (pmb != NULL) {
    pmb->ppar->Drift(pm->time, dth);
    pmb->ppar->Kick(pm->time, pm->dt);
    pmb->ppar->Drift(pm->time + dth, dth);
    pmb = pmb->next;
  }
}

//--------------------------------------------------------------------------------------
//! \fn Particles::FormattedTableOutput()
//  \brief outputs the particle data in tabulated format.

#include <fstream>
#include <iomanip>
#include <sstream>

void Particles::FormattedTableOutput(Mesh *pm, OutputParameters op)
{
  MeshBlock *pmb = pm->pblock;
  Particles *ppar;
  std::stringstream fname, msg;
  std::ofstream os;

  // Loop over MeshBlocks
  while (pmb != NULL) {
    ppar = pmb->ppar;

    // Create the filename.
    fname << op.file_basename
          << ".block" << pmb->gid << '.' << op.file_id
          << '.' << std::setw(5) << std::right << std::setfill('0') << op.file_number
          << '.' << "par.tab";

    // Open the file for write.
    os.open(fname.str());
    if (!os.is_open()) {
      msg << "### FATAL ERROR in function [Particles::FormattedTableOutput]"
          << std::endl << "Output file '" << fname.str() << "' could not be opened"
          << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }

    // Write the particle data in the meshblock.
    for (long k = 0; k < ppar->npar; ++k)
      os << ppar->id(k) << "  "
         << ppar->x1(k) << "  " << ppar->x2(k) << "  " << ppar->x3(k) << "  "
         << ppar->v1(k) << "  " << ppar->v2(k) << "  " << ppar->v3(k) << std::endl;

    // Close the file and get the next meshblock.
    os.close();
    pmb = pmb->next;
  }
}
