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
