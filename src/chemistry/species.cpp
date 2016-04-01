//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file species.cpp
//  \brief implementation of functions in class ChemSpecies 
//======================================================================================

// Athena++ headers
#include "../athena.hpp"                // array access, macros, Real
#include "../athena_arrays.hpp"         // AthenaArray
#include "../mesh.hpp"                  // MeshBlock, Mesh
#include "../parameter_input.hpp"       //ParameterInput

// this class header
#include "species.hpp"

//TODO: MPI & OMP header?

// constructor, initializes data structures and parameters

ChemSpecies::ChemSpecies(MeshBlock *pmb, ParameterInput *pin) {
  pmy_block = pmb;

  //read the number of species from input file
  nspec = pin->GetInteger("Chemistry, nspec");

  //allocate memory for fractional abundance of species
  int ncells1 = pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmy_block->block_size.nx2 > 1) ncells2 = pmy_block->block_size.nx2 + 2*(NGHOST);
  if (pmy_block->block_size.nx3 > 1) ncells3 = pmy_block->block_size.nx3 + 2*(NGHOST);

  s.NewAthenaArray(ncells3, ncells2, ncells1, nspec);

  //allocate memory for the copy of s at intermediate step
  s1.NewAthenaArray(ncells1, nspec);

  //construct ptrs to objects related to solving chemistry source term.
  pchemnet = new ChemNet(this, pin);
  podew = new ODEWrapper(this, pin);

}

ChemSpecies::~ChemSpecies() {
  s.DeleteAthenaArray();
  s1.DeleteAthenaArray();
  delete pchemnet;
  delete podew;
}
