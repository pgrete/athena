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

//c++ headers
#include <algorithm> //std:find()
#include <stdexcept>  // std::runtime_error()
#include <sstream>    // stringstream
#include <iostream> //std::cout

// Athena++ headers
#include "../athena.hpp"                // array access, macros, Real
#include "../athena_arrays.hpp"         // AthenaArray
#include "../mesh.hpp"                  // MeshBlock, Mesh
#include "../parameter_input.hpp"       //ParameterInput
#include "network/network.hpp"          //ChemNetwork

// this class header
#include "species.hpp"

//TODO: MPI & OMP header?

//function to split a string into a vector, be sure to have no repeat
static std::vector<std::string> SplitNoRep(std::string str, char delimiter);
//function to get rid of white space leading/trailing a string
static void Trim(std::string &s);

// constructor, initializes data structures and parameters
ChemSpecies::ChemSpecies(MeshBlock *pmb, ParameterInput *pin) {
  pmy_block = pmb;

  //read the number of species from input file
	std::string str_species = pin->GetString("chemistry", "species");
	//put different species into a vector container
	species_names = SplitNoRep(str_species, ',');
	//number of species
	nspec = species_names.size();

  //allocate memory for fractional abundance of species
  int ncells1 = pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmy_block->block_size.nx2 > 1) ncells2 = pmy_block->block_size.nx2 + 2*(NGHOST);
  if (pmy_block->block_size.nx3 > 1) ncells3 = pmy_block->block_size.nx3 + 2*(NGHOST);

  s.NewAthenaArray(nspec, ncells3, ncells2, ncells1);

  //allocate memory for the copy of s at intermediate step
  s1.NewAthenaArray(ncells1, nspec);

  //construct ptrs to objects related to solving chemistry source term.
  pchemnet = new ChemNetwork(this, pin);
  podew = new ODEWrapper(this, pin);
}

ChemSpecies::~ChemSpecies() {
  s.DeleteAthenaArray();
  s1.DeleteAthenaArray();
  delete pchemnet;
  delete podew;
}

int ChemSpecies::GetSpeciesIndex(std::string name) {
  std::stringstream msg; //error message
	int pos = std::find( species_names.begin(), species_names.end(), name )
		          - species_names.begin();
	if (pos < nspec) {
		return pos;
	} else {
			msg << "### FATAL ERROR in ChemSpecies::GetSpeciesIndex()" << std::endl
				<< "Sepcies " << name << " not known" << std::endl;
			throw std::runtime_error(msg.str().c_str());
	}
}

//======================================================================================
//! \fn static std::vector<std::string> SplitNoRep(std::string str, char delimiter)
//  \brief split a string, and store sub strings in a vector, make sure no
//  repetitive sub-strings.
//======================================================================================
static std::vector<std::string> SplitNoRep(std::string str, char delimiter) {
  std::vector<std::string> internal;
  std::stringstream ss(str); // Turn the string into a stream.
  std::string tok;
  std::stringstream msg; //error message
  
  while(getline(ss, tok, delimiter)) {
    Trim(tok);
		//check of the element still in vector
		if( std::find(internal.begin(), internal.end(), tok) == internal.end() ) {
			internal.push_back(tok);
		} else {
			msg << "### FATAL ERROR in ChemSpecies::ChemSpecies() [SplitNoRep]" << std::endl
				<< "Repetitive species: " << tok << std::endl;
			throw std::runtime_error(msg.str().c_str());
		}
  }
  
  return internal;
}

//======================================================================================
//! \fn static void Trim(std::string &s)
//  \brief get rid of white spaces leading and trailing a string
//======================================================================================
static void Trim(std::string &s)
{
  size_t p = s.find_first_not_of(" \t\n");
  s.erase(0, p);

  p = s.find_last_not_of(" \t\n");
  if (p != std::string::npos) {
    s.erase(p+1);
  }
}
