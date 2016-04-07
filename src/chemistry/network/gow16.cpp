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
//! \file gow16.cpp
//  \brief implementation of functions in class ChemNetwork, using the GOW16
//  network, see paper by Gong, Ostriker, Wolfire 2016 
//======================================================================================

// this class header
#include "network.hpp"

//athena++ header
#include "../species.hpp"
#include "../../parameter_input.hpp"       //ParameterInput

//c++ header
#include <vector>
#include <stdexcept>  // std::runtime_error()
#include <sstream>    // stringstream
#include <iostream>   // endl

//species names
std::string ChemNetwork::species_names[NSPECIES] = 
{"He+", "OHx", "CHx", "CO", "C+", "HCO+", "H2", "H+", "H3+", "H2+", "S+", "Si+", "E"};

//below are ghost species. The aboundances of ghost species are
// recalculated in RHS everytime by other species.
static const int ngs_ = 7;
static const std::string ghost_species_names_[ngs_] = 
{"*Si", "*S", "*C", "*O", "*He", "*e", "*H"};

//parameters of the netowork
static Real zdg_, xHe_, xC_std_, xO_std_, xS_std_, xSi_std_, cr_rate_, bCO_;
//index of species
static int iCO_,  iH2_, iH2plus_, iHCOplus_, iCHx_, iCplus_, 
		       iOHx_, iHeplus_, iH3plus_, iHplus_, iSplus_, iSiplus_, iE_,
					 igS_, igSi_, igO_, igH_, ige_, igC_, igHe_; 

//find the index of element in the array of strings.
//report error if find repetitive elements
static int FindStrIndex(const std::string *str_arr, const int len,
		                    const std::string name);


ChemNetwork::ChemNetwork(ChemSpecies *pspec, ParameterInput *pin) {
	//number of species and a list of name of species
  pmy_spec_ = pspec;
	pmy_mb_ = pspec->pmy_block;
	//set the parameters from input file
	zdg_ = pin->GetOrAddReal("chemistry", "Zdg", 1.);//dust and gas metallicity
	xHe_ = pin->GetOrAddReal("chemistry", "xHe", 0.1);//He aboundance per H
	//metal abundance at Z=1
	xC_std_ = pin->GetOrAddReal("chemistry", "xC", 1.6e-4); 
	xO_std_ = pin->GetOrAddReal("chemistry", "xO", 3.2e-4);
	xS_std_ = pin->GetOrAddReal("chemistry", "xS", 3.5e-6);
	xSi_std_ = pin->GetOrAddReal("chemistry", "xSi", 1.7e-6);
	//cosmic ray ionization rate per H
	cr_rate_ = pin->GetOrAddReal("chemistry", "CR", 2e-16);
	//Veolcity dispersion of CO in km/s for calculating NCOeff in CO cooling
	bCO_ = pin->GetOrAddReal("chemistry", "bCO", 1.);

	//assign index of species
	iHeplus_ = FindStrIndex(species_names, NSPECIES, "He+");
	iOHx_ = FindStrIndex(species_names, NSPECIES, "OHx");
	iCHx_ = FindStrIndex(species_names, NSPECIES, "CHx");
	iCO_ = FindStrIndex(species_names, NSPECIES, "CO");
	iCplus_ = FindStrIndex(species_names, NSPECIES, "C+");
	iHCOplus_ = FindStrIndex(species_names, NSPECIES, "HCO+");
	iH2_ = FindStrIndex(species_names, NSPECIES, "H2");
	iHplus_ = FindStrIndex(species_names, NSPECIES, "H+");
	iH3plus_ = FindStrIndex(species_names, NSPECIES, "H3+");
	iH2plus_ = FindStrIndex(species_names, NSPECIES, "H2+");
	iSplus_ = FindStrIndex(species_names, NSPECIES, "S+");
	iSiplus_ = FindStrIndex(species_names, NSPECIES, "Si+");
	iE_ = FindStrIndex(species_names, NSPECIES, "E");
	//assign index of ghost species
	igSi_ = FindStrIndex(ghost_species_names_, ngs_, "*Si");
	igS_ = FindStrIndex(ghost_species_names_, ngs_, "*S");
	igC_ = FindStrIndex(ghost_species_names_, ngs_, "*C");
	igO_ = FindStrIndex(ghost_species_names_, ngs_, "*O");
	igHe_ = FindStrIndex(ghost_species_names_, ngs_, "*He");
	ige_ = FindStrIndex(ghost_species_names_, ngs_, "*e");
	igH_ = FindStrIndex(ghost_species_names_, ngs_, "*H");

}

ChemNetwork::~ChemNetwork() {}

int ChemNetwork::RHS(const realtype t, const N_Vector y, N_Vector ydot) {
  return 0;
}

int ChemNetwork::Jacobian(const long int N, const realtype t,
    const N_Vector y, const N_Vector fy, 
    DlsMat J, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  return 0;
}

void ChemNetwork::Initialize() {
	return;
}

static int FindStrIndex(const std::string *str_arr, const int len,
		                    const std::string name) {
	std::vector<int> ifind;
  std::stringstream msg; //error message
	for (int i=0; i<len; i++) {
		if (str_arr[i] == name) {
			ifind.push_back(i);
		}
	}
	if (ifind.size() == 1) {
		return ifind[0];
	} else if (ifind.size() == 0) {
		msg <<  "### FATAL ERROR in ChemNetwork [FindStrIndex]" << std::endl
			<< name << " not found." << std::endl; 
      throw std::runtime_error(msg.str().c_str());
	} else {
		msg <<  "### FATAL ERROR in ChemNetwork [FindStrIndex]" << std::endl
			<< name << " found more than once (" << ifind.size() << ")."  << std::endl; 
      throw std::runtime_error(msg.str().c_str());
	}
}
