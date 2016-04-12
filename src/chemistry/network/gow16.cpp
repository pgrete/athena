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
#include <math.h> //a^x = pow(a,x)

//species names
std::string ChemNetwork::species_names[NSPECIES] = 
{"He+", "OHx", "CHx", "CO", "C+", "HCO+", "H2", "H+", "H3+", "H2+", "S+", "Si+", "E"};

//below are ghost species. The aboundances of ghost species are
// recalculated in RHS everytime by other species.
static const int ngs_ = 7;
static const std::string ghost_species_names_[ngs_] = 
{"*Si", "*S", "*C", "*O", "*He", "*e", "*H"};
static std::string species_names_all_[NSPECIES+ngs_];//all species

//find the index of element in the array of strings.
//report error if find repetitive elements
static int FindStrIndex(const std::string *str_arr, const int len,
		                    const std::string name);

//parameters of the netowork
static Real zdg_, xHe_, xC_std_, xO_std_, xS_std_, xSi_std_, cr_rate_, bCO_;
//index of species
static int iHeplus_ = FindStrIndex(ChemNetwork::species_names, NSPECIES, "He+");
static int iOHx_ = FindStrIndex(ChemNetwork::species_names, NSPECIES, "OHx");
static int iCHx_ = FindStrIndex(ChemNetwork::species_names, NSPECIES, "CHx");
static int iCO_ = FindStrIndex(ChemNetwork::species_names, NSPECIES, "CO");
static int iCplus_ = FindStrIndex(ChemNetwork::species_names, NSPECIES, "C+");
static int iHCOplus_ = FindStrIndex(ChemNetwork::species_names, NSPECIES, "HCO+");
static int iH2_ = FindStrIndex(ChemNetwork::species_names, NSPECIES, "H2");
static int iHplus_ = FindStrIndex(ChemNetwork::species_names, NSPECIES, "H+");
static int iH3plus_ = FindStrIndex(ChemNetwork::species_names, NSPECIES, "H3+");
static int iH2plus_ = FindStrIndex(ChemNetwork::species_names, NSPECIES, "H2+");
static int iSplus_ = FindStrIndex(ChemNetwork::species_names, NSPECIES, "S+");
static int iSiplus_ = FindStrIndex(ChemNetwork::species_names, NSPECIES, "Si+");
static int iE_ = FindStrIndex(ChemNetwork::species_names, NSPECIES, "E");
//index of ghost species
static int igSi_ = FindStrIndex(ghost_species_names_, ngs_, "*Si") + NSPECIES;
static int igS_ = FindStrIndex(ghost_species_names_, ngs_, "*S") + NSPECIES;
static int igC_ = FindStrIndex(ghost_species_names_, ngs_, "*C") + NSPECIES;
static int igO_ = FindStrIndex(ghost_species_names_, ngs_, "*O") + NSPECIES;
static int igHe_ = FindStrIndex(ghost_species_names_, ngs_, "*He") + NSPECIES;
static int ige_ = FindStrIndex(ghost_species_names_, ngs_, "*e") + NSPECIES;
static int igH_ = FindStrIndex(ghost_species_names_, ngs_, "*H") + NSPECIES;
//physical constants
static const Real temp_coll_ = 7.0e2;
static const Real mH_ = 1.67e-24;
static const Real mCO_ = 4.68e-23;


//-------------------chemical network---------------------
//number of different reactions
static const int n_cr_ = 8;
static const int n_2body_ = 29;
static const int n_ph_ = 8;
static const int n_gr_ = 6;
static const int nE_ = 15;//number of heating and cooling processes
static Real kcr_[n_cr_]; //rates for cosmic-ray reactions.
static Real k2body_[n_2body_]; //rates for 2 body reacrtions.
static Real kph_[n_ph_]; //rates for photo- reactions.
static Real kgr_[n_gr_]; //rates for grain assisted reactions.

//cosmic ray chemistry network
// (0) cr + H2 -> H2+ + *e
// (1) cr + *He -> He+ + *e 
// (2) cr + *H  -> H+ + *e 
// -----added as Clark + Glover 2015---- 
// (3) cr + *C -> C+ + *e     --including direct and cr induce photo reactions 
// (4) crphoto + CO -> *O + *C       
// (5) cr + CO -> HCO+ + e  --schematic for cr + CO -> CO+ + e
// -----S, CR induced photo ionization, experimenting----
// (6) cr + S -> S+ + e, simply use 2 times rate of C, as in UMIST12
// -----Si, CR induced photo ionization, experimenting----
// (7) cr + Si -> Si+ + e, UMIST12
static const int incr_[n_cr_] = {iH2_, igHe_, igH_, 
                                 igC_, iCO_, iCO_,
                                 igS_, igSi_};
static const int outcr_[n_cr_] = {iH2plus_, iHeplus_, iHplus_, 
                                  iCplus_, igO_, iHCOplus_,
                                  iSplus_, iSiplus_};
static const Real kcr_base_[n_cr_] = 
{2.0, 1.1, 1.0, 
1020., 10., 6.52,
2040., 4200.}; 

/*2 body reactions*/
/*NOTE: photons from recombination are ignored*/
/* Reactions are, in order.
 * -- are equations of special rate treatment in Glover, Federrath+ 2010:
 (0) H3+ + *C -> CH + H2
 (1) H3+ + *O -> OH + H2        
 (2) H3+ + CO -> HCO+ + H2
 (3) He+ + H2 -> H+ + *He + *H   --(89) exp(-35/T)
 (4) He+ + CO -> C+ + *O + *He   
 (5) C+ + H2 -> CH + *H         -- schematic reaction for C+ + H2 -> CH2+
 (6) C+ + OH -> HCO+             -- Schematic equation for C+ + OH -> CO+ + H.
 Use rates in KIDA website.
 (7) CH + *O -> CO + *H
 (8) OH + *C -> CO + *H          --exp(0.108/T)
 (9) He+ + *e -> *He             --(17) Case B
 (10) H3+ + *e -> H2 + *H
 (11) C+ + *e -> *C              -- Include RR and DR, Badnell2003, 2006.
 (12) HCO+ + *e -> CO + *H
 ----added in GO2012--------
 (13) H2+ + H2 -> H3+ + *H       --(54) exp(-T/46600)
 (14) H+ + *e -> *H              --(12) Case B
 ---collisional dissociation, only important at high temperature T>1e3---
 (15) H2 + *H -> 3 *H            --(9) Density dependent. See Glover+MacLow2007
 (16) H2 + H2 -> H2 + 2 *H       --(10) Density dependent. See Glover+MacLow2007
 (17) *H + *e -> H+ + 2 *e       --(11) Relates to Te
 ----added for H3+ destruction in addition to (10)----
 (18) H3+ + *e -> *3H            --(111)
 ----added He+ destruction in addtion to (3), from UMIST12----
 (19) He+ + H2 -> H2+ + *He
 ----added CH reaction to match for abundances of CH---
 (20) CH + *H -> H2 + *C         
 ----added to match the Meudon code ---
 (21) OH + *O -> *O + *O + *H
 ---branching of C+ + H2 ------
 (22) C+ + H2 + *e -> *C + *H + *H
 ---S , rate from UMIST12---
 (23) S+ + *e -> *S
 (24) C+ + *S -> S+ + *C
 ---Si , rate from UMIST12---
 (25) Si+ + *e -> *Si
 (26) C+ + *Si -> Si+ + *C
 --- H2O+ + e reaction ---
 (27) H3+ + *O + *e -> H2 + *O + *H
 --- OH destruction with He+
 (28) He+ + OH -> *H + *He + *O(O+)
 */
static const int in2body1_[n_2body_] = 
          {iH3plus_, iH3plus_, iH3plus_, iHeplus_, iHeplus_,    
           iCplus_, iCplus_, iCHx_, iOHx_, iHeplus_,
           iH3plus_, iCplus_, iHCOplus_,
					 iH2plus_, iHplus_, iH2_, iH2_, igH_,
					 iH3plus_, iHeplus_, 
           iCHx_, iOHx_, iCplus_, iSplus_, iCplus_,
           iSiplus_, iCplus_, iH3plus_, iHeplus_};
static const int in2body2_[n_2body_] = 
          {igC_, igO_, iCO_, iH2_, iCO_,   
           iH2_, iOHx_, igO_, igC_, ige_,   
           ige_, ige_, ige_,
					 iH2_, ige_, igH_, iH2_, ige_,
					 ige_, iH2_, 
           igH_, igO_, iH2_, ige_, igS_,
           ige_, igSi_, igO_, iOHx_};
/*Note: output to ghost species doesn't matter. The abundances of ghost species
 * are updated using the other species at every timestep*/
static const int out2body1_[n_2body_] = 
          {iCHx_, iOHx_, iHCOplus_, iHplus_, iCplus_,   
           iCHx_, iHCOplus_, iCO_, iCO_, igHe_,   
           iH2_, igC_, iCO_,
					 iH3plus_, igH_, igH_, iH2_, iHplus_,
					 igH_, iH2plus_, 
           iH2_, igO_, igC_, igS_, iSplus_,
           igSi_, iSiplus_, iH2_, igH_};
static const int out2body2_[n_2body_] = 
          {iH2_, iH2_, iH2_, igHe_, igO_,   
           igH_, igH_, igH_, igH_, igH_,   
           igH_, igH_, igH_,
					 igH_, igH_, igH_, igH_, ige_,
					 igH_, igHe_, 
           igC_, igH_, igH_, igH_, igC_,
           igH_, igC_, igO_, igHe_};
static const Real k2Texp_[n_2body_] = 
 {0.0, -0.190, 0.0, 0.0, 0.0, 
  -1.3, 0.0, 0.0, -0.339, -0.5, 
  -0.52, 0.0, -0.64,
  0.042, 0.0, 0.0, 0.0, 0.0,
  -0.52, 0.0,
  0.26, 0.0, -1.3, -0.59, 0.0,
  -0.62, 0.0, -0.190, 0.0};
static const Real k2body_base_[n_2body_] = 
                {2.0e-9, 1.99e-9, 1.7e-9, 3.7e-14, 1.6e-9, 
                 3.3e-13 * 0.7, 1.00, 7.0e-11, 7.95e-10, 1.0e-11, 
                 4.54e-7, 1.00, 1.15e-5,
								 2.84e-9, 2.753e-14, 1.00, 1.00, 1.00,
								 8.46e-7, 7.20e-15, 
                 2.81e-11, 3.5e-11, 3.3e-13 * 0.3, 1.6e-10, 5e-11,
                 1.46e-10, 2.1e-9, 1.99e-9, 1.00};

/* photo reactions.
 * Reaction rates in Drain 1978 field units.
 * Reactions are, in order:
 (0) h nu + *C -> C+ + *e
 (1) h nu + CH -> *C + *H
 (2) h nu + CO -> *C + *O            --self-shielding and shielding by H2
 (3) h nu + OH -> *O + *H
 (4) h nu + HCO+ -> CO + *H
 ----added in GO2012--------
 (5) h nu + H2 -> *H + *H            --self- and dust shielding
 ----S, from UMIST12
 (6) h nu + *S -> S+
 ----Si, from UMIST12
 (7) h nu + *Si -> Si+
 */
static const int inph_[n_ph_] = {
              igC_, iCHx_, iCO_,
              iOHx_, iHCOplus_, iH2_, igS_, igSi_};
static const int outph1_[n_ph_] = {
              iCplus_, igC_, igC_,
              igO_, iCO_, igH_, iSplus_, iSiplus_};
static const Real kph_base_[n_ph_] = {3.1e-10, 9.2e-10, 2.6e-10/*Visser2009*/,   
																			  3.9e-10, 5.40e-12, 5.6e-11, 
                                        6e-10, 3.1e-9}; 
static const Real kph_avfac_[n_ph_] = {3.33, 1.72, 3.53/*Visser2009*/,  
	                                       2.24, 3.32, 3.74/*Draine+Bertoldi1996*/,
                                         3.10, 2.3};

/* Grain assisted recombination of H, H2, C+ and H+
 (0) *H + *H + gr -> H2 + gr
 (1) H+ + *e + gr -> *H + gr
 (2) C+ + *e + gr -> *C + gr
 (3) He+ + *e + gr -> *He + gr
 ------S, from WD2001-----
 (4) S+ + *e + gr -> *S + gr
 ------Si, from WD2001-----
 (5) Si+ + *e + gr -> *Si + gr
 */
static const int ingr_[n_gr_] = {igH_, iHplus_, iCplus_, iHeplus_, 
                                 iSplus_, iSiplus_};
static const int outgr_[n_gr_] = {iH2_, igH_, igC_, igHe_, 
                                 igS_, igSi_};
static const Real cHp_[7] = {12.25, 8.074e-6, 1.378, 5.087e2, 1.586e-2,
															 0.4723, 1.102e-5}; 
static const Real cCp_[7] = {45.58, 6.089e-3, 1.128, 4.331e2, 4.845e-2,
                               0.8120, 1.333e-4};
static const Real cHep_[7] = {5.572, 3.185e-7, 1.512, 5.115e3, 3.903e-7,
                                0.4956, 5.494e-7};
static const Real cSp_[7] = {3.064, 7.769e-5, 1.319, 1.087e2, 3.475e-1,
                               0.4790, 4.689e-2};
static const Real cSip_[7] = {2.166, 5.678e-8, 1.874, 4.375e4, 1.635e-6,
                               0.8964, 7.538e-5};
//-----------------end of chemical network---------------------


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

  //initialize rates to zero
  for (int i=0; i<n_cr_; i++) {
    kcr_[i] = 0;
  }
  for (int i=0; i<n_2body_; i++) {
    k2body_[i] = 0;
  }
  for (int i=0; i<n_ph_; i++) {
    kph_[i] = 0;
  }
  for (int i=0; i<n_gr_; i++) {
    kgr_[i] = 0;
  }
  //copy species to a full list of species names
  for (int i=0; i<NSPECIES; i++) {
    species_names_all_[i] = species_names[i];
  }
  for (int i=NSPECIES; i<NSPECIES+ngs_; i++) {
    species_names_all_[i] = ghost_species_names_[i-NSPECIES];
  }

}

ChemNetwork::~ChemNetwork() {}

void ChemNetwork::RHS(const Real t, const Real y[NSPECIES], Real ydot[NSPECIES]) {
  for (int i=0; i<NSPECIES; i++) {
    ydot[i] = i;
  }
  return;
}

void ChemNetwork::Jacobian(const long int N, const Real t,
               const Real y[NSPECIES], const Real fy[NSPECIES], 
               Real J[NSPECIES][NSPECIES],
               Real tmp1[NSPECIES], Real tmp2[NSPECIES], Real tmp3[NSPECIES]) {
  for (int i=0; i<NSPECIES; i++) {
    for (int j=0; j<NSPECIES; j++) {
      J[i][j] = 0;
    }
  }
  return;
}

void ChemNetwork::Initialize() {
	return;
}

void ChemNetwork::OutputProperties(FILE *pf) const {
  //output the reactions and base rates
	for (int i=0; i<n_cr_; i++) {
		fprintf(pf, "cr  + %4s -> %4s,     kcr = %.2e ksi s-1 H-1\n", 
		 species_names_all_[incr_[i]].c_str(), species_names_all_[outcr_[i]].c_str(), kcr_base_[i]);
	}
	for (int i=0; i<n_2body_; i++) {
		fprintf(pf, "%4s  + %4s -> %4s  + %4s,     k2body = %.1e T^%.2f cm3 s-1 H-1\n", 
		 species_names_all_[in2body1_[i]].c_str(), species_names_all_[in2body2_[i]].c_str(),
		 species_names_all_[out2body1_[i]].c_str(), species_names_all_[out2body2_[i]].c_str(),
		 k2body_base_[i], k2Texp_[i]);
	}
	for (int i=0; i<n_ph_; i++) {
		fprintf(pf, "h nu  + %4s -> %4s,     kph = %.1e G0 exp(-%.1f Av) s-1 H-1\n", 
		 species_names_all_[inph_[i]].c_str(), species_names_all_[outph1_[i]].c_str(),
		 kph_base_[i], kph_avfac_[i]);
	}
	for (int i=0; i<n_gr_; i++) {
		fprintf(pf, "gr  + %4s -> %4s,     kgr = %.1e s-1 H-1\n", 
		 species_names_all_[ingr_[i]].c_str(), species_names_all_[outgr_[i]].c_str(),
		 kgr_[i]);
	}
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
