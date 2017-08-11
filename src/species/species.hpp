#ifndef SPECIES_HPP
#define SPECIES_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file species.hpp
//  \brief definitions for chemical species, network, and ode solver classes.
//======================================================================================

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#ifdef INCLUDE_CHEMISTRY
#include "../chemistry/network/network.hpp" 
#include "../chemistry/ode_wrapper.hpp"
#include  CHEMNETWORK_HEADER //ChemNetwork class
#endif //INCLUDE_CHEMISTRY

class MeshBlock;
class ParameterInput;

//! \class Species
//  \brief Species data and functions
class Species {
  friend class ODEWrapper;
public:
  //constructor: allocate memory. 
  Species(MeshBlock *pmb, ParameterInput *pin); 
  ~Species();

  MeshBlock *pmy_block; //ptr to a meshblock containing the chemical species

  AthenaArray<Real> s;  //fractional abundance of species s = n(s)/n(H)

#ifdef INCLUDE_CHEMISTRY
  //chemistry source term
  //s(ispec, k, j, i). read in s1(i, ispec), and loop over i,
  //maybe parallelize i with openmpi later.
  AthenaArray<Real> s1; //abundance of species copy at intermediate step
  ChemNetwork *pchemnet; //pointer to chemical network
  ODEWrapper *podew; //pointer to ode solver
#endif //INCLUDE_CHEMISTRY

  //Advection term
  void AddAdvectionTerm();
#ifdef INCLUDE_CHEMISTRY
private:
  AthenaArray<Real> h; //next stepsize in chemistry solver
#endif //INCLUDE_CHEMISTRY
};

#endif // SPECIES_HPP
