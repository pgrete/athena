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
#include "network/network.hpp" //ChemNetwork class

class MeshBlock;
class ParameterInput;
class ODEWrapper;

//! \class ChemSpecies
//  \brief Chemical species data and functions
class ChemSpecies {
public:
  //constructor: allocate memory. Need number of species in input file.
  ChemSpecies(MeshBlock *pmb, ParameterInput *pin); 
  ~ChemSpecies();

  MeshBlock *pmy_block; //ptr to a meshblock containing the chemical species

  //s(ispec, k, j, i). read in s1(i, ispec), and loop over i,
  //maybe parallelize i with openmpi later.
  AthenaArray<Real> s;  //fractional abundance of species s = n(s)/n(H)
  AthenaArray<Real> s1; //abundance of species copy at intermediate step

  ChemNetwork *pchemnet; //pointer to chemical network
  ODEWrapper *podew; //pointer to ode solver
};

//! \class ODEWrapper
//  \brief Wrapper for ODE solver, CVODE
class ODEWrapper {
public:
  //Constructor: Initialize CVODE, allocate memory for the ODE solver.
  ODEWrapper(ChemSpecies *pspec, ParameterInput *pin);
  ~ODEWrapper();
  //Update abundance in species over time dt.
  // For each cell:
  // Step 1: Set the radiation field strength in ChemNetwork.
  // Depends on the data structure of radiation field, this can be copying
  // the value from Radiation class to ChemNetwork class, or just pass a pointer.
  //
  // Step 2: re-initialize CVODE with starting time t, and starting abundance
  // y. If x(k, j, i, ispec), we can just pass a pointer to CVODE, otherwise,
  // we need to copy the abundance of species to an array.
  //
  // Step 3: Integration. Update the array of species abundance in that
  // cell over time dt.
  // 
  // Note that this will be not vectorizable(?).
  void Integrate();

  //solve the chemical abundance to equilibrium. Useful for post-processing.
  void SolveEq();

  //void SethInit(const Real h_init);
  //Get the last step size
  Real GetLastStep() const;
  //Get the next step size
  Real GetNextStep() const;
  //Get the number of steps between two reinits.
  long int GetNsteps() const;

private:
  ChemSpecies *pmy_spec_;
  Real reltol_;//relative tolerance
  Real abstol_[NSPECIES];//absolute tolerance
  void *cvode_mem_;
  N_Vector y_;
  Real *ydata_;

  //CVODE checkflag
  void CheckFlag(const void *flagvalue, const char *funcname, 
                 const int opt) const;
};


#endif // SPECIES_HPP
