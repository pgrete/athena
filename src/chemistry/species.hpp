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

//c++ headers
#include <string> //std::string
#include <vector>     // vector container

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

//CVODE headers. TODO: need to do this when Jim settle library with Kengo
#include <sundials/sundials_types.h> // realtype type
#include <nvector/nvector_serial.h> // N_Vector type
#include <sundials/sundials_dense.h> // definitions DlsMat DENSE_ELEM
//TODO: maybe move sundial.h here.
//#include "sundial.h" /*Ith IJth macro and CheckFlag function*/

class MeshBlock;
class ParameterInput;
class NetworkWrapper;
class ChemNetwork;
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
  
	//a list of species name
	std::vector<std::string> species_names;
  int nspec; //number of species

	//return the index of a species in the vector species_names
	int GetSpeciesIndex(std::string name);
private:
};

//! \class NetworkWrapper
//  \brief Wrapper of the chemical network class.
class NetworkWrapper {
public:
  NetworkWrapper();
  virtual ~NetworkWrapper();
  static int WrapJacobian(const long int N, const realtype t,
                          const N_Vector y, const N_Vector fy, 
                          DlsMat J, void *user_data,
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int WrapRHS(const realtype t, const N_Vector y,
                     N_Vector ydot, void *user_data);

  //------------All functions below has to be overloaded------------
  // Note that the RHS and Jac does NOT have user_data. All parameters should
  // be passed to the class as private variables.
  // right hand side of ode
  virtual int RHS(const realtype t, const N_Vector y,
                  N_Vector ydot) = 0;
  //Jacobian
  virtual int Jacobian(const long int N, const realtype t,
                       const N_Vector y, const N_Vector fy, 
                       DlsMat J, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) = 0;
};

//! \class ChemNetwork
//  \brief Chemical Network that defines the reaction rates between species.
class ChemNetwork : public NetworkWrapper {
public:
  ChemNetwork(ChemSpecies *pspec, ParameterInput *pin);
  ~ChemNetwork();
  //RHS: right-hand-side of ODE. dy/dt = ydot(t, y). Here y are the abundance
  //of species.
  //realtype is float/double, defined in CVODE header file.
  //N_Vector is a struct that used to represent vectors used in CVODE.
  //N_Vector contains information about the vector (e.g. dimension), and a
  //pointer to the array data. 
  //details see CVODE package documentation.
  int RHS(const realtype t, const N_Vector y, N_Vector ydot);
  //Jacobian: Jacobian of ODE. CVODE can also numerically calculate Jacobian if
  //this is not specified. Details see CVODE package documentation.
  int Jacobian(const long int N, const realtype t,
               const N_Vector y, const N_Vector fy, 
               DlsMat J, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

  //Set parameters in the chemistry network. This should be done at
  //construction, from the ParameterInput
  //For example, set total carbon abundance
  //void SetxCtot(const double xCtot);

  //Set radiation field strength.
  //TODO: need to think about the way to do this.
  //void SetRadField(/* Multi-frequency radiation field strength for
  //                    photo-ionization of different species*/);
  //... other parameters, such as metallicity.....

private:
  ChemSpecies *pmy_spec_;
};

//! \class ODEWrapper
//  \brief Wrapper for ODE solver, CVODE
class ODEWrapper {
public:
  //Constructor: Initialize CVODE, allocate memory for the ODE solver.
  ODEWrapper(ChemSpecies *pspec, ParameterInput *pin);
  ~ODEWrapper();
  //Setting the parameter for the solver.
  //void SetRelTol(const double reltol); This would be done from input block,
  //such as <chemistry> relative tolerance
  //... other parameters, such as absolute tolerance, maximum step...

  //Update abundance in species over time dt.
  //For post-processing, can design a function to solve to equilibrium. 
  //For time-dependent chemistry, dt = hydro timestep ?
  //
  /* The integration will look like:
   *
   * For each cell:
   * Step 1: Set the radiation field strength in ChemNetwork.
   * Depends on the data structure of radiation field, this can be copying
   * the value from Radiation class to ChemNetwork class, or just pass a pointer.
   *
   * Step 2: re-initialize CVODE with starting time t, and starting abundance
   * y. If x(k, j, i, ispec), we can just pass a pointer to CVODE, otherwise,
   * we need to copy the abundance of species to an array.
   *
   * Step 3: Integration. Update the array of species abundance in that
   * cell over time dt.
   * 
   * Note that this will be not vectorizable(?).
   */
  void Integrate(const double dt);

  //solve the chemical abundance to equilibrium. Useful for post-processing.
  void SolveEq();

private:
  ChemSpecies *pmy_spec_;
};


#endif // SPECIES_HPP
