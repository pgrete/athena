#ifndef NETWORK_HPP
#define NETWORK_HPP
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
#include <stdio.h> //FILE, fprintf()

// Athena++ classes headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"

//CVODE headers. 
#include <sundials/sundials_types.h> // realtype type
#include <nvector/nvector_serial.h> // N_Vector type
#include <sundials/sundials_dense.h> // definitions DlsMat DENSE_ELEM

class ChemSpecies;
class ParameterInput;
class NetworkWrapper;

//! \class NetworkWrapper
//  \brief Wrapper of the chemical network class.
class NetworkWrapper {
public:
  NetworkWrapper();
  virtual ~NetworkWrapper();
  static int WrapJacobian(const long int n, const realtype t,
                          const N_Vector y, const N_Vector fy, 
                          DlsMat jac, void *user_data,
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
  static int WrapRHS(const realtype t, const N_Vector y,
                     N_Vector ydot, void *user_data);

  //------------All functions below has to be overloaded------------
  // Note that the RHS and Jac does NOT have user_data. All parameters should
  // be passed to the class as private variables.
  // right hand side of ode
  virtual void RHS(const Real t, const Real y[NSPECIES], Real ydot[NSPECIES]) = 0;
  //Jacobian
  virtual void Jacobian(const Real t,
               const Real y[NSPECIES], const Real fy[NSPECIES], 
               Real jac[NSPECIES][NSPECIES], Real tmp1[NSPECIES],
               Real tmp2[NSPECIES], Real tmp3[NSPECIES]) = 0;
};

//! \class ChemNetwork
//  \brief Chemical Network that defines the reaction rates between species.
class ChemNetwork : public NetworkWrapper {
public:
  ChemNetwork(ChemSpecies *pspec, ParameterInput *pin);
  ~ChemNetwork();

	//a list of species name, used in output
	static std::string species_names[NSPECIES];

	//Set the rates of chemical reactions, eg. through density and radiation field.
  //k, j, i are the corresponding index of the grid
  void InitializeNextStep(const int k, const int j, const int i);
  //output properties of network. Can be used in eg. ProblemGenerator.
  void OutputProperties(FILE *pf) const;

  //RHS: right-hand-side of ODE. dy/dt = ydot(t, y). Here y are the abundance
  //of species.
  //realtype is float/double, defined in CVODE header file.
  //N_Vector is a struct that used to represent vectors used in CVODE.
  //N_Vector contains information about the vector (e.g. dimension), and a
  //pointer to the array data. 
  //details see CVODE package documentation.
  void RHS(const Real t, const Real y[NSPECIES], Real ydot[NSPECIES]);
  //Jacobian: Jacobian of ODE. CVODE can also numerically calculate Jacobian if
  //this is not specified. Details see CVODE package documentation.
  void Jacobian(const Real t,
               const Real y[NSPECIES], const Real fy[NSPECIES], 
               Real jac[NSPECIES][NSPECIES],
               Real tmp1[NSPECIES], Real tmp2[NSPECIES], Real tmp3[NSPECIES]);

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
	MeshBlock *pmy_mb_;

};

#endif // NETWORK_HPP
