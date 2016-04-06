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

// Athena++ classes headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"

//CVODE headers. TODO: need to do this when Jim settle library with Kengo
#include <sundials/sundials_types.h> // realtype type
#include <nvector/nvector_serial.h> // N_Vector type
#include <sundials/sundials_dense.h> // definitions DlsMat DENSE_ELEM
//TODO: maybe move sundial.h here.
//#include "sundial.h" /*Ith IJth macro and CheckFlag function*/

class ChemSpecies;
class ParameterInput;
class NetworkWrapper;

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

#endif // NETWORK_HPP
