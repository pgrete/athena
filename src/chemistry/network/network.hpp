#ifndef NETWORK_HPP
#define NETWORK_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file network.hpp
//  \brief definitions for chemical  network.
//======================================================================================

// Athena++ classes headers
#include "../../athena.hpp"

//c++ headers
#include <string>

//CVODE headers. 
#include <sundials/sundials_types.h> // realtype type
#include <nvector/nvector_serial.h> // N_Vector type
#include <sundials/sundials_dense.h> // definitions DlsMat DENSE_ELEM


class Species;
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

#endif // NETWORK_HPP
