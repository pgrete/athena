#ifndef CGK_UTILS_HPP
#define CGK_UTILS_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file cgk_utils.hpp
//  \brief prototypes of utility functions for Chang-Goo Kim's galactic disk
//  simulations.
//======================================================================================

#include "../athena.hpp"

namespace CGKUtility
{
  //calculate tempereature in Kelvin from t1=m_H P / (rho * k_B)
  Real get_temp_from_t1(const Real t1);
  //calculate temperature in Kelvin
  //pressure: in units of [dens][v]^2, can be found in phydro.w(IEN, k, j, i). 
  //dens: density in units of [dens] = mu_H m_H cm^-3,
  //in phydro.w(IDN, k, j, i) and phydro.u(IDN, k, j, i)
  //veolcity units: [v] = km/s
  Real get_temp(Real pressure, Real dens);
}
#endif // CGK_UTILS_HPP
