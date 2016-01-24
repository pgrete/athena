#ifndef RADIATION_HPP
#define RADIATION_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file radiation.hpp
//  \brief definitions for Radiation class
//======================================================================================

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

class MeshBlock;
class ParameterInput;

//! \class Radiation
//  \brief radiation data and functions

class Radiation {
public:
  Radiation(MeshBlock *pmb, ParameterInput *pin);
  ~Radiation();
    
  AthenaArray<Real> ir, ir1; // radiation specific intensity
  AthenaArray<Real> sigma_s, sigma_a; //   opacity
  AthenaArray<Real> mu, wmu; // angles and weight
  
  Real prat, crat; // prat=aT^4/P_0, crat=c/c_s
  Real reduced_c; // reduced speed of light
  
  int nang, nfreq;

private:


};
#endif // RADIATION_HPP
