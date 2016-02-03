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


// prototype for user-defined opacity function for radiative transfer
typedef void (*Opacity_t)(MeshBlock *pmb, AthenaArray<Real> &prim);

// Array indices for radiation moments
enum {IER=0, IFR1=1, IFR2=2, IFR3=3, IPR11=4, IPR12=5, IPR13=6, IPR21=7,
      IPR22=8, IPR23=9, IPR31=10, IPR32=11, IPR33=12};

class Radiation {
public:
  Radiation(MeshBlock *pmb, ParameterInput *pin);
  ~Radiation();
    
  AthenaArray<Real> ir, ir1; // radiation specific intensity
  AthenaArray<Real> rad_mom; // frequency integrated radiation moments
  AthenaArray<Real> sigma_s, sigma_a; //   opacity
  AthenaArray<Real> mu, wmu; // angles and weight
  AthenaArray<Real> wfreq; // weight in frequency space
  
  Real prat, crat; // prat=aT^4/P_0, crat=c/c_s
  Real reduced_c; // reduced speed of light
  
  int nang, nfreq, noct;

  MeshBlock* pmy_block;    // ptr to MeshBlock containing this Fluid

  
  
  //Function in problem generators to update opacity
  void EnrollOpacityFunction(Opacity_t MyOpacityFunction);
  
  // The function pointer for the opacity
  Opacity_t UpdateOpacity;
  
  //functin to calculate the radiation moments
  void CalculateMoment();
  
  void AngularGrid(int angle_flag, int nmu);

  void FrequencyGrid();


private:


};

#endif // RADIATION_HPP
