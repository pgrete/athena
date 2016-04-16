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
#include <string>



class MeshBlock;
class ParameterInput;
class RadIntegrator;

//! \class Radiation
//  \brief radiation data and functions


// prototype for user-defined opacity function for radiative transfer
typedef void (*Opacity_t)(MeshBlock *pmb, AthenaArray<Real> &prim);
typedef void (*OutInternal_t)(MeshBlock *pmb);

// Array indices for radiation moments
enum {IER=0, IFR1=1, IFR2=2, IFR3=3, IPR11=4, IPR12=5, IPR13=6, IPR21=7,
      IPR22=8, IPR23=9, IPR31=10, IPR32=11, IPR33=12};

class Radiation {
  friend class RadIntegrator;
public:
  Radiation(MeshBlock *pmb, ParameterInput *pin);
  ~Radiation();
    
  AthenaArray<Real> ir, ir1; // radiation specific intensity
  AthenaArray<Real> rad_mom; // frequency integrated radiation moments
  AthenaArray<Real> rad_mom_cm; // co-moving frame Er, Frx, Fry, Frz
  AthenaArray<Real> sigma_s, sigma_a, sigma_ae; //   opacity
                       //sigma_a T and sigma_ae I
  AthenaArray<Real> grey_sigma_s, grey_sigma_a; // frequency integrated opacity
  AthenaArray<Real> rad_ifov; // internal radiation variable
  AthenaArray<Real> mu, wmu; // angles and weight
  AthenaArray<Real> wfreq; // weight in frequency space
  
  AthenaArray<Real> flux[3]; // store transport flux, also need for refinement
  
  Real prat, crat; // prat=aT^4/P_0, crat=c/c_s
  Real vmax;
  Real reduced_c; // reduced speed of light
  Real tunit, telectron; // gas temperature cgs unit,
                         // effective electron scattering temperature
  
  int nang, nfreq, noct, n_fre_ang; // n_fre_ang=nang*nfreq

  MeshBlock* pmy_block;    // ptr to MeshBlock containing this Fluid
  
  RadIntegrator *pradintegrator;
  
  
  //Function in problem generators to update opacity
  void EnrollOpacityFunction(Opacity_t MyOpacityFunction);
  
    //Function in problem generators to update opacity
  void EnrollInternalVariableFunction(OutInternal_t MyOutputInternalFunction);
  
  // The function pointer for the opacity
  Opacity_t UpdateOpacity;
  
  // Function pointer to load internal variable output
  OutInternal_t LoadInternalVariable;
  
  //functin to calculate the radiation moments
  void CalculateMoment(AthenaArray<Real> &ir_in);
  void CalculateComMoment();

  
  void AngularGrid(int angle_flag, int nmu);

  void FrequencyGrid();


private:
  Real t_floor_; // temperature floor
  

};

#endif // RADIATION_HPP
