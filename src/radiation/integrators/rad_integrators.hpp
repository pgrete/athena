#ifndef RADINTEGRATORS_HPP
#define RADINTEGRATORS_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file rad_integrators.hpp
//  \brief definitions for RadIntegrator class
//======================================================================================

// Athena++ classes headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../radiation.hpp" // radiation
#ifdef INCLUDE_CHEMISTRY
#include "../../chemistry/species.hpp"
#endif

class MeshBlock;
class ParameterInput;
class Radiation;
class NeighborBlock;

//! \class RadIntegrator
//  \brief integrate algorithm for radiative transfer


class RadIntegrator {
  friend class Radiation;
public:
  RadIntegrator(Radiation *prad, ParameterInput *pin);
  ~RadIntegrator();
  
  Radiation *pmy_rad;
  MeshBlock *pmy_mb;

#ifdef INCLUDE_CHEMISTRY
  int ncol;
  AthenaArray<Real> col, col_avg;
  AthenaArray<Real> col_Htot, col_CO, col_H2;//TODO:for output
  ChemNetwork* pmy_chemnet;
  //calcuate column within each meshblock
  void GetColMB(int direction);
  //calcuate total column and update radiation
  void UpdateRadiation(int direction);
  //copy column density and average radiation field to output
  void CopyToOutput();
  void SetSixRayNeighbors();
#endif
private:
  Real rad_G0_; //unshielded radiation field strengh, uniform.
  Real unit_length_in_cm_;
  //index for direction of rays in six-ray
  static const int IXP = 0; //+x
  static const int IYP = 1; //+y
  static const int IZP = 2; //+z
  static const int IXM = 3; //-x
  static const int IYM = 4; //-y
  static const int IZM = 5; //-z
  NeighborBlock* pneighbors_[6]; //TODO:uniform mesh for now
};

#endif // RADINTEGRATORS_HPP
