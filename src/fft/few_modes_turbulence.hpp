#ifndef FFT_FEW_MODES_TURBULENCE_HPP_
#define FFT_FEW_MODES_TURBULENCE_HPP_

//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file turbulence.hpp
//  \brief defines Turbulence class

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

class Mesh;
class MeshBlock;
class ParameterInput;
class Coordinates;

//! \class FewModesTurbulenceDriver
//  \brief FewModesTurbulence Driver

class FewModesTurbulenceDriver {
 public:
  FewModesTurbulenceDriver(Mesh *pm, ParameterInput *pin, int64_t rseed_in = -1);
  ~FewModesTurbulenceDriver();
  void Driving(void);
  void Generate(Real dt);
  void Perturb(Real dt);
  void CopyAccelToOutputVars(AthenaArray<Real> &user_out_var);
  void SetPhases(void);
  void RestoreFromRestart(int64_t rseed_in);

 private:
  int64_t rseed;
  Real kpeak;
  int gis, gjs, gks; // global indices
  int num_modes;
  Real tcorr, sol_weight, accel_rms;
  AthenaArray<Real> accel, k_vec;
  AthenaArray<Complex> accel_hat, accel_hat_new;
  AthenaArray<Complex> phases_i, phases_j, phases_k;
  Mesh *pm;
};

#endif // FFT_FEW_MODES_TURBULENCE_HPP_
