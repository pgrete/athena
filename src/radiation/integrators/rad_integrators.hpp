#ifndef RADINTEGRATORS_HPP
#define RADINTEGRATORS_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file radiation.hpp
//  \brief definitions for Radiation class
//======================================================================================

// Athena++ classes headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../radiation.hpp" // radiation

class MeshBlock;
class ParameterInput;
class Radiation;

//! \class RadIntegrator
//  \brief integrate algorithm for radiative transfer


class RadIntegrator {
public:
  RadIntegrator(Radiation *prad, ParameterInput *pin);
  ~RadIntegrator();
  
  Radiation *pmy_rad;
  
  void FluxDivergence(MeshBlock *pmb, AthenaArray<Real> &ir, const int step);
  void CalculateFluxes(MeshBlock *pmb, AthenaArray<Real> &ir, const int step);

  void FirstOrderFluxX1(const int k, const int j, const int il,
        const int iu, const AthenaArray<Real> &q,
        const AthenaArray<Real> &vel, AthenaArray<Real> &flx);
  
  void FirstOrderFluxX2(const int k, const int j, const int il,
        const int iu, const AthenaArray<Real> &q,
        const AthenaArray<Real> &vel, AthenaArray<Real> &flx);
  
  void FirstOrderFluxX3(const int k, const int j, const int il,
        const int iu, const AthenaArray<Real> &q,
        const AthenaArray<Real> &vel, AthenaArray<Real> &flx);
 
  
  void SecondOrderFluxX1(const int k, const int j, const int il,
        const int iu, const AthenaArray<Real> &q,
        const AthenaArray<Real> &vel, AthenaArray<Real> &flx);
  
  void SecondOrderFluxX2(const int k, const int j, const int il,
        const int iu, const AthenaArray<Real> &q,
        const AthenaArray<Real> &vel, AthenaArray<Real> &flx);
  
  void SecondOrderFluxX3(const int k, const int j, const int il,
        const int iu, const AthenaArray<Real> &q,
        const AthenaArray<Real> &vel, AthenaArray<Real> &flx);
  
  void ThirdOrderFluxX1(const int k, const int j, const int il,
        const int iu, const AthenaArray<Real> &q,
        const AthenaArray<Real> &vel, AthenaArray<Real> &flx);
  
  void ThirdOrderFluxX2(const int k, const int j, const int il,
        const int iu, const AthenaArray<Real> &q,
        const AthenaArray<Real> &vel, AthenaArray<Real> &flx);
  
  void ThirdOrderFluxX3(const int k, const int j, const int il,
        const int iu, const AthenaArray<Real> &q,
        const AthenaArray<Real> &vel, AthenaArray<Real> &flx);

private:
  AthenaArray<Real> flx_, vel_; // temporary array to store the flux, velocity
  AthenaArray<Real> temp_i1_, temp_i2_; // temporary array to store Div q
                                    // also used to store co-moving frame intensity
  AthenaArray<Real> x1face_area_, x2face_area_, x3face_area_;
  AthenaArray<Real> x2face_area_p1_, x3face_area_p1_;
  AthenaArray<Real> cell_volume_;

};

#endif // RADINTEGRATORS_HPP
