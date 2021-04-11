//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file turb.cpp
//  \brief Problem generator for turbulence generator
//

// C++ headers
#include <array>
#include <cmath>
#include <ctime>
#include <random>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

namespace cooling {

Real T_low;    // lower temperature threshold for which cooling is active
Real tau_cool; // cooling timescale
Real max_diff; // max. fraction of change of internal energy density per step
Real T0;       // Target temperature (to be heated and cooled to)

void Cooling(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {
  const auto coeff = -dt / (tau_cool * (pmb->peos->GetGamma() - 1.0));
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
#pragma omp simd
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        const auto T = prim(IPR, k, j, i) / prim(IDN, k, j, i); // same "units" as T_low
        if (T > T_low) {
          cons(IEN, k, j, i) += coeff * prim(IDN, k, j, i) * (T - T0);
        }
      }
    }
  }
}

Real CoolDt(MeshBlock *pmb) {
  auto min_dt = std::numeric_limits<Real>::max();
  const AthenaArray<Real> &prim = pmb->phydro->w;
  const auto gm1 = pmb->peos->GetGamma() - 1.0;
  const auto coeff = 1.0 / (tau_cool * gm1);
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        const auto T = prim(IPR, k, j, i) / prim(IDN, k, j, i); // same "units" as T_low
        if (T > T_low) {
          const auto Edot = coeff * prim(IDN, k, j, i) * (T - T0);
          const auto E = prim(IPR, k, j, i) / gm1;
          const Real dt = std::abs(max_diff * E / (Edot + TINY_NUMBER));
          min_dt = std::min(min_dt, dt);
        }
      }
    }
  }
  return min_dt;
}
} // namespace cooling

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief
//========================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {

  const auto p0 = pin->GetReal("problem", "p0");
  const auto rho0 = pin->GetReal("problem", "rho0");
  const auto thres = pin->GetOrAddReal("problem", "cool_threshold_ratio", 1. / 30.);
  // this is a proxy for the temp, which we can use given that gamma and c_V are const
  cooling::T0 = p0 / rho0;
  cooling::T_low = cooling::T0 * thres;
  cooling::max_diff = pin->GetOrAddReal("problem", "cool_max_diff", 0.1);
  cooling::tau_cool = pin->GetOrAddReal("problem", "tau_cool", 0.00125);

  // EnrollUserExplicitSourceFunction(cooling::Cooling);
  // EnrollUserTimeStepFunction(cooling::CoolDt);
}

//========================================================================================
//! \fn void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in MeshBlock class.  Can also be
//  used to initialize variables which are global to other functions in this file.
//  Called in MeshBlock constructor before ProblemGenerator.
//========================================================================================

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) { return; }

//========================================================================================
//! \fn void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
//  \brief Function called before generating output files
//========================================================================================

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) { return; }

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real gamma = peos->GetGamma();
  Real gm1 = gamma - 1.0;
  Real p0 = pin->GetOrAddReal("problem", "p0", 0.3);
  Real rho0 = pin->GetOrAddReal("problem", "rho0", 1.0);
  Real Bx0 = pin->GetOrAddReal("problem", "Bx0", 0.056117);

  // All uniform so we can directly initialize (including face fields) with larger bounds
  for (int k = ks; k <= ke + 1; k++) {
    for (int j = js; j <= je + 1; j++) {
      for (int i = is; i <= ie + 1; i++) {
        phydro->u(IDN, k, j, i) = rho0;

        phydro->u(IM1, k, j, i) = 0.0;
        phydro->u(IM2, k, j, i) = 0.0;
        phydro->u(IM3, k, j, i) = 0.0;
        pfield->b.x1f(k, j, i) = Bx0;
        pfield->b.x2f(k, j, i) = 0.0;
        pfield->b.x3f(k, j, i) = 0.0;

        phydro->u(IEN, k, j, i) = p0 / gm1 + 0.5 * SQR(Bx0);
      }
    }
  }
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {}

namespace {

constexpr int num_blast = 30;
std::array<std::array<Real, 3>, num_blast> blasts = {{
    {7.825E-07, 1.32E-02, 7.56E-02},     {-5.413E-02, -4.672E-02, -7.810E-02},
    {-3.211E-02, 6.793E-02, 9.346E-02},  {-6.165E-02, 5.194E-02, -1.690E-02},
    {5.346E-03, 5.297E-02, 6.711E-02},   {7.698E-04, -6.165E-02, -9.331E-02},
    {4.174E-02, 6.867E-02, 5.889E-02},   {9.304E-02, -1.538E-02, 5.269E-02},
    {9.196E-03, -3.460E-02, -5.840E-02}, {7.011E-02, 9.103E-02, -2.378E-02},
    {-7.375E-02, 4.746E-03, -2.639E-02}, {3.653E-02, 2.470E-02, -1.745E-03},
    {7.268E-03, -3.683E-02, 8.847E-02},  {-7.272E-02, 4.364E-02, 7.664E-02},
    {4.777E-02, -7.622E-02, -7.250E-02}, {-1.023E-02, 9.08E-03, 6.06E-03},
    {-9.534E-03, -4.954E-02, 5.162E-02}, {-9.092E-02, -5.223E-03, 7.374E-03},
    {9.138E-02, 5.297E-02, -5.355E-02},  {9.409E-02, -9.499E-02, 7.615E-02},
    {7.702E-02, 8.278E-02, -8.746E-02},  {-7.306E-02, -5.846E-02, 5.373E-02},
    {4.679E-02, 2.872E-02, -8.216E-02},  {7.482E-02, 5.545E-02, 8.907E-02},
    {6.248E-02, -1.579E-02, -8.402E-02}, {-9.090E-02, 2.745E-02, -5.857E-02},
    {-1.130E-02, 6.520E-02, -8.496E-02}, {-3.186E-02, 3.858E-02, 3.877E-02},
    {4.997E-02, -8.524E-02, 5.871E-02},  {8.455E-02, -4.098E-02, -4.438E-02},
}};

} // namespace

//========================================================================================
//! \fn void MeshBlock::UserWorkInLoop()
//! \brief Function called once every time step for user-defined work.
//========================================================================================

void MeshBlock::UserWorkInLoop() {
  const Real dt_between_blasts = 0.00125;
  Real time_this_blast = -1.0;
  int blast_i = -1; // numer of blast to use. Negative -> no blast

  for (int i = 0; i < num_blast; i++) {
    time_this_blast = static_cast<Real>(i + 1) * dt_between_blasts;
    // this blast falls in intervall of current timestep
    if ((time_this_blast >= pmy_mesh->time) &&
        (time_this_blast < pmy_mesh->time + pmy_mesh->dt)) {
      blast_i = i;
      break;
    }
  }

  if (blast_i < 0) {
    return;
  }
  const auto gm1 = peos->GetGamma() - 1.0;

  for (int k = ks - NGHOST; k <= ke + NGHOST; k++) {
    for (int j = js - NGHOST; j <= je + NGHOST; j++) {
#pragma omp simd
      for (int i = is - NGHOST; i <= ie + NGHOST; i++) {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);
        Real dist = std::sqrt(SQR(x - blasts.at(blast_i).at(0)) +
                              SQR(y - blasts.at(blast_i).at(1)) +
                              SQR(z - blasts.at(blast_i).at(2)));

        if (dist < 0.005) {
          // Updating both primitive and converved var for consistency as this function is
          // called AFTER cons2prim at the very end of each timestep (not substage).
          // Therefore, we also update ghost cells for consistency.
          phydro->w(IPR,k,j,i) = 13649.6;
          phydro->u(IEN, k, j, i) =
              13649.6 / gm1 +
              0.5 * (SQR(0.5 * (pfield->b.x1f(k, j, i) + pfield->b.x1f(k, j, i + 1))) +
                     SQR(0.5 * (pfield->b.x2f(k, j, i) + pfield->b.x2f(k, j + 1, i))) +
                     SQR(0.5 * (pfield->b.x3f(k, j, i) + pfield->b.x3f(k + 1, j, i)))) +
              (0.5 / phydro->u(IDN, k, j, i)) *
                  (SQR(phydro->u(IM1, k, j, i)) + SQR(phydro->u(IM2, k, j, i)) +
                   SQR(phydro->u(IM3, k, j, i)));
        }
      }
    }
  }
}
