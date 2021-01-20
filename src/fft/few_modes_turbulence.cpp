//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file few_modes_turbulence.cpp
//  \brief implementation of functions in class Turbulence

// C/C++ headers
#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>   // sstream
#include <stdexcept> // runtime_error
#include <string>    // c_str()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../utils/utils.hpp"
#include "few_modes_turbulence.hpp"

//----------------------------------------------------------------------------------------
//! \fn FewModesTurbulenceDriver::FewModesTurbulenceDriver(Mesh *pm,
//                                                         ParameterInput *pin,
//                                                         int64_t rseed_in)
//  \brief FewModesTurbulenceDriver constructor. If rseed_in != -1 -> restarted sim.

FewModesTurbulenceDriver::FewModesTurbulenceDriver(Mesh *pm, ParameterInput *pin,
                                                   int64_t rseed_in) {
  rseed = pin->GetOrAddInteger("problem", "rseed", -1); // seed for random number.

  kpeak = pin->GetOrAddReal("problem", "kpeak", 0.0);  // peak of the forcing spec
  accel_rms = pin->GetReal("problem", "accel_rms");    // turbulence amplitude
  tcorr = pin->GetReal("problem", "corr_time");        // forcing autocorrelation time
  sol_weight = pin->GetReal("problem", "sol_weight");  // solenoidal weight
  num_modes = pin->GetInteger("problem", "num_modes"); // number of wavemodes

  if ((num_modes > 100) && (Globals::my_rank == 0)) {
    std::cout << "### WARNING using more than 100 explicit modes will significantly "
              << "increase the runtime." << std::endl
              << "If many modes are required in the acceleration field consider using "
              << "the driving mechanism based on full FFTs." << std::endl;
  }

  if (pm->fmturb_flag == 0) {
    std::stringstream msg;
    msg << "### FATAL ERROR in FewModesTurbulenceDriver::TurbulenceDriver" << std::endl
        << "Turbulence flag is set to zero! Shouldn't reach here!" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  if (pm->nblist[Globals::my_rank] != 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in FewModesTurbulenceDriver::TurbulenceDriver" << std::endl
        << "Currently only one meshblock per process is supported!" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  auto Lx = pin->GetReal("mesh", "x1max") - pin->GetReal("mesh", "x1min");
  auto Ly = pin->GetReal("mesh", "x2max") - pin->GetReal("mesh", "x2min");
  auto Lz = pin->GetReal("mesh", "x3max") - pin->GetReal("mesh", "x3min");
  if ((Lx != 1.0) || (Ly != 1.0) || (Lz != 1.0)) {
    std::stringstream msg;
    msg << "### FATAL ERROR in FewModesTurbulenceDriver::TurbulenceDriver" << std::endl
        << "Only domain sizes with edge lengths of 1 are supported." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  auto gnx1 = pm->mesh_size.nx1;
  auto gnx2 = pm->mesh_size.nx2;
  auto gnx3 = pm->mesh_size.nx3;
  if ((gnx1 != gnx2) || (gnx2 != gnx3)) {
    std::stringstream msg;
    msg << "### FATAL ERROR in FewModesTurbulenceDriver::TurbulenceDriver" << std::endl
        << "Only cubic mesh sizes are supported." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  int nx1 = pm->my_blocks(0)->block_size.nx1;
  int nx2 = pm->my_blocks(0)->block_size.nx2;
  int nx3 = pm->my_blocks(0)->block_size.nx3;

  gis = pm->my_blocks(0)->loc.lx1 * pm->my_blocks(0)->block_size.nx1;
  gjs = pm->my_blocks(0)->loc.lx2 * pm->my_blocks(0)->block_size.nx2;
  gks = pm->my_blocks(0)->loc.lx3 * pm->my_blocks(0)->block_size.nx3;

  this->pm = pm;

  // velocity array does not contain ghost zones
  accel.NewAthenaArray(3, nx3, nx2, nx1);

  // Acceleration field in Fourier space using complex to real transform.
  accel_hat.NewAthenaArray(3, num_modes);
  for (int n = 0; n < 3; n++)
    for (int m = 0; m < num_modes; m++)
      accel_hat(n,m) = Complex(0., 0.);

  // no need to init the following vars as they're set every cycle
  accel_hat_new.NewAthenaArray(3, num_modes);
  phases_i.NewAthenaArray(nx1, num_modes);
  phases_j.NewAthenaArray(nx2, num_modes);
  phases_k.NewAthenaArray(nx3, num_modes);

  // list of wavenumber vectors
  k_vec.NewAthenaArray(num_modes, 3);

  // temp storage of random numbers
  random_num.NewAthenaArray(3, num_modes, 2);

  int k_vec_in[3];
  for (int i = 1; i <= num_modes; i++) {
    pin->GetIntegerVector("modes", "k_" + std::to_string(i), k_vec_in);
    k_vec(i - 1, 0) = k_vec_in[0];
    k_vec(i - 1, 1) = k_vec_in[1];
    k_vec(i - 1, 2) = k_vec_in[2];
  }

  SetPhases();

  // if this is a restarted simulation recover acc field and current random seed
  if (rseed_in != -1) {
    RestoreFromRestart(rseed_in);
  }
}

// destructor
FewModesTurbulenceDriver::~FewModesTurbulenceDriver() {
  // TODO(pgrete) need cleanup of Athenaarrays?
}

//----------------------------------------------------------------------------------------
//! \fn void FewModesTurbulenceDriver::RestoreFromRestart(int64_t rseed_in)
//  \brief Restores the random seed and spectral acc field from restart dumps

void FewModesTurbulenceDriver::RestoreFromRestart(int64_t rseed_in) {
  // fast forward original seed to current seed
  // (this way we don't need to save the full ran2 state)
  while (rseed_in != rseed) {
    ran2(&rseed);
  }

  for (int n = 0; n < 3; n++) {
    for (int m = 0; m < num_modes; m++) {
      accel_hat(n, m) = Complex(pm->my_blocks(0)->ruser_meshblock_data[0](n, m, 0),
                                pm->my_blocks(0)->ruser_meshblock_data[0](n, m, 1));
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void FewModesTurbulenceDriver::SetPhases(void)
//  \brief Precalculate phase factors used in the inverse transform

void FewModesTurbulenceDriver::SetPhases(void) {
  Complex I(0.0, 1.0);
  auto gnx1 = pm->mesh_size.nx1;
  auto gnx2 = pm->mesh_size.nx2;
  auto gnx3 = pm->mesh_size.nx3;
  int nx1 = pm->my_blocks(0)->block_size.nx1;
  int nx2 = pm->my_blocks(0)->block_size.nx2;
  int nx3 = pm->my_blocks(0)->block_size.nx3;

  for (int i = 0; i < nx1; i++) {
    Real gi = static_cast<Real>(i + gis);
    Real w_kx;

    for (int m = 0; m < num_modes; m++) {
      w_kx = k_vec(m, 0) * 2. * PI / static_cast<Real>(gnx1);
      // adjust phase factor to Complex->Real IFT: u_hat*(k) = u_hat(-k)
      if (k_vec(m, 0) == 0.0) {
        phases_i(i, m) = 0.5 * std::exp(I * w_kx * gi);
      } else {
        phases_i(i, m) = std::exp(I * w_kx * gi);
      }
    }
  }

  for (int j = 0; j < nx2; j++) {
    Real gj = static_cast<Real>(j + gjs);
    Real w_ky;

    for (int m = 0; m < num_modes; m++) {
      w_ky = k_vec(m, 1) * 2. * PI / static_cast<Real>(gnx2);
      phases_j(j, m) = std::exp(I * w_ky * gj);
    }
  }

  for (int k = 0; k < nx3; k++) {
    Real gk = static_cast<Real>(k + gks);
    Real w_kz;

    for (int m = 0; m < num_modes; m++) {
      w_kz = k_vec(m, 2) * 2. * PI / static_cast<Real>(gnx3);
      phases_k(k, m) = std::exp(I * w_kz * gk);
    }
  }
}

//----------------------------------------------------------------------------------------
//! \fn void FewModesTurbulenceDriver::Driving(void)
//  \brief Generate and Perturb the velocity field

void FewModesTurbulenceDriver::Driving(void) {
  // evolve forcing
  Generate(pm->dt);

  switch (pm->fmturb_flag) {
  case 1: // fmturb_flag == 1 : continuously driven turbulence
    Perturb(pm->dt);
    break;
  default:
    std::stringstream msg;
    msg << "### FATAL ERROR in FewModesTurbulenceDriver::Driving" << std::endl
        << "Turbulence flag " << pm->fmturb_flag << " is not supported!" << std::endl;
    ATHENA_ERROR(msg);
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void FewModesTurbulenceDriver::Generate()
//  \brief Generate velocity pertubation.

void FewModesTurbulenceDriver::Generate(Real dt) {
  int nx1 = pm->my_blocks(0)->block_size.nx1;
  int nx2 = pm->my_blocks(0)->block_size.nx2;
  int nx3 = pm->my_blocks(0)->block_size.nx3;
  Complex I(0.0, 1.0);

  Real v1, v2, v_sqr;
  for (int n = 0; n < 3; n++)
    for (int m = 0; m < num_modes; m++) {
      do {
        v1 = 2.0 * ran2(&rseed) - 1.0;
        v2 = 2.0 * ran2(&rseed) - 1.0;
        v_sqr = v1 * v1 + v2 * v2;
      } while (v_sqr >= 1.0 || v_sqr == 0.0);

      random_num(n, m, 0) = v1;
      random_num(n, m, 1) = v2;
    }

  // generate new power spectrum (injection)
  for (int n = 0; n < 3; n++) {
    #pragma omp simd
    for (int m = 0; m < num_modes; m++) {

      Real kmag, tmp, norm, v_sqr;

      Real &kx = k_vec(m, 0);
      Real &ky = k_vec(m, 1);
      Real &kz = k_vec(m, 2);

      kmag = std::sqrt(kx * kx + ky * ky + kz * kz);

      accel_hat_new(n, m) = Complex(0., 0.);

      tmp = std::pow(kmag / kpeak, 2.) * (2. - std::pow(kmag / kpeak, 2.));
      if (tmp < 0.) tmp = 0.;
      v_sqr = SQR(random_num(n, m, 0)) + SQR(random_num(n, m, 1));
      norm = std::sqrt(-2.0 * std::log(v_sqr) / v_sqr);

      accel_hat_new(n, m) =
          Complex(tmp * norm * random_num(n, m, 0), tmp * norm * random_num(n, m, 1));
    }
  }

  // enforce symmetry of complex to real transform
  for (int n = 0; n < 3; n++) {
    for (int m = 0; m < num_modes; m++) {
      if (k_vec(m, 0) == 0.) {
        for (int m2 = 0; m2 < m; m2++) {
          if (k_vec(m, 1) == -k_vec(m2, 1) && k_vec(m, 2) == -k_vec(m2, 2))
            accel_hat_new(n, m) =
                Complex(accel_hat_new(n, m2).real(), -accel_hat_new(n, m2).imag());
        }
      }
    }
  }

  // project
  for (int m = 0; m < num_modes; m++) {
    Real kmag;

    Real kx = k_vec(m, 0);
    Real ky = k_vec(m, 1);
    Real kz = k_vec(m, 2);

    kmag = std::sqrt(kx * kx + ky * ky + kz * kz);

    // setting kmag to 1 as a "continue" doesn't work within the parallel_for construct
    // and it doesn't affect anything (there should never be power in the k=0 mode)
    if (kmag == 0.) kmag = 1.;

    // make it a unit vector
    kx /= kmag;
    ky /= kmag;
    kz /= kmag;

    Complex dot(accel_hat_new(0, m).real() * kx + accel_hat_new(1, m).real() * ky +
                    accel_hat_new(2, m).real() * kz,
                accel_hat_new(0, m).imag() * kx + accel_hat_new(1, m).imag() * ky +
                    accel_hat_new(2, m).imag() * kz);

    accel_hat_new(0, m) = Complex(accel_hat_new(0, m).real() * sol_weight +
                                      (1. - 2. * sol_weight) * dot.real() * kx,
                                  accel_hat_new(0, m).imag() * sol_weight +
                                      (1. - 2. * sol_weight) * dot.imag() * kx);
    accel_hat_new(1, m) = Complex(accel_hat_new(1, m).real() * sol_weight +
                                      (1. - 2. * sol_weight) * dot.real() * ky,
                                  accel_hat_new(1, m).imag() * sol_weight +
                                      (1. - 2. * sol_weight) * dot.imag() * ky);
    accel_hat_new(2, m) = Complex(accel_hat_new(2, m).real() * sol_weight +
                                      (1. - 2. * sol_weight) * dot.real() * kz,
                                  accel_hat_new(2, m).imag() * sol_weight +
                                      (1. - 2. * sol_weight) * dot.imag() * kz);
  }

  // evolve
  Real c_drift = std::exp(-dt / tcorr);
  Real c_diff = std::sqrt(1.0 - c_drift * c_drift);

  for (int n = 0; n < 3; n++) {
    for (int m = 0; m < num_modes; m++) {
      accel_hat(n, m) =
          Complex(accel_hat(n, m).real() * c_drift + accel_hat_new(n, m).real() * c_diff,
                  accel_hat(n, m).imag() * c_drift + accel_hat_new(n, m).imag() * c_diff);
    }
  }

  // inverse FT
  // implictly assuming cubic box of size L=1
  for (auto n = 0; n <= 2; n++)
    for (auto k = 0; k <= nx3 - 1; k++)
      for (auto j = 0; j <= nx2 - 1; j++)
#pragma omp simd
        for (auto i = 0; i <= nx1 - 1; i++) {

          Complex phase;
          accel(n, k, j, i) = 0.0;

          for (int m = 0; m < num_modes; m++) {
            phase = phases_i(i, m) * phases_j(j, m) * phases_k(k, m);
            accel(n, k, j, i) += 2. * (accel_hat(n, m).real() * phase.real() -
                                       accel_hat(n, m).imag() * phase.imag());
          }
        }
}

//----------------------------------------------------------------------------------------
//! \fn void FewModesTurbulenceDriver::Perturb(Real dt)
//  \brief Add velocity perturbation to the hydro variables

void FewModesTurbulenceDriver::Perturb(Real dt) {
  std::stringstream msg;

  int is = pm->my_blocks(0)->is, ie = pm->my_blocks(0)->ie;
  int js = pm->my_blocks(0)->js, je = pm->my_blocks(0)->je;
  int ks = pm->my_blocks(0)->ks, ke = pm->my_blocks(0)->ke;

  Real m[4] = {0};

  for (auto k = ks; k <= ke; k++)
    for (auto j = js; j <= je; j++)
      for (auto i = is; i <= ie; i++) {
        Real den = pm->my_blocks(0)->phydro->u(IDN, k, j, i);
        m[0] += den;
        m[1] += den * accel(0, k - ks, j - js, i - is);
        m[2] += den * accel(1, k - ks, j - js, i - is);
        m[3] += den * accel(2, k - ks, j - js, i - is);
      }

#ifdef MPI_PARALLEL
  Real gm[4];
  int mpierr;

  // Sum the perturbations over all processors
  mpierr = MPI_Allreduce(m, gm, 4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (mpierr) {
    msg << "[normalize]: MPI_Allreduce error = " << mpierr << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  // store result of reduction in local var
  for (int i = 0; i < 4; i++)
    m[i] = gm[i];
#endif // MPI_PARALLEL

  // remove mean momentum and calc normalization
  Real ampl_sum = 0.;
  for (auto n = 0; n <= 2; n++)
    for (auto k = 0; k <= ke - ks; k++)
      for (auto j = 0; j <= je - js; j++)
        for (auto i = 0; i <= ie - is; i++) {
          accel(n, k, j, i) -= m[n + 1] / m[0];
          ampl_sum += SQR(accel(n, k, j, i));
        }

#ifdef MPI_PARALLEL
  // Sum the perturbations over all processors
  mpierr = MPI_Allreduce(&ampl_sum, gm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (mpierr) {
    msg << "[normalize]: MPI_Allreduce error = " << mpierr << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  ampl_sum = gm[0];
#endif // MPI_PARALLEL

  Real gsize =
      (static_cast<Real>(pm->mesh_size.nx1) * static_cast<Real>(pm->mesh_size.nx2) *
       static_cast<Real>(pm->mesh_size.nx3));

  Real norm = accel_rms / std::sqrt(ampl_sum / gsize);
  // apply momentum perturb
  for (auto k = ks; k <= ke; k++)
    for (auto j = js; j <= je; j++)
#pragma omp simd
      for (auto i = is; i <= ie; i++) {

        // normalizing accel field here so that the actual values are used in the output
        accel(0, k - ks, j - js, i - is) *= norm;
        accel(1, k - ks, j - js, i - is) *= norm;
        accel(2, k - ks, j - js, i - is) *= norm;

        Real qa = dt * pm->my_blocks(0)->phydro->u(IDN, k, j, i);
        if (NON_BAROTROPIC_EOS) {
          pm->my_blocks(0)->phydro->u(IEN, k, j, i) +=
              (pm->my_blocks(0)->phydro->u(IM1, k, j, i) * dt *
                   accel(0, k - ks, j - js, i - is) +
               pm->my_blocks(0)->phydro->u(IM2, k, j, i) * dt *
                   accel(1, k - ks, j - js, i - is) +
               pm->my_blocks(0)->phydro->u(IM3, k, j, i) * dt *
                   accel(2, k - ks, j - js, i - is) +
               (SQR(accel(0, k - ks, j - js, i - is)) +
                SQR(accel(1, k - ks, j - js, i - is)) +
                SQR(accel(2, k - ks, j - js, i - is))) *
                   qa * qa / (2 * pm->my_blocks(0)->phydro->u(IDN, k, j, i)));
        }

        pm->my_blocks(0)->phydro->u(IM1, k, j, i) +=
            qa * accel(0, k - ks, j - js, i - is);
        pm->my_blocks(0)->phydro->u(IM2, k, j, i) +=
            qa * accel(1, k - ks, j - js, i - is);
        pm->my_blocks(0)->phydro->u(IM3, k, j, i) +=
            qa * accel(2, k - ks, j - js, i - is);
      }

  return;
}
//----------------------------------------------------------------------------------------
//! \fn void FewModesTurbulenceDriver::CopyAccelToOutputVars(
//        AthenaArray<Real> &user_out_var)
//  \brief Write acceleration field to output buffer

void FewModesTurbulenceDriver::CopyAccelToOutputVars(AthenaArray<Real> &user_out_var) {

  // ks,ke ... are member variables, and so need to be aliased to be
  // captured by the lambda
  auto is = pm->my_blocks(0)->is, ie = pm->my_blocks(0)->ie;
  auto js = pm->my_blocks(0)->js, je = pm->my_blocks(0)->je;
  auto ks = pm->my_blocks(0)->ks, ke = pm->my_blocks(0)->ke;
  // Copy acceleration field for analysis dump
  for (auto n = 0; n <= 2; n++)
    for (auto k = ks; k <= ke; k++)
      for (auto j = js; j <= je; j++)
#pragma omp simd
        for (auto i = is; i <= ie; i++) {
          user_out_var(n, k, j, i) = accel(n, k - ks, j - js, i - is);
        }

  // Copy accel_hat to ruser_meshblock_data
  for (auto n = 0; n < 3; n++)
    for (auto m = 0; m < num_modes; m++) {
      pm->my_blocks(0)->ruser_meshblock_data[0](n, m, 0) = accel_hat(n, m).real();
      pm->my_blocks(0)->ruser_meshblock_data[0](n, m, 1) = accel_hat(n, m).imag();
    }

  pm->rseed_rst = rseed;

  return;
}
