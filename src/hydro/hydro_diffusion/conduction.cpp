//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
// Anisotropic conduction implemented by Philipp Grete adapted from Michael Jennings
//========================================================================================
//! \file conduction.cpp
//! \brief

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../eos/eos.hpp"
#include "../../field/field.hpp"
#include "../../mesh/mesh.hpp"
#include "../hydro.hpp"
#include "hydro_diffusion.hpp"

//---------------------------------------------------------------------------------------
//! Calculate isotropic thermal conduction

void HydroDiffusion::ThermalFluxIso(const AthenaArray<Real> &prim,
                                    const AthenaArray<Real> &cons,
                                    AthenaArray<Real> *cndflx) {
  const bool f2 = pmb_->pmy_mesh->f2;
  const bool f3 = pmb_->pmy_mesh->f3;
  AthenaArray<Real> &x1flux = cndflx[X1DIR];
  int il, iu, jl, ju, kl, ku;
  int is = pmb_->is;
  int js = pmb_->js;
  int ks = pmb_->ks;
  int ie = pmb_->ie;
  int je = pmb_->je;
  int ke = pmb_->ke;
  Real kappaf, denf, dTdx, dTdy, dTdz;

  // i-direction
  jl = js, ju = je, kl = ks, ku = ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if (f2) {
      if (!f3) // 2D
        jl = js - 1, ju = je + 1, kl = ks, ku = ke;
      else // 3D
        jl = js - 1, ju = je + 1, kl = ks - 1, ku = ke + 1;
    }
  }
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
#pragma omp simd
      for (int i = is; i <= ie + 1; ++i) {
        kappaf = 0.5 * (kappa(DiffProcess::iso, k, j, i) +
                        kappa(DiffProcess::iso, k, j, i - 1));
        denf = 0.5 * (prim(IDN, k, j, i) + prim(IDN, k, j, i - 1));
        dTdx = (prim(IPR, k, j, i) / prim(IDN, k, j, i) -
                prim(IPR, k, j, i - 1) / prim(IDN, k, j, i - 1)) /
               pco_->dx1v(i - 1);
        x1flux(k, j, i) -= kappaf * denf * dTdx;
      }
    }
  }

  // j-direction
  il = is, iu = ie, kl = ks, ku = ke;
  if (MAGNETIC_FIELDS_ENABLED) {
    if (!f3) // 2D
      il = is - 1, iu = ie + 1, kl = ks, ku = ke;
    else // 3D
      il = is - 1, iu = ie + 1, kl = ks - 1, ku = ke + 1;
  }
  if (f2) { // 2D or 3D
    AthenaArray<Real> &x2flux = cndflx[X2DIR];
    for (int k = kl; k <= ku; ++k) {
      for (int j = js; j <= je + 1; ++j) {
#pragma omp simd
        for (int i = il; i <= iu; ++i) {
          kappaf = 0.5 * (kappa(DiffProcess::iso, k, j, i) +
                          kappa(DiffProcess::iso, k, j - 1, i));
          denf = 0.5 * (prim(IDN, k, j, i) + prim(IDN, k, j - 1, i));
          dTdy = (prim(IPR, k, j, i) / prim(IDN, k, j, i) -
                  prim(IPR, k, j - 1, i) / prim(IDN, k, j - 1, i)) /
                 pco_->h2v(i) / pco_->dx2v(j - 1);
          x2flux(k, j, i) -= kappaf * denf * dTdy;
        }
      }
    }
  } // zero flux for 1D

  // k-direction
  il = is, iu = ie, jl = js, ju = je;
  if (MAGNETIC_FIELDS_ENABLED) {
    if (f2) // 2D or 3D
      il = is - 1, iu = ie + 1, jl = js - 1, ju = je + 1;
    else // 1D
      il = is - 1, iu = ie + 1;
  }
  if (f3) { // 3D
    AthenaArray<Real> &x3flux = cndflx[X3DIR];
    for (int k = ks; k <= ke + 1; ++k) {
      for (int j = jl; j <= ju; ++j) {
#pragma omp simd
        for (int i = il; i <= iu; ++i) {
          kappaf = 0.5 * (kappa(DiffProcess::iso, k, j, i) +
                          kappa(DiffProcess::iso, k - 1, j, i));
          denf = 0.5 * (prim(IDN, k, j, i) + prim(IDN, k - 1, j, i));
          dTdz = (prim(IPR, k, j, i) / prim(IDN, k, j, i) -
                  prim(IPR, k - 1, j, i) / prim(IDN, k - 1, j, i)) /
                 pco_->dx3v(k - 1) / pco_->h31v(i) / pco_->h32v(j);
          x3flux(k, j, i) -= kappaf * denf * dTdz;
        }
      }
    }
  } // zero flux for 1D/2D
  return;
}

//---------------------------------------------------------------------------------------
//! Calculate anisotropic thermal conduction

void HydroDiffusion::ThermalFluxAniso(const AthenaArray<Real> &prim,
                                      const AthenaArray<Real> &cons,
                                      AthenaArray<Real> *cndflx) {
#if MAGNETIC_FIELDS_ENABLED
  //  Anisotropic flux only implemented (and tested) for 3D
  const bool f3 = pmb_->pmy_mesh->f3;
  if (!f3) {
    std::stringstream msg;
    msg << "Anisotropic thermal conduction only implemented and teste in 3D.";
    ATHENA_ERROR(msg);
  }
  Field *pf = pmb_->pfield;
  const auto &bcc = pf->bcc;
  const auto &b = pf->b;

  AthenaArray<Real> &x1flux = cndflx[X1DIR];
  AthenaArray<Real> &x2flux = cndflx[X2DIR];
  AthenaArray<Real> &x3flux = cndflx[X3DIR];
  int il, iu, jl, ju, kl, ku;
  int is = pmb_->is;
  int js = pmb_->js;
  int ks = pmb_->ks;
  int ie = pmb_->ie;
  int je = pmb_->je;
  int ke = pmb_->ke;
  Real Bx, By, Bz, B02, dTc, dTl, dTr, lim_slope, dTdx, dTdy, dTdz, bDotGradT, denf,
      kappaf;

  /* Compute heat fluxes in 1-direction  --------------------------------------*/
  // i-direction
  // 3D with B field limits (see iso flux above)
  jl = js - 1, ju = je + 1, kl = ks - 1, ku = ke + 1;

  for (int k = kl; k <= ku; k++) {
    for (int j = jl; j <= ju; j++) {
#pragma omp simd
      for (int i = is; i <= ie + 1; i++) {
        /* Monotonized temperature difference dT/dy */
        dTdy =
            limiters::lim4(prim(IPR, k, j + 1, i) / prim(IDN, k, j + 1, i) -
                               prim(IPR, k, j, i) / prim(IDN, k, j, i),
                           prim(IPR, k, j, i) / prim(IDN, k, j, i) -
                               prim(IPR, k, j - 1, i) / prim(IDN, k, j - 1, i),
                           prim(IPR, k, j + 1, i - 1) / prim(IDN, k, j + 1, i - 1) -
                               prim(IPR, k, j, i - 1) / prim(IDN, k, j, i - 1),
                           prim(IPR, k, j, i - 1) / prim(IDN, k, j, i - 1) -
                               prim(IPR, k, j - 1, i - 1) / prim(IDN, k, j - 1, i - 1));
        dTdy /= pco_->dx2v(j);

        /* Monotonized temperature difference dT/dz, 3D problem ONLY */
        dTdz =
            limiters::lim4(prim(IPR, k + 1, j, i) / prim(IDN, k + 1, j, i) -
                               prim(IPR, k, j, i) / prim(IDN, k, j, i),
                           prim(IPR, k, j, i) / prim(IDN, k, j, i) -
                               prim(IPR, k - 1, j, i) / prim(IDN, k - 1, j, i),
                           prim(IPR, k + 1, j, i - 1) / prim(IDN, k + 1, j, i - 1) -
                               prim(IPR, k, j, i - 1) / prim(IDN, k, j, i - 1),
                           prim(IPR, k, j, i - 1) / prim(IDN, k, j, i - 1) -
                               prim(IPR, k - 1, j, i - 1) / prim(IDN, k - 1, j, i - 1));
        dTdz /= pco_->dx3v(k);

        /* Add flux at x1-interface, 2D PROBLEM */

        By = 0.5 * (bcc(IB2, k, j, i - 1) + bcc(IB2, k, j, i));
        Bz = 0.5 * (bcc(IB3, k, j, i - 1) + bcc(IB3, k, j, i));
        B02 = SQR(b.x1f(k, j, i)) + SQR(By) + SQR(Bz);
        B02 = std::max(B02, TINY_NUMBER); /* limit in case B=0 */
        bDotGradT = b.x1f(k, j, i) *
                        (prim(IPR, k, j, i) / prim(IDN, k, j, i) -
                         prim(IPR, k, j, i - 1) / prim(IDN, k, j, i - 1)) /
                        pco_->dx1v(i) +
                    By * dTdy + Bz * dTdz;
        kappaf = 0.5 * (kappa(DiffProcess::aniso, k, j, i) +
                        kappa(DiffProcess::aniso, k, j, i - 1));
        denf = 0.5 * (prim(IDN, k, j, i) + prim(IDN, k, j, i - 1));
        x1flux(k, j, i) -= kappaf * denf * (b.x1f(k, j, i) * bDotGradT) / B02;
      } // i
    }   // j
  }     // k

  /* Compute heat fluxes in 2-direction  --------------------------------------*/
  // 3D with B field limits (see iso flux above)
  il = is - 1, iu = ie + 1, kl = ks - 1, ku = ke + 1;

  for (int k = kl; k <= ku; k++) {
    for (int j = js; j <= je + 1; j++) {
#pragma omp simd
      for (int i = il; i <= iu; i++) {
        /* Monotonized temperature difference dT/dx */
        dTdx =
            limiters::lim4(prim(IPR, k, j, i + 1) / prim(IDN, k, j, i + 1) -
                               prim(IPR, k, j, i) / prim(IDN, k, j, i),
                           prim(IPR, k, j, i) / prim(IDN, k, j, i) -
                               prim(IPR, k, j, i - 1) / prim(IDN, k, j, i - 1),
                           prim(IPR, k, j - 1, i + 1) / prim(IDN, k, j - 1, i + 1) -
                               prim(IPR, k, j - 1, i) / prim(IDN, k, j - 1, i),
                           prim(IPR, k, j - 1, i) / prim(IDN, k, j - 1, i) -
                               prim(IPR, k, j - 1, i - 1) / prim(IDN, k, j - 1, i - 1));
        dTdx /= pco_->dx1v(i);

        /* Monotonized temperature difference dT/dz, 3D problem ONLY */
        dTdz =
            limiters::lim4(prim(IPR, k + 1, j, i) / prim(IDN, k + 1, j, i) -
                               prim(IPR, k, j, i) / prim(IDN, k, j, i),
                           prim(IPR, k, j, i) / prim(IDN, k, j, i) -
                               prim(IPR, k - 1, j, i) / prim(IDN, k - 1, j, i),
                           prim(IPR, k + 1, j - 1, i) / prim(IDN, k + 1, j - 1, i) -
                               prim(IPR, k, j - 1, i) / prim(IDN, k, j - 1, i),
                           prim(IPR, k, j - 1, i) / prim(IDN, k, j - 1, i) -
                               prim(IPR, k - 1, j - 1, i) / prim(IDN, k - 1, j - 1, i));
        dTdz /= pco_->dx3v(k);

        /* Add flux at x2-interface, 3D PROBLEM */

        Bx = 0.5 * (bcc(IB1, k, j - 1, i) + bcc(IB1, k, j, i));
        Bz = 0.5 * (bcc(IB3, k, j - 1, i) + bcc(IB3, k, j, i));
        B02 = SQR(Bx) + SQR(b.x2f(k, j, i)) + SQR(Bz);
        B02 = std::max(B02, TINY_NUMBER); /* limit in case B=0 */
        bDotGradT = b.x2f(k, j, i) *
                        (prim(IPR, k, j, i) / prim(IDN, k, j, i) -
                         prim(IPR, k, j - 1, i) / prim(IDN, k, j - 1, i)) /
                        pco_->dx2v(j) +
                    Bx * dTdx + Bz * dTdz;
        kappaf = 0.5 * (kappa(DiffProcess::aniso, k, j, i) +
                        kappa(DiffProcess::aniso, k, j - 1, i));
        denf = 0.5 * (prim(IDN, k, j, i) + prim(IDN, k, j - 1, i));
        x2flux(k, j, i) -= kappaf * denf * (b.x2f(k, j, i) * bDotGradT) / B02;
      } // i
    }   // j
  }     // k

  /* Compute heat fluxes in 3-direction, 3D problem ONLY  ---------------------*/
  // 3D with B field limits (see iso flux above)
  il = is - 1, iu = ie + 1, jl = js - 1, ju = je + 1;

  for (int k = ks; k <= ke + 1; k++) {
    for (int j = jl; j <= ju; j++) {
#pragma omp simd
      for (int i = il; i <= iu; i++) {
        /* Monotonized temperature difference dT/dx */
        dTdx =
            limiters::lim4(prim(IPR, k, j, i + 1) / prim(IDN, k, j, i + 1) -
                               prim(IPR, k, j, i) / prim(IDN, k, j, i),
                           prim(IPR, k, j, i) / prim(IDN, k, j, i) -
                               prim(IPR, k, j, i - 1) / prim(IDN, k, j, i - 1),
                           prim(IPR, k - 1, j, i + 1) / prim(IDN, k - 1, j, i + 1) -
                               prim(IPR, k - 1, j, i) / prim(IDN, k - 1, j, i),
                           prim(IPR, k - 1, j, i) / prim(IDN, k - 1, j, i) -
                               prim(IPR, k - 1, j, i - 1) / prim(IDN, k - 1, j, i - 1));
        dTdx /= pco_->dx1v(i);

        /* Monotonized temperature difference dT/dy */
        dTdy =
            limiters::lim4(prim(IPR, k, j + 1, i) / prim(IDN, k, j + 1, i) -
                               prim(IPR, k, j, i) / prim(IDN, k, j, i),
                           prim(IPR, k, j, i) / prim(IDN, k, j, i) -
                               prim(IPR, k, j - 1, i) / prim(IDN, k, j - 1, i),
                           prim(IPR, k, j + 1, i) / prim(IDN, k, j + 1, i) -
                               prim(IPR, k - 1, j, i) / prim(IDN, k - 1, j, i),
                           prim(IPR, k - 1, j, i) / prim(IDN, k - 1, j, i) -
                               prim(IPR, k - 1, j - 1, i) / prim(IDN, k - 1, j - 1, i));
        dTdy /= pco_->dx2v(j);

        /* Add flux at x3-interface, 3D PROBLEM */

        Bx = 0.5 * (bcc(IB1, k - 1, j, i) + bcc(IB1, k, j, i));
        By = 0.5 * (bcc(IB2, k - 1, j, i) + bcc(IB2, k, j, i));
        B02 = SQR(Bx) + SQR(By) + SQR(b.x3f(k, j, i));
        B02 = std::max(B02, TINY_NUMBER); /* limit in case B=0 */
        bDotGradT = b.x3f(k, j, i) *
                        (prim(IPR, k, j, i) / prim(IDN, k, j, i) -
                         prim(IPR, k - 1, j, i) / prim(IDN, k - 1, j, i)) /
                        pco_->dx3v(k) +
                    Bx * dTdx + By * dTdy;
        kappaf = 0.5 * (kappa(DiffProcess::aniso, k, j, i) +
                        kappa(DiffProcess::aniso, k - 1, j, i));
        denf = 0.5 * (prim(IDN, k, j, i) + prim(IDN, k - 1, j, i));
        x3flux(k, j, i) -= kappaf * denf * (b.x3f(k, j, i) * bDotGradT) / B02;
      } // i
    }   // j
  }     // k
#endif
}

//----------------------------------------------------------------------------------------
//! constant viscosity

void ConstConduction(HydroDiffusion *phdif, MeshBlock *pmb, const AthenaArray<Real> &prim,
                     const AthenaArray<Real> &bcc, int is, int ie, int js, int je, int ks,
                     int ke) {
  if (phdif->kappa_iso > 0.0) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
#pragma omp simd
        for (int i = is; i <= ie; ++i)
          phdif->kappa(HydroDiffusion::DiffProcess::iso, k, j, i) = phdif->kappa_iso;
      }
    }
  }
  if (phdif->kappa_aniso > 0.0) {
    for (int k = ks; k <= ke; ++k) {
      for (int j = js; j <= je; ++j) {
#pragma omp simd
        for (int i = is; i <= ie; ++i)
          phdif->kappa(HydroDiffusion::DiffProcess::aniso, k, j, i) = phdif->kappa_aniso;
      }
    }
  }
  return;
}
