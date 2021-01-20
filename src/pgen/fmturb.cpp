//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file turb.cpp
//  \brief Problem generator for turbulence generator
//

// C++ headers
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
#include "../fft/few_modes_turbulence.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

static Real hst_turbulence(MeshBlock *pmb, int iout) {
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;

  Real gam, c_s;
  if (NON_BAROTROPIC_EOS) {
    gam = pmb->peos->GetGamma();
  } else {
    c_s =  pmb->peos->GetIsoSoundSpeed();
  }

  AthenaArray<Real> vol(pmb->ncells1);
  Real sum = 0;

  for (int k = ks; k <= ke; k++) {
    for (int j = js; j <= je; j++) {
      pmb->pcoord->CellVolume(k, j, is, ie, vol);
      for (int i = is; i <= ie; i++) {
        Real vel2 = (pmb->phydro->w(IVX, k, j, i) * pmb->phydro->w(IVX, k, j, i) +
                     pmb->phydro->w(IVY, k, j, i) * pmb->phydro->w(IVY, k, j, i) +
                     pmb->phydro->w(IVZ, k, j, i) * pmb->phydro->w(IVZ, k, j, i));

        if (NON_BAROTROPIC_EOS) {
          c_s = std::sqrt(gam * pmb->phydro->w(IPR, k, j, i) /
                          pmb->phydro->w(IDN, k, j, i)); // speed of sound
        }

        Real e_kin = 0.5 * pmb->phydro->w(IDN, k, j, i) * vel2;

        if (iout == 0) { // Ms
          sum += vol(i) * std::sqrt(vel2) / c_s;
        } else if (iout == 2) { // Ekin
          sum += vol(i) * e_kin;
        }

        if (MAGNETIC_FIELDS_ENABLED) {
          Real B2 = (pmb->pfield->bcc(IB1, k, j, i) * pmb->pfield->bcc(IB1, k, j, i) +
                     pmb->pfield->bcc(IB2, k, j, i) * pmb->pfield->bcc(IB2, k, j, i) +
                     pmb->pfield->bcc(IB3, k, j, i) * pmb->pfield->bcc(IB3, k, j, i));

          Real e_mag = 0.5 * B2;

          if (iout == 1) { // Ma
            sum += vol(i) * std::sqrt(e_kin / e_mag);
          } else if (iout == 3) { // Emag
            sum += vol(i) * e_mag;
          } else if (iout == 4) { // plasma beta
            if (NON_BAROTROPIC_EOS) {
              sum += vol(i) * pmb->phydro->w(IPR, k, j, i) / e_mag;
            } else {
              sum += vol(i) * pmb->phydro->w(IDN, k, j, i) * c_s * c_s / e_mag;
            }
          }
        }
      }
    }
  }
  return sum;
}

namespace cooling {

Real T_low;    // lower temperature threshold for which cooling is active
Real c_cool;   // arbitrary cooling coeff
Real max_diff; // max. fraction of change of internal energy density per step

void Cooling(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {
  const auto coeff = -dt * c_cool / std::sqrt(pmb->peos->GetGamma() - 1.0);
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      #pragma omp simd
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        const auto T = prim(IPR,k,j,i) / prim(IDN,k,j,i); // same "units" as T_low
        if (T > T_low) {
          cons(IEN,k,j,i) += coeff * prim(IDN,k,j,i) * prim(IDN,k,j,i) * std::sqrt(T);
        }
      }
    }
  }
}

Real CoolDt(MeshBlock *pmb) {
  const AthenaArray<Real> &prim = pmb->phydro->w;
  const auto gm1 = pmb->peos->GetGamma() - 1.0;
  const auto coeff = c_cool / std::sqrt(gm1);
  auto min_dt = std::numeric_limits<Real>::max();
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        const auto T = prim(IPR,k,j,i) / prim(IDN,k,j,i); // same "units" as T_low
        if (T > T_low) {
          const auto Edot = coeff * prim(IDN,k,j,i) * prim(IDN,k,j,i) * std::sqrt(T);
          const auto E = prim(IPR,k,j,i) / gm1;
          const Real dt = max_diff * E / Edot;
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
  fmturb_flag = pin->GetInteger("problem", "fmturb_flag");
  AllocateUserHistoryOutput(5);
  EnrollUserHistoryOutput(0, hst_turbulence, "Mean Ms");
  EnrollUserHistoryOutput(1, hst_turbulence, "Mean Ma");
  EnrollUserHistoryOutput(2, hst_turbulence, "Mean Ekin");
  EnrollUserHistoryOutput(3, hst_turbulence, "Mean Emag");
  EnrollUserHistoryOutput(4, hst_turbulence, "Mean p_b");

  auto cool_flag = pin->GetOrAddInteger("problem", "cool_flag", 0);
  if (cool_flag > 0) {
    const auto p0 = pin->GetReal("problem", "p0");
    const auto rho0 = pin->GetReal("problem", "rho0");
    const auto thres = pin->GetOrAddReal("problem", "cool_threshold_ratio", 0.0);
    // this is a proxy for the temp, which we can use given that gamma and c_V are const
    const auto T0 = p0/rho0;
    cooling::T_low = T0 * thres;
    cooling::c_cool = pin->GetOrAddReal("problem", "cool_coeff", 0.0);
    cooling::max_diff = pin->GetOrAddReal("problem", "cool_max_diff", 0.1);

    EnrollUserExplicitSourceFunction(cooling::Cooling);
    EnrollUserTimeStepFunction(cooling::CoolDt);
  }
}

//========================================================================================
//! \fn void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in MeshBlock class.  Can also be
//  used to initialize variables which are global to other functions in this file.
//  Called in MeshBlock constructor before ProblemGenerator.
//========================================================================================

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(3);
  SetUserOutputVariableName(0, "acceleration_x");
  SetUserOutputVariableName(1, "acceleration_y");
  SetUserOutputVariableName(2, "acceleration_z");

  auto num_modes = pin->GetInteger("problem", "num_modes"); // number of wavemodes
  AllocateRealUserMeshBlockDataField(1);
  ruser_meshblock_data[0].NewAthenaArray(3, num_modes, 2);
  return;
}

//========================================================================================
//! \fn void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin)
//  \brief Function called before generating output files
//========================================================================================

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  pmy_mesh->pfmtrbd->CopyAccelToOutputVars(user_out_var);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real gamma = peos->GetGamma();
  Real gm1 = gamma - 1.0;
  Real p0 = pin->GetReal("problem", "p0");
  Real rho0 = pin->GetReal("problem", "rho0");
  Real x3min = pmy_mesh->mesh_size.x3min;
  Real Lz = pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min;
  Real kz = 2.0 * PI / Lz;

  for (int k = ks; k <= ke; k++) {
    for (int j = js; j <= je; j++) {
      for (int i = is; i <= ie; i++) {
        phydro->u(IDN, k, j, i) = rho0;

        phydro->u(IM1, k, j, i) = 0.0;
        phydro->u(IM2, k, j, i) = 0.0;
        phydro->u(IM3, k, j, i) = 0.0;

        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN, k, j, i) = p0 / gm1;
        }
      }
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {
    int nx1 = (ie - is) + 1 + 2 * (NGHOST);
    int nx2 = (je - js) + 1 + 2 * (NGHOST);
    int nx3 = (ke - ks) + 1 + 2 * (NGHOST);

    AthenaArray<Real> ax, ay, az;
    ax.NewAthenaArray(nx3, nx2, nx1);
    ay.NewAthenaArray(nx3, nx2, nx1);
    az.NewAthenaArray(nx3, nx2, nx1);
    for (int k = ks; k <= ke + 1; k++) {
      for (int j = js; j <= je + 1; j++) {
#pragma omp simd
        for (int i = is; i <= ie + 1; i++) {
          ax(k, j, i) = 0.0;
          ay(k, j, i) = 0.0;
          az(k, j, i) = 0.0;
        }
      }
    }

    Real b0 = pin->GetReal("problem", "b0");
    auto b_config = pin->GetInteger("problem", "b_config");

    uint32_t seed_val = gid; // ensure different rng for each MB
    std::mt19937 rng;
    rng.seed(seed_val);
    std::normal_distribution<Real> normal_dist(0.0, b0);

    if (b_config == 3) { // random field
      // initial vector potential
      for (int k = ks; k <= ke + 1; k++) {
        for (int j = js; j <= je + 1; j++) {
          for (int i = is; i <= ie + 1; i++) {
            ax(k, j, i) = normal_dist(rng);
            ay(k, j, i) = normal_dist(rng);
            az(k, j, i) = normal_dist(rng);
          }
        }
      }
    }

    if (b_config == 4) { // field loop
      // the origin of the initial loop
      Real x0 = pin->GetOrAddReal("problem", "x0", 0.5);
      Real y0 = pin->GetOrAddReal("problem", "y0", 0.5);
      Real z0 = pin->GetOrAddReal("problem", "z0", 0.5);
      Real rad = pin->GetOrAddReal("problem", "loop_rad", 0.25);

      for (int k = ks; k <= ke + 1; k++) {
        for (int j = js; j <= je + 1; j++) {
          for (int i = is; i <= ie + 1; i++) {
            if ((SQR(pcoord->x1f(i) - x0) + SQR(pcoord->x2f(j) - y0)) < rad * rad) {
              az(k, j, i) =
                  (rad - std::sqrt(SQR(pcoord->x1f(i) - x0) + SQR(pcoord->x2f(j) - y0)));
            }
          }
        }
      }
    }

    Real local_mag_en = 0.0; // used for normalization

    for (int k = ks; k <= ke; k++) {
      for (int j = js; j <= je; j++) {
        for (int i = is; i <= ie + 1; i++) {
          pfield->b.x1f(k, j, i) = 0.0;

          if (b_config == 0) { // uniform field
            pfield->b.x1f(k, j, i) = b0;
          }
          if (b_config == 1) { // no net flux with uniform fieldi
            if (pcoord->x3v(k) < x3min + Lz / 2.0) {
              pfield->b.x1f(k, j, i) = b0;
            } else {
              pfield->b.x1f(k, j, i) = -b0;
            }
          }
          if (b_config == 2) { // no net flux with sin(z) shape
            // sqrt(0.5) is used so that resulting e_mag is approx b_0^2/2 similar to
            // other b_configs
            pfield->b.x1f(k, j, i) = b0 / std::sqrt(0.5) * sin(kz * pcoord->x3v(k));
          }

          pfield->b.x1f(k, j, i) += (az(k, j + 1, i) - az(k, j, i)) / pcoord->dx2f(j) -
                                    (ay(k + 1, j, i) - ay(k, j, i)) / pcoord->dx3f(k);
          if (i > is) {
            local_mag_en +=
                0.5 * SQR(0.5 * (pfield->b.x1f(k, j, i - 1) + pfield->b.x1f(k, j, i)));
          }
        }
      }
    }
    for (int k = ks; k <= ke; k++) {
      for (int j = js; j <= je + 1; j++) {
        for (int i = is; i <= ie; i++) {
          pfield->b.x2f(k, j, i) += (ax(k + 1, j, i) - ax(k, j, i)) / pcoord->dx3f(k) -
                                    (az(k, j, i + 1) - az(k, j, i)) / pcoord->dx1f(i);
          if (j > js) {
            local_mag_en +=
                0.5 * SQR(0.5 * (pfield->b.x2f(k, j - 1, i) + pfield->b.x2f(k, j, i)));
          }
        }
      }
    }
    for (int k = ks; k <= ke + 1; k++) {
      for (int j = js; j <= je; j++) {
        for (int i = is; i <= ie; i++) {
          pfield->b.x3f(k, j, i) += (ay(k, j, i + 1) - ay(k, j, i)) / pcoord->dx1f(i) -
                                    (ax(k, j + 1, i) - ax(k, j, i)) / pcoord->dx2f(j);
          if (k > ks) {
            local_mag_en +=
                0.5 * SQR(0.5 * (pfield->b.x3f(k - 1, j, i) + pfield->b.x3f(k, j, i)));
          }
        }
      }
    }

    Real global_mag_en = 0.0;
#ifdef MPI_PARALLEL
    int mpierr;
    mpierr = MPI_Allreduce(&local_mag_en, &global_mag_en, 1, MPI_DOUBLE, MPI_SUM,
                           MPI_COMM_WORLD);
    if (mpierr) {
      std::stringstream msg;
      msg << "[normalize]: MPI_Allreduce error = " << mpierr << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
#else
    global_mag_en = local_mag_en;
#endif // MPI_PARALLEL

    Real gsize = (static_cast<Real>(pmy_mesh->mesh_size.nx1) *
                  static_cast<Real>(pmy_mesh->mesh_size.nx2) *
                  static_cast<Real>(pmy_mesh->mesh_size.nx3));
    auto b_norm = std::sqrt(global_mag_en / gsize / (0.5 * b0 * b0));
    if (Globals::my_rank == 0) {
      std::cout << "Applying norm factor of " << b_norm << " to B field."
                << " Orig mean E_mag = " << (global_mag_en / gsize) << std::endl;
    }
    for (int k = ks; k <= ke; k++) {
      for (int j = js; j <= je; j++) {
        for (int i = is; i <= ie + 1; i++) {
          pfield->b.x1f(k, j, i) /= b_norm;
          // setting cell centered energy after the right hand side face has been set
          if (NON_BAROTROPIC_EOS && i > is) {
            phydro->u(IEN, k, j, i - 1) +=
                0.5 * SQR(0.5 * (pfield->b.x1f(k, j, i - 1) + pfield->b.x1f(k, j, i)));
          }
        }
      }
    }
    for (int k = ks; k <= ke; k++) {
      for (int j = js; j <= je + 1; j++) {
        for (int i = is; i <= ie; i++) {
          pfield->b.x2f(k, j, i) /= b_norm;
          // setting cell centered energy after the right hand side face has been set
          if (NON_BAROTROPIC_EOS && j > js) {
            phydro->u(IEN, k, j - 1, i) +=
                0.5 * SQR(0.5 * (pfield->b.x2f(k, j - 1, i) + pfield->b.x2f(k, j, i)));
          }
        }
      }
    }
    for (int k = ks; k <= ke + 1; k++) {
      for (int j = js; j <= je; j++) {
        for (int i = is; i <= ie; i++) {
          pfield->b.x3f(k, j, i) /= b_norm;
          // setting cell centered energy after the right hand side face has been set
          if (NON_BAROTROPIC_EOS && k > ks) {
            phydro->u(IEN, k - 1, j, i) +=
                0.5 * SQR(0.5 * (pfield->b.x3f(k - 1, j, i) + pfield->b.x3f(k, j, i)));
          }
        }
      }
    }
  }
}

//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {}
