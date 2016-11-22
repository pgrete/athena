//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================
//! \file loc_jeans.cpp
//  \brief implementation of radiation integrators: local jeans shielding
//======================================================================================

//c++ headers
#include <math.h>

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"
#include "../../utils/cgk_utils.hpp"
#include "../../hydro/hydro.hpp"
#include "../../eos/eos.hpp"
#include "../radiation.hpp"
#ifdef INCLUDE_CHEMISTRY
#include "../../chemistry/species.hpp"
#include "../../chemistry/shielding.hpp"
#include "../../chemistry/thermo.hpp"
#endif

// Class header
#include "rad_integrators.hpp"

static Real GetLJ(Real cs, Real dens) {
  const Real G = 6.67259e-8;//gravitational constant in cgs
  return cs * sqrt( PI / (G * dens)  );
}

RadIntegrator::RadIntegrator(Radiation *prad, ParameterInput *pin)
{
  pmy_rad = prad;
  rad_G0_ = pin->GetReal("problem", "G0");
}

RadIntegrator::~RadIntegrator() {}


#ifdef INCLUDE_CHEMISTRY
void RadIntegrator::UpdateRadiation() {
  const Real Tceiling = 40.; //temperature ceiling in Kelvin
  const Real Tfloor = 1.; //temperature floor in Kelvin
  const Real bH2 = 3.0e5; //H2 velocity dispersion
  const Real sigmaPE = Thermo::sigmaPE_;
  MeshBlock *pb = pmy_rad->pmy_block;
  const Real Zd = pb->pspec->pchemnet->zdg_;
  const int iH2 = pb->pspec->pchemnet->iH2_;
  const int iCO = pb->pspec->pchemnet->iCO_;
  const int iph_H2 = ChemNetwork::iph_H2_;
  const int iph_CO = ChemNetwork::iph_CO_;
  const int iPE = pb->pspec->pchemnet->index_gpe_;
  Real press, dens, dens_cgs, temp, cs, cs_sq;
  Real gm = pb->peos->GetGamma();
  Real Lshield; //in cm
  Real NH, NH2, NCO, AV, yH2, ye;
  Real xHe = pb->pspec->pchemnet->xHe_;
  Real dens_small_ = 1e-20;
  Real mu_i;


  //loop over each cell
  for (int k=pb->ks; k<=pb->ke; ++k) {
    for (int j=pb->js; j<=pb->je; ++j) {
      for (int i=pb->is; i<=pb->ie; ++i) {
        //calculate the shielding length
        press = pb->phydro->w(IEN, k, j, i);
        dens = pb->phydro->w(IDN, k, j, i);
        dens_cgs = dens *  CGKUtility::unitD;
        yH2 = pb->pspec->s(pb->pspec->pchemnet->iH2_, k, j, i);
        ye = 0.;
        temp = pb->pspec->s(pb->pspec->pchemnet->iE_, k, j, i) / 
          Thermo::CvCold(yH2, xHe, ye);
        //set T ceilling and floor
        if (temp > Tceiling) {
          temp = Tceiling;
        } else if (temp < Tfloor) {
          temp = Tfloor;
        }
        if (dens < dens_small_) {
          Lshield = 0.;
        } else {
          mu_i = 1.5 * CGKUtility::kB * (1. + 4.*xHe) / 
            Thermo::CvCold(yH2, xHe, ye);
          cs_sq = ( gm * temp * CGKUtility::kB ) / ( mu_i * CGKUtility::mH ); //cs in cgs
          cs = sqrt(cs_sq);
          Lshield = GetLJ(cs, dens_cgs);
        }
        //calculate NH, NH2, NCO
        if (dens < 0) {
          NH = 0;
        } else {
          NH = dens * Lshield;
        }
        AV = NH * Zd / 1.87e21;
        NH2 = pb->pspec->s(iH2, k, j, i) * NH;
        NCO = pb->pspec->s(iCO, k, j, i) * NH;
        //set CO column for calculating CO cooling
        pb->pspec->pchemnet->colCO_(k, j, i) = NCO;
        //dust shielding
        //photo-reactions
        for (int ifreq=0; ifreq < pmy_rad->nfreq-2; ++ifreq) {
          pmy_rad->ir(k, j, i, ifreq * pmy_rad->nang) = rad_G0_
            * exp( -ChemNetwork::kph_avfac_[ifreq] * AV );
        }
        //H2 and CO self-sheilding
        Real fs_CO = Shielding::fShield_CO_V09(NCO, NH2);
        Real fs_H2 = Shielding::fShield_H2(NH2, bH2);
        pmy_rad->ir(k, j, i, iph_H2 * pmy_rad->nang) *= 
          fs_H2;
        pmy_rad->ir(k, j, i, iph_CO * pmy_rad->nang) *= 
          fs_CO;
        //debug
#ifdef DEBUG
        if (dens > 5.) {
          printf("RadIntegrator: nH=%.2e, temp=%.2e, cs=%.2e, Lshield=%.2e, NCO=%.2e, NH2=%.2e",
                 dens, temp, cs, Lshield, NCO, NH2);
          printf(", mui=%.2e, yH2=%.2e, ye=%.2e\n",mu_i, yH2, ye);
        }
        if (isnan(fs_CO) ) {
          printf("RadIntegrator::UpdateRadJeans(): fs_CO=nan, Lshield=%.2e\n", Lshield);
          printf("press=%.2e, dens=%.2e, temp=%.2e, cs=%.2e, CO=%.2e, NCO=%.2e, NH2=%.2e\n",
                  press, dens, temp, cs, pb->pspec->s(iCO, k, j, i), NCO, NH2);
        }
#endif

        //GPE 
        pmy_rad->ir(k, j, i, iPE * pmy_rad->nang) = rad_G0_
          *  exp(-NH * sigmaPE * Zd);

        //copy to other angles
        for (int ifreq=0; ifreq < pmy_rad->nfreq; ++ifreq) {
          for (int iang=1; iang < pmy_rad->nang; ++iang) {
            pmy_rad->ir(k, j, i, ifreq * pmy_rad->nang + iang) = 
              pmy_rad->ir(k, j, i, ifreq * pmy_rad->nang);
          }
        }

      }
    }
  }
}
#endif

