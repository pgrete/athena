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
//! \file six_ray.cpp
//  \brief implementation of radiation integrators: six-ray
//======================================================================================

//c++ headers
#include <math.h>
#include <stdexcept>  // std::runtime_error()
#include <sstream>    // stringstream
#include <string>     // c_str()
#include <iostream>   // endl

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

RadIntegrator::RadIntegrator(Radiation *prad, ParameterInput *pin)
{
  std::stringstream msg; //error message
  pmy_mb = prad->pmy_block;
  pmy_rad = prad;
  rad_G0_ = pin->GetReal("problem", "G0");
  unit_length_in_cm_ = pin->GetReal("chemistry", "unit_length_in_cm");
  if (pmy_rad->nang != 6) {
    msg << "### FATAL ERROR in RadIntegrator constructor [RadIntegrator]" << std::endl
      << "Six-ray scheme with nang != 6." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  for (int i=0; i<6; i++) {
    pneighbors_[i] = NULL;
  }
#ifdef INCLUDE_CHEMISTRY
  pmy_chemnet = pmy_mb->pspec->pchemnet;
  ncol = pmy_chemnet->n_cols_;
  //allocate array for column density
  int ncells1 = pmy_mb->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmy_mb->block_size.nx2 > 1) ncells2 = pmy_mb->block_size.nx2 + 2*(NGHOST);
  if (pmy_mb->block_size.nx3 > 1) ncells3 = pmy_mb->block_size.nx3 + 2*(NGHOST);
  col.NewAthenaArray(6, ncells3, ncells2, ncells1, ncol);
  col_avg.NewAthenaArray(ncol, ncells3, ncells2, ncells1);
  //TODO:for output
  col_Htot.NewAthenaArray(6, ncells3, ncells2, ncells1);
  col_H2.NewAthenaArray(6, ncells3, ncells2, ncells1);
  col_CO.NewAthenaArray(6, ncells3, ncells2, ncells1);
#endif 
}

RadIntegrator::~RadIntegrator() {
#ifdef INCLUDE_CHEMISTRY
  col.DeleteAthenaArray();
  col_avg.DeleteAthenaArray();
  //TODO:for output
  col_Htot.DeleteAthenaArray();
  col_H2.DeleteAthenaArray();
  col_CO.DeleteAthenaArray();
#endif 
}


#ifdef INCLUDE_CHEMISTRY
//direction: 0:+x, 1:+y, 2:+z, 3:-x, 4:-y, 5:-z
void RadIntegrator::GetColMB(int direction) {
  const int iH2 = pmy_chemnet->iH2_;
  const int iCO = pmy_chemnet->iCO_;
  Real NHtot_cell, NHtot_cell_prev;
  std::stringstream msg; //error message
  int is = pmy_mb->is;
  int js = pmy_mb->js;
  int ks = pmy_mb->ks;
  int ie = pmy_mb->ie;
  int je = pmy_mb->je;
  int ke = pmy_mb->ke;
  if (direction == IXP) {
    //+x
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          NHtot_cell = pmy_mb->phydro->w(IDN, k, j, i) * pmy_mb->pcoord->dx1f(i) 
            * unit_length_in_cm_;
          if (i == is) {
            col(IXP, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2.;
            col(IXP, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i);
            col(IXP, k, j, i, pmy_chemnet->iNCO_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i);
          } else {
            NHtot_cell_prev = pmy_mb->phydro->w(IDN, k, j, i-1)
              * pmy_mb->pcoord->dx1f(i-1) * unit_length_in_cm_;
            col(IXP, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2. +  NHtot_cell_prev/2
              + col(IXP, k, j, i-1, pmy_chemnet->iNHtot_);
            col(IXP, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iH2, k, j, i-1)
              + col(IXP, k, j, i-1, pmy_chemnet->iNH2_);
            col(IXP, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iCO, k, j, i-1)
              + col(IXP, k, j, i-1, pmy_chemnet->iNCO_);
          }
        }
      }
    }
  } else if (direction == IXM) {
    //-x
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=ie; i>=is; --i) {
          NHtot_cell = pmy_mb->phydro->w(IDN, k, j, i) * pmy_mb->pcoord->dx1f(i)
            * unit_length_in_cm_;
          if (i == ie) {
            col(IXM, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2.;
            col(IXM, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i);
            col(IXM, k, j, i, pmy_chemnet->iNCO_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i);
          } else {
            NHtot_cell_prev = pmy_mb->phydro->w(IDN, k, j, i+1)
              * pmy_mb->pcoord->dx1f(i+1) * unit_length_in_cm_;
            col(IXM, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2. + NHtot_cell_prev/2.
              + col(IXM, k, j, i+1, pmy_chemnet->iNHtot_);
            col(IXM, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iH2, k, j, i+1)
              + col(IXM, k, j, i+1, pmy_chemnet->iNH2_);
            col(IXM, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i+1)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iCO, k, j, i+1)
              + col(IXM, k, j, i+1, pmy_chemnet->iNCO_);
          }
        }
      }
    }
  } else if (direction == IYP) {
    //+y
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          NHtot_cell = pmy_mb->phydro->w(IDN, k, j, i) * pmy_mb->pcoord->dx2f(j)
            * unit_length_in_cm_;
          if (j == js) {
            col(IYP, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2.;
            col(IYP, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i);
            col(IYP, k, j, i, pmy_chemnet->iNCO_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i);
          } else {
            NHtot_cell_prev = pmy_mb->phydro->w(IDN, k, j-1, i)
              * pmy_mb->pcoord->dx2f(j-1) * unit_length_in_cm_;
            col(IYP, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2. + NHtot_cell_prev/2.
              + col(IYP, k, j-1, i, pmy_chemnet->iNHtot_);
            col(IYP, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iH2, k, j-1, i)
              + col(IYP, k, j-1, i, pmy_chemnet->iNH2_);
            col(IYP, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iCO, k, j-1, i)
              + col(IYP, k, j-1, i, pmy_chemnet->iNCO_);
          }
        }
      }
    }
  } else if (direction == IYM) {
    //-y
    for (int k=ks; k<=ke; ++k) {
      for (int j=je; j>=js; --j) {
        for (int i=is; i<=ie; ++i) {
          NHtot_cell = pmy_mb->phydro->w(IDN, k, j, i) * pmy_mb->pcoord->dx2f(j)
            * unit_length_in_cm_;
          if (j == je) {
            col(IYM, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2.;
            col(IYM, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i);
            col(IYM, k, j, i, pmy_chemnet->iNCO_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i);
          } else {
            NHtot_cell_prev = pmy_mb->phydro->w(IDN, k, j+1, i)
              * pmy_mb->pcoord->dx2f(j+1) * unit_length_in_cm_;
            col(IYM, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2. + NHtot_cell_prev/2.
              + col(IYM, k, j+1, i, pmy_chemnet->iNHtot_);
            col(IYM, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iH2, k, j+1, i)
              + col(IYM, k, j+1, i, pmy_chemnet->iNH2_);
            col(IYM, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iCO, k, j+1, i)
              + col(IYM, k, j+1, i, pmy_chemnet->iNCO_);
          }
        }
      }
    }
  } else if (direction == IZP) {
    //+z
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          NHtot_cell = pmy_mb->phydro->w(IDN, k, j, i) * pmy_mb->pcoord->dx3f(k)
            * unit_length_in_cm_;
          if (k == ks) {
            col(IZP, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2.;
            col(IZP, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i);
            col(IZP, k, j, i, pmy_chemnet->iNCO_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i);
          } else {
            NHtot_cell_prev = pmy_mb->phydro->w(IDN, k-1, j, i) 
              * pmy_mb->pcoord->dx3f(k-1) * unit_length_in_cm_;
            col(IZP, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2. + NHtot_cell_prev/2.
              + col(IZP, k-1, j, i, pmy_chemnet->iNHtot_);
            col(IZP, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iH2, k-1, j, i)
              + col(IZP, k-1, j, i, pmy_chemnet->iNH2_);
            col(IZP, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iCO, k-1, j, i)
              + col(IZP, k-1, j, i, pmy_chemnet->iNCO_);
          }
        }
      }
    }
  } else if (direction == IZM) {
    //-z
    for (int k=ke; k>=ks; --k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          NHtot_cell = pmy_mb->phydro->w(IDN, k, j, i) * pmy_mb->pcoord->dx3f(k)
            * unit_length_in_cm_;
          if (k == ke) {
            col(IZM, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2.;
            col(IZM, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i);
            col(IZM, k, j, i, pmy_chemnet->iNCO_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i);
          } else {
            NHtot_cell_prev = pmy_mb->phydro->w(IDN, k+1, j, i) 
              * pmy_mb->pcoord->dx3f(k+1) * unit_length_in_cm_;
            col(IZM, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2. + NHtot_cell_prev/2.
              + col(IZM, k+1, j, i, pmy_chemnet->iNHtot_);
            col(IZM, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iH2, k+1, j, i)
              + col(IZM, k+1, j, i, pmy_chemnet->iNH2_);
            col(IZM, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iCO, k+1, j, i)
              + col(IZM, k+1, j, i, pmy_chemnet->iNCO_);
          }
        }
      }
    }
  } else {
    msg << "### FATAL ERROR in RadIntegrator six_ray [GetColMB]" << std::endl
      << "direction {0,1,2,3,4,5}:" << direction << " unknown." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  return;
}

void RadIntegrator::UpdateRadiation(int direction) {
  const int iang = direction;
  const Real Zd = pmy_chemnet->zdg_;
  const Real bH2 = 3.0e5; //H2 velocity dispersion
  const int iph_H2 = ChemNetwork::iph_H2_;
  const int iph_CO = ChemNetwork::iph_CO_;
  const int iPE = pmy_chemnet->index_gpe_;
  const Real sigmaPE = Thermo::sigmaPE_;
  Real NH, AV, NCO, NH2, fs_CO, fs_H2;
  for (int k=pmy_mb->ks; k<=pmy_mb->ke; ++k) {
    for (int j=pmy_mb->js; j<=pmy_mb->je; ++j) {
      for (int i=pmy_mb->is; i<=pmy_mb->ie; ++i) {
        NH = col(direction, k, j, i, pmy_chemnet->iNHtot_);
        NH2 = col(direction, k, j, i,  pmy_chemnet->iNH2_);
        NCO = col(direction, k, j, i,  pmy_chemnet->iNCO_);
        AV = NH * Zd / 1.87e21;
        //photo-reactions
        for (int ifreq=0; ifreq < pmy_rad->nfreq-1; ++ifreq) {
          pmy_rad->ir(k, j, i, ifreq * pmy_rad->nang+iang) = rad_G0_
            * exp( -ChemNetwork::kph_avfac_[ifreq] * AV );
        }
        //H2 and CO self-sheilding
        fs_CO = Shielding::fShield_CO_V09(NCO, NH2);
        fs_H2 = Shielding::fShield_H2(NH2, bH2);
        pmy_rad->ir(k, j, i, iph_H2 * pmy_rad->nang+iang) *= 
          fs_H2;
        pmy_rad->ir(k, j, i, iph_CO * pmy_rad->nang+iang) *= 
          fs_CO;
        //GPE 
        pmy_rad->ir(k, j, i, iPE * pmy_rad->nang+iang) = rad_G0_
          *  exp(-NH * sigmaPE * Zd);
      }
    }
  }
  return;
}

void RadIntegrator::CopyToOutput() {
  int is = pmy_mb->is;
  int js = pmy_mb->js;
  int ks = pmy_mb->ks;
  int ie = pmy_mb->ie;
  int je = pmy_mb->je;
  int ke = pmy_mb->ke;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        for (int iang=0; iang < 6; iang++) {
          //column densities
          col_Htot(iang, k, j, i) = col(iang, k, j, i, pmy_chemnet->iNHtot_); 
          col_H2(iang, k, j, i) = col(iang, k, j, i, pmy_chemnet->iNH2_); 
          col_CO(iang, k, j, i) = col(iang, k, j, i, pmy_chemnet->iNCO_); 
          //angel averaged column densities
          for (int icol=0; icol<ncol; icol++) {
            if (iang == 0) {
              col_avg(icol, k, j, i) = 0;
            }
            col_avg(icol, k, j, i) += col(iang, k, j, i, icol)/6.;
          }
          //radiation
          for (int ifreq=0; ifreq < pmy_rad->nfreq; ++ifreq) {
            if (iang == 0) {
              pmy_rad->ir_avg(ifreq, k, j, i) = 0;
            }
            pmy_rad->ir_avg(ifreq, k, j, i) += 
              pmy_rad->ir(k, j, i, ifreq * pmy_rad->nang+iang)/6.;
          }
        }
      }
    }
  }
}

void RadIntegrator::SetSixRayNeighbors() {
  //assign neighbors for six-ray
  NeighborBlock* nb;
  for(int n=0; n<pmy_mb->nneighbor; n++) {
    nb = &pmy_mb->neighbor[n];
    if (nb->fid == INNER_X1) {
      pneighbors_[0] = nb;
      std::cout << "INNER_X1" << std::endl;
    } else if (nb->fid == INNER_X2) {
      pneighbors_[1] = nb;
      std::cout << "INNER_X2" << std::endl;
    } else if (nb->fid == INNER_X3) {
      pneighbors_[2] = nb;
      std::cout << "INNER_X3" << std::endl;
    } else if (nb->fid == OUTER_X1) {
      pneighbors_[3] = nb;
      std::cout << "OUTER_X1" << std::endl;
    } else if (nb->fid == OUTER_X2) {
      pneighbors_[4] = nb;
      std::cout << "OUTER_X2" << std::endl;
    } else if (nb->fid == OUTER_X3) {
      pneighbors_[5] = nb;
      std::cout << "OUTER_X3" << std::endl;
    }
  }
  return;
}
#endif
