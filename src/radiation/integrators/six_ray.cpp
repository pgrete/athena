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
    pfacenb_[i] = NULL;
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
  if (direction == INNER_X1) {
    //+x
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          NHtot_cell = pmy_mb->phydro->w(IDN, k, j, i) * pmy_mb->pcoord->dx1f(i) 
            * unit_length_in_cm_;
          if (i == is) {
            col(INNER_X1, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2.;
            col(INNER_X1, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i);
            col(INNER_X1, k, j, i, pmy_chemnet->iNCO_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i);
          } else {
            NHtot_cell_prev = pmy_mb->phydro->w(IDN, k, j, i-1)
              * pmy_mb->pcoord->dx1f(i-1) * unit_length_in_cm_;
            col(INNER_X1, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2. +  NHtot_cell_prev/2
              + col(INNER_X1, k, j, i-1, pmy_chemnet->iNHtot_);
            col(INNER_X1, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iH2, k, j, i-1)
              + col(INNER_X1, k, j, i-1, pmy_chemnet->iNH2_);
            col(INNER_X1, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iCO, k, j, i-1)
              + col(INNER_X1, k, j, i-1, pmy_chemnet->iNCO_);
          }
        }
      }
    }
  } else if (direction == OUTER_X1) {
    //-x
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=ie; i>=is; --i) {
          NHtot_cell = pmy_mb->phydro->w(IDN, k, j, i) * pmy_mb->pcoord->dx1f(i)
            * unit_length_in_cm_;
          if (i == ie) {
            col(OUTER_X1, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2.;
            col(OUTER_X1, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i);
            col(OUTER_X1, k, j, i, pmy_chemnet->iNCO_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i);
          } else {
            NHtot_cell_prev = pmy_mb->phydro->w(IDN, k, j, i+1)
              * pmy_mb->pcoord->dx1f(i+1) * unit_length_in_cm_;
            col(OUTER_X1, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2. + NHtot_cell_prev/2.
              + col(OUTER_X1, k, j, i+1, pmy_chemnet->iNHtot_);
            col(OUTER_X1, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iH2, k, j, i+1)
              + col(OUTER_X1, k, j, i+1, pmy_chemnet->iNH2_);
            col(OUTER_X1, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iCO, k, j, i+1)
              + col(OUTER_X1, k, j, i+1, pmy_chemnet->iNCO_);
          }
        }
      }
    }
  } else if (direction == INNER_X2) {
    //+y
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          NHtot_cell = pmy_mb->phydro->w(IDN, k, j, i) * pmy_mb->pcoord->dx2f(j)
            * unit_length_in_cm_;
          if (j == js) {
            col(INNER_X2, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2.;
            col(INNER_X2, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i);
            col(INNER_X2, k, j, i, pmy_chemnet->iNCO_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i);
          } else {
            NHtot_cell_prev = pmy_mb->phydro->w(IDN, k, j-1, i)
              * pmy_mb->pcoord->dx2f(j-1) * unit_length_in_cm_;
            col(INNER_X2, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2. + NHtot_cell_prev/2.
              + col(INNER_X2, k, j-1, i, pmy_chemnet->iNHtot_);
            col(INNER_X2, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iH2, k, j-1, i)
              + col(INNER_X2, k, j-1, i, pmy_chemnet->iNH2_);
            col(INNER_X2, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iCO, k, j-1, i)
              + col(INNER_X2, k, j-1, i, pmy_chemnet->iNCO_);
          }
        }
      }
    }
  } else if (direction == OUTER_X2) {
    //-y
    for (int k=ks; k<=ke; ++k) {
      for (int j=je; j>=js; --j) {
        for (int i=is; i<=ie; ++i) {
          NHtot_cell = pmy_mb->phydro->w(IDN, k, j, i) * pmy_mb->pcoord->dx2f(j)
            * unit_length_in_cm_;
          if (j == je) {
            col(OUTER_X2, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2.;
            col(OUTER_X2, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i);
            col(OUTER_X2, k, j, i, pmy_chemnet->iNCO_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i);
          } else {
            NHtot_cell_prev = pmy_mb->phydro->w(IDN, k, j+1, i)
              * pmy_mb->pcoord->dx2f(j+1) * unit_length_in_cm_;
            col(OUTER_X2, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2. + NHtot_cell_prev/2.
              + col(OUTER_X2, k, j+1, i, pmy_chemnet->iNHtot_);
            col(OUTER_X2, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iH2, k, j+1, i)
              + col(OUTER_X2, k, j+1, i, pmy_chemnet->iNH2_);
            col(OUTER_X2, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iCO, k, j+1, i)
              + col(OUTER_X2, k, j+1, i, pmy_chemnet->iNCO_);
          }
        }
      }
    }
  } else if (direction == INNER_X3) {
    //+z
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          NHtot_cell = pmy_mb->phydro->w(IDN, k, j, i) * pmy_mb->pcoord->dx3f(k)
            * unit_length_in_cm_;
          if (k == ks) {
            col(INNER_X3, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2.;
            col(INNER_X3, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i);
            col(INNER_X3, k, j, i, pmy_chemnet->iNCO_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i);
          } else {
            NHtot_cell_prev = pmy_mb->phydro->w(IDN, k-1, j, i) 
              * pmy_mb->pcoord->dx3f(k-1) * unit_length_in_cm_;
            col(INNER_X3, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2. + NHtot_cell_prev/2.
              + col(INNER_X3, k-1, j, i, pmy_chemnet->iNHtot_);
            col(INNER_X3, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iH2, k-1, j, i)
              + col(INNER_X3, k-1, j, i, pmy_chemnet->iNH2_);
            col(INNER_X3, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iCO, k-1, j, i)
              + col(INNER_X3, k-1, j, i, pmy_chemnet->iNCO_);
          }
        }
      }
    }
  } else if (direction == OUTER_X3) {
    //-z
    for (int k=ke; k>=ks; --k) {
      for (int j=js; j<=je; ++j) {
        for (int i=is; i<=ie; ++i) {
          NHtot_cell = pmy_mb->phydro->w(IDN, k, j, i) * pmy_mb->pcoord->dx3f(k)
            * unit_length_in_cm_;
          if (k == ke) {
            col(OUTER_X3, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2.;
            col(OUTER_X3, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i);
            col(OUTER_X3, k, j, i, pmy_chemnet->iNCO_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i);
          } else {
            NHtot_cell_prev = pmy_mb->phydro->w(IDN, k+1, j, i) 
              * pmy_mb->pcoord->dx3f(k+1) * unit_length_in_cm_;
            col(OUTER_X3, k, j, i, pmy_chemnet->iNHtot_) = NHtot_cell/2. + NHtot_cell_prev/2.
              + col(OUTER_X3, k+1, j, i, pmy_chemnet->iNHtot_);
            col(OUTER_X3, k, j, i, pmy_chemnet->iNH2_) = 
              NHtot_cell/2. * pmy_mb->pspec->s(iH2, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iH2, k+1, j, i)
              + col(OUTER_X3, k+1, j, i, pmy_chemnet->iNH2_);
            col(OUTER_X3, k, j, i, pmy_chemnet->iNCO_) =
              NHtot_cell/2. * pmy_mb->pspec->s(iCO, k, j, i)
              + NHtot_cell_prev/2. * pmy_mb->pspec->s(iCO, k+1, j, i)
              + col(OUTER_X3, k+1, j, i, pmy_chemnet->iNCO_);
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

void RadIntegrator::UpdateCol(int direction) {
  const int iH2 = pmy_chemnet->iH2_;
  const int iCO = pmy_chemnet->iCO_;
  int is = pmy_mb->is;
  int js = pmy_mb->js;
  int ks = pmy_mb->ks;
  int ie = pmy_mb->ie;
  int je = pmy_mb->je;
  int ke = pmy_mb->ke;
  std::stringstream msg; //error message
  Real NH_ghostzone;
  Real NH_boundary, NH2_boundary, NCO_boundary;
  if (direction == INNER_X1) {
    //+x
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        NH_ghostzone = unit_length_in_cm_ * 
          pmy_mb->phydro->w(IDN, k, j, is-1) * pmy_mb->pcoord->dx1f(is-1) / 2.;
        NH_boundary = col(direction, k, j, is-1, pmy_chemnet->iNHtot_) + NH_ghostzone; 
        NH2_boundary = col(direction, k, j, is-1, pmy_chemnet->iNH2_)
          + pmy_mb->pspec->s(iH2, k, j, is-1) * NH_ghostzone; 
        NCO_boundary = col(direction, k, j, is-1, pmy_chemnet->iNCO_)
          + pmy_mb->pspec->s(iCO, k, j, is-1) * NH_ghostzone;
        for (int i=is; i<=ie; ++i) {
          col(direction, k, j, i, pmy_chemnet->iNHtot_) += NH_boundary;
          col(direction, k, j, i, pmy_chemnet->iNH2_) += NH2_boundary;
          col(direction, k, j, i, pmy_chemnet->iNCO_) += NCO_boundary;
        }
      }
    }
  } else if (direction == OUTER_X1) {
    //-x
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        NH_ghostzone = unit_length_in_cm_ *
          pmy_mb->phydro->w(IDN, k, j, ie+1) * pmy_mb->pcoord->dx1f(ie+1) / 2.;
        NH_boundary = col(direction, k, j, ie+1, pmy_chemnet->iNHtot_) + NH_ghostzone; 
        NH2_boundary = col(direction, k, j, ie+1, pmy_chemnet->iNH2_)
          + pmy_mb->pspec->s(iH2, k, j, ie+1) * NH_ghostzone; 
        NCO_boundary = col(direction, k, j, ie+1, pmy_chemnet->iNCO_)
          + pmy_mb->pspec->s(iCO, k, j, ie+1) * NH_ghostzone;
        for (int i=ie; i>=is; --i) {
          col(direction, k, j, i, pmy_chemnet->iNHtot_) += NH_boundary;
          col(direction, k, j, i, pmy_chemnet->iNH2_) += NH2_boundary;
          col(direction, k, j, i, pmy_chemnet->iNCO_) += NCO_boundary;
        }
      }
    }
  } else if (direction == INNER_X2) {
    //+y
    for (int k=ks; k<=ke; ++k) {
      for (int i=is; i<=ie; ++i) {
        NH_ghostzone = unit_length_in_cm_ *
          pmy_mb->phydro->w(IDN, k, js-1, i) * pmy_mb->pcoord->dx2f(js-1) / 2.;
        NH_boundary = col(direction, k, js-1, i, pmy_chemnet->iNHtot_) + NH_ghostzone; 
        NH2_boundary = col(direction, k, js-1, i, pmy_chemnet->iNH2_)
          + pmy_mb->pspec->s(iH2, k, js-1, i) * NH_ghostzone; 
        NCO_boundary = col(direction, k, js-1, i, pmy_chemnet->iNCO_)
          + pmy_mb->pspec->s(iCO, k, js-1, i) * NH_ghostzone;
        for (int j=js; j<=je; ++j) {
          col(direction, k, j, i, pmy_chemnet->iNHtot_) += NH_boundary;
          col(direction, k, j, i, pmy_chemnet->iNH2_) += NH2_boundary;
          col(direction, k, j, i, pmy_chemnet->iNCO_) += NCO_boundary;
        }
      }
    }
  } else if (direction == OUTER_X2) {
    //-y
    for (int k=ks; k<=ke; ++k) {
      for (int i=is; i<=ie; ++i) {
        NH_ghostzone = unit_length_in_cm_ *
          pmy_mb->phydro->w(IDN, k, je+1, i) * pmy_mb->pcoord->dx2f(je+1) / 2.;
        NH_boundary = col(direction, k, je+1, i, pmy_chemnet->iNHtot_) + NH_ghostzone; 
        NH2_boundary = col(direction, k, je+1, i, pmy_chemnet->iNH2_)
          + pmy_mb->pspec->s(iH2, k, je+1, i) * NH_ghostzone; 
        NCO_boundary = col(direction, k, je+1, i, pmy_chemnet->iNCO_)
          + pmy_mb->pspec->s(iCO, k, je+1, i) * NH_ghostzone;
        for (int j=je; j>=js; --j) {
          col(direction, k, j, i, pmy_chemnet->iNHtot_) += NH_boundary;
          col(direction, k, j, i, pmy_chemnet->iNH2_) += NH2_boundary;
          col(direction, k, j, i, pmy_chemnet->iNCO_) += NCO_boundary;
        }
      }
    }
  } else if (direction == INNER_X3) {
    //+z
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        NH_ghostzone = unit_length_in_cm_ *
          pmy_mb->phydro->w(IDN, ks-1, j, i) * pmy_mb->pcoord->dx3f(ks-1) / 2.;
        NH_boundary = col(direction, ks-1, j, i, pmy_chemnet->iNHtot_) + NH_ghostzone; 
        NH2_boundary = col(direction, ks-1, j, i, pmy_chemnet->iNH2_)
          + pmy_mb->pspec->s(iH2, ks-1, j, i) * NH_ghostzone; 
        NCO_boundary = col(direction, ks-1, j, i, pmy_chemnet->iNCO_)
          + pmy_mb->pspec->s(iCO, ks-1, j, i) * NH_ghostzone;
        for (int k=ks; k<=ke; ++k) {
          col(direction, k, j, i, pmy_chemnet->iNHtot_) += NH_boundary;
          col(direction, k, j, i, pmy_chemnet->iNH2_) += NH2_boundary;
          col(direction, k, j, i, pmy_chemnet->iNCO_) += NCO_boundary;
        }
      }
    }
  } else if (direction == OUTER_X3) {
    //-z
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        NH_ghostzone = unit_length_in_cm_ *
          pmy_mb->phydro->w(IDN, ke+1, j, i) * pmy_mb->pcoord->dx3f(ke+1) / 2.;
        NH_boundary = col(direction, ke+1, j, i, pmy_chemnet->iNHtot_) + NH_ghostzone; 
        NH2_boundary = col(direction, ke+1, j, i, pmy_chemnet->iNH2_)
          + pmy_mb->pspec->s(iH2, ke+1, j, i) * NH_ghostzone; 
        NCO_boundary = col(direction, ke+1, j, i, pmy_chemnet->iNCO_)
          + pmy_mb->pspec->s(iCO, ke+1, j, i) * NH_ghostzone;
        for (int k=ke; k>=ks; --k) {
          col(direction, k, j, i, pmy_chemnet->iNHtot_) += NH_boundary;
          col(direction, k, j, i, pmy_chemnet->iNH2_) += NH2_boundary;
          col(direction, k, j, i, pmy_chemnet->iNCO_) += NCO_boundary;
        }
      }
    }
  } else {
    msg << "### FATAL ERROR in RadIntegrator six_ray [UpdateCol]" << std::endl
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
  int iang_arr[6] = {INNER_X1, INNER_X2, INNER_X3, OUTER_X1, OUTER_X2, OUTER_X3};
  for (int k=ks-NGHOST; k<=ke+NGHOST; ++k) {
    for (int j=js-NGHOST; j<=je+NGHOST; ++j) {
      for (int i=is-NGHOST; i<=ie+NGHOST; ++i) {
        for (int iang=0; iang < 6; iang++) {
          //column densities
          col_Htot(iang, k, j, i) = col(iang_arr[iang], k, j, i, pmy_chemnet->iNHtot_);
          col_H2(iang, k, j, i) = col(iang_arr[iang], k, j, i, pmy_chemnet->iNH2_); 
          col_CO(iang, k, j, i) = col(iang_arr[iang], k, j, i, pmy_chemnet->iNCO_); 
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
      pfacenb_[INNER_X1] = nb;
    } else if (nb->fid == INNER_X2) {
      pfacenb_[INNER_X2] = nb;
    } else if (nb->fid == INNER_X3) {
      pfacenb_[INNER_X3] = nb;
    } else if (nb->fid == OUTER_X1) {
      pfacenb_[OUTER_X1] = nb;
    } else if (nb->fid == OUTER_X2) {
      pfacenb_[OUTER_X2] = nb;
    } else if (nb->fid == OUTER_X3) {
      pfacenb_[OUTER_X3] = nb;
    }
  }
  return;
}
#endif
