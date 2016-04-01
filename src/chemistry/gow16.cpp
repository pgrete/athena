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
//! \file gow16.cpp
//  \brief implementation of functions in class ChemNetwork, using the GOW16
//  network, see paper by Gong, Ostriker, Wolfire 2016 
//======================================================================================

// this class header
#include "species.hpp"

ChemNetwork::ChemNetwork(ChemSpecies *pspec, ParameterInput *pin) {
  pmy_spec_ = pspec;
}

ChemNetwork::~ChemNetwork() {}

int ChemNetwork::RHS(const realtype t, const N_Vector y, N_Vector ydot) {
  return 0;
}

int ChemNetwork::Jacobian(const long int N, const realtype t,
    const N_Vector y, const N_Vector fy, 
    DlsMat J, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  return 0;
}
