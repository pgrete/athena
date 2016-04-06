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
//! \file network_wrapper.cpp
//  \brief implementation of functions in class NetworkWrapper 
//======================================================================================

// this class header
#include "network.hpp"

NetworkWrapper::NetworkWrapper() {}

NetworkWrapper::~NetworkWrapper() {}

int NetworkWrapper::WrapJacobian(const long int N, const realtype t,
                          const N_Vector y, const N_Vector fy, 
                          DlsMat J, void *user_data,
                          N_Vector tmp1, N_Vector tmp2, N_Vector tmp3) {
  int r;
  NetworkWrapper *meptr = (NetworkWrapper*) user_data;
  r = meptr->Jacobian(N, t, y, fy, J, tmp1, tmp2, tmp3);
  return r;
}

int NetworkWrapper::WrapRHS(const realtype t, const N_Vector y,
                     N_Vector ydot, void *user_data) {
  int r;
  NetworkWrapper *meptr = (NetworkWrapper*) user_data;
  r = meptr->RHS(t, y, ydot);
  return r;
}
