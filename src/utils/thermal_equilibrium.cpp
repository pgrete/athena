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

// Athena headers
#include "../athena.hpp"

#include <iostream>
#include <sstream>  // msg
#include <stdexcept> // runtime erro


//======================================================================================
//! \file thermal_equilibrium.cpp
//======================================================================================

//--------------------------------------------------------------------------------------
//! \fn  Tequilibrium(Real temperature, Real coef1, Real coef2, Real coef3,
//      Real coef4, Real *fval, Real *dfval)

// \brief Function for thermal energy exchange


void Tequilibrium(Real temperature, Real coef1, Real coef2, Real coef3,
      Real coef4, Real *fval, Real *dfval)
{
  // function is
  //  coef1 * T^4 + coef2 * T + coef3 == 0
  Real temp3 = temperature * temperature * temperature;
  Real temp4 = temp3 * temperature;
 

  *fval = coef1 * temp4 + coef2 * temperature + coef3;
  *dfval = 4.0 * coef1 * temp3 + coef2;

  return;

}



// Function to calculate the equilibrium state due to Compton scattering */
void Tcompton(Real temperature, Real coef1, Real coef2, Real coef3, Real coef4,
      Real *fval, Real *dfval)
{

// function is
//  coef1 * T^8 + coef2 * T^5 + coef3 * T^4 + coef4 == 0 *
  Real temp3 = temperature * temperature * temperature;
  Real temp4 = temp3 * temperature;
  Real temp8 = temp4 * temp4;
  Real temp7 = temp3 * temp4;
  Real temp5 = temperature * temp4;

  *fval = coef1 * temp8 + coef2 * temp5 + coef3 * temp4 + coef4;
  *dfval = 8.0 * coef1 * temp7 + 5.0 * coef2 * temp4 + 4.0 * coef3 * temp3;

  return;
}
