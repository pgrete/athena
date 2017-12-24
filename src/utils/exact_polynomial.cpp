//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// Athena headers
#include "../athena.hpp"

#include <iostream>


//======================================================================================
//! \file exact_polynomial.cpp
//======================================================================================

//--------------------------------------------------------------------------------------
//! \fn  double Rtsafe(void (*funcd)(double, double, double, double, double, double *,
//        double *),  double x1, double x2, double xacc,
//      double coef1, double coef2, double coef3, double coef4)

// \brief laguer method to find root of polynomial, taken from numerical recipes
// for two special Polynomials to increase performance

// Exact solution for fourth order polynomical with the format
// coef4 * x^4 + coef * x + tconst == 0

int ExactPolynomial(const Real coef4, const Real coef, const Real tconst, Real &root)
{

  Real coef4_2 = coef4 * coef4;
  Real coef4_3 = coef4_2 * coef4;
  Real coef_2 = coef * coef;
  Real coef_4 = coef_2 * coef_2;
  Real tconst_3 = tconst * tconst * tconst;
  
  Real del1 = 27.0 * coef4_2 * coef_4 - 256.0 * coef4_3 * tconst_3;
  Real sqrt3 = sqrt(3.0);
  
  if(del1 < 0.0){
    return -1;
  }
  Real del2 = 9.0 * coef4 * coef_2 + sqrt3 * sqrt(del1);
  del2 = pow(del2, 1.0/3.0);
  // 4*(2/3)^1/3 = 4*0.8735804647362989
  Real top0 = 3.4943218589451956 * tconst;
  // 2^1/3 * 3^2/3 = 2.620741394208897
  Real top1 = 2.620741394208897 * coef4;
  Real del3 = top0/del2 + del2/top1;
  if(del3 <= 0.0){
      return -1;
  }
  
  Real del4 = sqrt(del3);
  root = -0.5 * del4 + 0.5 * sqrt(-del3+2.0*coef/(coef4*del4));

  return 0;
}
