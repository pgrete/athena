//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================

// Athena headers
#include "../athena.hpp"

#include <iostream>
#include <sstream>  // msg
#include <stdexcept> // runtime erro
#include <cmath>

#define MR 8
#define MT 10
#define MAXIT (MT*MR)

//======================================================================================
//! \file laguer.cpp
//======================================================================================

//--------------------------------------------------------------------------------------
//! \fn  double Rtsafe(void (*funcd)(double, double, double, double, double, double *,
//        double *),  double x1, double x2, double xacc,
//      double coef1, double coef2, double coef3, double coef4)

// \brief laguer method to find root of polynomial, taken from numerical recipes
// for two special Polynomials to increase performance


void Laguer(Real *coef, int m, Real &root)
{
  const Real eps = 1.e-7;
  int mr=8, mt=10;
  int max_it= mr * mt;
  
  std::stringstream msg;


  Real err, abp, abm;
  Real b, d, f, g, g2, h, gp, sq, gm, dx, x1;
  static Real frac[MR+1] = {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
  
  for(int iter=1; iter<=MAXIT; ++iter){
    
    d=0.0;
    f=0.0;
    b=coef[m];
    err=std::abs(b);
    for(int i=m-1; i>=0; --i){
      f=(root)*f + d;
      d=(root)*d + b;
      b=(root)*b + coef[i];
      err=std::abs(b)+root*err;
    }
    err *= eps;
    if (std::abs(b) <= err) return;
    g=d/b;
    g2=g*g;
    h=g2-2.0*f/b;
    Real diff=double(m-1)*(double(m)*h-g2);
    if (diff < 0.0 || root < 0.0) {
      std::cout << "[laguer]: Need better estimate " << m << "h: " << h <<
            "g2: " << g2 << "coefm" << coef[m] << "root" << root <<
            "coef0" << coef[0] << "coef1" << coef[1] << "\n" << std::endl;
      sq = 0.0;
    }else{
      sq=sqrt(diff);
    }
    gp=g+sq;
    gm=g-sq;
    abp=std::abs(gp);
    abm=std::abs(gm);
    if (abp < abm) gp=gm;
    dx=std::max(abp,abm) > 0.0 ? double(m)/gp : 0.0;
    x1=root-dx;
    if (root==x1) return;
    if (iter % MT != 0) root =x1;
    else root -= std::abs(frac[iter/MT]*dx);
  
  }// End itr



  std::cout << "[laguer]:Maximum number of iterations exceeded in laguer: root:"
            << root <<  "coef0: " << coef[0] << "coef1: "
            << coef[1] <<  "coef4: " << coef[4] << "\n" << std::endl;


  msg << "### FATAL ERROR in function [laguer]" << std::endl;
    throw std::runtime_error(msg.str().c_str());

  return;
}
