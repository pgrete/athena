#ifndef UTILS_HPP
#define UTILS_HPP
//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
// See LICENSE file for full public license information.
//======================================================================================
//! \file utils.hpp
//  \brief prototypes of utility functions in utils/*.cpp
//======================================================================================

void ChangeRunDir(const char *pdir);
double ran2(long int *idum);
void ShowConfig();
int Permutation(int i, int j, int k, int np, AthenaArray<int> &pl);
void MatrixMult(int m, int n, AthenaArray<Real> &a,
                AthenaArray<Real> &b, AthenaArray<Real> &c);
void Gauleg(int n, Real x1, Real x2,  AthenaArray<Real> &x,
            AthenaArray<Real> &w);
void Ludcmp_nr(int n, AthenaArray<Real> &a, AthenaArray<int> &indx,
               Real *d);
void Lubksb_nr(int n, AthenaArray<Real> &a, AthenaArray<int> &indx,
               AthenaArray<Real> &b);
void InverseMatrix(int n, AthenaArray<Real> &a, AthenaArray<Real> &b);

double Rtsafe(void (*funcd)(double, double, double, double, double, double *, double *),
      double x1, double x2, double xacc,
      double coef1, double coef2, double coef3, double coef4);

void Tequilibrium(Real temperature, Real coef1, Real coef2, Real coef3,
      Real coef4, Real *fval, Real *dfval);

void Tcompton(Real temperature, Real coef1, Real coef2, Real coef3, Real coef4,
      Real *fval, Real *dfval);

void Laguer(Real *coef, int m, Real &root);

void wtime(double *t);

int ExactPolynomial(const Real coef4, const Real coef, const Real tconst, Real &root);
int FouthPolyRoot(const Real coef4, const Real tconst, Real &root);

namespace WallTimeLimit {
  void InitWTLimit(void);
  void SendWTLimit(int nwtlimit);
  bool TestWTLimit(int &nwtlimit);
  void FinalizeWTLimit(int wtflag);
}

#endif // UTILS_HPP
