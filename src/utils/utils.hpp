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


#endif // UTILS_HPP
