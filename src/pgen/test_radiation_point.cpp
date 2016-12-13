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
//! \file test_radiation.cpp
//  \brief problem generator, initalize mesh by read in vtk files.
//======================================================================================

// C++ headers
#include <string>     // c_str()
#include <iostream>   // endl
#include <vector>     // vector container
#include <sstream>    // stringstream
#include <stdio.h>    // c style file
#include <string.h>   // strcmp()
#include <algorithm>  // std::find()
#include <stdexcept>  // std::runtime_error()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../chemistry/species.hpp"
#include "../chemistry/thermo.hpp"
#include "../radiation/radiation.hpp"
#include "../field/field.hpp"
#include "../eos/eos.hpp"
#include "../coordinates/coordinates.hpp"

void SixRayBoundaryInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void SixRayBoundaryOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void SixRayBoundaryInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void SixRayBoundaryOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void SixRayBoundaryInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);
void SixRayBoundaryOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke);


//======================================================================================
//! \fn void Mesh::TerminateUserMeshProperties(void)
//  \brief Clean up the Mesh properties
//======================================================================================
void Mesh::UserWorkAfterLoop(ParameterInput *pin)
{
  FILE *pf = fopen("chem_network.dat", "w");
  pblock->pspec->pchemnet->OutputProperties(pf);
  fclose(pf);
  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief initialize problem by reading in vtk file.
//======================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  //dimensions of meshblock
  const int Nx = ie - is + 1;
  const int Ny = je - js + 1;
  const int Nz = ke - ks + 1;
	//read density and radiation field strength
	const Real nH = pin->GetReal("problem", "nH");
	const Real G0 = pin->GetReal("problem", "G0");
	const Real s_init = pin->GetReal("problem", "s_init");
	const Real xH2_init = pin->GetReal("problem", "xH2_init");
	const Real xCO_init = pin->GetReal("problem", "xCO_init");
	const int ipx = pin->GetInteger("problem", "ipx");
	const int ipy = pin->GetInteger("problem", "ipy");
	const int ipz = pin->GetInteger("problem", "ipz");
	//set density
  int gis = loc.lx1 * Nx;
  int gjs = loc.lx2 * Ny;
  int gks = loc.lx3 * Nz;
  
	for (int k=ks; k<=ke; ++k) {
		for (int j=js; j<=je; ++j) {
			for (int i=is; i<=ie; ++i) {
        if (i-is+gis == ipx and j-js+gjs == ipy and k-ks+gks == ipz) {
          phydro->u(IDN, k, j, i) = nH;//nH;
        } else {
          phydro->u(IDN, k, j, i) = 0;//nH;
        }
			}
		}
	}

	//intialize radiation field
	if (RADIATION_ENABLED) {
      for (int k=ks; k<=ke; ++k) {
        for (int j=js; j<=je; ++j) {
          for (int i=is; i<=ie; ++i) {
						for (int ifreq=0; ifreq < prad->nfreq; ++ifreq) {
							for (int iang=0; iang < prad->nang; ++iang) {
								prad->ir(k, j, i, ifreq * prad->nang + iang) = G0;
							}
						}
          }
        }
      }
	}

	//intialize chemical species
#ifdef INCLUDE_CHEMISTRY
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      for (int i=is; i<=ie; ++i) {
        for (int ispec=0; ispec < NSPECIES; ++ispec) {
          pspec->s(ispec, k, j, i) = s_init;
        }
        pspec->s(pspec->pchemnet->iH2_, k, j, i) = xH2_init;
        pspec->s(pspec->pchemnet->iCO_, k, j, i) = xCO_init;
        //temperature
        pspec->s(pspec->pchemnet->iE_, k, j, i) = 20. *
          Thermo::CvCold(0., 0.1, 0.);
      }
    }
  }
#endif

  return;
}


void SixRayBoundaryInnerX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
    FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          pmb->pspec->s(n,k,j,is-i) = 0;
        }
      }
    }
  }
  //set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          prim(n,k,j,is-i) = 0;
        }
      }
    }
  }
  return;
}

void SixRayBoundaryInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
    FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          pmb->pspec->s(n,k,js-j,i) = 0;
        }
      }
    }
  }
  //set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          prim(n,k,js-j,i) = 0;
        }
      }
    }
  }
  return;
}

void SixRayBoundaryInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
    FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=1; k<=(NGHOST); ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          pmb->pspec->s(n,ks-k,j,i) = 0;
        }
      }
    }
  }
  //set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=1; k<=(NGHOST); ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          prim(n,ks-k,j,i) = 0;
        }
      }
    }
  }
  return;
}

void SixRayBoundaryOuterX1(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
    FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          pmb->pspec->s(n,k,j,ie+i) = 0;
        }
      }
    }
  }
  //set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=1; i<=(NGHOST); ++i) {
          prim(n,k,j,ie+i) = 0;
        }
      }
    }
  }
  return;
}

void SixRayBoundaryOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
    FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          pmb->pspec->s(n,k,je+j,i) = 0;
        }
      }
    }
  }
  //set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=1; j<=(NGHOST); ++j) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          prim(n,k,je+j,i) = 0;
        }
      }
    }
  }
  return;
}

void SixRayBoundaryOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
    FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke) {
  // set species and column boundary to zero
  for (int n=0; n<(NSPECIES); ++n) {
    for (int k=1; k<=(NGHOST); ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          pmb->pspec->s(n,ke+k,j,i) = 0;
        }
      }
    }
  }
  //set hydro variables to zero
  for (int n=0; n<(NHYDRO); ++n) {
    for (int k=1; k<=(NGHOST); ++k) {
      for (int j=js; j<=je; ++j) {
#pragma simd
        for (int i=is; i<=ie; ++i) {
          prim(n,ke+k,j,i) = 0;
        }
      }
    }
  }
  return;
}

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Function to initialize problem-specific data in mesh class.  Can also be used
//  to initialize variables which are global to (and therefore can be passed to) other
//  functions in this file.  Called in Mesh constructor.
//========================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  EnrollUserBoundaryFunction(INNER_X1, SixRayBoundaryInnerX1);
  EnrollUserBoundaryFunction(OUTER_X1, SixRayBoundaryOuterX1);
  EnrollUserBoundaryFunction(INNER_X2, SixRayBoundaryInnerX2);
  EnrollUserBoundaryFunction(OUTER_X2, SixRayBoundaryOuterX2);
  EnrollUserBoundaryFunction(INNER_X3, SixRayBoundaryInnerX3);
  EnrollUserBoundaryFunction(OUTER_X3, SixRayBoundaryOuterX3);
  return;
}
