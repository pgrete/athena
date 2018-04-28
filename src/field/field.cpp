//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file field.cpp
//  \brief implementation of functions in class Field

// C++ headers
#include <string>

// Athena++ headers
#include "field.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"

// constructor, initializes data structures and parameters

Field::Field(MeshBlock *pmb, ParameterInput *pin) {
  pmy_block = pmb;

  // Allocate memory for interface fields, but only when needed.
  if (MAGNETIC_FIELDS_ENABLED) {
    int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
    int ncells2 = 1, ncells3 = 1;
    if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
    if (pmb->block_size.nx3 > 1) ncells3 = pmb->block_size.nx3 + 2*(NGHOST);

    //  Note the extra cell in each longitudinal dirn for interface fields
    b.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
    b.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
    b.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );

    b1.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
    b1.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
    b1.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );
    // Allocated 3rd registers if user-requested time integrator is type 3N or 3S*
    std::string integrator = pin->GetOrAddString("time","integrator","vl2");
    if (integrator == "ssprk5_4") {
      // future extension may add "int nregister" to Hydro class
      b2.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
      b2.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
      b2.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );
    }

    bcc.NewAthenaArray(NFIELD, ncells3, ncells2, ncells1);

    //-------- begin allocations for fourth-order MHD
    // TODO(kfelker): add xorder==4 switch for array allocation
    b_fc.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
    b_fc.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
    b_fc.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );
    bcc_center.NewAthenaArray(NFIELD, ncells3, ncells2, ncells1);

    // fourth-order UCT reconstructions at corners
    // TODO(kfelker): cut down on these temporary arrays
    // TODO(kfelker): check array limits-- add +1 for upper face?
    by_W.NewAthenaArray(ncells3, ncells2, ncells1);
    by_E.NewAthenaArray(ncells3, ncells2, ncells1);
    bx_S.NewAthenaArray(ncells3, ncells2, ncells1);
    bx_N.NewAthenaArray(ncells3, ncells2, ncells1);

    // TODO(kfelker): expand to 3D
    v_NE.NewAthenaArray(2, ncells3, ncells2, ncells1);
    v_SE.NewAthenaArray(2, ncells3, ncells2, ncells1);
    v_NW.NewAthenaArray(2, ncells3, ncells2, ncells1);
    v_SW.NewAthenaArray(2, ncells3, ncells2, ncells1);
    vl_temp_.NewAthenaArray(2, ncells3, ncells2, ncells1);
    vr_temp_.NewAthenaArray(2, ncells3, ncells2, ncells1);
    //-------- end allocations for fourth-order MHD

    e.x1e.NewAthenaArray((ncells3+1),(ncells2+1), ncells1   );
    e.x2e.NewAthenaArray((ncells3+1), ncells2   ,(ncells1+1));
    e.x3e.NewAthenaArray( ncells3   ,(ncells2+1),(ncells1+1));

    wght.x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
    wght.x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
    wght.x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );

    e2_x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
    e3_x1f.NewAthenaArray( ncells3   , ncells2   ,(ncells1+1));
    e1_x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
    e3_x2f.NewAthenaArray( ncells3   ,(ncells2+1), ncells1   );
    e1_x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );
    e2_x3f.NewAthenaArray((ncells3+1), ncells2   , ncells1   );

    // Allocate memory for scratch vectors
    cc_e_.NewAthenaArray(ncells3,ncells2,ncells1);

    face_area_.NewAthenaArray(ncells1);
    edge_length_.NewAthenaArray(ncells1);
    edge_length_p1_.NewAthenaArray(ncells1);
    if (GENERAL_RELATIVITY) {
      g_.NewAthenaArray(NMETRIC,ncells1);
      gi_.NewAthenaArray(NMETRIC,ncells1);
    }

    // fourth-order MHD
    // 4D scratch arrays
    scr2_nkji_.NewAthenaArray(NFIELD, ncells3, ncells2, ncells1);
  }
}

// destructor

Field::~Field() {
  b.x1f.DeleteAthenaArray();
  b.x2f.DeleteAthenaArray();
  b.x3f.DeleteAthenaArray();
  b1.x1f.DeleteAthenaArray();
  b1.x2f.DeleteAthenaArray();
  b1.x3f.DeleteAthenaArray();
  // b2 only allocated if integrator was 3S* integrator
  b2.x1f.DeleteAthenaArray();
  b2.x2f.DeleteAthenaArray();
  b2.x3f.DeleteAthenaArray();
  bcc.DeleteAthenaArray();

  //-------- begin deallocations for fourth-order MHD
  b_fc.x1f.DeleteAthenaArray();
  b_fc.x2f.DeleteAthenaArray();
  b_fc.x3f.DeleteAthenaArray();
  bcc_center.DeleteAthenaArray();

  by_W.DeleteAthenaArray();
  by_E.DeleteAthenaArray();
  bx_N.DeleteAthenaArray();
  bx_S.DeleteAthenaArray();

  v_NE.DeleteAthenaArray();
  v_SE.DeleteAthenaArray();
  v_NW.DeleteAthenaArray();
  v_SW.DeleteAthenaArray();
  vl_temp_.DeleteAthenaArray();
  vr_temp_.DeleteAthenaArray();
  //----------end deallocations for fourth-order MHD

  e.x1e.DeleteAthenaArray();
  e.x2e.DeleteAthenaArray();
  e.x3e.DeleteAthenaArray();
  wght.x1f.DeleteAthenaArray();
  wght.x2f.DeleteAthenaArray();
  wght.x3f.DeleteAthenaArray();
  e2_x1f.DeleteAthenaArray();
  e3_x1f.DeleteAthenaArray();
  e1_x2f.DeleteAthenaArray();
  e3_x2f.DeleteAthenaArray();
  e1_x3f.DeleteAthenaArray();
  e2_x3f.DeleteAthenaArray();

  cc_e_.DeleteAthenaArray();
  face_area_.DeleteAthenaArray();
  edge_length_.DeleteAthenaArray();
  edge_length_p1_.DeleteAthenaArray();
  if (GENERAL_RELATIVITY) {
    g_.DeleteAthenaArray();
    gi_.DeleteAthenaArray();
  }
}

//----------------------------------------------------------------------------------------
// \! fn
// \! brief

void Field::CalculateCellCenteredField(const FaceField &bf, AthenaArray<Real> &bc,
            Coordinates *pco, int is, int ie, int js, int je, int ks, int ke) {
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
    // calc cell centered fields first
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        const Real& b1_i   = bf.x1f(k,j,i  );
        const Real& b1_ip1 = bf.x1f(k,j,i+1);
        const Real& b2_j   = bf.x2f(k,j  ,i);
        const Real& b2_jp1 = bf.x2f(k,j+1,i);
        const Real& b3_k   = bf.x3f(k  ,j,i);
        const Real& b3_kp1 = bf.x3f(k+1,j,i);

        Real& bcc1 = bc(IB1,k,j,i);
        Real& bcc2 = bc(IB2,k,j,i);
        Real& bcc3 = bc(IB3,k,j,i);

        // cell center B-fields are defined as spatial interpolation at the volume center
        const Real& x1f_i  = pco->x1f(i);
        const Real& x1f_ip = pco->x1f(i+1);
        const Real& x1v_i  = pco->x1v(i);
        const Real& dx1_i  = pco->dx1f(i);
        Real lw=(x1f_ip-x1v_i)/dx1_i;
        Real rw=(x1v_i -x1f_i)/dx1_i;
        bcc1 = lw*b1_i + rw*b1_ip1;
        const Real& x2f_j  = pco->x2f(j);
        const Real& x2f_jp = pco->x2f(j+1);
        const Real& x2v_j  = pco->x2v(j);
        const Real& dx2_j  = pco->dx2f(j);
        lw=(x2f_jp-x2v_j)/dx2_j;
        rw=(x2v_j -x2f_j)/dx2_j;
        bcc2 = lw*b2_j + rw*b2_jp1;
        const Real& x3f_k  = pco->x3f(k);
        const Real& x3f_kp = pco->x3f(k+1);
        const Real& x3v_k  = pco->x3v(k);
        const Real& dx3_k  = pco->dx3f(k);
        lw=(x3f_kp-x3v_k)/dx3_k;
        rw=(x3v_k -x3f_k)/dx3_k;
        bcc3 = lw*b3_k + rw*b3_kp1;
      }
    }
  }
  return;
}
