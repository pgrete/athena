//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_cc.cpp
//  \brief functions that apply BCs for CELL_CENTERED variables

// C++ headers
#include <iostream>   // endl
#include <iomanip>
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <cstring>    // memcpy
#include <cstdlib>
#include <cmath>

// Athena++ classes headers
#include "bvals.hpp"
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../coordinates/coordinates.hpp"
#include "../parameter_input.hpp"
#include "../utils/buffer_utils.hpp"
#include "../utils/cgk_utils.hpp"
#include "../radiation/radiation.hpp"
#include "../radiation/integrators/rad_integrators.hpp"
#ifdef INCLUDE_CHEMISTRY
#include "../chemistry/species.hpp"
#endif

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif



//----------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadCellCenteredBoundaryBufferSameLevel(AthenaArray<Real> &src,
//                                  int ns, int ne, Real *buf, const NeighborBlock& nb)
//  \brief Set hydro boundary buffers for sending to a block on the same level

int BoundaryValues::LoadCellCenteredBoundaryBufferSameLevel(AthenaArray<Real> &src,
                                    int ns, int ne, Real *buf, const NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_block_;
  int si, sj, sk, ei, ej, ek;

  si=(nb.ox1>0)?(pmb->ie-NGHOST+1):pmb->is;
  ei=(nb.ox1<0)?(pmb->is+NGHOST-1):pmb->ie;
  sj=(nb.ox2>0)?(pmb->je-NGHOST+1):pmb->js;
  ej=(nb.ox2<0)?(pmb->js+NGHOST-1):pmb->je;
  sk=(nb.ox3>0)?(pmb->ke-NGHOST+1):pmb->ks;
  ek=(nb.ox3<0)?(pmb->ks+NGHOST-1):pmb->ke;
  int p=0;
  BufferUtility::Pack4DData(src, buf, ns, ne, si, ei, sj, ej, sk, ek, p);
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadCellCenteredBoundaryBufferToCoarser(AthenaArray<Real> &src,
//          int ns, int ne, AthenaArray<Real> &cbuf, Real *buf, const NeighborBlock& nb)
//  \brief Set hydro boundary buffers for sending to a block on the coarser level

int BoundaryValues::LoadCellCenteredBoundaryBufferToCoarser(AthenaArray<Real> &src,
    int ns, int ne, Real *buf, AthenaArray<Real> &cbuf, const NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_block_;
  MeshRefinement *pmr=pmb->pmr;
  int si, sj, sk, ei, ej, ek;
  int cn=pmb->cnghost-1;

  si=(nb.ox1>0)?(pmb->cie-cn):pmb->cis;
  ei=(nb.ox1<0)?(pmb->cis+cn):pmb->cie;
  sj=(nb.ox2>0)?(pmb->cje-cn):pmb->cjs;
  ej=(nb.ox2<0)?(pmb->cjs+cn):pmb->cje;
  sk=(nb.ox3>0)?(pmb->cke-cn):pmb->cks;
  ek=(nb.ox3<0)?(pmb->cks+cn):pmb->cke;

  int p=0;
  pmr->RestrictCellCenteredValues(src, cbuf, ns, ne,
                                  si, ei, sj, ej, sk, ek);
  BufferUtility::Pack4DData(cbuf, buf, ns, ne,
                            si, ei, sj, ej, sk, ek, p);
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadCellCenteredBoundaryBufferToFiner(AthenaArray<Real> &src,
//                                  int ns, int ne, Real *buf, const NeighborBlock& nb)
//  \brief Set hydro boundary buffers for sending to a block on the finer level

int BoundaryValues::LoadCellCenteredBoundaryBufferToFiner(AthenaArray<Real> &src,
                                    int ns, int ne, Real *buf, const NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_block_;
  int si, sj, sk, ei, ej, ek;
  int cn=pmb->cnghost-1;

  si=(nb.ox1>0)?(pmb->ie-cn):pmb->is;
  ei=(nb.ox1<0)?(pmb->is+cn):pmb->ie;
  sj=(nb.ox2>0)?(pmb->je-cn):pmb->js;
  ej=(nb.ox2<0)?(pmb->js+cn):pmb->je;
  sk=(nb.ox3>0)?(pmb->ke-cn):pmb->ks;
  ek=(nb.ox3<0)?(pmb->ks+cn):pmb->ke;

  // send the data first and later prolongate on the target block
  // need to add edges for faces, add corners for edges
  if(nb.ox1==0) {
    if(nb.fi1==1)   si+=pmb->block_size.nx1/2-pmb->cnghost;
    else            ei-=pmb->block_size.nx1/2-pmb->cnghost;
  }
  if(nb.ox2==0 && pmb->block_size.nx2 > 1) {
    if(nb.ox1!=0) {
      if(nb.fi1==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
      else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
    }
    else {
      if(nb.fi2==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
      else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
    }
  }
  if(nb.ox3==0 && pmb->block_size.nx3 > 1) {
    if(nb.ox1!=0 && nb.ox2!=0) {
      if(nb.fi1==1) sk+=pmb->block_size.nx3/2-pmb->cnghost;
      else          ek-=pmb->block_size.nx3/2-pmb->cnghost;
    }
    else {
      if(nb.fi2==1) sk+=pmb->block_size.nx3/2-pmb->cnghost;
      else          ek-=pmb->block_size.nx3/2-pmb->cnghost;
    }
  }

  int p=0;
  BufferUtility::Pack4DData(src, buf, ns, ne, si, ei, sj, ej, sk, ek, p);
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendCellCenteredBoundaryBuffers(AthenaArray<Real> &src,
//                                                           enum CCBoundaryType type)
//  \brief Send boundary buffers of cell-centered variables
void BoundaryValues::SendCellCenteredBoundaryBuffers(AthenaArray<Real> &src,
                                                     enum CCBoundaryType type)
{
  MeshBlock *pmb=pmy_block_, *pbl;
  int mylevel=pmb->loc.level;
  int ns, ne;
  Real *sbuf, *rbuf;
  AthenaArray<Real> cbuf;
  enum BoundaryStatus *flag;
#ifdef MPI_PARALLEL
  MPI_Request *req;
#endif

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    if(nb.rank == Globals::my_rank) // on the same process
      pbl=pmb->pmy_mesh->FindMeshBlock(nb.gid);
    switch(type) {
      case HYDRO_CONS:
        sbuf=hydro_send_[nb.bufid];
        if(nb.level<mylevel) cbuf.InitWithShallowCopy(pmb->pmr->coarse_cons_);
        ns=0, ne=NHYDRO-1;
        if(nb.rank == Globals::my_rank) {
          rbuf=pbl->pbval->hydro_recv_[nb.targetid];
          flag=&(pbl->pbval->hydro_flag_[nb.targetid]);
        }
#ifdef MPI_PARALLEL
        else
          req=&(req_hydro_send_[nb.bufid]);
#endif
        break;
      case HYDRO_PRIM:
        sbuf=hydro_send_[nb.bufid];
        if(nb.level<mylevel) cbuf.InitWithShallowCopy(pmb->pmr->coarse_prim_);
        ns=0, ne=NHYDRO-1;
        if(nb.rank == Globals::my_rank) {
          rbuf=pbl->pbval->hydro_recv_[nb.targetid];
          flag=&(pbl->pbval->hydro_flag_[nb.targetid]);
        }
#ifdef MPI_PARALLEL
        else
          req=&(req_hydro_send_[nb.bufid]);
#endif
        break;
      defualt:
        std::stringstream msg;
        msg << "### FATAL ERROR in SendCellCenteredBoundaryBuffers" << std::endl
            << "invalid boundary buffer type." << std::endl;
        throw std::runtime_error(msg.str().c_str());
        break;
    }
    int ssize;
    if(nb.level==mylevel)
      ssize=LoadCellCenteredBoundaryBufferSameLevel(src, ns, ne, sbuf, nb);
    else if(nb.level<mylevel)
      ssize=LoadCellCenteredBoundaryBufferToCoarser(src, ns, ne, sbuf, cbuf, nb);
    else
      ssize=LoadCellCenteredBoundaryBufferToFiner(src, ns, ne, sbuf, nb);
    if(nb.rank == Globals::my_rank) {
      std::memcpy(rbuf, sbuf, ssize*sizeof(Real));
      *flag=BNDRY_ARRIVED;
    }
#ifdef MPI_PARALLEL
    else // MPI
      MPI_Start(req);
#endif
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetCellCenteredBoundarySameLevel(AthenaArray<Real> &dst,
//                         int ns, int ne, Real *buf, const NeighborBlock& nb, bool *flip)
//  \brief Set hydro boundary received from a block on the same level

void BoundaryValues::SetCellCenteredBoundarySameLevel(AthenaArray<Real> &dst,
                     int ns, int ne, Real *buf, const NeighborBlock& nb, bool *flip)
{
  MeshBlock *pmb=pmy_block_;
  int si, sj, sk, ei, ej, ek;

  if(nb.ox1==0)     si=pmb->is,        ei=pmb->ie;
  else if(nb.ox1>0) si=pmb->ie+1,      ei=pmb->ie+NGHOST;
  else              si=pmb->is-NGHOST, ei=pmb->is-1;
  if(nb.ox2==0)     sj=pmb->js,        ej=pmb->je;
  else if(nb.ox2>0) sj=pmb->je+1,      ej=pmb->je+NGHOST;
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  if(nb.ox3==0)     sk=pmb->ks,        ek=pmb->ke;
  else if(nb.ox3>0) sk=pmb->ke+1,      ek=pmb->ke+NGHOST;
  else              sk=pmb->ks-NGHOST, ek=pmb->ks-1;

  int p=0;
  if (nb.polar) {
    for (int n=ns; n<=ne; ++n) {
      Real sign = 1.0;
      if(flip!=NULL) sign =flip[n] ? -1.0 : 1.0;
      for (int k=sk; k<=ek; ++k) {
        for (int j=ej; j>=sj; --j) {
#pragma simd
          for (int i=si; i<=ei; ++i)
            dst(n,k,j,i) = sign * buf[p++];
        }
      }
    }
  }
  else
    BufferUtility::Unpack4DData(buf, dst, ns, ne, si, ei, sj, ej, sk, ek, p);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetCellCenteredBoundaryFromCoarser(int ns, int ne,
//           Real *buf, AthenaArray<Real> &cbuf, const NeighborBlock& nb, bool *flip)
//  \brief Set hydro prolongation buffer received from a block on a coarser level

void BoundaryValues::SetCellCenteredBoundaryFromCoarser(int ns, int ne,
             Real *buf, AthenaArray<Real> &cbuf, const NeighborBlock& nb, bool *flip)
{
  MeshBlock *pmb=pmy_block_;
  MeshRefinement *pmr=pmb->pmr;

  int si, sj, sk, ei, ej, ek;
  int cng=pmb->cnghost;

  if(nb.ox1==0) {
    si=pmb->cis, ei=pmb->cie;
    if((pmb->loc.lx1&1L)==0L) ei+=cng;
    else             si-=cng; 
  }
  else if(nb.ox1>0)  si=pmb->cie+1,   ei=pmb->cie+cng;
  else               si=pmb->cis-cng, ei=pmb->cis-1;
  if(nb.ox2==0) {
    sj=pmb->cjs, ej=pmb->cje;
    if(pmb->block_size.nx2 > 1) {
      if((pmb->loc.lx2&1L)==0L) ej+=cng;
      else             sj-=cng; 
    }
  }
  else if(nb.ox2>0)  sj=pmb->cje+1,   ej=pmb->cje+cng;
  else               sj=pmb->cjs-cng, ej=pmb->cjs-1;
  if(nb.ox3==0) {
    sk=pmb->cks, ek=pmb->cke;
    if(pmb->block_size.nx3 > 1) {
      if((pmb->loc.lx3&1L)==0L) ek+=cng;
      else             sk-=cng; 
    }
  }
  else if(nb.ox3>0)  sk=pmb->cke+1,   ek=pmb->cke+cng;
  else               sk=pmb->cks-cng, ek=pmb->cks-1;

  int p=0;
  if (nb.polar) {
    for (int n=ns; n<=ne; ++n) {
      Real sign = 1.0;
      if(flip!=NULL) sign = flip[n] ? -1.0 : 1.0;
      for (int k=sk; k<=ek; ++k) {
        for (int j=ej; j>=sj; --j) {
#pragma simd
          for (int i=si; i<=ei; ++i)
            cbuf(n,k,j,i) = sign * buf[p++];
        }
      }
    }
  }
  else
    BufferUtility::Unpack4DData(buf, cbuf, ns, ne, si, ei, sj, ej, sk, ek, p);
  return;
}


//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetCellCenteredBoundaryFromFiner(AthenaArray<Real> &dst,
//                   int ns, int ne, Real *buf, const NeighborBlock& nb, bool *flip)
//  \brief Set hydro boundary received from a block on a finer level

void BoundaryValues::SetCellCenteredBoundaryFromFiner(AthenaArray<Real> &dst,
                     int ns, int ne, Real *buf, const NeighborBlock& nb, bool *flip)
{
  MeshBlock *pmb=pmy_block_;
  // receive already restricted data
  int si, sj, sk, ei, ej, ek;

  if(nb.ox1==0) {
    si=pmb->is, ei=pmb->ie;
    if(nb.fi1==1)   si+=pmb->block_size.nx1/2;
    else            ei-=pmb->block_size.nx1/2;
  }
  else if(nb.ox1>0) si=pmb->ie+1,      ei=pmb->ie+NGHOST;
  else              si=pmb->is-NGHOST, ei=pmb->is-1;
  if(nb.ox2==0) {
    sj=pmb->js, ej=pmb->je;
    if(pmb->block_size.nx2 > 1) {
      if(nb.ox1!=0) {
        if(nb.fi1==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      }
      else {
        if(nb.fi2==1) sj+=pmb->block_size.nx2/2;
        else          ej-=pmb->block_size.nx2/2;
      }
    }
  }
  else if(nb.ox2>0) sj=pmb->je+1,      ej=pmb->je+NGHOST;
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  if(nb.ox3==0) {
    sk=pmb->ks, ek=pmb->ke;
    if(pmb->block_size.nx3 > 1) {
      if(nb.ox1!=0 && nb.ox2!=0) {
        if(nb.fi1==1) sk+=pmb->block_size.nx3/2;
        else          ek-=pmb->block_size.nx3/2;
      }
      else {
        if(nb.fi2==1) sk+=pmb->block_size.nx3/2;
        else          ek-=pmb->block_size.nx3/2;
      }
    }
  }
  else if(nb.ox3>0) sk=pmb->ke+1,      ek=pmb->ke+NGHOST;
  else              sk=pmb->ks-NGHOST, ek=pmb->ks-1;

  int p=0;
  if (nb.polar) {
    for (int n=ns; n<=ne; ++n) {
      Real sign=1.0;
      if(flip!=NULL) sign = flip[n] ? -1.0 : 1.0;
      for (int k=sk; k<=ek; ++k) {
        for (int j=sj; j<=ej; ++j) {
#pragma simd
          for (int i=si; i<=ei; ++i)
            dst(n,k,j,i) = sign * buf[p++];
        }
      }
    }
  }
  else 
    BufferUtility::Unpack4DData(buf, dst, ns, ne, si, ei, sj, ej, sk, ek, p);
  return;
}

//----------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveCellCenteredBoundaryBuffers(AthenaArray<Real> &dst,
//                                                              enum CCBoundaryType type)
//  \brief receive the boundary data

bool BoundaryValues::ReceiveCellCenteredBoundaryBuffers(AthenaArray<Real> &dst,
                                                        enum CCBoundaryType type)
{
  MeshBlock *pmb=pmy_block_;
  bool bflag=true;
  bool *flip=NULL;
  Real *rbuf;
  AthenaArray<Real> cbuf;
  int ns, ne;
  enum BoundaryStatus *flag;
#ifdef MPI_PARALLEL
  MPI_Request *req;
#endif

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    switch(type) {
      case HYDRO_CONS:
        ns=0, ne=NHYDRO-1;
        rbuf=hydro_recv_[nb.bufid];
        if(nb.level<pmb->loc.level) cbuf.InitWithShallowCopy(pmb->pmr->coarse_cons_);
        flag=&(hydro_flag_[nb.bufid]);
        flip=flip_across_pole_hydro;
#ifdef MPI_PARALLEL
        if(nb.rank!=Globals::my_rank)
          req=&(req_hydro_recv_[nb.bufid]);
#endif
        break;
      case HYDRO_PRIM:
        ns=0, ne=NHYDRO-1;
        rbuf=hydro_recv_[nb.bufid];
        if(nb.level<pmb->loc.level) cbuf.InitWithShallowCopy(pmb->pmr->coarse_prim_);
        flag=&(hydro_flag_[nb.bufid]);
        flip=flip_across_pole_hydro;
#ifdef MPI_PARALLEL
        if(nb.rank!=Globals::my_rank)
          req=&(req_hydro_recv_[nb.bufid]);
#endif
        break;
      defualt:
        std::stringstream msg;
        msg << "### FATAL ERROR in ReceiveCellCenteredBoundaryBuffers" << std::endl
            << "invalid boundary buffer type." << std::endl;
        throw std::runtime_error(msg.str().c_str());
        break;
    }
    if(*flag==BNDRY_COMPLETED) continue;
    if(*flag==BNDRY_WAITING) {
      if(nb.rank==Globals::my_rank) {// on the same process
        bflag=false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // MPI boundary
        int test;
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
        MPI_Test(req,&test,MPI_STATUS_IGNORE);
        if(test==false) {
          bflag=false;
          continue;
        }
        *flag = BNDRY_ARRIVED;
      }
#endif
    }
    if(nb.level==pmb->loc.level)
      SetCellCenteredBoundarySameLevel(dst, ns, ne, rbuf, nb, flip);
    else if(nb.level<pmb->loc.level) // this set only the prolongation buffer
      SetCellCenteredBoundaryFromCoarser(ns, ne, rbuf, cbuf, nb, flip);
    else
      SetCellCenteredBoundaryFromFiner(dst, ns, ne, rbuf, nb, flip);
    *flag = BNDRY_COMPLETED; // completed
  }

  if(bflag && (pmb->block_bcs[INNER_X2]==POLAR_BNDRY
           ||  pmb->block_bcs[OUTER_X2]==POLAR_BNDRY))
     PolarSingleCellCentered(dst, ns, ne);
  return bflag;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ReceiveCellCenteredBoundaryBuffersWithWait
//                                      (AthenaArray<Real> &dst, enum CCBoundaryType type)
//  \brief receive the boundary data for initialization

void BoundaryValues::ReceiveCellCenteredBoundaryBuffersWithWait(AthenaArray<Real> &dst,
                                                                enum CCBoundaryType type)
{
  MeshBlock *pmb=pmy_block_;
  bool *flip=NULL;
  Real *rbuf;
  AthenaArray<Real> cbuf;
  int ns, ne;
  enum BoundaryStatus *flag;
#ifdef MPI_PARALLEL
  MPI_Request *req;
#endif

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    switch(type) {
      case HYDRO_CONS:
        ns=0, ne=NHYDRO-1;
        rbuf=hydro_recv_[nb.bufid];
        if(nb.level<pmb->loc.level) cbuf.InitWithShallowCopy(pmb->pmr->coarse_cons_);
        flag=&(hydro_flag_[nb.bufid]);
        flip=flip_across_pole_hydro;
#ifdef MPI_PARALLEL
        if(nb.rank!=Globals::my_rank)
          req=&(req_hydro_recv_[nb.bufid]);
#endif
        break;
      case HYDRO_PRIM:
        ns=0, ne=NHYDRO-1;
        rbuf=hydro_recv_[nb.bufid];
        if(nb.level<pmb->loc.level) cbuf.InitWithShallowCopy(pmb->pmr->coarse_prim_);
        flag=&(hydro_flag_[nb.bufid]);
        flip=flip_across_pole_hydro;
#ifdef MPI_PARALLEL
        if(nb.rank!=Globals::my_rank)
          req=&(req_hydro_recv_[nb.bufid]);
#endif
        break;
      defualt:
        std::stringstream msg;
        msg << "### FATAL ERROR in ReceiveCellCenteredBoundaryBuffers" << std::endl
            << "invalid boundary buffer type." << std::endl;
        throw std::runtime_error(msg.str().c_str());
        break;
    }
#ifdef MPI_PARALLEL
    if(nb.rank!=Globals::my_rank)
      MPI_Wait(req,MPI_STATUS_IGNORE);
#endif
    if(nb.level==pmb->loc.level)
      SetCellCenteredBoundarySameLevel(dst, ns, ne, rbuf, nb, flip);
    else if(nb.level<pmb->loc.level)
      SetCellCenteredBoundaryFromCoarser(ns, ne, rbuf, cbuf, nb, flip);
    else
      SetCellCenteredBoundaryFromFiner(dst, ns, ne, rbuf, nb, flip);
    *flag = BNDRY_COMPLETED; // completed
  }
 
  if (pmb->block_bcs[INNER_X2]==POLAR_BNDRY || pmb->block_bcs[OUTER_X2]==POLAR_BNDRY)
    PolarSingleCellCentered(dst, ns, ne);

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::PolarSingleCellCentered(AthenaArray<Real> &dst,
//                                                   int ns, int ne)
//
// \brief  single CPU in the azimuthal direction for the polar boundary 

void BoundaryValues::PolarSingleCellCentered(AthenaArray<Real> &dst, int ns, int ne)
{
  MeshBlock *pmb=pmy_block_;
  if(pmb->loc.level == pmb->pmy_mesh->root_level && pmb->pmy_mesh->nrbx3 == 1){

    if(pmb->block_bcs[INNER_X2]==POLAR_BNDRY){
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int n=ns; n<=ne; ++n) {
        for (int j=pmb->js-NGHOST; j<=pmb->js-1; ++j) {
         for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i){
           for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
             exc_(k)=dst(n,k,j,i);
           }
           for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
             int k_shift = k;
             k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
             dst(n,k,j,i)=exc_(k_shift);
           }
         }
        }
      }
    }

    if(pmb->block_bcs[OUTER_X2]==POLAR_BNDRY){
      int nx3_half = (pmb->ke - pmb->ks + 1) / 2;
      for (int n=ns; n<=ne; ++n) {
        for (int j=pmb->je+1; j<=pmb->je+NGHOST; ++j) {
         for (int i=pmb->is-NGHOST; i<=pmb->ie+NGHOST; ++i){
           for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
             exc_(k)=dst(n,k,j,i);
           }
           for (int k=pmb->ks-NGHOST; k<=pmb->ke+NGHOST; ++k) {
             int k_shift = k;
             k_shift += (k < (nx3_half+NGHOST) ? 1 : -1) * nx3_half;
             dst(n,k,j,i)=exc_(k_shift);
           }
         }
        }
      }
    }
  }
  return;
}

#ifdef INCLUDE_CHEMISTRY
//----------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadSpeciesBoundaryBufferSameLevel(AthenaArray<Real> &src,
//                                                 Real *buf, const NeighborBlock& nb)
//  \brief Set species boundary buffers for sending to a block on the same level

int BoundaryValues::LoadSpeciesBoundaryBufferSameLevel(AthenaArray<Real> &src,
    Real *buf, const NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_block_;
  int si, sj, sk, ei, ej, ek;

  si=(nb.ox1>0)?(pmb->ie-NGHOST+1):pmb->is;
  ei=(nb.ox1<0)?(pmb->is+NGHOST-1):pmb->ie;
  sj=(nb.ox2>0)?(pmb->je-NGHOST+1):pmb->js;
  ej=(nb.ox2<0)?(pmb->js+NGHOST-1):pmb->je;
  sk=(nb.ox3>0)?(pmb->ke-NGHOST+1):pmb->ks;
  ek=(nb.ox3<0)?(pmb->ks+NGHOST-1):pmb->ke;
  int p=0;
  BufferUtility::Pack4DData(src, buf, 0, NSPECIES-1, si, ei, sj, ej, sk, ek, p);
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendSpeciesBoundaryBuffers(AthenaArray<Real> &src)
//  \brief Send boundary buffers

void BoundaryValues::SendSpeciesBoundaryBuffers(AthenaArray<Real> &src)
{
  std::stringstream msg;
  MeshBlock *pmb=pmy_block_;
  int mylevel=pmb->loc.level;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    int ssize;
    if(nb.level==mylevel) {
      ssize=LoadSpeciesBoundaryBufferSameLevel(src, species_send_[nb.bufid],nb);
    } else {
      //need to add AMR here
      msg << "### FATAL ERROR in SendSpeciesBoundaryBuffers()" << std::endl
          << "Chemistry BC not yet working with AMR yet." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }

    if(nb.rank == Globals::my_rank) { // on the same process
      MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(nb.gid);
      std::memcpy(pbl->pbval->species_recv_[nb.targetid],
                  species_send_[nb.bufid], ssize*sizeof(Real));
      pbl->pbval->species_flag_[nb.targetid]=BNDRY_ARRIVED;
    }
#ifdef MPI_PARALLEL
    else // MPI
      MPI_Start(&req_species_send_[nb.bufid]);
#endif
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetSpeciesBoundarySameLevel(AthenaArray<Real> &dst,
//                                           Real *buf, const NeighborBlock& nb)
//  \brief Set species boundary received from a block on the same level

void BoundaryValues::SetSpeciesBoundarySameLevel(AthenaArray<Real> &dst, Real *buf,
                                               const NeighborBlock& nb)
{
  std::stringstream msg;
  MeshBlock *pmb=pmy_block_;
  int si, sj, sk, ei, ej, ek;

  if(nb.ox1==0)     si=pmb->is,        ei=pmb->ie;
  else if(nb.ox1>0) si=pmb->ie+1,      ei=pmb->ie+NGHOST;
  else              si=pmb->is-NGHOST, ei=pmb->is-1;
  if(nb.ox2==0)     sj=pmb->js,        ej=pmb->je;
  else if(nb.ox2>0) sj=pmb->je+1,      ej=pmb->je+NGHOST;
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  if(nb.ox3==0)     sk=pmb->ks,        ek=pmb->ke;
  else if(nb.ox3>0) sk=pmb->ke+1,      ek=pmb->ke+NGHOST;
  else              sk=pmb->ks-NGHOST, ek=pmb->ks-1;

  int p=0;
  if (nb.polar) {
      //need to add polar cooridnates here
      msg << "### FATAL ERROR in SetSpeciesBoundarySameLevel()" << std::endl
          << "Chemistry BC not yet working with polar coordinates yet." << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  else {
    BufferUtility::Unpack4DData(buf, dst, 0, NSPECIES-1, si, ei, sj, ej, sk, ek, p);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveSpeciesBoundaryBuffers(AthenaArray<Real> &dst)
//  \brief receive the boundary data

bool BoundaryValues::ReceiveSpeciesBoundaryBuffers(AthenaArray<Real> &dst)
{
  std::stringstream msg;
  MeshBlock *pmb=pmy_block_;
  bool flag=true;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    if(species_flag_[nb.bufid]==BNDRY_COMPLETED) continue;
    if(species_flag_[nb.bufid]==BNDRY_WAITING) {
      if(nb.rank==Globals::my_rank) {// on the same process
        flag=false;
        continue;
      }
#ifdef MPI_PARALLEL
      else { // MPI boundary
        int test;
        MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
        MPI_Test(&req_species_recv_[nb.bufid],&test,MPI_STATUS_IGNORE);
        if(test==false) {
          flag=false;
          continue;
        }
        species_flag_[nb.bufid] = BNDRY_ARRIVED;
      }
#endif
    }
    if(nb.level==pmb->loc.level) {
      SetSpeciesBoundarySameLevel(dst, species_recv_[nb.bufid], nb);
    } else {
      //need to add AMR here
      msg << "### FATAL ERROR in ReceiveSpeciesBoundaryBuffers()" << std::endl
          << "Chemistry BC not yet working with AMR yet." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    species_flag_[nb.bufid] = BNDRY_COMPLETED; // completed
  }

  if(flag&& (pmb->block_bcs[INNER_X2]==POLAR_BNDRY
         ||  pmb->block_bcs[OUTER_X2]==POLAR_BNDRY)) {
      //need to add polar cooridnates here
      msg << "### FATAL ERROR in ReceiveSpeciesBoundaryBuffers()" << std::endl
          << "Chemistry BC not yet working with polar coordinates yet." << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  return flag;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ReceiveSpeciesBoundaryBuffersWithWait(AthenaArray<Real> &dst)
//  \brief receive the boundary data for initialization

void BoundaryValues::ReceiveSpeciesBoundaryBuffersWithWait(AthenaArray<Real> &dst)
{
  std::stringstream msg;
  MeshBlock *pmb=pmy_block_;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
#ifdef MPI_PARALLEL
    if(nb.rank!=Globals::my_rank)
      MPI_Wait(&req_species_recv_[nb.bufid],MPI_STATUS_IGNORE);
#endif
    if(nb.level==pmb->loc.level) {
      SetSpeciesBoundarySameLevel(dst, species_recv_[nb.bufid], nb);
    } else {
      //need to add AMR here
      msg << "### FATAL ERROR in ReceiveSpeciesBoundaryBuffers()" << std::endl
          << "Chemistry BC not yet working with AMR yet." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    species_flag_[nb.bufid] = BNDRY_COMPLETED; // completed
  }
 
  if (pmb->block_bcs[INNER_X2]==POLAR_BNDRY||pmb->block_bcs[OUTER_X2]==POLAR_BNDRY) {
      //need to add polar cooridnates here
      msg << "### FATAL ERROR in ReceiveSpeciesBoundaryBuffers()" << std::endl
          << "Chemistry BC not yet working with polar coordinates yet." << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }

  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::StartReceivingSpecies(void)
//  \brief initiate MPI_Irecv for species

void BoundaryValues::StartReceivingSpecies(void)
{
  firsttime_=true;
#ifdef MPI_PARALLEL
  MeshBlock *pmb=pmy_block_;
  int mylevel=pmb->loc.level;
  for(int n=0;n<pmb->nneighbor;n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    if(nb.rank!=Globals::my_rank) { 
      MPI_Start(&req_species_recv_[nb.bufid]);
    }
  }
#endif
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ClearBoundarySpecies(void)
//  \brief clean up the boundary flags after each loop for species

void BoundaryValues::ClearBoundarySpecies(void)
{
  MeshBlock *pmb=pmy_block_;

  // Clear non-polar boundary communications
  for(int n=0;n<pmb->nneighbor;n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    species_flag_[nb.bufid] = BNDRY_WAITING;
#ifdef MPI_PARALLEL
    if(nb.rank!=Globals::my_rank) {
      MPI_Wait(&req_species_send_[nb.bufid],MPI_STATUS_IGNORE); // Wait for Isend
    }
#endif
  }
  return;
}

//six-ray
//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::StartReceivingSixray(void)
//  \brief initiate MPI_Irecv for six-ray column

void BoundaryValues::StartReceivingSixray(void)
{
  firsttime_=true;
#ifdef MPI_PARALLEL
  MeshBlock *pmb=pmy_block_;
  int mylevel=pmb->loc.level;
  for(int n=0;n<6;n++) {
    NeighborBlock *nb = pmb->prad->pradintegrator->pfacenb_[n];
    if(nb != NULL && nb->rank!=Globals::my_rank) { 
      MPI_Start(&req_sixray_recv_[n]);
    }
  }
#endif
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ClearBoundarySixray(void)
//  \brief clean up the boundary flags after each loop for sixray

void BoundaryValues::ClearBoundarySixray(void)
{
  MeshBlock *pmb=pmy_block_;

  // Clear non-polar boundary communications
  for(int n=0;n<6;n++) {
    NeighborBlock *nb = pmb->prad->pradintegrator->pfacenb_[n];
    sixray_flag_[n] = BNDRY_WAITING;
#ifdef MPI_PARALLEL
    if(nb != NULL && nb->rank!=Globals::my_rank) {
      MPI_Wait(&req_sixray_send_[n],MPI_STATUS_IGNORE); // Wait for Isend
    }
#endif
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadSixrayBoundaryBufferSameLevel(AthenaArray<Real> &src,
//                                                 Real *buf, const NeighborBlock& nb,
//                                                 const int direction)
//  \brief Set six-ray boundary buffers for sending to a block on the same level

int BoundaryValues::LoadSixrayBoundaryBufferSameLevel(AthenaArray<Real> &src,
    Real *buf, const NeighborBlock& nb, const int direction)
{
  MeshBlock *pmb=pmy_block_;
  int si, sj, sk, ei, ej, ek;

  si=(nb.ox1>0)?(pmb->ie):pmb->is;
  ei=(nb.ox1<0)?(pmb->is):pmb->ie;
  sj=(nb.ox2>0)?(pmb->je):pmb->js;
  ej=(nb.ox2<0)?(pmb->js):pmb->je;
  sk=(nb.ox3>0)?(pmb->ke):pmb->ks;
  ek=(nb.ox3<0)?(pmb->ks):pmb->ke;
  int p=0;
  BufferUtility::Pack5DData(src, buf, direction, direction, 0, pmb->pspec->pchemnet->n_cols_-1,
                            si, ei, sj, ej, sk, ek, p);
  return p;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetSixrayBoundarySameLevel(AthenaArray<Real> &dst,
//                                           Real *buf, const NeighborBlock& nb,
//                                           const int direction)
//  \brief Set six-ray boundary received from a block on the same level

void BoundaryValues::SetSixrayBoundarySameLevel(AthenaArray<Real> &dst, Real *buf,
                                               const NeighborBlock& nb,
                                               const int direction)
{
  std::stringstream msg;
  MeshBlock *pmb=pmy_block_;
  int si, sj, sk, ei, ej, ek;

  if(nb.ox1==0)     si=pmb->is,        ei=pmb->ie;
  else if(nb.ox1>0) si=pmb->ie+1,      ei=pmb->ie+1;
  else              si=pmb->is-1,      ei=pmb->is-1;
  if(nb.ox2==0)     sj=pmb->js,        ej=pmb->je;
  else if(nb.ox2>0) sj=pmb->je+1,      ej=pmb->je+1;
  else              sj=pmb->js-1,      ej=pmb->js-1;
  if(nb.ox3==0)     sk=pmb->ks,        ek=pmb->ke;
  else if(nb.ox3>0) sk=pmb->ke+1,      ek=pmb->ke+1;
  else              sk=pmb->ks-1,      ek=pmb->ks-1;

  int p=0;
  if (nb.polar) {
      //need to add polar cooridnates here
      msg << "### FATAL ERROR in SetSixrayBoundarySameLevel()" << std::endl
          << "Chemistry BC not yet working with polar coordinates yet." << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  else {
    BufferUtility::Unpack5DData(buf, dst, direction, direction, 0, pmb->pspec->pchemnet->n_cols_-1, 
                                si, ei, sj, ej, sk, ek, p);
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendSixrayBoundaryBuffers(AthenaArray<Real> &src,
//const int direction)
//  \brief Send boundary buffers

void BoundaryValues::SendSixrayBoundaryBuffers(AthenaArray<Real> &src, 
                                               const int direction)
{
  std::stringstream msg;
  MeshBlock *pmb=pmy_block_;
  int mylevel=pmb->loc.level;

  NeighborBlock *nb = pmb->prad->pradintegrator->pfacenb_[direction];

  int recv_direction = CGKUtility::GetOppositeDirection(direction);

  int ssize;
  if(nb->level==mylevel) {
    ssize=LoadSixrayBoundaryBufferSameLevel(src, col_sixray_send_[direction],*nb,
                                            recv_direction);
  } else {
    //need to add AMR here
    msg << "### FATAL ERROR in SendSixrayBoundaryBuffers()" << std::endl
      << "Chemistry BC not yet working with AMR yet." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  if(nb->rank == Globals::my_rank) { // on the same process
    MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(nb->gid);
    std::memcpy(pbl->pbval->col_sixray_recv_[recv_direction],
        col_sixray_send_[direction], ssize*sizeof(Real));
    pbl->pbval->sixray_flag_[recv_direction]=BNDRY_ARRIVED;
  }
#ifdef MPI_PARALLEL
  else // MPI
    MPI_Start(&req_sixray_send_[direction]);
#endif

  return;
}

//----------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveSixrayBoundaryBuffers(AthenaArray<Real> &dst,
//const int direction)
//  \brief receive the boundary data
//  return true if boundary already completed, or just received and completed. 
//  Return fasle if boundary has not arrived yet, or is at the same processor.

bool BoundaryValues::ReceiveSixrayBoundaryBuffers(AthenaArray<Real> &dst,
                                                  const int direction)
{
  std::stringstream msg;
  MeshBlock *pmb=pmy_block_;
  bool flag=true;

  NeighborBlock *nb = pmb->prad->pradintegrator->pfacenb_[direction];
  if(sixray_flag_[direction]==BNDRY_COMPLETED) {
    flag = true;  
  } else if (sixray_flag_[direction]==BNDRY_WAITING) {
    if(nb->rank==Globals::my_rank) {// on the same process
      flag=false;
    }
#ifdef MPI_PARALLEL
    else { // MPI boundary
      int test;
      MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,MPI_STATUS_IGNORE);
      MPI_Test(&req_sixray_recv_[direction],&test,MPI_STATUS_IGNORE);
      if(test==false) {
        flag=false;
      } else {
        sixray_flag_[direction] = BNDRY_ARRIVED;
      }
    }
#endif
  } 

  if (sixray_flag_[direction] == BNDRY_ARRIVED) {
    if(nb->level==pmb->loc.level) {
      SetSixrayBoundarySameLevel(dst, col_sixray_recv_[direction], *nb, direction);
    } else {
      //need to add AMR here
      msg << "### FATAL ERROR in ReceiveSixrayBoundaryBuffers()" << std::endl
        << "Chemistry BC not yet working with AMR yet." << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    sixray_flag_[direction] = BNDRY_COMPLETED; // completed
  }

  if(flag&& (pmb->block_bcs[INNER_X2]==POLAR_BNDRY
         ||  pmb->block_bcs[OUTER_X2]==POLAR_BNDRY)) {
      //need to add polar cooridnates here
      msg << "### FATAL ERROR in ReceiveSixrayBoundaryBuffers()" << std::endl
          << "Chemistry BC not yet working with polar coordinates yet." << std::endl;
      throw std::runtime_error(msg.str().c_str());
  }
  return flag;
}
#endif //INCLUDE_CHEMISTRY
