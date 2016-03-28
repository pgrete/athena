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
//! \file cell_center_bvals.cpp
//  \brief implements functions that initialize/apply BCs on for cell-center variables
//======================================================================================

// C++ headers
#include <string>     // c_str()
#include <cstring>    // memcpy
#include <cstdlib>
#include <cmath>

// Athena++ classes headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../athena_arrays.hpp"
#include "../mesh_refinement/mesh_refinement.hpp"
#include "../mesh.hpp"

#include "../radiation/radiation.hpp"
#include "../parameter_input.hpp"
#include "../utils/buffer_utils.hpp"

// this class header
#include "bvals.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadCenterBoundaryBufferSameLevel(AthenaArray<Real> &src,
//                                        Real *buf, const NeighborBlock& nb, int phys)
//  \brief Set hydro boundary buffers for sending to a block on the same level
int BoundaryValues::LoadCenterBoundaryBufferSameLevel(AthenaArray<Real> &src,
                             Real *buf, const NeighborBlock& nb, int phys)
{
  MeshBlock *pmb=pmy_mblock_;
  
  int si, sj, sk, ei, ej, ek;

  // calculate the index
  si=(nb.ox1>0)?(pmb->ie-NGHOST+1):pmb->is;
  ei=(nb.ox1<0)?(pmb->is+NGHOST-1):pmb->ie;
  sj=(nb.ox2>0)?(pmb->je-NGHOST+1):pmb->js;
  ej=(nb.ox2<0)?(pmb->js+NGHOST-1):pmb->je;
  sk=(nb.ox3>0)?(pmb->ke-NGHOST+1):pmb->ks;
  ek=(nb.ox3<0)?(pmb->ks+NGHOST-1):pmb->ke;
  int p=0;
  if(phys == HYDRO){
    BufferUtility::Pack4DData(src,buf, 0, NHYDRO-1,
                             sk, ek, sj, ej, si, ei, p);
  }else if(phys == RAD){
    BufferUtility::Pack4DData(src,buf, sk, ek, sj, ej,
                             si, ei, 0, pmb->prad->n_fre_ang-1, p);
  }
  
  return p;
}


//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadCenterBoundaryBufferToCoarser(AthenaArray<Real> &src,
//                                                 Real *buf, const NeighborBlock& nb,
//                                                 bool conserved_values)
//  \brief Set hydro boundary buffers for sending to a block on the coarser level
int BoundaryValues::LoadCenterBoundaryBufferToCoarser(AthenaArray<Real> &src, Real *buf,
                      const NeighborBlock& nb, int phys, bool conserved_values)
{
  MeshBlock *pmb=pmy_mblock_;
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
  // restrict the data before sending
  if(phys==HYDRO){
    if(conserved_values){
      pmr->RestrictCellCenteredValues(src, pmr->coarse_cons_, HYDRO,
                                    0, NHYDRO-1,sk, ek, sj, ej, si, ei);
      BufferUtility::Pack4DData(pmr->coarse_cons_, buf,
                            0, NHYDRO-1, sk, ek, sj, ej, si, ei, p);
    }else{
    // must be initialization; need to restrict but not send primitives
      pmr->RestrictCellCenteredValues(src, pmr->coarse_prim_, HYDRO,
                                    0, NHYDRO-1,sk, ek, sj, ej, si, ei);
    }
  }else if(phys==RAD){
     pmr->RestrictCellCenteredValues(src, pmr->coarse_ir_, RAD,
                            sk, ek, sj, ej, si, ei, 0, pmb->prad->n_fre_ang-1);
     BufferUtility::Pack4DData(pmr->coarse_ir_, buf,
                sk, ek, sj, ej, si, ei, 0, pmb->prad->n_fre_ang-1, p);
  }

  return p;
}


//--------------------------------------------------------------------------------------
//! \fn int BoundaryValues::LoadCenterBoundaryBufferToFiner(AthenaArray<Real> &src,
//                                                 Real *buf, const NeighborBlock& nb)
//  \brief Set hydro boundary buffers for sending to a block on the finer level
int BoundaryValues::LoadCenterBoundaryBufferToFiner(AthenaArray<Real> &src, Real *buf,
                                    const NeighborBlock& nb, int phys)
{
  MeshBlock *pmb=pmy_mblock_;
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
    } else {
      if(nb.fi2==1) sj+=pmb->block_size.nx2/2-pmb->cnghost;
      else          ej-=pmb->block_size.nx2/2-pmb->cnghost;
    }
  }
  if(nb.ox3==0 && pmb->block_size.nx3 > 1) {
    if(nb.ox1!=0 && nb.ox2!=0) {
      if(nb.fi1==1) sk+=pmb->block_size.nx3/2-pmb->cnghost;
         else          ek-=pmb->block_size.nx3/2-pmb->cnghost;
      } else {
        if(nb.fi2==1) sk+=pmb->block_size.nx3/2-pmb->cnghost;
        else          ek-=pmb->block_size.nx3/2-pmb->cnghost;
      }
  }
  int p=0;
  if(phys==HYDRO){
      BufferUtility::Pack4DData(src, buf, 0,
                    NHYDRO-1, sk, ek, sj, ej, si, ei, p);
      
  }else if(phys==RAD){
      BufferUtility::Pack4DData(src, buf, sk, ek,
                   sj, ej, si, ei, 0, pmb->prad->n_fre_ang-1, p);
  }

  return p;
}





//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SendHydroBoundaryBuffers(AthenaArray<Real> &src, int step,
//                                                    bool conserved_values)
//  \brief Send boundary buffers
void BoundaryValues::SendCenterBoundaryBuffers(AthenaArray<Real> &src, int phys,
                                                int step, bool conserved_values)
{
  MeshBlock *pmb=pmy_mblock_;
  int mylevel=pmb->loc.level;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    int ssize;
    if(nb.level==mylevel){
      if(phys == HYDRO)
       ssize=LoadCenterBoundaryBufferSameLevel(src, hydro_send_[step][nb.bufid],
                                               nb, phys);
      else if(phys == RAD)
       ssize=LoadCenterBoundaryBufferSameLevel(src, rad_send_[step][nb.bufid],
                                               nb, phys);
    }
    else if(nb.level<mylevel){
      if(phys == HYDRO)
        ssize=LoadCenterBoundaryBufferToCoarser(src, hydro_send_[step][nb.bufid],
                                               nb, phys, conserved_values);
      else if(phys == RAD)
        ssize=LoadCenterBoundaryBufferToCoarser(src, rad_send_[step][nb.bufid],
                                               nb, phys, conserved_values);
    }else{
      if(phys==HYDRO)
        ssize=LoadCenterBoundaryBufferToFiner(src, hydro_send_[step][nb.bufid],
                                               nb, phys);
      else if(phys==RAD)
        ssize=LoadCenterBoundaryBufferToFiner(src, rad_send_[step][nb.bufid],
                                               nb, phys);
    
    }
    if(nb.rank == Globals::my_rank) { // on the same process
      MeshBlock *pbl=pmb->pmy_mesh->FindMeshBlock(nb.gid);
      if(phys==HYDRO){
        std::memcpy(pbl->pbval->hydro_recv_[step][nb.targetid],
                  hydro_send_[step][nb.bufid], ssize*sizeof(Real));
        pbl->pbval->hydro_flag_[step][nb.targetid]=boundary_arrived;
      }else if(phys==RAD){
        std::memcpy(pbl->pbval->rad_recv_[step][nb.targetid],
                  rad_send_[step][nb.bufid], ssize*sizeof(Real));
        pbl->pbval->rad_flag_[step][nb.targetid]=boundary_arrived;
      }
    }
#ifdef MPI_PARALLEL
    else{ // MPI
      if(phys==HYDRO)
        MPI_Start(&req_hydro_send_[step][nb.bufid]);
      else if(phys==RAD)
        MPI_Start(&req_rad_send_[step][nb.bufid]);
    }
#endif
  }

  return;
}




//--------------------------------------------------------------------------------------
//! \fn bool BoundaryValues::ReceiveCenterBoundaryBuffers(AthenaArray<Real> &dst,
//                        int phys, int step)
//  \brief receive the boundary data
bool BoundaryValues::ReceiveCenterBoundaryBuffers(AthenaArray<Real> &dst, int phys,
                                                  int step)
{
  MeshBlock *pmb=pmy_mblock_;
  enum boundary_status *status_flag;
  Real *recv_buf;
  
  
  bool flag=true;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    if(phys==HYDRO) {
       status_flag = &(hydro_flag_[step][nb.bufid]);
       recv_buf = hydro_recv_[step][nb.bufid];
    }else if(phys==RAD){
       status_flag = &(rad_flag_[step][nb.bufid]);
       recv_buf = rad_recv_[step][nb.bufid];
    }

    if(*status_flag==boundary_completed) continue;
    if(*status_flag==boundary_waiting) {
       if(nb.rank==Globals::my_rank) {// on the same process
         flag=false;
         continue;
       }
#ifdef MPI_PARALLEL
       else { // MPI boundary
         int test;
         MPI_Iprobe(MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&test,
                     MPI_STATUS_IGNORE);
         if(phys==HYDRO)
           MPI_Test(&req_hydro_recv_[step][nb.bufid],&test,MPI_STATUS_IGNORE);
         else if(phys==RAD)
           MPI_Test(&req_rad_recv_[step][nb.bufid],&test,MPI_STATUS_IGNORE);
         
         if(test==false) {
           flag=false;
           continue;
         }
         *status_flag = boundary_arrived;
       }
#endif
    }
    
    
    
    if(nb.level==pmb->loc.level)
      SetCenterBoundarySameLevel(dst,recv_buf,phys,nb);
    else if(nb.level<pmb->loc.level) // this set only the prolongation buffer
      SetCenterBoundaryFromCoarser(recv_buf,phys,nb, true);
    else
      SetCenterBoundaryFromFiner(dst, recv_buf, phys, nb);
      
      
    *status_flag = boundary_completed; // completed
    
  }

  if((flag&&(phys==HYDRO))&&(pmb->block_bcs[INNER_X2]==POLAR_BNDRY||
       pmb->block_bcs[OUTER_X2]==POLAR_BNDRY))
       PolarSingleHydro(dst);
  

  if((flag&&(phys==RAD))&&(pmb->block_rad_bcs[INNER_X2]==POLAR_BNDRY||
       pmb->block_rad_bcs[OUTER_X2]==POLAR_BNDRY))
       PolarSingleRad(dst);
  
  
  return flag;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::ReceiveCenterBoundaryBuffersWithWait(
//                              AthenaArray<Real> &dst,int phys, int step)
//  \brief receive the boundary data for initialization
void BoundaryValues::ReceiveCenterBoundaryBuffersWithWait(AthenaArray<Real> &dst,
                                               int phys, int step)
{
  MeshBlock *pmb=pmy_mblock_;
  enum boundary_status *status_flag;
  Real *recv_buf;

  for(int n=0; n<pmb->nneighbor; n++) {
    NeighborBlock& nb = pmb->neighbor[n];
    if(phys==HYDRO) {
      status_flag = &(hydro_flag_[step][nb.bufid]);
      recv_buf = hydro_recv_[step][nb.bufid];
    }else if(phys==RAD){
      status_flag = &(rad_flag_[step][nb.bufid]);
      recv_buf = rad_recv_[step][nb.bufid];
    }
#ifdef MPI_PARALLEL
    if(nb.rank!=Globals::my_rank){
      if(phys == HYDRO)
        MPI_Wait(&req_hydro_recv_[step][nb.bufid],MPI_STATUS_IGNORE);
      else if(phys==RAD)
        MPI_Wait(&req_rad_recv_[step][nb.bufid],MPI_STATUS_IGNORE);
    }
#endif
    if(nb.level==pmb->loc.level)
      SetCenterBoundarySameLevel(dst, recv_buf, phys, nb);
    else if(nb.level<pmb->loc.level)
      SetCenterBoundaryFromCoarser(recv_buf, phys, nb, step == 0);
    else
      SetCenterBoundaryFromFiner(dst, recv_buf, phys, nb);
    
    *status_flag = boundary_completed; // completed
  }
 
  if ((phys==HYDRO)&&(pmb->block_bcs[INNER_X2]==POLAR_BNDRY
        ||pmb->block_bcs[OUTER_X2]==POLAR_BNDRY))
        PolarSingleHydro(dst);
  
  if((phys==RAD)&&(pmb->block_rad_bcs[INNER_X2]==POLAR_BNDRY||
       pmb->block_rad_bcs[OUTER_X2]==POLAR_BNDRY))
       PolarSingleRad(dst);


  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetCenterBoundarySameLevel(AthenaArray<Real> &dst,
//                            Real *buf, int phys, const NeighborBlock& nb)
//  \brief Set cell center boundary received from a block on the same level
void BoundaryValues::SetCenterBoundarySameLevel(AthenaArray<Real> &dst, Real *buf,
                                    int phys, const NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
  int si, sj, sk, ei, ej, ek;
  int s1, e1, s2, e2, s3, e3, s4, e4;

  if(nb.ox1==0)     si=pmb->is,        ei=pmb->ie;
  else if(nb.ox1>0) si=pmb->ie+1,      ei=pmb->ie+NGHOST;
  else              si=pmb->is-NGHOST, ei=pmb->is-1;
  if(nb.ox2==0)     sj=pmb->js,        ej=pmb->je;
  else if(nb.ox2>0) sj=pmb->je+1,      ej=pmb->je+NGHOST;
  else              sj=pmb->js-NGHOST, ej=pmb->js-1;
  if(nb.ox3==0)     sk=pmb->ks,        ek=pmb->ke;
  else if(nb.ox3>0) sk=pmb->ke+1,      ek=pmb->ke+NGHOST;
  else              sk=pmb->ks-NGHOST, ek=pmb->ks-1;
  
  if(phys==HYDRO){
    s1=si; e1=ei; s2=sj; e2=ej; s3=sk; e3=ek; s4=0; e4=(NHYDRO-1);
  }else if(phys==RAD){
    s1=0; e1=(pmb->prad->n_fre_ang-1); s2=si; e2=ei; s3=sj; e3=ej; s4=sk; e4=ek;
  }

  int p=0;
  if (nb.polar) {
    if(phys==HYDRO){
      for (int n=s4; n<=e4; ++n) {
       Real sign = flip_across_pole_hydro[n] ? -1.0 : 1.0;
       for (int k=s3; k<=e3; ++k) {
         for (int j=e2; j>=s2; --j) {
#pragma simd
          for (int i=s1; i<=e1; ++i)
            dst(n,k,j,i) = sign * buf[p++];
         }
       }
     }
    }else if(phys==RAD){
      //Need need to flip radiation
      // The miu_theta, miu_phi flip automatically based on the angles
      for (int k=s4; k<=e4; ++k) {
       for (int j=e3; j>=s3; --j) {
         for (int i=s2; i<=e2; ++i) {
#pragma simd
           for(int n=s1; n<=e1; ++n){
               dst(k,j,i,n) = buf[p++];

           }// end ifr
         }
       }
      }
    
    }// end phys=Rad
  }
  else
    BufferUtility::Unpack4DData(buf, dst, s4, e4, s3, e3, s2, e2, s1, e1, p);
  return;
}


//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetCenterBoundaryFromCoarser(Real *buf,
//                                                       const NeighborBlock& nb)
//  \brief Set cell center prolongation buffer received from a block on a coarser level
void BoundaryValues::SetCenterBoundaryFromCoarser(Real *buf, int phys,
                      const NeighborBlock& nb, bool conserved_values)
{
  MeshBlock *pmb=pmy_mblock_;
  MeshRefinement *pmr=pmb->pmr;

  int si, sj, sk, ei, ej, ek;
  int s1, e1, s2, e2, s3, e3, s4, e4;
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
  
  if(phys==HYDRO){
    s1=si; e1=ei; s2=sj; e2=ej; s3=sk; e3=ek; s4=0; e4=(NHYDRO-1);
  }else if(phys==RAD){
    s1=0; e1=(pmb->prad->n_fre_ang-1); s2=si; e2=ei; s3=sj; e3=ej; s4=sk; e4=ek;
  }


  int p=0;
  if (nb.polar) {
    if(phys==HYDRO){
     for (int n=s4; n<=e4; ++n) {
      Real sign = flip_across_pole_hydro[n] ? -1.0 : 1.0;
      for (int k=s3; k<=e3; ++k) {
        for (int j=e2; j>=s2; --j) {
#pragma simd
          for (int i=s1; i<=e1; ++i){
            if (conserved_values)
              pmr->coarse_cons_(n,k,j,i) = sign * buf[p++];
            else
              pmr->coarse_prim_(n,k,j,i) = sign * buf[p++];
          }
        }
      }
     }
    }else if(phys==RAD){
      for (int k=s4; k<=e4; ++k) {
       for (int j=e3; j>=s3; --j) {
         for (int i=s2; i<=e2; ++i) {
           for(int n=s1; n<=e1; ++n){
              pmr->coarse_ir_(k,j,i,n) = buf[p++];
           }// end nang
         }
       }
      }
    }// End phys==RAD
  }
  else{
    if(phys==HYDRO){
      if (conserved_values)
        BufferUtility::Unpack4DData(buf, pmr->coarse_cons_, s4, e4, s3,
                   e3, s2, e2, s1, e1, p);
      else
        BufferUtility::Unpack4DData(buf, pmr->coarse_prim_, s4, e4, s3,
                   e3, s2, e2, s1, e1, p);
    }
    else if(phys==RAD)
      BufferUtility::Unpack4DData(buf, pmr->coarse_ir_, s4, e4, s3,
                   e3, s2, e2, s1, e1, p);
  }
  return;
}



//--------------------------------------------------------------------------------------
//! \fn void BoundaryValues::SetCenterBoundaryFromFiner(AthenaArray<Real> &dst,
//                               Real *buf, int phys, const NeighborBlock& nb)
//  \brief Set hydro boundary received from a block on a finer level
void BoundaryValues::SetCenterBoundaryFromFiner(AthenaArray<Real> &dst, Real *buf,
                                      int phys, const NeighborBlock& nb)
{
  MeshBlock *pmb=pmy_mblock_;
  // receive already restricted data
  int si, sj, sk, ei, ej, ek;
  int s1, e1, s2, e2, s3, e3, s4, e4;
  
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
  
  if(phys==HYDRO){
    s1=si; e1=ei; s2=sj; e2=ej; s3=sk; e3=ek; s4=0; e4=(NHYDRO-1);
  }else if(phys==RAD){
    s1=0; e1=(pmb->prad->n_fre_ang-1); s2=si; e2=ei; s3=sj; e3=ej; s4=sk; e4=ek;
  }


  int p=0;
  if (nb.polar) {
   if(phys==HYDRO){
    for (int n=s4; n<=e4; ++n) {
      Real sign = flip_across_pole_hydro[n] ? -1.0 : 1.0;
      for (int k=s3; k<=e3; ++k) {
        for (int j=s2; j<=e2; ++j) {
#pragma simd
          for (int i=s1; i<=e1; ++i)
            dst(n,k,j,i) = sign * buf[p++];
        }
      }
    }
   }else if(phys==RAD){
      for (int k=s4; k<=e4; ++k) {
       for (int j=s3; j<=e3; ++j) {
         for (int i=s2; i<=e2; ++i) {
#pragma simd
           for(int n=s1; n<=e1; ++n){
             dst(k,j,i,n) = buf[p++];
           }// end ifr
         }
       }
      }
   }//end phys==RAD
  }
  else 
    BufferUtility::Unpack4DData(buf, dst, s4, e4, s3, e3, s2, e2, s1, e1, p);
  return;
}


