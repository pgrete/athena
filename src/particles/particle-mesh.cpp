//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file particle-mesh.cpp
//  \brief implements ParticleMesh class used for operations involved in particle-mesh 
//         methods.

// Standard library
#include <cstring>
#include <sstream>

// Athena++ classes headers
#include "../athena.hpp"
#include "../globals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../utils/buffer_utils.hpp"
#include "particles.hpp"

// Class variable initialization
bool ParticleMesh::initialized_ = false;
int ParticleMesh::nmeshaux = 0;
int ParticleMesh::iweight = -1;
#ifdef MPI_PARALLEL
MPI_Comm ParticleMesh::my_comm = MPI_COMM_NULL;
#endif

// Local function prototypes.
static Real _WeightFunction(Real dxi);

// Local constants.
const int OFFSET = NGHOST - NGPM;  // offset between meshblock and meshaux

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::Initialize(ParameterInput *pin)
//  \brief initiates the ParticleMesh class.

void ParticleMesh::Initialize(ParameterInput *pin)
{
  if (initialized_) return;

  // Add weight in meshaux.
  iweight = AddMeshAux();

#ifdef MPI_PARALLEL
  // Get my MPI communicator.
  MPI_Comm_dup(MPI_COMM_WORLD, &my_comm);
#endif

  initialized_ = true;
}

//--------------------------------------------------------------------------------------
//! \fn int ParticleMesh::AddMeshAux()
//  \brief adds one auxiliary to the mesh and returns the index.

int ParticleMesh::AddMeshAux()
{
  return nmeshaux++;
}

//--------------------------------------------------------------------------------------
//! \fn ParticleMesh::ParticleMesh(Particles *ppar, int nmeshaux)
//  \brief constructs a new ParticleMesh instance.

ParticleMesh::ParticleMesh(Particles *ppar)
{
  // Save some inputs.
  ppar_ = ppar;
  pmb_ = ppar->pmy_block;
  pmesh_ = pmb_->pmy_mesh;
  pbval_ = pmb_->pbval;

  // Determine active dimensions.
  RegionSize& block_size = pmb_->block_size;
  active1_ = block_size.nx1 > 1;
  active2_ = block_size.nx2 > 1;
  active3_ = block_size.nx3 > 1;

  // Determine the range of a particle cloud.
  dxi1_ = active1_ ? RINF : 0;
  dxi2_ = active2_ ? RINF : 0;
  dxi3_ = active3_ ? RINF : 0;

  // Determine the dimensions of the block for boundary communication.
  int dim = 0;
  nx1_ = nx2_ = nx3_ = 1;

  if (active1_) {
    ++dim;
    is = NGPM;
    ie = NGPM + block_size.nx1 - 1;
    nx1_ = block_size.nx1 + 2 * NGPM;
  } else
    is = ie = 0;

  if (active2_) {
    ++dim;
    js = NGPM;
    je = NGPM + block_size.nx2 - 1;
    nx2_ = block_size.nx2 + 2 * NGPM;
  } else
    js = je = 0;

  if (active3_) {
    ++dim;
    ks = NGPM;
    ke = NGPM + block_size.nx3 - 1;
    nx3_ = block_size.nx3 + 2 * NGPM;
  } else
    ks = ke = 0;

  // Allocate the block for particle-mesh.
  meshaux.NewAthenaArray(nmeshaux, nx3_, nx2_, nx1_);
  ncells_ = nx1_ * nx2_ * nx3_;

  // Get a shorthand to weights.
  weight.InitWithShallowSlice(meshaux, 4, iweight, 1);

  // Initialize boundary data.
  bd_.nbmax = 56;
  for (int n = 0; n < bd_.nbmax; n++) {
    bd_.flag[n] = BNDRY_WAITING;
    bd_.send[n] = NULL;
    bd_.recv[n] = NULL;
#ifdef MPI_PARALLEL
    bd_.req_recv[n] = MPI_REQUEST_NULL;
    bd_.req_send[n] = MPI_REQUEST_NULL;
#endif
  }
}

//--------------------------------------------------------------------------------------
//! \fn ParticleMesh::~ParticleMesh()
//  \brief destructs a ParticleMesh instance.

ParticleMesh::~ParticleMesh()
{
  // Destroy the particle meshblock.
  weight.DeleteAthenaArray();
  meshaux.DeleteAthenaArray();

  // Destroy boundary data.
  for (int n = 0; n < bd_.nbmax; n++) {
    if (bd_.recv[n] != NULL) delete [] bd_.recv[n];
#ifdef MPI_PARALLEL
    if (bd_.send[n] != NULL) delete [] bd_.send[n];
    if (bd_.req_recv[n] != MPI_REQUEST_NULL) MPI_Request_free(&bd_.req_recv[n]);
    if (bd_.req_send[n] != MPI_REQUEST_NULL) MPI_Request_free(&bd_.req_send[n]);
#endif
  }
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::InterpolateMeshToParticles(
//               const AthenaArray<Real>& meshsrc, int ms1, int ms2,
//               AthenaArray<Real>& par, int p1)
//  \brief interpolates meshsrc from property index ms1 to ms2 onto particle array par
//         (realprop, auxprop, or work in Particles class) from property index p1 and
//         up.

void ParticleMesh::InterpolateMeshToParticles(
         const AthenaArray<Real>& meshsrc, int ms1, int ms2,
         AthenaArray<Real>& par, int p1)
{
  // Zero out the particle arrays.
  int nprop = ms2 - ms1 + 1;
  Real *pp[nprop];
  for (int n = 0; n < nprop; ++n) {
    Real *p = pp[n] = &par(p1+n,0);
    for (int k = 0; k < ppar_->npar; ++k)
      *p++ = 0.0;
  }

  // Loop over each particle.
  for (int k = 0; k < ppar_->npar; ++k) {
    // Find the domain the particle influences.
    Real xi1 = ppar_->xi1(k), xi2 = ppar_->xi2(k), xi3 = ppar_->xi3(k);
    int ix1s = int(xi1 - dxi1_), ix1e = int(xi1 + dxi1_);
    int ix2s = int(xi2 - dxi2_), ix2e = int(xi2 + dxi2_);
    int ix3s = int(xi3 - dxi3_), ix3e = int(xi3 + dxi3_);

    // Weight each cell and accumulate the mesh properties onto the particles.
    for (int ix3 = ix3s; ix3 <= ix3e; ++ix3) {
      Real w3 = active3_ ? _WeightFunction(ix3 + 0.5 - xi3) : 1.0;

      for (int ix2 = ix2s; ix2 <= ix2e; ++ix2) {
        Real w23 = w3 * (active2_ ? _WeightFunction(ix2 + 0.5 - xi2) : 1.0);

        for (int ix1 = ix1s; ix1 <= ix1e; ++ix1) {
          Real weight = w23 * (active1_ ?
                                    _WeightFunction(ix1 + 0.5 - xi1) : 1.0);

          for (int n = 0; n < nprop; ++n)
            *pp[n] += weight * meshsrc(ms1+n,ix3,ix2,ix1);
        }
      }
    }
    for (int n = 0; n < nprop; ++n)
      ++pp[n];
  }
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::AssignParticlesToMeshAux(
//               const AthenaArray<Real>& par, int p1, int p2, int ma1)
//  \brief assigns par (realprop, auxprop, or work in Particles class) from property
//         index p1 to p2 onto meshaux from property index ma1 and up.

void ParticleMesh::AssignParticlesToMeshAux(
         const AthenaArray<Real>& par, int p1, int p2, int ma1)
{
  // Zero out meshaux.
  int nprop = p2 - p1 + 1;
  Real *pmw0 = &weight(0,0,0), *pmw = pmw0;
  for (int i = 0; i < ncells_; ++i)
    *pmw++ = 0.0;

  Real *pm0[nprop];
  for (int n = 0; n < nprop; ++n) {
    Real *p = pm0[n] = &meshaux(ma1+n,0,0,0);
    for (int i = 0; i < ncells_; ++i)
      *p++ = 0.0;
  }

  // Loop over each particle.
  for (int k = 0; k < ppar_->npar; ++k) {
    // Find the domain the particle influences.
    Real xi1 = ppar_->xi1(k) - (active1_ ? OFFSET : 0),
         xi2 = ppar_->xi2(k) - (active2_ ? OFFSET : 0),
         xi3 = ppar_->xi3(k) - (active3_ ? OFFSET : 0);
    int ix1s = int(xi1 - dxi1_), ix1e = int(xi1 + dxi1_);
    int ix2s = int(xi2 - dxi2_), ix2e = int(xi2 + dxi2_);
    int ix3s = int(xi3 - dxi3_), ix3e = int(xi3 + dxi3_);

    // Prepare for pointer operations.
    int dpm1 = ix1s + nx1_ * (ix2s + nx2_ * ix3s),
        dpm2 = nx1_ - ix1e + ix1s - 1,
        dpm3 = nx1_ * (nx2_ - ix2e + ix2s - 1);
    pmw = pmw0 + dpm1;

    Real p[nprop], *pm[nprop];
    for (int n = 0; n < nprop; ++n) {
      p[n] = par(p1+n,k);
      pm[n] = pm0[n] + dpm1;
    }

    // Weight each cell and accumulate particle property onto meshaux.
    for (int ix3 = ix3s; ix3 <= ix3e; ++ix3) {
      Real w3 = active3_ ? _WeightFunction(ix3 + 0.5 - xi3) : 1.0;

      for (int ix2 = ix2s; ix2 <= ix2e; ++ix2) {
        Real w23 = w3 * (active2_ ? _WeightFunction(ix2 + 0.5 - xi2) : 1.0);

        for (int ix1 = ix1s; ix1 <= ix1e; ++ix1) {
          Real w = w23 * (active1_ ?  _WeightFunction(ix1 + 0.5 - xi1) : 1.0);
          *pmw++ += w;

          for (int n = 0; n < nprop; ++n)
            *pm[n]++ += w * p[n];
        }
        pmw += dpm2;
        for (int n = 0; n < nprop; ++n)
          pm[n] += dpm2;
      }
      pmw += dpm3;
      for (int n = 0; n < nprop; ++n)
        pm[n] += dpm3;
    }
  }

  // Treat neighbors of different levels.
  if (pmesh_->multilevel)
    AssignParticlesToDifferentLevels(par, p1, p2, ma1);
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::InterpolateMeshAndAssignParticles(
//               const AthenaArray<Real>& meshsrc, int ms1, int ms2,
//               AthenaArray<Real>& pardst, int pd1,
//               const AthenaArray<Real>& parsrc, int ps1, int ps2, int ma1)
//  \brief interpolates meshsrc from property index ms1 to ms2 onto particle array
//      pardst from index pd1 and up, and assigns parsrc from property index ps1 to ps2
//      onto meshaux from ma1 and up.  The arrays parsrc and pardst can be realprop,
//      auxprop, or work in Particles class.

void ParticleMesh::InterpolateMeshAndAssignParticles(
         const AthenaArray<Real>& meshsrc, int ms1, int ms2,
         AthenaArray<Real>& pardst, int pd1,
         const AthenaArray<Real>& parsrc, int ps1, int ps2, int ma1)
{
  // Zero out destination particle arrays.
  int ni = ms2 - ms1 + 1;
  for (int n = 0; n < ni; ++n) {
    Real *p = &pardst(pd1+n,0);
    for (int k = 0; k < ppar_->npar; ++k)
      *p++ = 0.0;
  }

  // Zero out meshaux.
  Real *p = &weight(0,0,0);
  for (int i = 0; i < ncells_; ++i)
    *p++ = 0.0;

  int na = ps2 - ps1 + 1;
  for (int n = 0; n < na; ++n) {
    Real *p = &meshaux(ma1+n,0,0,0);
    for (int i = 0; i < ncells_; ++i)
      *p++ = 0.0;
  }

  // Transpose meshsrc.
  AthenaArray<Real> u;
  u.NewAthenaArray(meshsrc.GetDim3(), meshsrc.GetDim2(), meshsrc.GetDim1(), ni);
  for (int n = 0; n < ni; ++n)
    for (int k = 0; k < meshsrc.GetDim3(); ++k)
      for (int j = 0; j < meshsrc.GetDim2(); ++j)
        for (int i = 0; i < meshsrc.GetDim1(); ++i)
          u(k,j,i,n) = meshsrc(ms1+n,k,j,i);

  // Get the dimensions of each particle cloud.
  int npc1 = active1_ ? 2 * NGPM + 1 : 1,
      npc2 = active2_ ? 2 * NGPM + 1 : 1,
      npc3 = active3_ ? 2 * NGPM + 1 : 1;

  // Loop over each particle.
  for (int k = 0; k < ppar_->npar; ++k) {
    // Find the domain the particle influences.
    Real xi1 = ppar_->xi1(k), xi2 = ppar_->xi2(k), xi3 = ppar_->xi3(k);
    int imb1 = int(xi1 - dxi1_), imb2 = int(xi2 - dxi2_), imb3 = int(xi3 - dxi3_);
    int ima1 = imb1 - (active1_ ? OFFSET : 0),
        ima2 = imb2 - (active2_ ? OFFSET : 0),
        ima3 = imb3 - (active3_ ? OFFSET : 0);
    xi1 = imb1 + 0.5 - xi1;
    xi2 = imb2 + 0.5 - xi2;
    xi3 = imb3 + 0.5 - xi3;

    // Initialize interpolated properties and fetch those of the particle for assignment.
    Real pd[ni], ps[na];
    for (int n = 0; n < ni; ++n)
      pd[n] = 0.0;
    for (int n = 0; n < na; ++n)
      ps[n] = parsrc(ps1+n,k);

    // Weigh each cell.
    for (int ipc3 = 0; ipc3 < npc3; ++ipc3) {
      Real w3 = active3_ ? _WeightFunction(xi3 + ipc3) : 1.0;

      for (int ipc2 = 0; ipc2 < npc2; ++ipc2) {
        Real w23 = w3 * (active2_ ? _WeightFunction(xi2 + ipc2) : 1.0);

        for (int ipc1 = 0; ipc1 < npc1; ++ipc1) {
          Real w = w23 * (active1_ ? _WeightFunction(xi1 + ipc1) : 1.0);

          // Record the weights.
          weight(ima3+ipc3,ima2+ipc2,ima1+ipc1) += w;

          // Interpolate meshsrc to particles.
          for (int n = 0; n < ni; ++n)
            pd[n] += w * u(ipc3,ipc2,ipc1,n);

          // Assign particles to meshaux.
          for (int n = 0; n < na; ++n)
            meshaux(ma1+n,ima3+ipc3,ima2+ipc2,ima1+ipc1) += w * ps[n];
        }
      }
    }

    // Record the final interpolated properties.
    for (int n = 0; n < ni; ++n)
      pardst(pd1+n,k) = pd[n];
  }

  // Release working array.
  u.DeleteAthenaArray();

  // Treat neighbors of different levels.
  if (pmesh_->multilevel)
    AssignParticlesToDifferentLevels(parsrc, ps1, ps2, ma1);
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::DepositMeshAux(
//               AthenaArray<Real>& u, int ma1, int ma2, int mb1)
//  \brief deposits data in meshaux from property index ma1 to ma2 to meshblock data u
//         from property index mb1 and up, divided by cell volume.

void ParticleMesh::DepositMeshAux(AthenaArray<Real>& u, int ma1, int ma2, int mb1)
{
  Coordinates *pc = pmb_->pcoord;
  int nprop = ma2 - ma1 + 1;
  for (int n = 0; n < nprop; ++n)
    for (int ka = ks, kb = pmb_->ks; ka <= ke; ++ka, ++kb)
      for (int ja = js, jb = pmb_->js; ja <= je; ++ja, ++jb)
        for (int ia = is, ib = pmb_->is; ia <= ie; ++ia, ++ib)
          u(mb1+n,kb,jb,ib) += meshaux(ma1+n,ka,ja,ia) / pc->GetCellVolume(kb,jb,ib);
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::InitiateBoundaryData()
//  \brief allocates space for boundary data.

void ParticleMesh::InitiateBoundaryData()
{
  for (int n = 0; n < pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    BoundaryAttributes& ba = ba_[n];

    int nrecv = (ba.ire - ba.irs + 1) *
                (ba.jre - ba.jrs + 1) *
                (ba.kre - ba.krs + 1) * nmeshaux;
    bd_.recv[nb.bufid] = new Real [nrecv];

#ifdef MPI_PARALLEL
    if (nb.rank != Globals::my_rank) {
      int nsend = ba.ngtot * nmeshaux;
      bd_.send[nb.bufid] = new Real [nsend];
      MPI_Recv_init(bd_.recv[nb.bufid], nrecv, MPI_ATHENA_REAL, nb.rank,
                    (pmb_->lid<<6) | nb.bufid, my_comm, &bd_.req_recv[nb.bufid]);
      MPI_Send_init(bd_.send[nb.bufid], nsend, MPI_ATHENA_REAL, nb.rank,
                    (nb.lid<<6) | nb.targetid, my_comm, &bd_.req_send[nb.bufid]);
    }
#endif
  }
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::SetBoundaryAttributes()
//  \brief initializes or reinitializes attributes for each boundary.

void ParticleMesh::SetBoundaryAttributes()
{
  const RegionSize& block_size = pmb_->block_size;
  const Real xi1mid = (pmb_->is + pmb_->ie + 1) / 2,
             xi2mid = (pmb_->js + pmb_->je + 1) / 2,
             xi3mid = (pmb_->ks + pmb_->ke + 1) / 2;
  const int mylevel = pmb_->loc.level;
  const int myfx1 = int(pbval_->loc.lx1 & 1L),
            myfx2 = int(pbval_->loc.lx2 & 1L),
            myfx3 = int(pbval_->loc.lx3 & 1L);
  const int nx1h = active1_ ? block_size.nx1 / 2 + NGPM : 1,
            nx2h = active2_ ? block_size.nx2 / 2 + NGPM : 1,
            nx3h = active3_ ? block_size.nx3 / 2 + NGPM : 1;

  // Loop over each neighbor block.
  for (int n = 0; n < pbval_->nneighbor; ++n) {
    NeighborBlock& nb = pbval_->neighbor[n];

    // Find the index domain.
    Real xi1min = pmb_->is, xi1max = pmb_->ie + 1,
         xi2min = pmb_->js, xi2max = pmb_->je + 1,
         xi3min = pmb_->ks, xi3max = pmb_->ke + 1;
    Real xi1_0 = xi1min, xi2_0 = xi2min, xi3_0 = xi3min;
    int irs = is, ire = ie, jrs = js, jre = je, krs = ks, kre = ke;
    int iss = is, ise = ie, jss = js, jse = je, kss = ks, kse = ke;

    // Find the radius of influence needed from the neighbor block.
    Real dxi;
    if (nb.level > mylevel)
      dxi = 0.5 * RINF;
    else if (nb.level < mylevel)
      dxi = 2 * RINF;
    else
      dxi = RINF;

    // Consider the normal directions.
    if (nb.ox1 > 0) {
      xi1min = xi1max - dxi;
      xi1_0 = xi1max;
      irs = ie - NGPM + 1;
      iss = ie + 1;
      ise += NGPM;
    } else if (nb.ox1 < 0) {
      xi1max = xi1min + dxi;
      xi1_0 = xi1min - dxi;
      ire = is + NGPM - 1;
      iss -= NGPM;
      ise = is - 1;
    }

    if (nb.ox2 > 0) {
      xi2min = xi2max - dxi;
      xi2_0 = xi2max;
      jrs = je - NGPM + 1;
      jss = je + 1;
      jse += NGPM;
    } else if (nb.ox2 < 0) {
      xi2max = xi2min + dxi;
      xi2_0 = xi2min - dxi;
      jre = js + NGPM - 1;
      jss -= NGPM;
      jse = js - 1;
    }

    if (nb.ox3 > 0) {
      xi3min = xi3max - dxi;
      xi3_0 = xi3max;
      krs = ke - NGPM + 1;
      kss = ke + 1;
      kse += NGPM;
    } else if (nb.ox3 < 0) {
      xi3max = xi3min + dxi;
      xi3_0 = xi3min - dxi;
      kre = ks + NGPM - 1;
      kss -= NGPM;
      kse = ks - 1;
    }

    // Consider the transverse directions.
    if (nb.level > mylevel) {  // Neighbor block is at a finer level.
      if (nb.type == NEIGHBOR_FACE) {
        if (nb.ox1 != 0) {
          if (active2_) {
            if (nb.fi1) {
              xi2min = xi2mid - dxi;
              xi2_0 = xi2mid;
              jrs = je - nx2h + 1;
            } else {
              xi2max = xi2mid + dxi;
              jre = js + nx2h - 1;
            }
          }
          if (active3_) {
            if (nb.fi2) {
              xi3min = xi3mid - dxi; 
              xi3_0 = xi3mid;
              krs = ke - nx3h + 1;
            } else {
              xi3max = xi3mid + dxi;
              kre = ks + nx3h - 1;
            }
          }
        } else if (nb.ox2 != 0) {
          if (active1_) {
            if (nb.fi1) {
              xi1min = xi1mid - dxi; 
              xi1_0 = xi1mid;
              irs = ie - nx1h + 1;
            } else {
              xi1max = xi1mid + dxi;
              ire = is + nx1h - 1;
            }
          }
          if (active3_) {
            if (nb.fi2) {
              xi3min = xi3mid - dxi;
              xi3_0 = xi3mid;
              krs = ke - nx3h + 1;
            } else {
              xi3max = xi3mid + dxi;
              kre = ks + nx3h - 1;
            }
          }
        } else {
          if (active1_) {
            if (nb.fi1) {
              xi1min = xi1mid - dxi;
              xi1_0 = xi1mid;
              irs = ie - nx1h + 1;
            } else {
              xi1max = xi1mid + dxi;
              ire = is + nx1h - 1;
            }
          }
          if (active2_) {
            if (nb.fi2) {
              xi2min = xi2mid - dxi;
              xi2_0 = xi2mid;
              jrs = je - nx2h + 1;
            } else {
              xi2max = xi2mid + dxi;
              jre = js + nx2h - 1;
            }
          }
        }
      } else if (nb.type == NEIGHBOR_EDGE) {
        if (nb.ox1 == 0) {
          if (active1_) {
            if (nb.fi1) {
              xi1min = xi1mid - dxi;
              xi1_0 = xi1mid;
              irs = ie - nx1h + 1;
            } else {
              xi1max = xi1mid + dxi;
              ire = is + nx1h - 1;
            }
          }
        } else if (nb.ox2 == 0) {
          if (active2_) {
            if (nb.fi1) {
              xi2min = xi2mid - dxi;
              xi2_0 = xi2mid;
              jrs = je - nx2h + 1;
            } else {
              xi2max = xi2mid + dxi;
              jre = js + nx2h - 1;
            }
          }
        } else
          if (active3_) {
            if (nb.fi1) {
              xi3min = xi3mid - dxi;
              xi3_0 = xi3mid;
              krs = ke - nx3h + 1;
            } else {
              xi3max = xi3mid + dxi;
              kre = ks + nx3h - 1;
            }
          }
      }
    } else if (nb.level < mylevel) {  // Neighbor block is at a coarser level.
      if (nb.type == NEIGHBOR_FACE) {
        if (nb.ox1 != 0) {
          if (active2_ && myfx2) xi2_0 = xi2min - dxi;
          if (active3_ && myfx3) xi3_0 = xi3min - dxi;
        } else if (nb.ox2 != 0) {
          if (active1_ && myfx1) xi1_0 = xi1min - dxi;
          if (active3_ && myfx3) xi3_0 = xi3min - dxi;
        } else {
          if (active1_ && myfx1) xi1_0 = xi1min - dxi;
          if (active2_ && myfx2) xi2_0 = xi2min - dxi;
        }
      } else if (nb.type == NEIGHBOR_EDGE) {
        if (nb.ox1 == 0) {
          if (active1_ && myfx1) xi1_0 = xi1min - dxi;
        } else if (nb.ox2 == 0) {
          if (active2_ && myfx2) xi2_0 = xi2min - dxi;
        } else
          if (active3_ && myfx3) xi3_0 = xi3min - dxi;
      }
    }

    // Set the domain that influences the ghost block.
    BoundaryAttributes& ba = ba_[n];
    ba.xi1min = xi1min;  ba.xi1max = xi1max;
    ba.xi2min = xi2min;  ba.xi2max = xi2max;
    ba.xi3min = xi3min;  ba.xi3max = xi3max;

    // Set the origin of the ghost block.
    ba.xi1_0 = xi1_0;
    ba.xi2_0 = xi2_0;
    ba.xi3_0 = xi3_0;

    // Set the dimensions of the ghost block.
    if (nb.level < mylevel) {
      ba.ngx1 = (nb.ox1 == 0) ? nx1h : NGPM;
      ba.ngx2 = (nb.ox2 == 0) ? nx2h : NGPM;
      ba.ngx3 = (nb.ox3 == 0) ? nx3h : NGPM;
    } else {
      ba.ngx1 = (nb.ox1 == 0) ? block_size.nx1 : NGPM;
      ba.ngx2 = (nb.ox2 == 0) ? block_size.nx2 : NGPM;
      ba.ngx3 = (nb.ox3 == 0) ? block_size.nx3 : NGPM;
    }
    ba.ngtot = ba.ngx1 * ba.ngx2 * ba.ngx3;

    // Set the indices in meshaux to send and receive.
    ba.irs = irs;  ba.ire = ire;
    ba.jrs = jrs;  ba.jre = jre;
    ba.krs = krs;  ba.kre = kre;

    ba.iss = iss;  ba.ise = ise;
    ba.jss = jss;  ba.jse = jse;
    ba.kss = kss;  ba.kse = kse;
  }
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::AssignParticlesToDifferentLevels(
//               const AthenaArray<Real>& par, int p1, int p2, int ma1)
//  \brief assigns particle array par from property index p1 to p2 to meshaux from
//         property index ma1 and up in neighbors of different levels.

void ParticleMesh::AssignParticlesToDifferentLevels(
         const AthenaArray<Real>& par, int p1, int p2, int ma1)
{
  const int mylevel = pmb_->loc.level;

  // Find neighbor blocks that are on a different level.
  int nprop = p2 - p1 + 1;
  for (int i = 0; i < pbval_->nneighbor; ++i) {
    NeighborBlock& nb = pbval_->neighbor[i];
    if (nb.level == mylevel) continue;

    // Identify the buffer for assignment.
    Real *pbuf0 = NULL;
    BoundaryData *pnbd = NULL;
    if (nb.rank == Globals::my_rank) {
      pnbd = &(pmesh_->FindMeshBlock(nb.gid)->ppar->ppm->bd_);
      pbuf0 = pnbd->recv[nb.targetid];
    }
#ifdef MPI_PARALLEL
    else
      pbuf0 = bd_.send[nb.bufid];
#endif

    // Zero out the buffer.
    BoundaryAttributes& ba = ba_[i];
    Real *pbufw = pbuf0 + iweight * ba.ngtot, *bufw = pbufw;
    for (int j = 0; j < ba.ngtot; ++j)
      *bufw++ = 0.0;

    Real *pbuf[nprop], *buf[nprop];
    for (int n = 0; n < nprop; ++n) {
      buf[n] = pbuf[n] = pbuf0 + (ma1 + n) * ba.ngtot;
      for (int j = 0; j < ba.ngtot; ++j)
        *buf[n]++ = 0.0;
    }

    // Find particles that influences the neighbor block.
    for (int k = 0; k < ppar_->npar; ++k) {
      Real xi1 = ppar_->xi1(k), xi2 = ppar_->xi2(k), xi3 = ppar_->xi3(k);
      if ((active1_ && (xi1 <= ba.xi1min || xi1 >= ba.xi1max)) ||
          (active2_ && (xi2 <= ba.xi2min || xi2 >= ba.xi2max)) ||
          (active3_ && (xi3 <= ba.xi3min || xi3 >= ba.xi3max))) continue;

      // Shift and scale the position index of the particle.
      xi1 -= ba.xi1_0;
      xi2 -= ba.xi2_0;
      xi3 -= ba.xi3_0;
      if (nb.level > mylevel) {
        xi1 *= 2;
        xi2 *= 2;
        xi3 *= 2;
      } else {
        xi1 /= 2;
        xi2 /= 2;
        xi3 /= 2;
      }

      // Find the region of the ghost block to assign the particle to.
      int ix1s = std::max(int(xi1 - dxi1_), 0),
          ix1e = std::min(int(xi1 + dxi1_), ba.ngx1-1),
          ix2s = std::max(int(xi2 - dxi2_), 0),
          ix2e = std::min(int(xi2 + dxi2_), ba.ngx2-1),
          ix3s = std::max(int(xi3 - dxi3_), 0),
          ix3e = std::min(int(xi3 + dxi3_), ba.ngx3-1);

      // Stack the properties of the particle and set the pointers to the buffer.
      Real prop[nprop];
      int dbuf1 = ix1s + ba.ngx1 * (ix2s + ba.ngx2 * ix3s),
           dbuf2 = ba.ngx1 - ix1e + ix1s - 1,
           dbuf3 = ba.ngx1 * (ba.ngx2 - ix2e + ix2s - 1);

      bufw = pbufw + dbuf1;
      for (int n = 0; n < nprop; ++n) {
        prop[n] = par(p1+n,k);
        buf[n] = pbuf[n] + dbuf1;
      }

      // Assign the particle.
      for (int ix3 = ix3s; ix3 <= ix3e; ++ix3) {
        Real w3 = active3_ ? _WeightFunction(ix3 + 0.5 - xi3) : 1.0;

        for (int ix2 = ix2s; ix2 <= ix2e; ++ix2) {
          Real w23 = w3 * (active2_ ?
                             _WeightFunction(ix2 + 0.5 - xi2) : 1.0);

          for (int ix1 = ix1s; ix1 <= ix1e; ++ix1) {
            Real w = w23 * (active1_ ?
                             _WeightFunction(ix1 + 0.5 - xi1) : 1.0);
            *bufw++ += w;

            for (int n = 0; n < nprop; ++n)
              *buf[n]++ += w * prop[n];
          }

          bufw += dbuf2;
          for (int n = 0; n < nprop; ++n)
            buf[n] += dbuf2;
        }

        bufw += dbuf3;
        for (int n = 0; n < nprop; ++n)
          buf[n] += dbuf3;
      }
    }

    // Set the boundary flag.
    if (nb.rank == Globals::my_rank)
      pnbd->flag[nb.targetid] = BNDRY_ARRIVED;
#ifdef MPI_PARALLEL
    else
      MPI_Start(&bd_.req_send[nb.bufid]);
#endif
  }
}

//--------------------------------------------------------------------------------------
//! \fn int ParticleMesh::LoadBoundaryBufferSameLevel(
//                            Real *buf, const BoundaryAttributes& ba)
//  \brief Fill boundary buffers for sending to a block on the same level

int ParticleMesh::LoadBoundaryBufferSameLevel(Real *buf, const BoundaryAttributes& ba)
{
  int p = 0;
  BufferUtility::Pack4DData(meshaux, buf, 0, nmeshaux - 1,
                            ba.iss, ba.ise, ba.jss, ba.jse, ba.kss, ba.kse, p);
  return p;
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::AddBoundaryBuffer(Real *buf, const BoundaryAttributes& ba)
//  \brief Add boundary buffer from a neighbor block to meshaux.

void ParticleMesh::AddBoundaryBuffer(Real *buf, const BoundaryAttributes& ba)
{
  // Add the data to the meshaux.
  for (int n = 0; n < nmeshaux; ++n)
    for (int k = ba.krs; k <= ba.kre; ++k)
      for (int j = ba.jrs; j <= ba.jre; ++j)
        for (int i = ba.irs; i <= ba.ire; ++i)
          meshaux(n,k,j,i) += *buf++;
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::ClearBoundary()
//  \brief clears boundary data to neighboring blocks.

void ParticleMesh::ClearBoundary()
{
  for (int n = 0; n < pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    bd_.flag[nb.bufid] = BNDRY_WAITING;
  }
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::SendBoundary()
//  \brief Send boundary values to neighboring blocks.

void ParticleMesh::SendBoundary()
{
  const int mylevel = pmb_->loc.level;

  for (int n = 0; n < pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];

    // Load boundary values.
    if (nb.level == mylevel) {
      if (nb.rank == Globals::my_rank) {
        BoundaryData *pnbd = &(pmesh_->FindMeshBlock(nb.gid)->ppar->ppm->bd_);
        LoadBoundaryBufferSameLevel(pnbd->recv[nb.targetid], ba_[n]);
        pnbd->flag[nb.targetid] = BNDRY_ARRIVED;
      }
#ifdef MPI_PARALLEL
      else {
        LoadBoundaryBufferSameLevel(bd_.send[nb.bufid], ba_[n]);
        MPI_Start(&bd_.req_send[nb.bufid]);
      }
#endif
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::StartReceiving()
//  \brief starts receiving meshaux near boundary from neighbor processes.

void ParticleMesh::StartReceiving()
{
#ifdef MPI_PARALLEL
  for (int n = 0; n < pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (nb.rank == Globals::my_rank) continue;
    MPI_Start(&bd_.req_recv[nb.bufid]);
  }
#endif
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::ReceiveBoundary()
//  \brief receives boundary values from neighboring blocks and add to my block and
//         returns a flag indicating if all receives are completed.

#include <iomanip>

bool ParticleMesh::ReceiveBoundary()
{
  bool completed = true;

  for (int n = 0; n < pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    enum BoundaryStatus& bstatus = bd_.flag[nb.bufid];

    int flag = 0;
    switch (bstatus) {

    case BNDRY_WAITING:
#ifdef MPI_PARALLEL
      if (nb.rank != Globals::my_rank) {
        MPI_Test(&bd_.req_recv[nb.bufid], &flag, MPI_STATUS_IGNORE);
        if (flag) bstatus = BNDRY_ARRIVED;
      }
#endif
      if (!flag) {
        completed = false;
        break;
      }

    case BNDRY_ARRIVED:
      AddBoundaryBuffer(bd_.recv[nb.bufid], ba_[n]);
      bstatus = BNDRY_COMPLETED;
      break;

    case BNDRY_COMPLETED:
      break;
    }
  }

  return completed;
}

//--------------------------------------------------------------------------------------
//! \fn Real _WeightFunction(Real dxi)
//  \brief evaluates the weight function given index distance.

Real _WeightFunction(Real dxi)
{
  dxi = std::min(std::abs(dxi), 1.5);
  return dxi < 0.5 ? 0.75 - dxi * dxi : 0.5 * ((1.5 - dxi) * (1.5 - dxi));
}
