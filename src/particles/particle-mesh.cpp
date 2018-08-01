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

// Local function prototypes.
Real _ParticleMeshWeightFunction(Real dxi);

// Local constants.
const int OFFSET = NGHOST - NGPM;  // offset between meshblock and meshaux

//--------------------------------------------------------------------------------------
//! \fn ParticleMesh::ParticleMesh(Particles *ppar, int nmeshaux)
//  \brief constructs a new ParticleMesh instance.

ParticleMesh::ParticleMesh(Particles *ppar, int nmeshaux_)
{
  // Save some inputs.
  ppar_ = ppar;
  pmb_ = ppar->pmy_block;
  pmesh_ = pmb_->pmy_mesh;
  pbval_ = pmb_->pbval;
  nmeshaux = nmeshaux_;

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

  // Find the maximum number of neighbors.
  bd_.nbmax = BoundaryBase::BufferID(dim, pmesh_->multilevel);

  // Initialize boundary data.
  for (int n = 0; n < bd_.nbmax; n++) {
    bd_.flag[n] = BNDRY_WAITING;
    bd_.send[n] = NULL;
    bd_.recv[n] = NULL;

    int size = ((pbval_->ni[n].ox1 == 0) ? block_size.nx1 : NGPM) *
               ((pbval_->ni[n].ox2 == 0) ? block_size.nx2 : NGPM) *
               ((pbval_->ni[n].ox3 == 0) ? block_size.nx3 : NGPM) * nmeshaux;
    bd_.send[n] = new Real [size];
    bd_.recv[n] = new Real [size];
  }

  // Set the boundary attributes.
  SetBoundaryAttributes();
}

//--------------------------------------------------------------------------------------
//! \fn ParticleMesh::~ParticleMesh()
//  \brief destructs a ParticleMesh instance.

ParticleMesh::~ParticleMesh()
{
  // Destroy the particle meshblock.
  meshaux.DeleteAthenaArray();

  // Destroy boundary data.
  for (int n = 0; n < bd_.nbmax; n++) {
    delete [] bd_.send[n];
    delete [] bd_.recv[n];
  }
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::InterpolateMeshToParticles(
//               const AthenaArray<Real>& meshsrc, const AthenaArray<int>& imeshsrc,
//               AthenaArray<Real>& par, const AthenaArray<int>& ipar)
//  \brief interpolates meshsrc at specified indices imeshsrc onto particle array par
//         (realprop, auxprop, or work in Particles class) at the corresponding indices
//         ipar.

void ParticleMesh::InterpolateMeshToParticles(
         const AthenaArray<Real>& meshsrc, const AthenaArray<int>& imeshsrc,
         AthenaArray<Real>& par, const AthenaArray<int>& ipar)
{
  // Check the index mapping.
  int nprop = imeshsrc.GetSize();
  if (nprop <= 0 || ipar.GetSize() != nprop) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [Particles::InterpolateMeshToParticles]"
        << std::endl
        << "index arrays imeshsrc and ipar do not match." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  // Zero out the particle arrays.
  Real *pp[nprop];
  for (int n = 0; n < nprop; ++n) {
    Real *p = pp[n] = &par(ipar(n),0);
    for (long k = 0; k < ppar_->npar; ++k)
      *p++ = 0.0;
  }

  // Loop over each particle.
  for (long k = 0; k < ppar_->npar; ++k) {
    // Find the domain the particle influences.
    Real xi1 = ppar_->xi1(k), xi2 = ppar_->xi2(k), xi3 = ppar_->xi3(k);
    int ix1s = int(xi1 - dxi1_), ix1e = int(xi1 + dxi1_);
    int ix2s = int(xi2 - dxi2_), ix2e = int(xi2 + dxi2_);
    int ix3s = int(xi3 - dxi3_), ix3e = int(xi3 + dxi3_);

    // Weight each cell and accumulate the mesh properties onto the particles.
    for (int ix3 = ix3s; ix3 <= ix3e; ++ix3) {
      Real w3 = active3_ ? _ParticleMeshWeightFunction(ix3 + 0.5 - xi3) : 1.0;

      for (int ix2 = ix2s; ix2 <= ix2e; ++ix2) {
        Real w23 = w3 * (active2_ ? _ParticleMeshWeightFunction(ix2 + 0.5 - xi2) : 1.0);

        for (int ix1 = ix1s; ix1 <= ix1e; ++ix1) {
          Real weight = w23 * (active1_ ?
                                    _ParticleMeshWeightFunction(ix1 + 0.5 - xi1) : 1.0);

          for (int n = 0; n < nprop; ++n)
            *pp[n] += weight * meshsrc(imeshsrc(n),ix3,ix2,ix1);
        }
      }
    }
    for (int n = 0; n < nprop; ++n)
      ++pp[n];
  }
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::AssignParticlesToMeshAux(
//               const AthenaArray<Real>& par, const AthenaArray<int>& ipar,
//               const AthenaArray<int>& imeshaux)
//  \brief assigns par (realprop, auxprop, or work in Particles class) at specified
//         indices ipar onto meshaux at the corresponding indices imeshaux.

void ParticleMesh::AssignParticlesToMeshAux(
         const AthenaArray<Real>& par, const AthenaArray<int>& ipar,
         const AthenaArray<int>& imeshaux)
{
  // Check the index mapping.
  int nprop = ipar.GetSize();
  if (nprop <= 0 || imeshaux.GetSize() != nprop) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [Particles::AssignParticlesToMeshAux]"
        << std::endl
        << "index arrays ipar and imeshaux does not match." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  // Zero out meshaux.
  for (int n = 0; n < nprop; ++n) {
    Real* pdata = &meshaux(imeshaux(n),0,0,0);
    for (int i = 0; i < ncells_; ++i)
      *pdata++ = 0.0;
  }

  // Allocate working array.
  AthenaArray<Real> p;
  p.NewAthenaArray(nprop);

  for (long k = 0; k < ppar_->npar; ++k) {
    // Find the domain the particle influences.
    Real xi1 = ppar_->xi1(k) - (active1_ ? OFFSET : 0),
         xi2 = ppar_->xi2(k) - (active2_ ? OFFSET : 0),
         xi3 = ppar_->xi3(k) - (active3_ ? OFFSET : 0);
    int ix1s = int(xi1 - dxi1_), ix1e = int(xi1 + dxi1_);
    int ix2s = int(xi2 - dxi2_), ix2e = int(xi2 + dxi2_);
    int ix3s = int(xi3 - dxi3_), ix3e = int(xi3 + dxi3_);

    // Copy the particle properties.
    for (int n = 0; n < nprop; ++n)
      p(n) = par(ipar(n),k);

    // Weight each cell and accumulate particle property onto meshaux.
    for (int ix3 = ix3s; ix3 <= ix3e; ++ix3) {
      Real w3 = active3_ ? _ParticleMeshWeightFunction(ix3 + 0.5 - xi3) : 1.0;

      for (int ix2 = ix2s; ix2 <= ix2e; ++ix2) {
        Real w23 = w3 * (active2_ ? _ParticleMeshWeightFunction(ix2 + 0.5 - xi2) : 1.0);

        for (int ix1 = ix1s; ix1 <= ix1e; ++ix1) {
          Real weight = w23 * (active1_ ?
                                    _ParticleMeshWeightFunction(ix1 + 0.5 - xi1) : 1.0);

          for (int n = 0; n < nprop; ++n)
            meshaux(imeshaux(n),ix3,ix2,ix1) += weight * p(n);
        }
      }
    }
  }

  // Release working array.
  p.DeleteAthenaArray();

  // Treat neighbors of different levels.
  if (pmesh_->multilevel)
    AssignParticlesToDifferentLevels(par, ipar, imeshaux);
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::InterpolateMeshAndAssignParticles(
//               const AthenaArray<Real>& meshsrc, const AthenaArray<int>& imeshsrc,
//               AthenaArray<Real>& pardst, const AthenaArray<int>& ipardst,
//               const AthenaArray<Real>& parsrc, const AthenaArray<int>& iparsrc,
//               const AthenaArray<int>& imeshaux)
//  \brief interpolates meshsrc at specified indices imeshsrc onto particle array pardst
//         the corresponding indices ipardst, and assigns parsrc at specified indices
//         iparsrc onto meshaux at the corresponding indices imeshaux.  The arrays
//         parsrc and pardst can be realprop, auxprop, or work in Particles class.

void ParticleMesh::InterpolateMeshAndAssignParticles(
         const AthenaArray<Real>& meshsrc, const AthenaArray<int>& imeshsrc,
         AthenaArray<Real>& pardst, const AthenaArray<int>& ipardst,
         const AthenaArray<Real>& parsrc, const AthenaArray<int>& iparsrc,
         const AthenaArray<int>& imeshaux)
{
  // Check the index mapping.
  int nmeshsrc = imeshsrc.GetSize();
  if (nmeshsrc <= 0 || ipardst.GetSize() != nmeshsrc) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [Particles::InterpolateMeshAndAssignParticles]"
        << std::endl
        << "index arrays imeshsrc and ipardst do not match." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  int nmeshdst = imeshaux.GetSize();
  if (nmeshdst <= 0 || iparsrc.GetSize() != nmeshdst) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [Particles::InterpolateMeshAndAssignParticles]"
        << std::endl
        << "index arrays imeshaux and iparsrc do not match." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  // Zero out the destination particle arrays.
  for (int n = 0; n < nmeshsrc; ++n) {
    Real* pdata = &pardst(ipardst(n),0);
    for (long k = 0; k < ppar_->npar; ++k)
      *pdata++ = 0.0;
  }

  // Zero out meshaux.
  for (int n = 0; n < nmeshdst; ++n) {
    Real* pdata = &meshaux(imeshaux(n),0,0,0);
    for (int i = 0; i < ncells_; ++i)
      *pdata++ = 0.0;
  }

  // Allocate working array.
  AthenaArray<Real> p;
  p.NewAthenaArray(nmeshdst);

  // Loop over each particle.
  for (long k = 0; k < ppar_->npar; ++k) {
    // Find the domain the particle influences.
    Real xi1 = ppar_->xi1(k), xi2 = ppar_->xi2(k), xi3 = ppar_->xi3(k);
    int imb1s = int(xi1 - dxi1_), imb1e = int(xi1 + dxi1_);
    int imb2s = int(xi2 - dxi2_), imb2e = int(xi2 + dxi2_);
    int imb3s = int(xi3 - dxi3_), imb3e = int(xi3 + dxi3_);
    int ima1s = imb1s - (active1_ ? OFFSET : 0);
    int ima2s = imb2s - (active2_ ? OFFSET : 0);
    int ima3s = imb3s - (active3_ ? OFFSET : 0);

    // Copy the particle properties.
    for (int n = 0; n < nmeshdst; ++n)
      p(n) = parsrc(iparsrc(n),k);

    // Weigh each cell.
    for (int imb3 = imb3s, ima3 = ima3s; imb3 <= imb3e; ++imb3, ++ima3) {
      Real w3 = active3_ ? _ParticleMeshWeightFunction(imb3 + 0.5 - xi3) : 1.0;

      for (int imb2 = imb2s, ima2 = ima2s; imb2 <= imb2e; ++imb2, ++ima2) {
        Real w23 = w3 * (active2_ ? _ParticleMeshWeightFunction(imb2 + 0.5 - xi2) : 1.0);

        for (int imb1 = imb1s, ima1 = ima1s; imb1 <= imb1e; ++imb1, ++ima1) {
          Real weight = w23 * (active1_ ?
                                    _ParticleMeshWeightFunction(imb1 + 0.5 - xi1) : 1.0);

          // Interpolate mesh to particles.
          for (int n = 0; n < nmeshsrc; ++n)
            pardst(ipardst(n),k) += weight * meshsrc(imeshsrc(n),imb3,imb2,imb1);

          // Assign particles to meshaux.
          for (int n = 0; n < nmeshdst; ++n)
            meshaux(imeshaux(n),ima3,ima2,ima1) += weight * p(n);
        }
      }
    }
  }

  // Release working array.
  p.DeleteAthenaArray();

  // Treat neighbors of different levels.
  if (pmesh_->multilevel)
    AssignParticlesToDifferentLevels(parsrc, iparsrc, imeshaux);
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::DepositMeshAux(AthenaArray<Real>& u,
//               const AthenaArray<int>& imeshaux, const AthenaArray<int>& imeshblock)
//  \brief deposits data in meshaux at specified indices imeshaux to meshblock u at the
//         corresponding indices imeshblock, divided by cell volume.

void ParticleMesh::DepositMeshAux(AthenaArray<Real>& u,
         const AthenaArray<int>& imeshaux, const AthenaArray<int>& imeshblock)
{
  // Check the index mapping.
  const int nprop = imeshaux.GetSize();
  if (nprop <= 0 || imeshblock.GetSize() != nprop) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [Particles::AssignParticlesToMeshAux]"
        << std::endl
        << "index arrays ipar and imeshaux does not match." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  // Add meshaux to meshblock.
  Coordinates *pc = pmb_->pcoord;
  for (int n = 0; n < nprop; ++n) {
    int ima = imeshaux(n), imb = imeshblock(n);

    for (int ka = ks, kb = pmb_->ks; ka <= ke; ++ka, ++kb)
      for (int ja = js, jb = pmb_->js; ja <= je; ++ja, ++jb)
        for (int ia = is, ib = pmb_->is; ia <= ie; ++ia, ++ib)
          u(imb,kb,jb,ib) += meshaux(ima,ka,ja,ia) / pc->GetCellVolume(kb,jb,ib);
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
    ba.ngx12 = ba.ngx1 * ba.ngx2;
    ba.ngtot = ba.ngx12 * ba.ngx3;

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
//               const AthenaArray<Real>& par, const AthenaArray<int>& ipar,
//               const AthenaArray<int>& imeshaux)
//  \brief assigns particles to neighbors of different levels.  The parameters are the
//         same as those in ParticleMesh::AssignParticlesToMeshAux().

void ParticleMesh::AssignParticlesToDifferentLevels(
         const AthenaArray<Real>& par, const AthenaArray<int>& ipar,
         const AthenaArray<int>& imeshaux)
{
  const int mylevel = pmb_->loc.level;

  // Check the index mapping.
  int nprop = ipar.GetSize();
  if (nprop <= 0 || imeshaux.GetSize() != nprop) {
    std::stringstream msg;
    msg << "### FATAL ERROR in function [Particles::AssignParticlesToDifferentLevels]"
        << std::endl
        << "index arrays ipar and imeshaux does not match." << std::endl;
    throw std::runtime_error(msg.str().c_str());
    return;
  }

  // Find neighbor blocks that are on a different level.
  for (int i = 0; i < pbval_->nneighbor; ++i) {
    NeighborBlock& nb = pbval_->neighbor[i];
    if (nb.level == mylevel) continue;

    // Identify the buffer for assignment.
    Real *pbuf0 = NULL;
    BoundaryData *pnbd = NULL;
    if (nb.rank == Globals::my_rank) {
      pnbd = &(pmesh_->FindMeshBlock(nb.gid)->ppar->ppm->bd_);
      pbuf0 = pnbd->recv[nb.targetid];
    } else
      pbuf0 = bd_.send[nb.bufid];

    // Zero out the buffer.
    Real *pbuf[nprop], *buf[nprop];
    BoundaryAttributes& ba = ba_[i];
    for (int n = 0; n < nprop; ++n) {
      buf[n] = pbuf[n] = pbuf0 + imeshaux(n) * ba.ngtot;
      for (int j = 0; j < ba.ngtot; ++j)
        *buf[n]++ = 0.0;
    }

    // Find particles that influences the neighbor block.
    for (long k = 0; k < ppar_->npar; ++k) {
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
      long dbuf1 = ix1s + ba.ngx1 * (ix2s + ba.ngx2 * ix3s),
           dbuf2 = ba.ngx1 - ix1e + ix1s - 1,
           dbuf3 = ba.ngx1 * (ba.ngx2 - ix2e + ix2s - 1);
      for (int n = 0; n < nprop; ++n) {
        prop[n] = par(ipar(n),k);
        buf[n] = pbuf[n] + dbuf1;
      }

      // Assign the particle.
      for (int ix3 = ix3s; ix3 <= ix3e; ++ix3) {
        Real w3 = active3_ ? _ParticleMeshWeightFunction(ix3 + 0.5 - xi3) : 1.0;

        for (int ix2 = ix2s; ix2 <= ix2e; ++ix2) {
          Real w23 = w3 * (active2_ ?
                         _ParticleMeshWeightFunction(ix2 + 0.5 - xi2) : 1.0);

          for (int ix1 = ix1s; ix1 <= ix1e; ++ix1) {
            Real weight = w23 * (active1_ ?
                         _ParticleMeshWeightFunction(ix1 + 0.5 - xi1) : 1.0);

            for (int n = 0; n < nprop; ++n)
              *buf[n]++ += weight * prop[n];
          }

          for (int n = 0; n < nprop; ++n)
            buf[n] += dbuf2;
        }

        for (int n = 0; n < nprop; ++n)
          buf[n] += dbuf3;
      }
    }

    // Set the boundary flag.
    if (nb.rank == Globals::my_rank)
      pnbd->flag[nb.targetid] = BNDRY_ARRIVED;
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
      } else
        LoadBoundaryBufferSameLevel(bd_.send[nb.bufid], ba_[n]);
    }
  }
}

//--------------------------------------------------------------------------------------
//! \fn void ParticleMesh::ReceiveBoundary()
//  \brief Receive boundary values from neighboring blocks and add to my block.

void ParticleMesh::ReceiveBoundary()
{
  for (int n = 0; n < pbval_->nneighbor; n++) {
    NeighborBlock& nb = pbval_->neighbor[n];
    if (bd_.flag[nb.bufid] == BNDRY_COMPLETED) continue;

    AddBoundaryBuffer(bd_.recv[nb.bufid], ba_[n]);
    bd_.flag[nb.bufid] = BNDRY_COMPLETED;
  }
}

//--------------------------------------------------------------------------------------
//! \fn Real _ParticleMeshWeightFunction(Real dxi)
//  \brief evaluates the weight function given index distance.

Real _ParticleMeshWeightFunction(Real dxi)
{
  if (dxi < 0) dxi = -dxi;

  if (dxi < 0.5)
    return 0.75 - dxi * dxi;

  if (dxi < 1.5) {
    dxi = 1.5 - dxi;
    return 0.5 * (dxi * dxi);
  }

  return 0;
}
