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
//! \file rad_transport.cpp
//  \brief implementation of radiation integrators
//======================================================================================


// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../parameter_input.hpp"
#include "../../mesh/mesh.hpp"
#include "../radiation.hpp"
#include "../../coordinates/coordinates.hpp" //


// class header
#include "rad_integrators.hpp"

// MPI/OpenMP header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif


// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif


void RadIntegrator::CalculateFluxes(MeshBlock *pmb, AthenaArray<Real> &w,
                     AthenaArray<Real> &ir, const int step)
{
  Radiation *prad=pmy_rad;
  Coordinates *pco = pmb->pcoord;
  
  int nang=pmy_rad->nang;
  int nfreq=pmy_rad->nfreq;
  Real invcrat=1.0/pmy_rad->crat;
  
  Real tau_fact;
  
  AthenaArray<Real> &x1flux=prad->flux[X1DIR];
  AthenaArray<Real> &x2flux=prad->flux[X2DIR];
  AthenaArray<Real> &x3flux=prad->flux[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  int il, iu, jl, ju, kl, ku;


  int tid=0;
  int nthreads = pmb->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) private(tid) num_threads(nthreads)
{
#ifdef OPENMP_PARALLEL
    tid=omp_get_thread_num();
#endif
   
    // First, prepare the Div(q), separate advection part and co-moving part
    int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
    int ncells2 = pmb->block_size.nx2;
    if(ncells2 > 1) ncells2 += 2*(NGHOST);
    int ncells3 = pmb->block_size.nx3;
    if(ncells3 > 1) ncells3 += 2*(NGHOST);

    temp_i1_ = ir;

    if(adv_flag_ > 0){
      for(int k=0; k<ncells3; ++k)
        for(int j=0; j<ncells2; ++j){
            if (ncells2 > 1) pco->CenterWidth2(k,j,0,ncells1-1,cwidth2_);
            if (ncells3 > 1) pco->CenterWidth3(k,j,0,ncells1-1,cwidth3_);
          for(int i=0; i<ncells1; ++i){
             Real vx = w(IVX,k,j,i);
             Real vy = w(IVY,k,j,i);
             Real vz = w(IVZ,k,j,i);
             for(int ifr=0; ifr<nfreq; ++ifr){
               Real ds = pco->dx1v(i);
               if(ncells2 > 1) ds += cwidth2_(i);
               if(ncells3 > 1) ds += cwidth3_(i);
               
               GetTaufactor(vx, vy, vz, ds,
                        prad->sigma_a(k,j,i,ifr)+prad->sigma_s(k,j,i,ifr), &tau_fact);
#pragma omp simd
               for(int n=0; n<nang; ++n){
               
                 Real vdotn = vx*prad->mu(0,k,j,i,n)+vy*prad->mu(1,k,j,i,n)
                             + vz*prad->mu(2,k,j,i,n);
                 
                 vdotn *= invcrat;
                 
                 Real adv_coef = tau_fact * vdotn * (3.0 + vdotn * vdotn);
                 Real q1 = ir(k,j,i,n+ifr*nang) * (1.0 - adv_coef);
                 temp_i1_(k,j,i,n+ifr*nang) = q1;
                 temp_i2_(k,j,i,n+ifr*nang) = adv_coef;
               
               }
            }
          }
        }
    } // End adv_flag
    
    
    
//--------------------------------------------------------------------------------------
// i-direction
    // set the loop limits
    jl=js, ju=je, kl=ks, ku=ke;
    for (int k=kl; k<=ku; ++k){
#pragma omp for schedule(static)
      for (int j=jl; j<=ju; ++j){
        // get the velocity
        for(int i=is; i<=ie+1; ++i){
          Real dxl = pco->x1f(i)-pco->x1v(i-1);
          Real dxr = pco->x1v(i) - pco->x1f(i);
          Real factl = dxr/(dxl+dxr);
          Real factr = dxl/(dxl+dxr);
          for(int ifr=0; ifr<nfreq; ++ifr){
#pragma omp simd
            for(int n=0; n<nang; ++n){
            // linear intepolation between x1v(i-1), x1f(i), x1v(i)
              vel_(k,j,i,n+ifr*nang) = prad->reduced_c *
                                       (factl * prad->mu(0,k,j,i-1,n)
                                      + factr * prad->mu(0,k,j,i,n));
              x1flux(k,j,i,n+ifr*nang) = 0.0;                         
            }
            if(adv_flag_ > 0){
#pragma omp simd
              for(int n=0; n<nang; ++n){
              // linear intepolation between x1v(i-1), x1f(i), x1v(i)
                
                vel2_(k,j,i,n+ifr*nang) = prad->reduced_c *
                    (factl * prad->mu(0,k,j,i-1,n) * temp_i2_(k,j,i-1,n+ifr*nang)
                     + factr * prad->mu(0,k,j,i,n) * temp_i2_(k,j,i,n+ifr*nang));
              }

            }// end adv_flag_
          }
        }
      }
    }
    // calculate the flux
    if(step == 1){
      FirstOrderFluxX1(pco, kl, ku, jl, ju, is, ie+1, temp_i1_, vel_, x1flux);
    }else{
      SecondOrderFluxX1(pco, kl, ku, jl, ju, is, ie+1, temp_i1_, vel_, x1flux);
//      ThirdOrderFluxX1Uniform(pco, kl, ku, jl, ju, is, ie+1, temp_i1_, vel_, x1flux);
    }

    if(adv_flag_ > 0){
      if(step == 1){
        FirstOrderFluxX1(pco, kl, ku, jl, ju, is, ie+1, ir, vel2_, x1flux);
      }else{
        SecondOrderFluxX1(pco, kl, ku, jl, ju, is, ie+1, ir, vel2_, x1flux);
//        ThirdOrderFluxX1Uniform(pco, kl, ku, jl, ju, is, ie+1, ir, vel2_, x1flux);
      }

    }
    //Add second order flux correction at certain regions
    if(step == 2 && flux_correct_flag_ > 0){
      for(int k=kl; k<=ku; ++k){
        for(int j=jl; j<=ju; ++j){
          for(int i=is; i<=ie+1; ++i){
            Real tauc=0.0;
            Real ds = pco->dx1v(i);
            for(int ifr=0; ifr<nfreq; ++ifr){
              tauc += ds * (prad->sigma_a(k,j,i,ifr)
                              +prad->sigma_s(k,j,i,ifr));
            }
            if(tauc < tau_limit_){
              // first, initialize the flux to be zero
              for(int n=0; n<nfreq*nang; ++n)
                x1flux(k,j,i,n) = 0.0; 
              SecondOrderCorrectFluxX1(pco, k, j, i, temp_i1_, vel_, x1flux);
              if(adv_flag_ > 0){
                SecondOrderCorrectFluxX1(pco, k, j, i, ir, vel2_, x1flux);               
              }
            }// end tau_lim
          }// end i
        }// end j
      }// end k
    }// end flag

    
// j-direction
   if(pmb->block_size.nx2 > 1){
     il=is; iu=ie; kl=ks; ku=ke;
     for (int k=kl; k<=ku; ++k){
#pragma omp for schedule(static)
       for (int j=js; j<=je+1; ++j){
        // get the velocity
         for(int i=il; i<=iu; ++i){
           Real dxl = pco->x2f(j)-pco->x2v(j-1);
           Real dxr = pco->x2v(j) - pco->x2f(j);
           Real factl = dxr/(dxl+dxr);
           Real factr = dxl/(dxl+dxr);
           for(int ifr=0; ifr<nfreq; ++ifr){
#pragma omp simd
             for(int n=0; n<nang; ++n){
            // linear intepolation between x2v(j-1), x2f(j), x2v(j)
               vel_(k,j,i,n+ifr*nang) = prad->reduced_c *
                                (factl * prad->mu(1,k,j-1,i,n)
                               + factr * prad->mu(1,k,j,i,n));
               x2flux(k,j,i,n+ifr*nang) = 0.0; 
             }
             if(adv_flag_ > 0){
#pragma omp simd
               for(int n=0; n<nang; ++n){
              // linear intepolation between x2v(j-1), x2f(j), x2v(j)                 
                 vel2_(k,j,i,n+ifr*nang) = prad->reduced_c *
                    (factl * prad->mu(1,k,j-1,i,n) * temp_i2_(k,j-1,i,n+ifr*nang)
                     + factr * prad->mu(1,k,j,i,n) * temp_i2_(k,j,i,n+ifr*nang));
               }

             }

           }
         }
       }// end j
     }// end k

    if(step == 1){
      FirstOrderFluxX2(pco, kl, ku, js, je+1, il, iu, temp_i1_, vel_, x2flux);
    }else{
      SecondOrderFluxX2(pco, kl, ku, js, je+1, il, iu, temp_i1_, vel_, x2flux);
//      ThirdOrderFluxX2Uniform(pco, kl, ku, js, je+1, il, iu, temp_i1_, vel_, x2flux);
    }

    if(adv_flag_){
  // calculate the flux
      if(step == 1){
        FirstOrderFluxX2(pco, kl, ku, js, je+1, il, iu, ir, vel2_, x2flux);
      }else{
        SecondOrderFluxX2(pco, kl, ku, js, je+1, il, iu, ir, vel2_, x2flux);
//        ThirdOrderFluxX2Uniform(pco, kl, ku, js, je+1, il, iu, ir, vel2_, x2flux);
      }

    }
      //Add second order flux correction at certain regions
      if(step == 2 && flux_correct_flag_ > 0){
        for(int k=kl; k<=ku; ++k){
          for(int j=js; j<=je+1; ++j){
            pco->CenterWidth2(k,j,il,iu,cwidth2_);
            for(int i=il; i<=iu; ++i){
              Real tauc=0.0;
              Real ds = cwidth2_(i);
              for(int ifr=0; ifr<nfreq; ++ifr){
                tauc += ds * (prad->sigma_a(k,j,i,ifr)
                             +prad->sigma_s(k,j,i,ifr));
              }
              if(tauc < tau_limit_){
                // first, initialize the flux to be zero
                for(int n=0; n<nfreq*nang; ++n)
                  x2flux(k,j,i,n) = 0.0; 
                SecondOrderCorrectFluxX2(pco, k, j, i, temp_i1_, vel_, x2flux);
                if(adv_flag_ > 0){
                  SecondOrderCorrectFluxX2(pco, k, j, i, ir, vel2_, x2flux);               
                }
              }// end tau_lim
            }// end i
          }// end j
        }// end k
      }// end flag
       
   }// finish ncell2 > 1
  // k-direction
   if(pmb->block_size.nx3 > 1){
     il=is; iu=ie; jl=js; ju=je;
     for (int k=ks; k<=ke+1; ++k){
#pragma omp for schedule(static)
       for (int j=jl; j<=ju; ++j){
        // get the velocity
         for(int i=il; i<=iu; ++i){
           Real dxl = pco->x3f(k)-pco->x3v(k-1);
           Real dxr = pco->x3v(k)-pco->x3f(k);
           Real factl = dxr/(dxl+dxr);
           Real factr = dxl/(dxl+dxr);
           for(int ifr=0; ifr<prad->nfreq; ++ifr){
#pragma omp simd
             for(int n=0; n<prad->nang; ++n){
            // linear intepolation between x3v(k-1), x3f(k), x3v(k)
               vel_(k,j,i,n+ifr*prad->nang) = prad->reduced_c *
                          (factl * prad->mu(2,k-1,j,i,n)
                        + factr * prad->mu(2,k,j,i,n));
               x3flux(k,j,i,n+ifr*nang) = 0.0;           
               
             }
             if(adv_flag_){
#pragma omp simd
               for(int n=0; n<prad->nang; ++n){
              // linear intepolation between x3v(k-1), x3f(k), x3v(k)                 
                 vel2_(k,j,i,n+ifr*nang) = prad->reduced_c *
                    (factl * prad->mu(2,k-1,j,i,n) * temp_i2_(k-1,j,i,n+ifr*nang)
                     + factr * prad->mu(2,k,j,i,n) * temp_i2_(k,j,i,n+ifr*nang));
               }

             }// End adv_flag

           }
         }
       }}// end j, k
        // calculate the flux
    if(step == 1){
      FirstOrderFluxX3(pco, ks, ke+1, jl, ju, il, iu, temp_i1_, vel_, x3flux);
    }else{
      SecondOrderFluxX3(pco, ks, ke+1, jl, ju, il, iu, temp_i1_, vel_, x3flux);
//      ThirdOrderFluxX3Uniform(pco, ks, ke+1, jl, ju, il, iu, temp_i1_, vel_, x3flux);
    }

    if(adv_flag_){
  // calculate the flux
      if(step == 1){
        FirstOrderFluxX3(pco, ks, ke+1, jl, ju, il, iu, ir, vel2_, x3flux);
      }else{
        SecondOrderFluxX3(pco, ks, ke+1, jl, ju, il, iu, ir, vel2_, x3flux);
 //       ThirdOrderFluxX3Uniform(pco, ks, ke+1, jl, ju, il, iu, ir, vel2_, x3flux);
      }

    }
    //Add second order flux correction at certain regions
    if(step == 2 && flux_correct_flag_ > 0){
      for(int k=ks; k<=ke+1; ++k){
        for(int j=jl; j<=ju; ++j){
          pco->CenterWidth3(k,j,il,iu,cwidth3_);
          for(int i=il; i<=iu; ++i){
            Real tauc=0.0;
            Real ds = cwidth3_(i);
            for(int ifr=0; ifr<nfreq; ++ifr){
              tauc += ds * (prad->sigma_a(k,j,i,ifr)
                           +prad->sigma_s(k,j,i,ifr));
            }
            if(tauc < tau_limit_){
              // first, initialize the flux to be zero
              for(int n=0; n<nfreq*nang; ++n)
                x3flux(k,j,i,n) = 0.0; 
              SecondOrderCorrectFluxX3(pco, k, j, i, temp_i1_, vel_, x3flux);
              if(adv_flag_ > 0){
                SecondOrderCorrectFluxX3(pco, k, j, i, ir, vel2_, x3flux);               
              }
            }// end tau_lim
          }// end i
        }// end j
      }// end k
    }// end flag
     
   }// finish k direction
  
}// end of omp parallel region
  

  
}


void RadIntegrator::FluxDivergence(MeshBlock *pmb, AthenaArray<Real> &ir_in1,
      AthenaArray<Real> &ir_in2, const Real wght[3], AthenaArray<Real> &ir_out)
{
  Radiation *prad=pmb->prad;

  AthenaArray<Real> &x1flux=prad->flux[X1DIR];
  AthenaArray<Real> &x2flux=prad->flux[X2DIR];
  AthenaArray<Real> &x3flux=prad->flux[X3DIR];
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;
  
  Real dt = pmb->pmy_mesh->dt;

    
  int tid=0;
  int nthreads = pmb->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) private(tid) num_threads(nthreads)
{
#ifdef OPENMP_PARALLEL
  tid=omp_get_thread_num();
#endif
  AthenaArray<Real> x1area, x2area, x2area_p1, x3area, x3area_p1, vol, dflx;
  x1area.InitWithShallowSlice(x1face_area_,2,tid,1);
  x2area.InitWithShallowSlice(x2face_area_,2,tid,1);
  x2area_p1.InitWithShallowSlice(x2face_area_p1_,2,tid,1);
  x3area.InitWithShallowSlice(x3face_area_,2,tid,1);
  x3area_p1.InitWithShallowSlice(x3face_area_p1_,2,tid,1);
  vol.InitWithShallowSlice(cell_volume_,2,tid,1);
  dflx.InitWithShallowSlice(flx_,3,tid,1);
    

#pragma omp for schedule(static)
  for (int k=ks; k<=ke; ++k) { 
    for (int j=js; j<=je; ++j) {

      // calculate x1-flux divergence 
      pmb->pcoord->Face1Area(k,j,is,ie+1,x1area);
      for(int i=is; i<=ie; ++i){
#pragma omp simd
        for(int n=0; n<prad->n_fre_ang; ++n){
          dflx(i,n) = (x1area(i+1) *x1flux(k,j,i+1,n) - x1area(i)*x1flux(k,j,i,n));
        }// end n
      }// End i

     // calculate x2-flux
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Face2Area(k,j  ,is,ie,x2area   );
        pmb->pcoord->Face2Area(k,j+1,is,ie,x2area_p1);
      for(int i=is; i<=ie; ++i){
#pragma omp simd
        for(int n=0; n<prad->n_fre_ang; ++n){
            dflx(i,n) += (x2area_p1(i)*x2flux(k,j+1,i,n) - x2area(i)*x2flux(k,j,i,n));
          }
        }
      }// end nx2

      // calculate x3-flux divergence
      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Face3Area(k  ,j,is,ie,x3area   );
        pmb->pcoord->Face3Area(k+1,j,is,ie,x3area_p1);
      for(int i=is; i<=ie; ++i){
#pragma omp simd
        for(int n=0; n<prad->n_fre_ang; ++n){
            dflx(i,n) += (x3area_p1(i)*x3flux(k+1,j,i,n) - x3area(i)*x3flux(k,j,i,n));
          }
        }
      }// end nx3
      // update variable with flux divergence
      pmb->pcoord->CellVolume(k,j,is,ie,vol);
      for(int i=is; i<=ie; ++i){
#pragma omp simd
        for(int n=0; n<prad->n_fre_ang; ++n){
          ir_out(k,j,i,n) = wght[0]*ir_in1(k,j,i,n) + wght[1]*ir_in2(k,j,i,n)
                         - wght[2]*dt*dflx(i,n)/vol(i);
        }
      }


    }// end j
  }// End k
  


}// end omp parallel region

}

void RadIntegrator::GetTaufactor(const Real vx, const Real vy, const Real vz,
                                 const Real ds, const Real sigma, Real *factor)
{

   Real invcrat = 1.0/pmy_rad->crat;
   Real vel = vx*vx+vy*vy+vz*vz;
   Real taucell = taufact_ * ds * sigma;
   Real tausq = taucell * taucell * (vel*invcrat*invcrat);
   if(tausq > 1.e-3) tausq = 1.0 - exp(-tausq);
   (*factor) = tausq;

}
