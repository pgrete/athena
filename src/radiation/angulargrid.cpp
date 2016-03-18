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
//! \file angulargrid.cpp
//  \brief implementation of functions in class Radiation
//======================================================================================

#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "radiation.hpp"
#include "../utils/utils.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh.hpp"


//--------------------------------------------------------------------------------------
// \!fn void Radiation::AngularGrid(int angle_flag, int nmu)

// \brief function to create the angular grid

void Radiation::AngularGrid(int angle_flag, int nmu)
{
  std::stringstream msg;
  // allocate some temporaray arrays
  AthenaArray<Real> mu2tmp, mutmp, wtmp2, wtmp;
  AthenaArray<Real> pmat, pinv, wpf;
  AthenaArray<int> plab, pl;
  
  int n_ang = nang/noct;
  
  mu2tmp.NewAthenaArray(nmu);
  mutmp.NewAthenaArray(n_ang,3);
  wtmp2.NewAthenaArray(nmu-1);
  wtmp.NewAthenaArray(nmu);
  
  pmat.NewAthenaArray(nmu,nmu);
  pinv.NewAthenaArray(nmu-1,nmu-1);
  plab.NewAthenaArray(n_ang);
  pl.NewAthenaArray(nmu,3);
  wpf.NewAthenaArray(nmu-1);
  
  
  // initialize coordinate direction
  int axisx=0, axisy=1, axisz=2;
  pmy_block->pcoord->AxisDirection(&axisx, &axisy, &axisz);
  
  // check the dimension of the problem
  
  int ndim = 1;
  int n1z = pmy_block->block_size.nx1 + 2*(NGHOST);
  int n2z = 1;
  int n3z = 1;
  
  if(pmy_block->block_size.nx2 > 1) {
    ndim = 2;
    n2z = pmy_block->block_size.nx2 + 2*(NGHOST);
  }
  
  if(pmy_block->block_size.nx3 > 1) {
    ndim = 3;
    n3z = pmy_block->block_size.nx3 + 2*(NGHOST);
  }
  
  
  if(angle_flag == 0){
    Real deltamu = 2.0 / (2 * nmu - 1);
    mu2tmp(0) = 1.0 / (3.0 * (2 * nmu - 1));
    for(int i=1; i<nmu; i++) {
      mu2tmp(i) = mu2tmp(i-1) + deltamu;
    }
    
    Real w2 = 4.0 * mu2tmp(0);
    Real wsum2 = sqrt(w2);
    wtmp2(0) = wsum2;
    
    for(int i=1; i<nmu-2; i++) {
      w2 += deltamu;
      wtmp2(i) = sqrt(w2);
      wsum2 += wtmp2(i);
    }
    
    if(nmu > 2)
      wtmp2(nmu-2) = 2.0*(nmu-1)/3.0 - wsum2;
    
    wtmp(0) = wtmp2(0);
    Real wsum = wtmp(0);
    for(int i=1; i<nmu-1; ++i) {
      wtmp(i) = wtmp2(i) - wtmp2(i-1);
      wsum += wtmp(i);
    }
    wtmp(nmu-1) = 1.0 - wsum;
    
    int np = 0;
    int iang = 0;
    
    // initialize to be zero
    for(int i=0; i<nmu; ++i)
      for(int j=0; j<nmu; ++j)
        pmat(i,j) = 0.0;
    
    for(int i=0; i<nmu; i++) {
      for(int j=0; j<nmu; j++) {
        for (int k=0; k<nmu; k++) {
          if (i + j + k == nmu - 1) {
            // assign cosines to temporary array grid
            mutmp(iang,0) = sqrt(mu2tmp(j));
            mutmp(iang,1) = sqrt(mu2tmp(k));
            mutmp(iang,2) = sqrt(mu2tmp(i));
            
            int ip=Permutation(i,j,k,np,pl);
            if (ip == -1) {
              pl(np,0) = i;
              pl(np,1) = j;
              pl(np,2) = k;
              pmat(i,np) += 1.0;
              plab(iang) = np;
              np++;
            } else {
              pmat(i,ip) += 1.0;
              plab(iang) = ip;
            }
            
            iang++;
          }// end if i+j+k
        }// end k nmu
      }// end j nmu
    }// end i nmu
    
    
    
    
    if (nmu > 1) {
      //  Invert matrix of Permutations families */
      InverseMatrix(nmu-1, pmat,pinv);
      // Solve for and assign weights for each Permutation family
      MatrixMult(nmu-1,nmu-1,pinv,wtmp,wpf);

      for(int l=0; l<noct; ++l){
        for (int i=0; i<n_ang; ++i){
            wmu(l*n_ang+i) = wpf(plab(i));
        }// end nang
      }// end noct

    } else {
      for(int l=0; l<noct; ++l){
          wmu(l) = 1.0;
      }// end l
    }// end nmu > 1
    

    
    

    
    if(ndim == 1){
      AthenaArray<Real> mutmp1d, wtmp1d;
      mutmp1d.NewAthenaArray(2*nmu);
      wtmp1d.NewAthenaArray(2*nmu);
      
      Gauleg(2*nmu, -1.0, 1.0, mutmp1d, wtmp1d);
      for(int i=nmu; i<2*nmu; ++i){
        wmu(i-nmu) = 0.5 * wtmp1d(i);
        wmu(nmu) = 0.5 * wtmp1d(i);
      }
      
      for(int n1=0; n1<n1z; ++n1){
        for(int i=nmu; i<2*nmu; ++i){
          mu(0,0,0,n1,i-nmu) = mutmp1d(i);
          mu(0,0,0,n1,nmu) = -mutmp1d(i);
          
        }
      }
      
      mutmp1d.DeleteAthenaArray();
      wtmp1d.DeleteAthenaArray();
      
      
    }else if (ndim == 2) {
      // for spherical coordinate system, it should be r-theta-phi
      for(int n2=0; n2<n2z; ++n2){
        for(int n1=0; n1<n1z; ++n1){
          for(int j=0; j<2; ++j) {
            for(int k=0; k<2; ++k) {
              int l=2*j+k;
              for (int i=0; i<n_ang; ++i) {
                int mi = l*n_ang + i;
                if (k == 0){
                  mu(0,0,n2,n1,mi) =  mutmp(i,axisx);
                }
                else{
                  mu(0,0,n2,n1,mi) = -mutmp(i,axisx);
                }
                if (j == 0){
                  mu(1,0,n2,n1,mi) =  mutmp(i,axisy);
                }
                else{
                  mu(1,0,n2,n1,mi) = -mutmp(i,axisy);
                }
                
              }// end nang
            }// end k
          }// end j
          
        }// end n1
      }// end n2
      
      // for the angular weight
      for (int i=0; i<noct * n_ang; ++i) {
          wmu(i) *= 0.25;
      }// end nang
      
    } else if(ndim == 3) {
      
      for(int n3=0; n3<n3z; ++n3){
        for(int n2=0; n2<n2z; ++n2){
          for(int n1=0; n1<n1z; ++n1){
            for(int j=0; j<2; ++j) {
              for(int k=0; k<2; ++k) {
                for(int l=0; l<2; ++l) {
                  int m=4*j+2*k+l;

                  for(int i=0; i<n_ang; ++i){
                    int mi = m*n_ang + i;
                    
                    if (l == 0){
                      mu(0,n3,n2,n1,mi) =  mutmp(i,axisx);
                    }
                    else{
                      mu(0,n3,n2,n1,mi) = -mutmp(i,axisx);
                    }
                    
                    if (k == 0){
                      mu(1,n3,n2,n1,mi) =  mutmp(i,axisy);
                    }
                    else{
                      mu(1,n3,n2,n1,mi) = -mutmp(i,axisy);
                    }
                    
                    if (j == 0){
                      mu(2,n3,n2,n1,mi) =  mutmp(i,axisz);
                    }
                    else{
                      mu(2,n3,n2,n1,mi) = -mutmp(i,axisz);
                    }
                    
                  }// end i
                }// end l
              }// end k
            }// end j
          }// end n1
        }// end n2
      }// end n3
      
      for(int i=0; i<noct * n_ang; ++i){
        wmu(i) *= 0.125;
      }
      
    }// end nDIM = 3
  }else if(angle_flag == 10){
    
    if ((ndim == 1) || (ndim == 3)) {
      msg <<"[Warning]: ang_quad = 10 should be used" <<
      "only for 2D problems. \n" << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    
    Real delmu=1.0/sqrt((Real)3.0);
    Real sintheta = sqrt((Real)2.0)/sqrt((Real)3.0);
    Real phi = 0.5 * PI / (Real) (2*nmu);

    for(int i=0; i<nmu; ++i){
      mutmp(i,0) = sintheta * cos(phi*(Real)(2*i+1));
      mutmp(i,1) = sintheta * sin(phi*(Real)(2*i+1));
      mutmp(i,2) = delmu;
    }
    for(int n2=0; n2<n2z; ++n2){
      for(int n1=0; n1<n1z; ++n1){
        for(int j=0; j<2; ++j) {
          for(int k=0; k<2; ++k) {
            int l=2*j+k;

            for(int i=0; i<n_ang; ++i){
              int mi = l*n_ang + i;
              
              if (k == 0){
                mu(0,0,n2,n1,mi) =  mutmp(i,axisx);
              }
              else{
                mu(0,0,n2,n1,mi) = -mutmp(i,axisx);
              }
              
              if (j == 0){
                mu(1,0,n2,n1,mi) =  mutmp(i,axisy);
              }
              else{
                mu(1,0,n2,n1,mi) = -mutmp(i,axisy);
              }
              
            }// end i
          }// end k
        }// end j
      }// end n1
    }// end n2
    
    for(int i=0; i<noct*n_ang; ++i)      wmu(i) = 0.25/((Real)nmu);
    
  }else{
    msg << "### FATAL ERROR in function [InitialAngle]" << std::endl
    << "Type of angular discretization unknow: "<< angle_flag << "\n ";
    throw std::runtime_error(msg.str().c_str());
  }
  
  // Now change the angle cosines to different coordinate systems
  pmy_block->pcoord->ConvertAngle(pmy_block,nang, mu);
  
  
  // free the temporary arrays

  
  
  mu2tmp.DeleteAthenaArray();
  mutmp.DeleteAthenaArray();
  wtmp2.DeleteAthenaArray();
  wtmp.DeleteAthenaArray();
  pmat.DeleteAthenaArray();
  pinv.DeleteAthenaArray();
  plab.DeleteAthenaArray();
  pl.DeleteAthenaArray();
  wpf.DeleteAthenaArray();
  
  
  
  return;
}


