import matplotlib
matplotlib.use('Agg')

import yt
import sys
sys.settrace
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pylab import *
import struct
import array
import os
import glob
import h5py
from scipy.interpolate import griddata

from mpi4py import MPI



# setup latex fonts
#rcParams['text.usetex']=True
#rcParams['text.latex.unicode']=True


#matplotlib.rc('font', family='serif', serif='cm10')
#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']



def Convertdata(filename,Frflag=0):

    file=h5py.File(filename,'r')


    x1v=file['x1v'].value # in unit of r_s
    x1v=x1v*2
    x2v=file['x2v'].value
    x1f=file['x1f'].value  # in unit of r_s
    x1f=x1f*2
    x2f=file['x2f'].value


    gamma=5.0/3.0
    #
    nr=len(x1v)
    ntheta=len(x2v)


    if Frflag > 0:
      Fr1=file['Fr1new'].value
      Fr2=file['Fr2new'].value
    else:
      Fr1=file['Fr1'].value
      Fr2=file['Fr2'].value

    Fr01=file['Fr01'].value
    Fr02=file['Fr02'].value

    vel1=file['vel1'].value
    vel2=file['vel2'].value

    rho=file['rho'].value

    kappa=(file['Sigma_a'].value+file['Sigma_s'].value)/rho
    Fr01kappa=(file['Fr01Sigma'].value)/rho
    Fr02kappa=(file['Fr02Sigma'].value)/rho


    Ekin1=file['Ekin1'].value
    Ekin2=file['Ekin2'].value
    Ekin3=file['Ekin3'].value

    Ekin=Ekin1+Ekin2+Ekin3


    rhovr=file['RhoVr'].value
    rhovtheta=file['RhoVtheta'].value

    Er=file['Er'].value

    vr_rho=rhovr/rho
    vtheta_rho=rhovtheta/rho

    vr=file['vel1'].value
    vtheta=file['vel2'].value
    tgas=file['pgas'].value/rho

    rho_1D=np.zeros(nr*ntheta)
    Er_1D=np.zeros(nr*ntheta)
    vr_1D=np.zeros(nr*ntheta)
    vtheta_1D=np.zeros(nr*ntheta)
    vr_rho_1D=np.zeros(nr*ntheta)
    vtheta_rho_1D=np.zeros(nr*ntheta)
    tgas_1D=np.zeros(nr*ntheta)
    Ekin_1D=np.zeros(nr*ntheta)
    kappa_1D=np.zeros(nr*ntheta)
    Fr01kappa_1D=np.zeros(nr*ntheta)
    Fr02kappa_1D=np.zeros(nr*ntheta)
    Fr1_1D=np.zeros(nr*ntheta)
    Fr2_1D=np.zeros(nr*ntheta)
    Fr01_1D=np.zeros(nr*ntheta)
    Fr02_1D=np.zeros(nr*ntheta)
    radius=np.zeros(ntheta*nr)
    theta=np.zeros(ntheta*nr)
    rhovr_1D=np.zeros(ntheta*nr)
    rhovtheta_1D=np.zeros(ntheta*nr)


    for j in range(ntheta):
      for i in range(nr):
        rho_1D[j*nr+i]=rho[j,i]
        Er_1D[j*nr+i]=Er[j,i]
        vr_1D[j*nr+i]=vr[j,i]
        vtheta_1D[j*nr+i]=vtheta[j,i]
        vr_rho_1D[j*nr+i]=vr_rho[j,i]
        vtheta_rho_1D[j*nr+i]=vtheta_rho[j,i]
        tgas_1D[j*nr+i]=tgas[j,i]
        Ekin_1D[j*nr+i]=Ekin[j,i]
        kappa_1D[j*nr+i]=kappa[j,i]
        Fr01kappa_1D[j*nr+i]=Fr01kappa[j,i]
        Fr02kappa_1D[j*nr+i]=Fr02kappa[j,i]
        Fr1_1D[j*nr+i]=Fr1[j,i]
        Fr2_1D[j*nr+i]=Fr2[j,i]
        Fr01_1D[j*nr+i]=Fr01[j,i]
        Fr02_1D[j*nr+i]=Fr02[j,i]
        rhovr_1D[j*nr+i]=rhovr[j,i]
        rhovtheta_1D[j*nr+i]=rhovtheta[j,i]
        radius[j*nr+i]=x1v[i]
        theta[j*nr+i]=x2v[j]


    # convert r-theta to r-z
    vx=vr_1D*np.sin(theta)+vtheta_1D*np.cos(theta)
    vy=vr_1D*np.cos(theta)-vtheta_1D*np.sin(theta)

    vx_rho=vr_rho_1D*np.sin(theta)+vtheta_rho_1D*np.cos(theta)
    vy_rho=vr_rho_1D*np.cos(theta)-vtheta_rho_1D*np.sin(theta)

    rhovx=rhovr_1D*np.sin(theta)+rhovtheta_1D*np.cos(theta)
    rhovy=rhovr_1D*np.cos(theta)-rhovtheta_1D*np.sin(theta)

    Frx=Fr1_1D*np.sin(theta)+Fr2_1D*np.cos(theta)
    Fry=Fr1_1D*np.cos(theta)-Fr2_1D*np.sin(theta)

    Fr0x=Fr01_1D*np.sin(theta)+Fr02_1D*np.cos(theta)
    Fr0y=Fr01_1D*np.cos(theta)-Fr02_1D*np.sin(theta)

    Fr0xkappa=Fr01kappa_1D*np.sin(theta)+Fr02kappa_1D*np.cos(theta)
    Fr0ykappa=Fr01kappa_1D*np.cos(theta)-Fr02kappa_1D*np.sin(theta)

    # convert to cartesian coordinate
    xcoord=radius*np.sin(theta)
    ycoord=radius*np.cos(theta)


    width=100
    height=800
    nx=512
    ny=int(2*height*nx/width)

    xmin=4

    xmax=width
    ymin=-height
    ymax=height

    xgrid=np.linspace(xmin, xmax, nx)
    ygrid=np.linspace(ymin, ymax, ny)

    xmesh,ymesh=np.meshgrid(xgrid,ygrid)

    rho_cart=griddata(np.c_[xcoord,ycoord],rho_1D,(xmesh,ymesh),method='nearest')
    Er_cart=griddata(np.c_[xcoord,ycoord],Er_1D,(xmesh,ymesh),method='nearest')
    vx_cart=griddata(np.c_[xcoord,ycoord],vx,(xmesh,ymesh),method='nearest')
    vy_cart=griddata(np.c_[xcoord,ycoord],vy,(xmesh,ymesh),method='nearest')
    vx_rho_cart=griddata(np.c_[xcoord,ycoord],vx_rho,(xmesh,ymesh),method='nearest')
    vy_rho_cart=griddata(np.c_[xcoord,ycoord],vy_rho,(xmesh,ymesh),method='nearest')
    Frx_cart=griddata(np.c_[xcoord,ycoord],Frx,(xmesh,ymesh),method='nearest')
    Fry_cart=griddata(np.c_[xcoord,ycoord],Fry,(xmesh,ymesh),method='nearest')
    Fr0x_cart=griddata(np.c_[xcoord,ycoord],Fr0x,(xmesh,ymesh),method='nearest')
    Fr0y_cart=griddata(np.c_[xcoord,ycoord],Fr0y,(xmesh,ymesh),method='nearest')
    Fr0xkappa_cart=griddata(np.c_[xcoord,ycoord],Fr0xkappa,(xmesh,ymesh),method='nearest')
    Fr0ykappa_cart=griddata(np.c_[xcoord,ycoord],Fr0ykappa,(xmesh,ymesh),method='nearest')
    Ekin_cart=griddata(np.c_[xcoord,ycoord],Ekin_1D,(xmesh,ymesh),method='nearest')
    tgas_cart=griddata(np.c_[xcoord,ycoord],tgas_1D,(xmesh,ymesh),method='nearest')
    rhovx_cart=griddata(np.c_[xcoord,ycoord],rhovx,(xmesh,ymesh),method='nearest')
    rhovy_cart=griddata(np.c_[xcoord,ycoord],rhovy,(xmesh,ymesh),method='nearest')

    data=[xgrid,ygrid,rho_cart,Er_cart,vx_cart,vy_cart,vx_rho_cart,vy_rho_cart, \
           + Frx_cart,Fry_cart,Fr0x_cart,Fr0y_cart,Ekin_cart,tgas_cart,rhovx_cart,rhovy_cart]

    return data




Crat=5694.76
Prat=43.6907
kappaes=5.04655e5
#rs=2GM/c^2==1
Medd = 20*np.pi*Crat/kappaes
Ledd = 2* np.pi * Crat * Crat * Crat/kappaes

vol_func = lambda rm,rp,thetam,thetap: \
           (1.0/3.0)*(rp**3-rm**3) * abs(np.cos(thetam)-np.cos(thetap)) * 2.0*np.pi

rarea_func = lambda r,thetam,thetap: \
           r**2.0 * abs(np.cos(thetam)-np.cos(thetap)) * 2.0*np.pi

#calculate Mdot at 5, 10, 15, 20 r_s

files=sorted(glob.glob('ave*athdf'))

num_file=len(files)


#######################


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

for i in range(10175,num_file,nprocs):
  fi=i+rank
  filename=files[fi]
  print fi, rank, filename
  out_name='Lum_'+'{:05d}'.format(fi)+'.txt'

#data=[[0]xgrid,[1]ygrid,[2]rho_cart,[3]Er_cart,[4]vx_cart,[5]vy_cart,[6]vx_rho_cart,
#[7]vy_rho_cart, [8]Frx_cart,[9]Fry_cart,[10]Fr0x_cart,[11]Fr0y_cart,
#[12]Ekin_cart,[13]tgas_cart,[14]rhovx_cart,[15]rhovy_cart]


  data=Convertdata(filename,Frflag=0)
  rho_cart=data[2]
  Er_cart=data[3]
  vx_cart=data[4]
  vy_cart=data[5]
  vx_rho_cart=data[6]
  vy_rho_cart=data[7]
  Fry_cart=data[9]
  Fry0_cart=data[11]
  tgas_cart=data[13]
  Ekin_cart=data[12]
  rhovx_cart=data[14]
  rhovy_cart=data[15]

  radius=data[0]
  height=data[1]

  nh=len(height)
  nr=len(radius)

  Lum=np.zeros(nh)
  LEkin=np.zeros(nh)
  Mout=np.zeros(nh)

  area=np.zeros(nr)

  outputarray=np.zeros((nh,12))


  rloc1=15
  rloc2=20
  rloc3=30

  nr1=np.abs(radius-rloc1).argmin() 
  nr2=np.abs(radius-rloc2).argmin()
  nr3=np.abs(radius-rloc3).argmin()

  dr=radius[1]-radius[0]
  # convert to unit rs^2
  for i in range(nr-1):
    area[i]=2.0*np.pi*radius[i]*(radius[i+1]-radius[i])*0.25

  for j in range(nh):
    Lum[j]=0.0
    LEkin[j]=0.0
    Mout[j]=0.0
    for i in range(nr1):
      Lum[j]=Lum[j]+np.sign(height[j])*Fry_cart[j,i]*area[i]*Prat*Crat
      LEkin[j]=LEkin[j]+np.sign(height[j])*vy_rho_cart[j,i]*area[i]*Ekin_cart[j,i]
      Mout[j]=Mout[j]+np.sign(height[j])*rhovy_cart[j,i]*area[i]
    Lum[j]=Lum[j]/Ledd
    LEkin[j]=LEkin[j]/Ledd
    Mout[j]=Mout[j]/Medd


  outputarray[:,0]=Lum
  outputarray[:,1]=LEkin
  outputarray[:,2]=Mout
  outputarray[:,3]=rho_cart[:,nr1]

  for j in range(nh):
    Lum[j]=0.0
    LEkin[j]=0.0
    Mout[j]=0.0
    for i in range(nr2):
      Lum[j]=Lum[j]+np.sign(height[j])*Fry_cart[j,i]*area[i]*Prat*Crat
      LEkin[j]=LEkin[j]+np.sign(height[j])*vy_rho_cart[j,i]*area[i]*Ekin_cart[j,i]
      Mout[j]=Mout[j]+np.sign(height[j])*rhovy_cart[j,i]*area[i]
    Lum[j]=Lum[j]/Ledd
    LEkin[j]=LEkin[j]/Ledd
    Mout[j]=Mout[j]/Medd

  outputarray[:,4]=Lum
  outputarray[:,5]=LEkin
  outputarray[:,6]=Mout
  outputarray[:,7]=rho_cart[:,nr2]

  for j in range(nh):
    Lum[j]=0.0
    LEkin[j]=0.0
    Mout[j]=0.0
    for i in range(nr3):
      Lum[j]=Lum[j]+np.sign(height[j])*Fry_cart[j,i]*area[i]*Prat*Crat
      LEkin[j]=LEkin[j]+np.sign(height[j])*vy_rho_cart[j,i]*area[i]*Ekin_cart[j,i]
      Mout[j]=Mout[j]+np.sign(height[j])*rhovy_cart[j,i]*area[i]
    Lum[j]=Lum[j]/Ledd
    LEkin[j]=LEkin[j]/Ledd
    Mout[j]=Mout[j]/Medd

  outputarray[:,8]=Lum
  outputarray[:,9]=LEkin
  outputarray[:,10]=Mout
  outputarray[:,11]=rho_cart[:,nr3]


  np.savetxt(out_name,outputarray,fmt='%4.5e')
  np.savetxt('Lum_height.txt',height,fmt='%4.5e')
  #############################################################
  #############################################################

