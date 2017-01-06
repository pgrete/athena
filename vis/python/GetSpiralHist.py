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




# setup latex fonts
#rcParams['text.usetex']=True
#rcParams['text.latex.unicode']=True


#matplotlib.rc('font', family='serif', serif='cm10')
#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

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


nrloc=4
rloc=np.zeros((nrloc,2))
rloc[0,0]=5.0
rloc[1,0]=10.0
rloc[2,0]=15.0
rloc[3,0]=20.0



#time, (Mdot, Mdotin, Mdotout)*3, meanrho, meanEr,

spiralhist=np.zeros((1+4*nrloc,num_file))
count=0
vol={}
x1area={}
nr=1
ntheta=1
tot_vol=0.0

rho_1d={}
B3_1d={}
Er_1d={}
gasT_1d={}
radius={}
theta={}

xcoord={}
ycoord={}

#########################
# The cartesian grid
width=40
height=400
nx=128
ny=int(height*nx/width)

xmin=0
# xmax=np.max(rcoord)
xmax=width
ymin=-height/2
ymax=height/2
xgrid=np.linspace(xmin, xmax, nx)
ygrid=np.linspace(ymin, ymax, ny)

xmesh,ymesh=np.meshgrid(xgrid,ygrid)

# find the location in the cartesian grid
#rloc_cart=np.zeros((nrloc,2))
#rloc_cart[:,0]=rloc[:,0]

#for i in range(0,nrloc):
#    rloc_cart[i,1]=np.abs(xgrid - rloc_cart[i,0]).argmin()


#######################

first=0

for filename in files:
  print filename
  f=h5py.File(filename, 'r')
  x1f=f['x1f'].value
  x2f=f['x2f'].value
  x1v=f['x1v'].value
  x2v=f['x2v'].value
  
  # time

  spiralhist[0,count]=(f['Time'].value)*Crat
  
  quantities=f.keys()
#####################################
# Only need to calculate once
  if first==0:
    # get the position of rloc
    for i in range(0,nrloc):
      rloc[i,1]=np.abs(x1v - rloc[i,0]).argmin()
    # get the volume, data size
    nr=len(x1v)
    ntheta=len(x2v)
    vol=np.zeros((ntheta,nr))
    x1area=np.zeros((ntheta,nr+1))

    radius=np.zeros(nr*ntheta)
    theta=np.zeros(nr*ntheta)
    
    for j in range(0,ntheta):
       for i in range(0,nr):
         cellsize=vol_func(x1f[i],x1f[i+1],x2f[j],x2f[j+1])
         vol[j,i]=cellsize
         tot_vol=tot_vol+cellsize
         radius[j*nr+i]=x1v[i]
         theta[j*nr+i]=x2v[j]
         
    ST_curvgrads=np.zeros((nr,num_file))
    ST_sigma_curvgrads=np.zeros((nr,num_file))
    ST_curvtrho=np.zeros((nr,num_file))
    ST_reynolds=np.zeros((nr,num_file))


    for j in range(0,ntheta):
      for i in range(0,nr+1):
         x1area[j,i]=rarea_func(x1f[i],x2f[j],x2f[j+1])
    first=1
  
  #only calculate Mdot near the midplane
  ntlim1=ntheta/4
  ntlim2=ntheta*3/4
#####################################

  gasT=np.divide(f['pgas'].value,f['rho'].value)
  rho=f['rho'].value
  B3=f['B3'].value
  Er=f['Er'].value
#  meanrhovrvphi=f['meanrhovrvphi'].value
  rhovrvphi=f['rhovrvphi'].value

  rhovr=f['RhoVr'].value
  rhovphi=f['RhoVphi'].value
  reynolds=rhovrvphi-rhovr*rhovphi/rho
#  reynolds=rhovrvphi-meanrhovrvphi
  curvgrads=f['curvgrads'].value
  curvgradssq=f['curvgradssq'].value
  sigma_curvgrads=(curvgradssq-curvgrads*curvgrads)**0.5
  
  curvtrho=f['curvtrho'].value
  
  mid_curvtrho=curvtrho[ntheta/2,:]
  mid_reynolds=reynolds[ntheta/2,:]
  mid_curvgrads=curvgrads[ntheta/2,:]
  mid_sigma_curvgrads=sigma_curvgrads[ntheta/2,:]
  
  ST_curvtrho[:,count]=mid_curvtrho
  ST_reynolds[:,count]=mid_reynolds
  ST_curvgrads[:,count]=mid_curvgrads
  ST_sigma_curvgrads[:,count]=mid_sigma_curvgrads
  
 
  for i in range(0,nrloc):
    rindx=np.int(rloc[i,1])
    spiralhist[1+i*4,count]=mid_reynolds[rindx]
    spiralhist[2+i*4,count]=mid_curvgrads[rindx]
    spiralhist[3+i*4,count]=mid_sigma_curvgrads[rindx]
    spiralhist[4+i*4,count]=mid_curvtrho[rindx]
  


  count=count+1
  f.close()


first_file=files[0][8:13]
end_file=files[count-1][8:13]

col_names3='time    Rey1    PV1    SigmaPV1    curv1    Rey2    PV2    SigmaPV2    curv2    Rey3    PV3    SigmaPV3    curv3'

np.savetxt('His_Spiral_'+first_file+'_'+end_file+'.txt',np.transpose(spiralhist),fmt='%4.3e',header=col_names3)

np.savetxt('ST_curvgrads.txt',ST_curvgrads,fmt='%4.5e')
np.savetxt('ST_sigma_curvgrads.txt',ST_sigma_curvgrads,fmt='%4.5e')
np.savetxt('ST_reynolds.txt',ST_reynolds,fmt='%4.5e')
np.savetxt('ST_curvtrho.txt',ST_curvtrho,fmt='%4.5e')

np.savetxt('radius.txt',x1v,fmt='%4.3e')
np.savetxt('time.txt',spiralhist[0,:],fmt='%4.4e')
