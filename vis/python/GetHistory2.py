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
rloc[0,0]=10.0
rloc[1,0]=20.0
rloc[2,0]=30.0
rloc[3,0]=40.0



#time, (Mdot, Mdotin, Mdotout)*3, meanrho, meanEr,
mdothist=np.zeros((1+3*nrloc+2,num_file))
mdothist2=np.zeros((1+3*nrloc,num_file))
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
  mdothist[0,count]=(f['Time'].value)*Crat
  mdothist2[0,count]=mdothist[0,count]
  
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
         
    # convert to cartesian coordinate
#    xcoord=radius*np.sin(theta)
#    ycoord=radius*np.cos(theta)
    ST_rho_r1=np.zeros((ntheta,num_file))
    ST_rho_r2=np.zeros((ntheta,num_file))
    ST_rho_r3=np.zeros((ntheta,num_file))
    ST_rho_r4=np.zeros((ntheta,num_file))

    ST_B3_r1=np.zeros((ntheta,num_file))
    ST_B3_r2=np.zeros((ntheta,num_file))
    ST_B3_r3=np.zeros((ntheta,num_file))
    ST_B3_r4=np.zeros((ntheta,num_file))

    ST_Er_r1=np.zeros((ntheta,num_file))
    ST_Er_r2=np.zeros((ntheta,num_file))
    ST_Er_r3=np.zeros((ntheta,num_file))
    ST_Er_r4=np.zeros((ntheta,num_file))

    ST_T_r1=np.zeros((ntheta,num_file))
    ST_T_r2=np.zeros((ntheta,num_file))
    ST_T_r3=np.zeros((ntheta,num_file))
    ST_T_r4=np.zeros((ntheta,num_file))

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
  
#  rho_1d=rho.reshape(nr*ntheta)
#  B3_1d=B3.reshape(nr*ntheta)
#  Er_1d=Er.reshape(nr*ntheta)
#  gasT_1d=gasT.reshape(nr*ntheta)
  
# get the cartesian data
#  for j in range(ntheta):
#    for i in range(nr):
    
#      rho_1d[j*nr+i]=f['rho'][j,i]
#      B3_1d[j*nr+i]=f['B3'][j,i]
#      Er_1d[j*nr+i]=f['Er'][j,i]
#      gasT_1d[j*nr+i]=gasT[j,i]

#  rho_cart=griddata(np.c_[xcoord,ycoord],rho_1d,(xmesh,ymesh),method='nearest')
#  B3_cart=griddata(np.c_[xcoord,ycoord],B3_1d,(xmesh,ymesh),method='nearest')
#  Er_cart=griddata(np.c_[xcoord,ycoord],Er_1d,(xmesh,ymesh),method='nearest')
#  gasT_cart=griddata(np.c_[xcoord,ycoord],gasT_1d,(xmesh,ymesh),method='nearest')


  ST_rho_r1[:,count]=rho[:,rloc[0,1]]
  ST_rho_r2[:,count]=rho[:,rloc[1,1]]
  ST_rho_r3[:,count]=rho[:,rloc[2,1]]
  ST_rho_r4[:,count]=rho[:,rloc[3,1]]

  ST_B3_r1[:,count]=B3[:,rloc[0,1]]
  ST_B3_r2[:,count]=B3[:,rloc[1,1]]
  ST_B3_r3[:,count]=B3[:,rloc[2,1]]
  ST_B3_r4[:,count]=B3[:,rloc[3,1]]


  ST_Er_r1[:,count]=Er[:,rloc[0,1]]
  ST_Er_r2[:,count]=Er[:,rloc[1,1]]
  ST_Er_r3[:,count]=Er[:,rloc[2,1]]
  ST_Er_r4[:,count]=Er[:,rloc[3,1]]

  ST_T_r1[:,count]=gasT[:,rloc[0,1]]
  ST_T_r2[:,count]=gasT[:,rloc[1,1]]
  ST_T_r3[:,count]=gasT[:,rloc[2,1]]
  ST_T_r4[:,count]=gasT[:,rloc[3,1]]

#####################################

# do this for each file
  # calculate the Mdot through each raii
  RhoVr=f['RhoVr'].value
  RhoVin=f['RhoVin'].value
  RhoVout=f['RhoVout'].value
  
  for i in range(0,nrloc):
    rindx=rloc[i,1]
    massflux=np.sum(RhoVr[:,rindx]*x1area[:,rindx])
    massfluxin=np.sum(RhoVin[:,rindx]*x1area[:,rindx])
    massfluxout=np.sum(RhoVout[:,rindx]*x1area[:,rindx])
    mdothist[1+i*3,count]=massflux/Medd
    mdothist[2+i*3,count]=massfluxin/Medd
    mdothist[3+i*3,count]=massfluxout/Medd
    massflux2=np.sum(RhoVr[ntlim1:ntlim2,rindx]*x1area[ntlim1:ntlim2,rindx])
    massfluxin2=np.sum(RhoVin[ntlim1:ntlim2,rindx]*x1area[ntlim1:ntlim2,rindx])
    massfluxout2=np.sum(RhoVout[ntlim1:ntlim2,rindx]*x1area[ntlim1:ntlim2,rindx])
    mdothist2[1+i*3,count]=massflux2/Medd
    mdothist2[2+i*3,count]=massfluxin2/Medd
    mdothist2[3+i*3,count]=massfluxout2/Medd
  
  rhorange=f['rho'][:,int(rloc[0,1]):int(rloc[nrloc-1,1])]
  Errange=f['Er'][:,int(rloc[0,1]):int(rloc[nrloc-1,1])]
  volrange=vol[:,int(rloc[0,1]):int(rloc[nrloc-1,1])]
  meanrho=np.average(rhorange,weights=volrange)
  meanEr=np.average(Errange,weights=volrange)

  mdothist[1+3*nrloc,count]=meanrho
  mdothist[1+3*nrloc+1,count]=meanEr

  count=count+1
  f.close()


first_file=files[0][9:13]
end_file=files[count-1][9:13]
col_names='time    Mdot1    Mdotin1    Mdotout1    Mdot2    Mdotin2    Mdotout2    Mdot3    Mdotin3    Mdotout3    Mdot4    Mdotin4    Mdotout4    meanRho    meanEr'

np.savetxt('History_'+first_file+'_'+end_file+'.txt',np.transpose(mdothist),fmt='%4.3e',header=col_names)

col_names2='time    Mdot1    Mdotin1    Mdotout1    Mdot2    Mdotin2    Mdotout2    Mdot3    Mdotin3    Mdotout3    Mdot4    Mdotin4    Mdotout4'

np.savetxt('History2_'+first_file+'_'+end_file+'.txt',np.transpose(mdothist2),fmt='%4.3e',header=col_names)


np.savetxt('ST_rho_'+'{:.0f}'.format(rloc[0,0])+'.txt',ST_rho_r1,fmt='%4.5e')
np.savetxt('ST_rho_'+'{:.0f}'.format(rloc[1,0])+'.txt',ST_rho_r2,fmt='%4.5e')
np.savetxt('ST_rho_'+'{:.0f}'.format(rloc[2,0])+'.txt',ST_rho_r3,fmt='%4.5e')
np.savetxt('ST_rho_'+'{:.0f}'.format(rloc[3,0])+'.txt',ST_rho_r4,fmt='%4.5e')

np.savetxt('ST_B3_'+'{:.0f}'.format(rloc[0,0])+'.txt',ST_B3_r1,fmt='%4.5e')
np.savetxt('ST_B3_'+'{:.0f}'.format(rloc[1,0])+'.txt',ST_B3_r2,fmt='%4.5e')
np.savetxt('ST_B3_'+'{:.0f}'.format(rloc[2,0])+'.txt',ST_B3_r3,fmt='%4.5e')
np.savetxt('ST_B3_'+'{:.0f}'.format(rloc[3,0])+'.txt',ST_B3_r4,fmt='%4.5e')

np.savetxt('ST_Er_'+'{:.0f}'.format(rloc[0,0])+'.txt',ST_Er_r1,fmt='%4.5e')
np.savetxt('ST_Er_'+'{:.0f}'.format(rloc[1,0])+'.txt',ST_Er_r2,fmt='%4.5e')
np.savetxt('ST_Er_'+'{:.0f}'.format(rloc[2,0])+'.txt',ST_Er_r3,fmt='%4.5e')
np.savetxt('ST_Er_'+'{:.0f}'.format(rloc[3,0])+'.txt',ST_Er_r4,fmt='%4.5e')

np.savetxt('ST_T_'+'{:.0f}'.format(rloc[0,0])+'.txt',ST_T_r1,fmt='%4.5e')
np.savetxt('ST_T_'+'{:.0f}'.format(rloc[1,0])+'.txt',ST_T_r2,fmt='%4.5e')
np.savetxt('ST_T_'+'{:.0f}'.format(rloc[2,0])+'.txt',ST_T_r3,fmt='%4.5e')
np.savetxt('ST_T_'+'{:.0f}'.format(rloc[3,0])+'.txt',ST_T_r4,fmt='%4.5e')

np.savetxt('height.txt',x2v,fmt='%4.3e')
np.savetxt('time.txt',mdothist[0,:],fmt='%4.4e')
