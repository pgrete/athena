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
from scipy.interpolate import griddata


from MidPlaneAthenaData import *

from mpi4py import MPI

# setup latex fonts
rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True


matplotlib.rc('font', family='serif', serif='cm10')
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

crat=17178.9
gm1=1.02737e5
time0=2.0*np.pi/gm1**0.5

rhomin=1.e-6
rhomax=0.001

Ermin=1.0
Ermax=5.e3

def MakeRhoVSlice(data_cart, vx_cart, vy_cart, minval, maxval, vlim1, vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,time,outputname,vel=0):
    plots, axes = plt.subplots(figsize=(10.,8.4),dpi=300)
    plt.xlabel('$ x/r_0$', size = 30)
    plt.ylabel('$ y/r_0$', size = 30)
    plt.subplots_adjust(left=0.16,right=0.78,top=0.85,bottom=0.1)
    plt.title(time,size=25,y=1.12)
    
    if vel>0:
      speed=np.sqrt(vx_cart**2+vy_cart**2)
      if vlim2 < vlim1:
         vlim2=speed.max()
      speed=np.clip(speed,vlim1,vlim2)
      logspeed=np.log10(speed)
      vcolor=axes.imshow(speed,cmap='jet',norm=LogNorm(vmin=vlim1,vmax=vlim2),origin='lower', extent=[xmin,xmax,ymin,ymax])
      velaxes = plots.add_axes([0.24, 0.89, 0.5, 0.03])
      velcbar = plots.colorbar(vcolor, cax = velaxes, orientation = 'horizontal')
      velcbar.ax.tick_params(labelsize=20)
      if speed.max() > vlim1:
        velfield=axes.streamplot(xgrid,ygrid,vx_cart,vy_cart,cmap='jet',density=1,color=logspeed)
      else:
        velfield=axes.streamplot(xgrid,ygrid,vx_cart,vy_cart,cmap='jet',density=1)
        
      axes.set_ylim([ymin,ymax])

       
    im = axes.imshow(data_cart,cmap='RdGy_r', norm=LogNorm(vmin = minval, vmax=maxval), \
                     origin='lower', extent=[xmin,xmax,ymin,ymax])
                     
    cbaxes = plots.add_axes([0.8,0.15,0.03,0.6])
    cbar=plots.colorbar(im,cax=cbaxes)
    cbar.set_label(label1, size=30)
    cbar.ax.tick_params(labelsize=25)
#    axes.set_xticks([-50,-30.,-10.,10.,30.,50.])
    axes.yaxis.set_tick_params(labelsize=20)
    axes.xaxis.set_tick_params(labelsize=20)
    
    axes.set_aspect('auto')
    plt.savefig(outputname)
    plt.close(plots)


def makemovies(ni,no,minval,maxval,vlim1,vlim2,width=23.44,vel=0,var1='rho',var2='vel1',var3='vel3',label1='$\\rho/\\rho_0$',label2='$v$'):
  
    for i in range(ni,no+1):
        print i
        filename='disk.out1.'+'{:05d}'.format(i)+'.athdf'
        data=MidPlaneData(filename)

        rhoslice=data[var1]

        x1v=data['x1v']
        x3v=data['x3v']

        nr=x1v.size
        ntheta=x3v.size
        

     
        rcoord=np.zeros(ntheta*nr)
        phi=np.zeros(ntheta*nr)
        
        for n in range(0,nr):
          rcoord[n*ntheta:(n+1)*ntheta]=x1v[n]
        
        for n in range(0,nr):
          phi[n*ntheta:(n+1)*ntheta]=x3v[0:ntheta]
    
        vel1=data[var2]
        vel2=data[var3]
        
        rhoslice=np.reshape(rhoslice,nr*ntheta,order='F')
        vel1=np.reshape(vel1,nr*ntheta,order='F')
        vel2=np.reshape(vel2,nr*ntheta,order='F')
        
        

        vx=vel1*np.cos(phi)-vel2*np.sin(phi)
        vy=vel1*np.sin(phi)+vel2*np.cos(phi)

# convert to cartesian coordinate
        xcoord=rcoord*np.cos(phi)
        ycoord=rcoord*np.sin(phi)
        

#create grid for plotting

        nx=2*nr
        ny=2*nr
        xmin=-width
# xmax=np.max(rcoord)
        xmax=width
        ymin=-xmax
        ymax=xmax
        xgrid=np.linspace(xmin, xmax, nx)
        ygrid=np.linspace(ymin, ymax, ny)

        xmesh,ymesh=np.meshgrid(xgrid,ygrid)
       
        rmesh=(xmesh**2.0+ymesh**2.0)**0.5

        rindix=rmesh < 1.0
        rindix2=rmesh > width
        

        rho_cart=griddata(np.c_[xcoord,ycoord],rhoslice,(xmesh,ymesh),method='nearest')

        vx_cart=griddata(np.c_[xcoord,ycoord],vx,(xmesh,ymesh),method='nearest')
        vy_cart=griddata(np.c_[xcoord,ycoord],vy,(xmesh,ymesh),method='nearest')
        
        rho_cart[rindix] = minval
        vx_cart[rindix] = 0.0
        vy_cart[rindix] = 0.0

        rho_cart[rindix2] = minval
        vx_cart[rindix2] = 0.0
        vy_cart[rindix2] = 0.0
        
# only scale to speed of light for velocity
#        if var2=='vel1':
#           vx_cart=vx_cart/crat
#           vy_cart=vy_cart/crat

        outputname='disk.'+'{:05d}'.format(i)+'_'+var1+'.png'

        labelname='$\\rho/\\rho_0$'

        time='${\\rm time}='+"%4.2f"%(data['Time']/time0)+'{\ t_0}$'

        MakeRhoVSlice(rho_cart, vx_cart, vy_cart, minval,maxval, vlim1,vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,time,outputname,vel)

#################################################


# The main program

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

ni=461
no=461

for i in range(ni,no+1,nprocs):
  fi=i+rank
  print fi, rank
  makemovies(fi,fi,rhomin,rhomax,1.0,1.e4,vel=1)
#  makemovies(fi,fi,Ermin,Ermax,1.e-3,10.0,vel=0,var1='Er',label1='$E_r/a_rT_0^4$')


