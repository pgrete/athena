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

import h5py

from mpi4py import MPI

# setup latex fonts
#rcParams['text.usetex']=True
#rcParams['text.latex.unicode']=True


#matplotlib.rc('font', family='serif', serif='cm10')
#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

crat=5694.76

rhomin=1.e-5
rhomax=1.0

Ermin=0.1
Ermax=5.e2




def MakeRhoVSlice(data_cart, vx_cart, vy_cart, minval, maxval, vlim1, vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,time,outputname,vel=0):
    plots, axes = plt.subplots(figsize=(6,10),dpi=300)
    plt.xlabel('$ r/r_g$', size = 30)
    plt.ylabel('$ z/r_g$', size = 30)
    plt.subplots_adjust(left=0.22,right=0.70,top=0.85,bottom=0.1)
    plt.title(time,size=25,y=1.12)
    
    if vel>0:
      speed=np.sqrt(vx_cart**2+vy_cart**2)
      if vlim2 < vlim1:
         vlim2=speed.max()
      speed=np.clip(speed,vlim1,vlim2)
      logspeed=np.log10(speed)
      vcolor=axes.imshow(speed,cmap='jet',norm=LogNorm(vmin=vlim1,vmax=vlim2),origin='lower', extent=[xmin,xmax,ymin,ymax])
      velaxes = plots.add_axes([0.22, 0.89, 0.5, 0.03])
      velcbar = plots.colorbar(vcolor, cax = velaxes, orientation = 'horizontal')
      velcbar.ax.tick_params(labelsize=20)
      if speed.max() > vlim1:
        velfield=axes.streamplot(xgrid,ygrid,vx_cart,vy_cart,cmap='jet',density=1,color=logspeed,arrowsize=3.0)
      else:
        velfield=axes.streamplot(xgrid,ygrid,vx_cart,vy_cart,cmap='jet',density=1,arrowsize=3.0)
        
      axes.set_ylim([ymin,ymax])

       
    im = axes.imshow(data_cart,cmap='RdYlBu_r', norm=LogNorm(vmin = minval, vmax=maxval), \
                     origin='lower', extent=[xmin,xmax,ymin,ymax])
                     
    cbaxes = plots.add_axes([0.72,0.15,0.03,0.6])
    cbar=plots.colorbar(im,cax=cbaxes)
    cbar.set_label(label1, size=30)
    cbar.ax.tick_params(labelsize=25)
    axes.set_xticks([0,20,40,60,80])
    axes.yaxis.set_tick_params(labelsize=20)
    axes.xaxis.set_tick_params(labelsize=20)
    
    axes.set_aspect('auto')
    plt.savefig(outputname)
    plt.close(plots)


def makemovies(ni,no,minval,maxval,vlim1,vlim2,width=80,ywidth=60,res=4,vel=0,var1='rho',var2='vel1',var3='vel2',label1='$\\rho/\\rho_0$',label2='$v/c$'):
  
    for nf in range(ni,no+1):
        print nf
        filename='average_'+'{:05d}'.format(nf)+'.athdf'
        f=h5py.File(filename, 'r')
        rho=f[var1].value
        vel1=f[var2].value
        vel2=f[var3].value

        radius=f['x1v'].value
        #change from rs to rg
        radius=2*radius
        theta=f['x2v'].value

        nr=radius.size
        ntheta=theta.size


#create grid for plotting
        nx=nr
        ny=ywidth*2*nx/width
        xmin=0
# xmax=np.max(rcoord)
        xmax=width
        ymin=-ywidth
        ymax=ywidth
        xgrid=np.linspace(xmin, xmax, nx)
        ygrid=np.linspace(ymin, ymax, ny)

        xmesh,ymesh=np.meshgrid(xgrid,ygrid)

        rmesh=(xmesh**2.0+ymesh**2.0)**0.5

        rindix=rmesh < 4.0

        # convert to 1D array
        rho_1D=np.zeros(nr*ntheta)
        vel1_1D=np.zeros(nr*ntheta)
        vel2_1D=np.zeros(nr*ntheta)
        radius_1D=np.zeros(nr*ntheta)
        theta_1D=np.zeros(nr*ntheta)

        for j in range(ntheta):
          for i in range(nr):
            rho_1D[j*nr+i]=rho[j,i]
            vel1_1D[j*nr+i]=vel1[j,i]
            vel2_1D[j*nr+i]=vel2[j,i]
            radius_1D[j*nr+i]=radius[i]
            theta_1D[j*nr+i]=theta[j]


        vx=vel1_1D*np.sin(theta_1D)+vel2_1D*np.cos(theta_1D)
        vy=vel1_1D*np.cos(theta_1D)-vel2_1D*np.sin(theta_1D)

# convert to cartesian coordinate
        xcoord=radius_1D*np.sin(theta_1D)
        ycoord=radius_1D*np.cos(theta_1D)



        rho_cart=griddata(np.c_[xcoord,ycoord],rho_1D,(xmesh,ymesh),method='nearest')

        rho_cart[rindix] = minval

        vx_cart=griddata(np.c_[xcoord,ycoord],vx,(xmesh,ymesh),method='nearest')
        vy_cart=griddata(np.c_[xcoord,ycoord],vy,(xmesh,ymesh),method='nearest')

        vx_cart[rindix] = 0.0
        vy_cart[rindix] = 0.0


# only scale to speed of light for velocity
        if var2=='vel1':
           vx_cart=vx_cart/crat
           vy_cart=vy_cart/crat

        outputname='disk_ave.'+'{:05d}'.format(nf)+'_'+var1+'.png'


        time='${\\rm time}='+"%4.2f"%(2.0*f['Time'].value*crat)+'{\ r_g/c}$'


        MakeRhoVSlice(rho_cart, vx_cart, vy_cart, minval,maxval, vlim1,vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,time,outputname,vel)

#################################################


# The main program

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

ni=1055
no=1055

for i in range(ni,no+1,nprocs):
  fi=i+rank
  print fi, rank
  makemovies(fi,fi,rhomin,rhomax,1.e-2,1.0,vel=1)
  makemovies(fi,fi,Ermin,Ermax,1.e-3,10.0,vel=1,var1='Er',var2='Bcc1',var3='Bcc2',label1='$E_r/a_rT_0^4$',label2='B')


