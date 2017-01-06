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

# setup latex fonts
rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True


matplotlib.rc('font', family='serif', serif='cm10')
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']




def MakeRhoVSlice(data_cart, vx_cart, vy_cart, minval, maxval, vlim1, vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname,vel=0,logscale=1):
    plots, axes = plt.subplots(figsize=(8,10),dpi=300)
    plt.xlabel('$ \\log(r/r_s)$', size = 30)
    plt.ylabel('$ z/r_s$', size = 30)
    plt.subplots_adjust(left=0.15,right=0.78,top=0.85,bottom=0.1)
#    plt.title(time,size=25,y=1.12)
    
    if vel>0:
      speed=np.sqrt(vx_cart**2+vy_cart**2)
      if vlim2 < vlim1:
         vlim2=speed.max()
      speed=np.clip(speed,vlim1,vlim2)
      logspeed=np.log10(speed)
      vcolor=axes.imshow(speed,cmap='jet',norm=LogNorm(vmin=vlim1,vmax=vlim2),origin='lower', extent=[xmin,xmax,ymin,ymax])
      velaxes = plots.add_axes([0.24, 0.89, 0.5, 0.03])
      velcbar = plots.colorbar(vcolor, cax = velaxes, orientation = 'horizontal')
      velcbar.set_label(label2, size=30, labelpad=-90)
      velcbar.ax.tick_params(labelsize=25)

      velfield=axes.streamplot(xgrid,ygrid,vx_cart,vy_cart,cmap='jet',density=1,color=logspeed)
        
      axes.set_ylim([ymin,ymax])

    if logscale > 0:
      im = axes.imshow(data_cart,cmap='RdGy_r', norm=LogNorm(vmin = minval, vmax=maxval), \
                       origin='lower', extent=[xmin,xmax,ymin,ymax])
    else:
      im = axes.imshow(data_cart,cmap='RdGy_r', vmin=minval, vmax=maxval, \
                       origin='lower', extent=[xmin,xmax,ymin,ymax])
    
    
    
    cbaxes = plots.add_axes([0.8,0.15,0.03,0.6])
    cbar=plots.colorbar(im,cax=cbaxes)
    cbar.set_label(label1, size=30)
    cbar.ax.tick_params(labelsize=25)
    axes.set_xticks([0,width/6,width/3,width/2,width*2/3,width*5/6])
    axes.yaxis.set_tick_params(labelsize=20)
    axes.xaxis.set_tick_params(labelsize=20)
    
    axes.set_aspect('auto')
    plt.savefig(outputname)
    plt.close(plots)


#################################################

# read in 2D average data
filename='Average/timeavedata_04000_05279.athdf'

crat=5694.76

rhomin=1.e-5
rhomax=20.0

kappaes=5.04655e5

Ermin=0.1
Ermax=3.e3


f=h5py.File(filename, 'r')



rho=f['rho'].value
x1v=f['x1v'].value
x2v=f['x2v'].value
vel1=f['vel1'].value
vel2=f['vel2'].value
Fr1=f['Fr1'].value
Fr2=f['Fr2'].value
Er=f['Er'].value
Fr01=f['Fr01'].value
Fr02=f['Fr02'].value
B1=f['B1'].value
B2=f['B2'].value
PB=f['PB'].value
kappa=f['kappa_s'].value+f['kappa_a'].value
kappa=kappa/kappaes

# convert 2D data to 1D array
nr=len(x1v)
ntheta=len(x2v)

radius=np.zeros(nr*ntheta)
theta=np.zeros(nr*ntheta)


for j in range(ntheta):
  for i in range(nr):
    radius[j*nr+i]=x1v[i]
    theta[j*nr+i]=x2v[j]

rho_1D=rho.reshape(nr*ntheta)
vel1_1D=vel1.reshape(nr*ntheta)
vel2_1D=vel2.reshape(nr*ntheta)
Fr1_1D=Fr1.reshape(nr*ntheta)
Fr2_1D=Fr2.reshape(nr*ntheta)
Fr01_1D=Fr01.reshape(nr*ntheta)
Fr02_1D=Fr02.reshape(nr*ntheta)
Er_1D=Er.reshape(nr*ntheta)
B1_1D=B1.reshape(nr*ntheta)
B2_1D=B2.reshape(nr*ntheta)
PB_1D=PB.reshape(nr*ntheta)
kappa_1D=kappa.reshape(nr*ntheta)

vx=vel1_1D*np.sin(theta)+vel2_1D*np.cos(theta)
vy=vel1_1D*np.cos(theta)-vel2_1D*np.sin(theta)

Bx=B1_1D*np.sin(theta)+B2_1D*np.cos(theta)
By=B1_1D*np.cos(theta)-B2_1D*np.sin(theta)


Frx=Fr1_1D*np.sin(theta)+Fr2_1D*np.cos(theta)
Fry=Fr1_1D*np.cos(theta)-Fr2_1D*np.sin(theta)

Fr0x=Fr01_1D*np.sin(theta)+Fr02_1D*np.cos(theta)
Fr0y=Fr01_1D*np.cos(theta)-Fr02_1D*np.sin(theta)

# convert to cartesian coordinate
xcoord=radius*np.sin(theta)
ycoord=radius*np.cos(theta)

width=60
height=200
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

rho_cart=griddata(np.c_[xcoord,ycoord],rho_1D,(xmesh,ymesh),method='nearest')

vx_cart=griddata(np.c_[xcoord,ycoord],vx,(xmesh,ymesh),method='nearest')
vy_cart=griddata(np.c_[xcoord,ycoord],vy,(xmesh,ymesh),method='nearest')

Bx_cart=griddata(np.c_[xcoord,ycoord],Bx,(xmesh,ymesh),method='nearest')
By_cart=griddata(np.c_[xcoord,ycoord],By,(xmesh,ymesh),method='nearest')

Er_cart=griddata(np.c_[xcoord,ycoord],Er_1D,(xmesh,ymesh),method='nearest')

Frx_cart=griddata(np.c_[xcoord,ycoord],Frx,(xmesh,ymesh),method='nearest')
Fry_cart=griddata(np.c_[xcoord,ycoord],Fry,(xmesh,ymesh),method='nearest')

Fr0x_cart=griddata(np.c_[xcoord,ycoord],Fr0x,(xmesh,ymesh),method='nearest')
Fr0y_cart=griddata(np.c_[xcoord,ycoord],Fr0y,(xmesh,ymesh),method='nearest')

PB_cart=griddata(np.c_[xcoord,ycoord],PB_1D,(xmesh,ymesh),method='nearest')

kappa_cart=griddata(np.c_[xcoord,ycoord],kappa_1D,(xmesh,ymesh),method='nearest')


vx_cart=vx_cart/crat
vy_cart=vy_cart/crat

outputname='ave_rhov.png'

label1='$\\rho/\\rho_0$'
label2='$v/c$'



vlim1=1.e-5
vlim2=1.0

MakeRhoVSlice(rho_cart, vx_cart, vy_cart, rhomin, rhomax, vlim1,vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname,1)

outputname='ave_ErFr.png'

label1='$E_r/a_rT_0^4$'
label2='$F_r$'


vlim1=1.e-3
vlim2=10.0

MakeRhoVSlice(Er_cart, Frx_cart, Fry_cart, Ermin,Ermax, vlim1,vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname,1)


outputname='ave_ErFr0.png'

label1='$E_r/a_rT_0^4$'
label2='$F_{r,0}$'


vlim1=1.e-5
vlim2=0.1

MakeRhoVSlice(Er_cart, Fr0x_cart, Fr0y_cart, Ermin,Ermax, vlim1,vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname,1)



outputname='ave_B.png'

label1='$P_B$'
label2='$B$'


MakeRhoVSlice(PB_cart, Bx_cart, By_cart, 1.e-3, 1.e4, 1.e-2, 5.e2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname,1)



outputname='ave_kappa.png'

label1='$\kappa/\kappa_{\\rm es}$'
label2=''


MakeRhoVSlice(kappa_cart,vx_cart, vy_cart, 5.e-1, 2.0, 1.e-2, 5.e2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname,0,logscale=0)


f.close()




