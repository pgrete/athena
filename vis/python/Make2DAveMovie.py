import matplotlib
matplotlib.use('Agg')

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


prat=0.000155307
crat=17178.9
kappaes=9.38e3
#binary frequency Omega
omega0=1.79925
gm=1.02737e5

time0=3.49211


vol_func = lambda rm,rp,thetam,thetap: \
            1.0/3.0*(rp**3-rm**3) * abs(np.cos(thetam)-np.cos(thetap)) * 2.0*np.pi

rarea_func = lambda r,thetam,thetap: \
           r**2.0 * abs(np.cos(thetam)-np.cos(thetap)) * 2.0*np.pi

#the tabulated color
tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
             (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
             (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
             (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
             (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]   

for i in range(len(tableau20)):    
    r, g, b = tableau20[i]    
    tableau20[i] = (r / 255., g / 255., b / 255.) 

def MakeRhoVSlice(data_cart, vx_cart, vy_cart, minval, maxval, vlim1, vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname,vel=0,logscale=1,xticks=None,taux=None,tauy=None,title=None):
    plots, axes = plt.subplots(figsize=(8,10),dpi=300)
    plt.xlabel('$r\\sin\\theta/r_0$', size = 30)
    plt.ylabel('$r\\cos\\theta/r_0$', size = 30)
    plt.subplots_adjust(left=0.15,right=0.8,top=0.85,bottom=0.1)

    if title is not None:
      plt.title(title,size=25,y=1.12)
    
    if vel>0:
      speed=np.sqrt(vx_cart**2+vy_cart**2)
      if vlim2 < vlim1:
         vlim2=speed.max()
      speed=np.clip(speed,vlim1,vlim2)
      logspeed=np.log10(speed)
      vcolor=axes.imshow(speed,cmap='jet',norm=LogNorm(vmin=vlim1,vmax=vlim2),origin='lower', extent=[xmin,xmax,ymin,ymax])
      velaxes = plots.add_axes([0.24, 0.9, 0.5, 0.03])
      velcbar = plots.colorbar(vcolor, cax = velaxes, orientation = 'horizontal')
      if title is None:
        velcbar.set_label(label2, size=30, labelpad=-90)
      velcbar.ax.tick_params(labelsize=25)

      velfield=axes.streamplot(xgrid,ygrid,vx_cart,vy_cart,cmap='jet',arrowsize=3,density=0.6,linewidth=1.5,color=logspeed)
        
      axes.set_ylim([ymin,ymax])

    if taux is not None:
      axes.plot(taux,tauy,color='black',linestyle='dashed',linewidth=3.0)

    if logscale > 0:
      im = axes.imshow(data_cart,cmap='RdGy_r', norm=LogNorm(vmin = minval, vmax=maxval), \
                       origin='lower', extent=[xmin,xmax,ymin,ymax])
    else:
      im = axes.imshow(data_cart,cmap='RdGy_r', vmin=minval, vmax=maxval, \
                       origin='lower', extent=[xmin,xmax,ymin,ymax])
    
    
    
    cbaxes = plots.add_axes([0.81,0.15,0.03,0.6])
    cbar=plots.colorbar(im,cax=cbaxes)
    cbar.set_label(label1, size=30)
    cbar.ax.tick_params(labelsize=25)
    if xticks is not None:
      axes.set_xticks(xticks)
    axes.yaxis.set_tick_params(labelsize=20)
    axes.xaxis.set_tick_params(labelsize=20)
  

    axes.set_aspect('auto')

    plt.savefig(outputname)
    plt.close(plots)



def PlotAveStructure(filename,rhomin,rhomax,Ermin,Ermax,PBmin,PBmax,Bin,Bmax,Frmin,Frmax,outputname1,outputname2,outputname3=None,outputname4=None,vecname1='Fr1',vecname2='Fr2',vecname3='RhoVr',vecname4='RhoVtheta'):
    f=h5py.File(filename, 'r')



    rho=f['rho'].value
    x1v=f['x1v'].value
    x2v=f['x2v'].value
    x1f=f['x1f'].value
    x2f=f['x2f'].value
    vel1=f[vecname3].value
    vel2=f[vecname4].value
    Fr1=f[vecname1].value
    Fr2=f[vecname2].value
    Er=f['Er'].value
    B1=f['Bcc1'].value
    B2=f['Bcc2'].value
    PB=f['PB'].value
    sigma_s=f['Sigma_s'].value
    sigma_a=f['Sigma_a'].value

  
    time=f['Time'].value

    sigma=sigma_s+sigma_a

    kappa=sigma/(rho*kappaes)

    if vecname3=='RhoVr':
      vel1=vel1/rho
    
    if vecname4=='RhoVtheta':
      vel2=vel2/rho

    #Get location of tau=1 from rotation axis

    #convert to 1D
    nr=len(x1v)
    ntheta=len(x2v)
    dx1=np.zeros(nr)
    cellvol=np.zeros((ntheta,nr))
    rho_1D=np.zeros(nr*ntheta)
    vel1_1D=np.zeros(nr*ntheta)
    vel2_1D=np.zeros(nr*ntheta)
    Fr1_1D=np.zeros(nr*ntheta)
    Fr2_1D=np.zeros(nr*ntheta)
    Er_1D=np.zeros(nr*ntheta)
    radius=np.zeros(nr*ntheta)
    theta=np.zeros(nr*ntheta)
    B1_1D=np.zeros(nr*ntheta)
    B2_1D=np.zeros(nr*ntheta)
    PB_1D=np.zeros(nr*ntheta)
    sigma_1D=np.zeros(nr*ntheta)
    kappa_1D=np.zeros(nr*ntheta)
    dx_1D=np.zeros(nr*ntheta)

    for i in range(nr):
       dx1[i]=x1f[i+1]-x1f[i]

    for j in range(ntheta):
      for i in range(nr):
        rho_1D[j*nr+i]=rho[j,i]
        vel1_1D[j*nr+i]=vel1[j,i]
        vel2_1D[j*nr+i]=vel2[j,i]
        Fr1_1D[j*nr+i]=Fr1[j,i]
        Fr2_1D[j*nr+i]=Fr2[j,i]
        Er_1D[j*nr+i]=Er[j,i]
        B1_1D[j*nr+i]=B1[j,i]
        B2_1D[j*nr+i]=B2[j,i]
        PB_1D[j*nr+i]=PB[j,i]
        sigma_1D[j*nr+i]=sigma[j,i]
        kappa_1D[j*nr+i]=kappa[j,i]
        cellvol[j,i]=vol_func(x1f[i],x1f[i+1],x2f[j],x2f[j+1])
        radius[j*nr+i]=2*x1v[i]
        theta[j*nr+i]=x2v[j]
        dx_1D[j*nr+i]=dx1[i]*np.sin(x2v[j])




    vx=vel1_1D*np.sin(theta)+vel2_1D*np.cos(theta)
    vy=vel1_1D*np.cos(theta)-vel2_1D*np.sin(theta)

    Bx=B1_1D*np.sin(theta)+B2_1D*np.cos(theta)
    By=B1_1D*np.cos(theta)-B2_1D*np.sin(theta)


    Frx=Fr1_1D*np.sin(theta)+Fr2_1D*np.cos(theta)
    Fry=Fr1_1D*np.cos(theta)-Fr2_1D*np.sin(theta)

    # convert to cartesian coordinate
    xcoord=radius*np.sin(theta)
    ycoord=radius*np.cos(theta)


    width=23
    height=23
    nx=256
    ny=int(2*height*nx/width)

    xmin=0
    rmin=1
    xmax=width
    ymin=-height
    ymax=height

    xgrid=np.linspace(xmin, xmax, nx)
    ygrid=np.linspace(ymin, ymax, ny)

    xmesh,ymesh=np.meshgrid(xgrid,ygrid)

    rho_cart=griddata(np.c_[xcoord,ycoord],rho_1D,(xmesh,ymesh),method='nearest')
    sigma_cart=griddata(np.c_[xcoord,ycoord],sigma_1D,(xmesh,ymesh),method='nearest')
    kappa_cart=griddata(np.c_[xcoord,ycoord],kappa_1D,(xmesh,ymesh),method='nearest')   
    dx_cart=griddata(np.c_[xcoord,ycoord],dx_1D,(xmesh,ymesh),method='nearest')

    vx_cart=griddata(np.c_[xcoord,ycoord],vx,(xmesh,ymesh),method='nearest')
    vy_cart=griddata(np.c_[xcoord,ycoord],vy,(xmesh,ymesh),method='nearest')

    vx_cart=vx_cart/crat
    vy_cart=vy_cart/crat

    PB_cart=griddata(np.c_[xcoord,ycoord],PB_1D,(xmesh,ymesh),method='nearest')

    Bx_cart=griddata(np.c_[xcoord,ycoord],Bx,(xmesh,ymesh),method='nearest')
    By_cart=griddata(np.c_[xcoord,ycoord],By,(xmesh,ymesh),method='nearest')

    Er_cart=griddata(np.c_[xcoord,ycoord],Er_1D,(xmesh,ymesh),method='nearest')

    Frx_cart=griddata(np.c_[xcoord,ycoord],Frx,(xmesh,ymesh),method='nearest')
    Fry_cart=griddata(np.c_[xcoord,ycoord],Fry,(xmesh,ymesh),method='nearest')

#calculate the optical depth from rotation axis
    tau=np.zeros((ny,nx))
    xtaupos=np.zeros(ny)
    for j in range(ny):
       xtaupos[j]=xmin
       tau[j,0] = 0.5*dx_cart[j,0]*sigma_cart[j,0]
       if tau[j,0] < 1:
          xtaupos[j] = xgrid[0]
       for i in range(1,nx):
           tau[j,i] = tau[j,i-1] + 0.5*dx_cart[j,i]*sigma_cart[j,i]
           if tau[j,i] < 1:
              xtaupos[j] = xgrid[i]


    rmesh=(xmesh**2.0+ymesh**2.0)**0.5

    rindix=rmesh < rmin

    rho_cart[rindix]=rhomin
    vx_cart[rindix]=0.0
    vy_cart[rindix]=0.0
    Er_cart[rindix]=Ermin
    Bx_cart[rindix]=0.0
    By_cart[rindix]=0.0
    Frx_cart[rindix]=0.0
    Fry_cart[rindix]=0.0
    kappa_cart[rindix]=0.0

    vlim1=1.e-3
    vlim2=1

    label1='$\\rho/\\rho_0$'
    label2='$v$'

    title='$t='+"%4.2f"%(time/time0)+'{\ t_0}$'

    MakeRhoVSlice(rho_cart, vx_cart, vy_cart, rhomin, rhomax, vlim1,vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname1,1,taux=xtaupos,tauy=ygrid,title=title)



    label1='$E_r/a_rT_0^4$'
    label2='$B/B_0$'


    vlim1=Bmin
    vlim2=Bmax

    MakeRhoVSlice(Er_cart, Bx_cart, By_cart, Ermin,Ermax, vlim1,vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname2,1,title=title)

    if outputname3 is not None:

      label1='$E_r/a_rT_0^4$'
      label2='$F_r/ca_rT_0^4$'


      vlim1=Frmin
      vlim2=Frmax

      MakeRhoVSlice(Er_cart, Frx_cart, Fry_cart, Ermin,Ermax, vlim1,vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname3,1,title=title)
      
    if outputname4 is not None:

      label1='$\\kappa/\kappa_0$'
      label2='$v$'
      
      MakeRhoVSlice(kappa_cart, vx_cart, vy_cart, 0.1,20, 1.e-3, 1, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,outputname4,1,logscale=1,title=title)




###################################################################
# The constants

prat=0.000155307
crat=17178.9
kappaes=9.38e3
#binary frequency Omega
omega0=1.79925
gm=1.02737e5


rhomin=1.e-6
rhomax=1.e-1
Ermin=1.e-3
Ermax=1.e2
PBmin=1.e-7
PBmax=1.e2
Bmin=1.e-4
Bmax=1.0
Frmin=5.e-3
Frmax=5.e-1






comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

ni=10577
no=14660

for i in range(ni,no+1,nprocs):
  fi=i+rank
  print fi, rank

  avefile='average_'+'{:05d}'.format(fi)+'.athdf'

  outputname1='Ave2D.'+'{:05d}'.format(fi)+'_rhov.png'
  outputname2='Ave2D.'+'{:05d}'.format(fi)+'_ErB.png'
  outputname3='Ave2D.'+'{:05d}'.format(fi)+'_ErFr.png'



  PlotAveStructure(avefile,rhomin,rhomax,Ermin,Ermax,PBmin,PBmax,Bmin,Bmax,Frmin,Frmax,outputname1,outputname2,outputname3=outputname3,vecname1='Fr1',vecname2='Fr2',vecname3='vel1',vecname4='vel2')
  

