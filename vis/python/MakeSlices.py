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

crat=6616.36
prat=56.8151
gm=7182.03
rpeak=25.0
grav=gm/rpeak**2.0
time0=(prat/3.0)**0.5/grav
kappaes=111.28




def PlotScatter(datax,datay,color,xmin,xmax,ymin,ymax,xname,yname,colorname,filename,xlog=0,ylog=0,clog=0,title=None):
    plots, axes = plt.subplots(figsize=(10,6),dpi=300)
    plt.xlabel(xname, size = 25)
    plt.ylabel(yname, size = 25)
    plt.subplots_adjust(left=0.12,right=0.83,top=0.9,bottom=0.15)
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    if title is not None:
      plt.title(title,size=25)
    if clog > 0:
      im=axes.scatter(datax,datay,c=color,s=40,norm=LogNorm())
    else:
      im=axes.scatter(datax,datay,c=color,s=40)     
    if xlog > 0:
      axes.set_xscale('log')
    if ylog > 0:
      axes.set_yscale('log')
    cbaxes=plots.add_axes([0.85, 0.2,0.03,0.6])
    cbar=plots.colorbar(im,cax=cbaxes)
    cbar.set_label(colorname,size=20)

    axes.yaxis.set_tick_params(labelsize=20)
    axes.xaxis.set_tick_params(labelsize=20)
    
    axes.set_aspect('auto')

    
    plt.savefig(filename)
    plt.close()

def MakeRhoVSlice(data_cart, vx_cart, vy_cart, minval, maxval, vlim1, vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,time,outputname,vel=0,xlabel='$ x/r_0$',ylabel='$ y/r_0$'):
    plots, axes = plt.subplots(figsize=(10.,6.0),dpi=300)
    plt.xlabel(xlabel, size = 30)
    plt.ylabel(ylabel, size = 30)
    plt.subplots_adjust(left=0.16,right=0.78,top=0.8,bottom=0.15)
    plt.title(time,size=25,y=1.2)
    
    if vel>0:
      speed=np.sqrt(vx_cart**2+vy_cart**2)
      if vlim2 < vlim1:
         vlim2=speed.max()
      speed=np.clip(speed,vlim1,vlim2)
      logspeed=np.log10(speed)
      vcolor=axes.imshow(speed,cmap='jet',norm=LogNorm(vmin=vlim1,vmax=vlim2),origin='lower', extent=[xmin,xmax,ymin,ymax])
      velaxes = plots.add_axes([0.24, 0.89, 0.5, 0.03])
      velcbar = plots.colorbar(vcolor, cax = velaxes, orientation = 'horizontal')
      velcbar.ax.tick_params(labelsize=25)
      if speed.max() > vlim1:
        velfield=axes.streamplot(xgrid,ygrid,vx_cart,vy_cart,cmap='jet',density=1,color=logspeed)
      else:
        velfield=axes.streamplot(xgrid,ygrid,vx_cart,vy_cart,cmap='jet',density=1)
        
      axes.set_ylim([ymin,ymax])

       
    im = axes.imshow(data_cart,cmap='RdGy_r', norm=LogNorm(vmin = minval, vmax=maxval), \
                     origin='lower', extent=[xmin,xmax,ymin,ymax])
                     
    cbaxes = plots.add_axes([0.82,0.15,0.03,0.6])
    cbar=plots.colorbar(im,cax=cbaxes)
    cbar.set_label(label1, size=30)
    cbar.ax.tick_params(labelsize=25)
#    axes.set_xticks([-50,-30.,-10.,10.,30.,50.])
    axes.yaxis.set_tick_params(labelsize=20)
    axes.xaxis.set_tick_params(labelsize=20)
    
    axes.set_aspect('auto')
    plt.savefig(outputname)
    plt.close(plots)

def makerphislices(rho,Er,kappa,vel1,vel2,vel3,Fr01,Fr02,Fr03,x1v,x2v,x3v,time,outwidth=50,width=14.23):

    rhomin=1.e-5
    rhomax=5.e3
    vlim1=1.e-2
    vlim2=5.e1

    Ermin=1.e-5
    Ermax=5.e4
    Frmin=1.e-6
    Frmax=1.e-2

    kappamin=1
    kappamax=30

    #first, plot theta slice through 90
    ntheta=np.abs(x2v - np.pi*0.5).argmin()

    rhoslice=rho[:,ntheta,:]
    vrslice=vel1[:,ntheta,:]
    vphislice=vel3[:,ntheta,:]
    Erslice=Er[:,ntheta,:]
    Frslice=Fr01[:,ntheta,:]
    Fphislice=Fr03[:,ntheta,:]

    kappaslice=kappa[:,ntheta,:]


    nr=x1v.size
    nphi=x3v.size
    
    rcoord=np.zeros(nphi*nr)
    phi=np.zeros(nphi*nr)
    
    for n in range(0,nr):
      rcoord[n*nphi:(n+1)*nphi]=x1v[n]
    
    for n in range(0,nr):
      phi[n*nphi:(n+1)*nphi]=x3v[0:nphi]

    
    rhoslice=np.reshape(rhoslice,nr*nphi,order='F')
    vrslice=np.reshape(vrslice,nr*nphi,order='F')
    vphislice=np.reshape(vphislice,nr*nphi,order='F')

    Erslice=np.reshape(Erslice,nr*nphi,order='F')
    Frslice=np.reshape(Frslice,nr*nphi,order='F')
    Fphislice=np.reshape(Fphislice,nr*nphi,order='F')

    vx=vrslice*np.cos(phi)-vphislice*np.sin(phi)
    vy=vrslice*np.sin(phi)+vphislice*np.cos(phi)

    Fx=Frslice*np.cos(phi)-Fphislice*np.sin(phi)
    Fy=Frslice*np.sin(phi)+Fphislice*np.cos(phi)

    kappaslice=np.reshape(kappaslice,nr*nphi,order='F')   

# convert to cartesian coordinate
    xcoord=rcoord*np.cos(phi)
    ycoord=rcoord*np.sin(phi)

#create grid for plotting

    nx=2*nr
    ny=nr
    xmin=-outwidth
# xmax=np.max(rcoord)
    xmax=outwidth
    ymin=width/2
    ymax=outwidth
    xgrid=np.linspace(xmin, xmax, nx)
    ygrid=np.linspace(ymin, ymax, ny)

    xmesh,ymesh=np.meshgrid(xgrid,ygrid)
   
    rmesh=(xmesh**2.0+ymesh**2.0)**0.5

    rindix=rmesh < width
    rindix2=rmesh > outwidth 
    

    rho_cart=griddata(np.c_[xcoord,ycoord],rhoslice,(xmesh,ymesh),method='nearest')

    vx_cart=griddata(np.c_[xcoord,ycoord],vx,(xmesh,ymesh),method='nearest')
    vy_cart=griddata(np.c_[xcoord,ycoord],vy,(xmesh,ymesh),method='nearest')
    
    rho_cart[rindix] = rhomin
    vx_cart[rindix] = 0.0
    vy_cart[rindix] = 0.0

    rho_cart[rindix2] = rhomin
    vx_cart[rindix2] = 0.0
    vy_cart[rindix2] = 0.0

    Er_cart=griddata(np.c_[xcoord,ycoord],Erslice,(xmesh,ymesh),method='nearest')

    Fx_cart=griddata(np.c_[xcoord,ycoord],Fx,(xmesh,ymesh),method='nearest')
    Fy_cart=griddata(np.c_[xcoord,ycoord],Fy,(xmesh,ymesh),method='nearest')
    
    Er_cart[rindix] = Ermin
    Fx_cart[rindix] = 0.0
    Fy_cart[rindix] = 0.0

    Er_cart[rindix2] = Ermin
    Fx_cart[rindix2] = 0.0
    Fy_cart[rindix2] = 0.0

    kappa_cart=griddata(np.c_[xcoord,ycoord],kappaslice,(xmesh,ymesh),method='nearest')

    kappa_cart[rindix] = kappamin
    kappa_cart[rindix2] = kappamin
    
# only scale to speed of light for velocity
#        if var2=='vel1':
#           vx_cart=vx_cart/crat
#           vy_cart=vy_cart/crat

    outputname='star.'+"%4.2f"%(time)+'_rhov.pdf'
    label1='$\\rho/\\rho_0$'
    label2='$v$'

    timelabel='${\\rm time}='+"%4.2f"%(time)+'{\ t_0}$'

    MakeRhoVSlice(rho_cart, vx_cart, vy_cart, rhomin,rhomax, vlim1,vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,timelabel,outputname,vel=1)

    outputname='star.'+"%4.2f"%(time)+'_ErFr.pdf'
    label1='$E_r/a_rT_0^4$'
    label2='$F_{r,0}$'


    MakeRhoVSlice(Er_cart, Fx_cart, Fy_cart, Ermin,Ermax, Frmin,Frmax, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,timelabel,outputname,vel=1)

    outputname='star.'+"%4.2f"%(time)+'_kappav.pdf'
    label1='$\kappa/\kappa_{\\rm es}$'
    label2='$v$'


    MakeRhoVSlice(kappa_cart, vx_cart, vy_cart, kappamin,kappamax, vlim1,vlim2, width, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,timelabel,outputname,vel=1)


#################################################


def makethetaphislices(rho,Er,kappa,vel1,vel2,vel3,Fr01,Fr02,Fr03,x1v,x2v,x3v,time,rloc,rhomin=1.e-5,rhomax=1.e3,Ermin=1.e-4,Ermax=1.e1,vlim1=1.e-1,vlim2=1.e2,Frmin=1.e-6,Frmax=1.e-2,aradmin=0.01,aradmax=10):


    nr1=np.abs(x1v - rloc).argmin()

    tgas=pgas/rho

    grav1=gm/rloc**2.0
    arad1=kappaes*kappa*Fr01*prat/grav1

    rho1=rho[:,:,nr1].transpose()
    vr1=vel1[:,:,nr1].transpose()
    vt1=vel2[:,:,nr1].transpose()
    vp1=vel3[:,:,nr1].transpose()
    kappa1=kappa[:,:,nr1].transpose()
    tgas1=tgas[:,:,nr1].transpose()
    arad1=arad1[:,:,nr1].transpose()
    Er1=Er[:,:,nr1].transpose()
    Frad1=Fr01[:,:,nr1].transpose()
    Frad2=Fr02[:,:,nr1].transpose()
    Frad3=Fr03[:,:,nr1].transpose()


#create grid for plotting
    ymin=np.min(x2v)
    ymax=np.max(x2v)
    xmin=np.min(x3v)
    xmax=np.max(x3v)

    ygrid=x2v
    xgrid=x3v


    outputname='star.'+"%4.2f"%(time)+'_'+"%3.1f"%(rloc)+'_rhov_thetaphi.pdf'
    label1='$\\rho/\\rho_0$'
    label2='$v$'

    timelabel='$r='+"%4.2f"%(rloc)+'r_{\\odot},\ '+'{\\rm time}='+"%4.2f"%(time)+'{\ t_0}$'


    MakeRhoVSlice(rho1, vp1, vt1, rhomin,rhomax, vlim1,vlim2, 0.0, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,timelabel,outputname,vel=1,xlabel='$\\phi$',ylabel='$\\theta$')

    outputname='star.'+"%4.2f"%(time)+'_'+"%3.1f"%(rloc)+'_ErFr_thetaphi.pdf'
    label1='$E_r/a_rT_0^4$'
    label2='$F_{r,0}$'


    MakeRhoVSlice(Er1, Frad3, Frad2, Ermin,Ermax, Frmin,Frmax, 0.0, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,timelabel,outputname,vel=1,xlabel='$\\phi$',ylabel='$\\theta$')

    outputname='star.'+"%4.2f"%(time)+'_'+"%3.1f"%(rloc)+'_arad_thetaphi.pdf'
    label1='$a_r/a_g$'
    label2='$v$'


    MakeRhoVSlice(arad1, vp1,vt1, aradmin,aradmax, vlim1,vlim2, 0.0, xmin, xmax, ymin, ymax,xgrid,ygrid,label1,label2,timelabel,outputname,vel=0,xlabel='$\\phi$',ylabel='$\\theta$')


#################################################

def makephaseplots(rho,pgas,kappa,vel1,Fr01,x1v,x2v,x3v,time,rloc,rhomin=1.e-5,rhomax=30.0,vlim1=1.e-2,vlim2=5.e1,Ermin=1.e-5,Ermax=5.e4,Frmin=1.e-6,Frmax=1.e-2,kappamin=1,kappamax=120,aradmin=0.05,aradmax=3,tmin=0.1,tmax=10.0):


    #first, plot theta slice through 90
    nr1=np.abs(x1v - rloc).argmin()

    tgas=pgas/rho

    grav1=gm/rloc**2.0
    arad1=kappaes*kappa*Fr01*prat/grav1

    rho1=rho[:,:,nr1]
    vr1=vel1[:,:,nr1]
    kappa1=kappa[:,:,nr1]
    tgas1=tgas[:,:,nr1]
    arad1=arad1[:,:,nr1]


    nr=x1v.size
    ntheta=x2v.size
    nphi=x3v.size

    rho1=np.reshape(rho1,ntheta*nphi,order='F')
    vr1=np.reshape(vr1,ntheta*nphi,order='F')
    kappa1=np.reshape(kappa1,ntheta*nphi,order='F')
    tgas1=np.reshape(tgas1,ntheta*nphi,order='F')

    
# only scale to speed of light for velocity
#        if var2=='vel1':
#           vx_cart=vx_cart/crat
#           vy_cart=vy_cart/crat

    outputname='star.'+"%4.2f"%(time)+'_'+"%3.1f"%(rloc)+'_kappa_T_rho.png'
 

    title='$r='+"%4.2f"%(rloc)+'r_{\\odot},\ '+'{\\rm time}='+"%4.2f"%(time)+'{\ t_0}$'


    PlotScatter(kappa1,tgas1,rho1,kappamin,kappamax,tmin,tmax,'$\\kappa/\\kappa_{\\rm es}$','$T_g/T_0$','$\\rho/\\rho_0$',outputname,xlog=1,ylog=1,clog=1,title=title)

    outputname='star.'+"%4.2f"%(time)+'_'+"%3.1f"%(rloc)+'_arad_rho_v.png'
 

    PlotScatter(arad1,rho1,vr1,aradmin,aradmax,rhomin,rhomax,'$a_r/a_g$','$\\rho/\\rho_0$','$v_r$',outputname,xlog=1,ylog=1,title=title)


 #   PlotScatter(kappa2,tgas2,rho2,kappamin, kappamax,0.2,2.0,'$\\kappa/\\kappa_{\\rm es}$','$T_g/T_0$','$\\rho/\\rho_0$',outputname,xlog=1,ylog=1,clog=1,title=title)

 

 #   PlotScatter(arad2,rho2,vr2,0.01, 10.0,1.e-4,100,'$a_r/a_g$','$\\rho/\\rho_0$','$v_r$',outputname,xlog=1,ylog=1,title=title)

#################################################





fi=1138

filename='vis.'+'{:05d}'.format(fi)+'.athdf'

f=h5py.File(filename, 'r')
#quantities=['rho','Er','Sigma_s','Sigma_a','pgas','vel1','vel2','vel3','Fr01','Fr02','Fr03']
x1f=f['x1f'].value
x1f=x1f[0]
x2f=f['x2f'].value
x2f=x2f[0]
x3f=f['x3f'].value
x3f=x3f[0]

time=f.attrs['Time']
time=time/time0

rho=f['quantities'][0][0]
Er=f['quantities'][1][0]
Sigma_s=f['quantities'][2][0]
Sigma_a=f['quantities'][3][0]
pgas=f['quantities'][4][0]
vel1=f['quantities'][5][0]
vel2=f['quantities'][6][0]
vel3=f['quantities'][7][0]
Fr01=f['quantities'][8][0]
Fr02=f['quantities'][9][0]
Fr03=f['quantities'][10][0]

kappa=(Sigma_s+Sigma_a)/(rho*kappaes)

size=rho.shape

nr=size[0]
ntheta=size[1]
nphi=size[2]

#get cell centered positions
x1v=np.zeros(nr)
x2v=np.zeros(ntheta)
x3v=np.zeros(nphi)

for ni in range(nr):
   x1v[ni]=0.75*(x1f[ni+1]**4.0 - x1f[ni]**4.0)/(x1f[ni+1]**3.0-x1f[ni]**3.0)

for nj in range(ntheta):
  x2v[nj]=((np.sin(x2f[nj+1]) - x2f[nj+1] * np.cos(x2f[nj+1])) \
      -(np.sin(x2f[nj]) - x2f[nj] * np.cos(x2f[nj]))) \
        / (np.cos(x2f[nj]) - np.cos(x2f[nj+1]))

for nk in range(nphi):
  x3v[nk] = 0.5 * (x3f[nk+1]+x3f[nk])


#plot slice images

# // makerphislices(rho,Er,kappa,vel1,vel2,vel3,Fr01,Fr02,Fr03,x1v,x2v,x3v,time)

#make phase plots between kappa, rho, T, vr

#// makephaseplots(rho,pgas,kappa,vel1,x1v,x2v,x3v,time)

# // makephaseplots(rho,pgas,kappa,vel1,Fr01,x1v,x2v,x3v,time)

makephaseplots(rho,pgas,kappa,vel1,Fr01,x1v,x2v,x3v,time,31.5,tmin=0.1,tmax=1.0)

makephaseplots(rho,pgas,kappa,vel1,Fr01,x1v,x2v,x3v,time,28.0,tmin=0.2,tmax=2.0,aradmin=0.01,aradmax=10.0,rhomin=1.e-4,rhomax=1.e2)

makethetaphislices(rho,Er,kappa,vel1,vel2,vel3,Fr01,Fr02,Fr03,x1v,x2v,x3v,time,31.5)

makethetaphislices(rho,Er,kappa,vel1,vel2,vel3,Fr01,Fr02,Fr03,x1v,x2v,x3v,time,28,rhomin=1.e-3,rhomax=5.e3,Ermin=1.e-4,Ermax=10.,vlim1=1.e-1,vlim2=1.e2,Frmin=1.e-6,Frmax=1.e-2,aradmin=0.01,aradmax=10)

