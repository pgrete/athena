
from athena_read import *

import sys
sys.settrace
import matplotlib
#matplotlib.use('Agg')
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pylab import *
import struct
import array
import os

import h5py


# setup latex fonts
rcParams['text.usetex']=True
rcParams['text.latex.unicode']=True


matplotlib.rc('font', family='serif', serif='cm10')
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

#Function to make line plot

def PloTwoLines(xpos1, xpos2,var1, var2, xmin,xmax,xname,ymin, ymax, yname, outname, logflag):
    plots, axes = plt.subplots(figsize=(12,9),dpi=300)
    plt.xlabel(xname, size = 20)
    plt.ylabel(yname, size = 20)
    plt.subplots_adjust(left=0.1,right=0.95,top=0.9,bottom=0.15)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlim([xmin,xmax])
    plt.ylim([ymin,ymax])
    if logflag > 0:
      axes.set_yscale('log')
      axes.plot(xpos1,math.log10(var1),linewidth=2.0, label='var1', linestyle='solid',color='black')
      axes.plot(xpos2,math.log10(var2),linewidth=2.0, label='var2', linestyle='dashed',color='red')
    else:
      axes.plot(xpos1,var1,linewidth=2.0, label='var1', linestyle='solid',color='black')
      axes.plot(xpos2,var2,linewidth=2.0, label='var2', linestyle='dashed',color='red')

    axes.set_aspect('auto')
    axes.legend(frameon=False)
    axes.spines["top"].set_visible(False)
    axes.spines["right"].set_visible(False)
    axes.get_xaxis().tick_bottom()
    axes.get_yaxis().tick_left()


    plt.savefig(outname,bbox_inches="tight")
    plt.close(plots)


#################################################

# function to make 2D image for scalar
def MakeImage(data, minval, maxval, xmin, xmax, ymin, ymax,logrticks,rticks,labelname,filename):
    plots, axes = plt.subplots(figsize=(12,9),dpi=300)
    plt.xlabel('$ \\log(r/r_s)$', size = 20)
    plt.ylabel('$ \\theta$', size = 20)
    plt.subplots_adjust(left=0.1,right=0.8,top=0.9,bottom=0.1)
    plt.xticks(logrticks,rticks)
    
    im = axes.imshow(data,cmap='RdGy_r', norm=LogNorm(vmin = minval, vmax=maxval), \
                     origin='lower', extent=[xmin,xmax,ymin,ymax])
    cbaxes = plots.add_axes([0.85,0.15,0.03,0.6])
    cbar=plots.colorbar(im,cax=cbaxes)
    cbar.set_label(labelname, size=20)
    cbar.ax.tick_params(labelsize=15)
    
    axes.yaxis.set_tick_params(labelsize=20)
    axes.xaxis.set_tick_params(labelsize=20)
    axes.set_aspect('auto')
    plt.savefig(filename)
    plt.close(plots)




vol = lambda x1m,x1p,x2m,x2p,x3m,x3p: 1.0/3.0 * (x1p**3-x1m**3) * abs(np.cos(x2m)-np.cos(x2p)) * (x3p-x3m)



filename='disk.out1.00155.athdf'


# first use python hdf5 reader to get basic information
df=h5py.File(filename,'r')
time=df.attrs[u'Time']
maxlevel=df.attrs[u'MaxLevel']

df.close()

data=athdf(filename,quantities=['rho'],level=maxlevel,vol_func=vol)

x1f=data['x1f']
x2f=data['x2f']
x3f=data['x3f']

rho=data['rho']

n1=x1f.size-1
n2=x2f.size-1
n3=x3f.size-1


x1v=np.zeros(n1)
x2v=np.zeros(n2)
x3v=np.zeros(n3)

# calculate cell centered coordinate in spherical polar
for i in range(0,n1):
  x1v[i] = 0.75*(np.power(x1f[i+1],4.0) - np.power(x1f[i],4.0))/(np.power(x1f[i+1],3.0) - np.power(x1f[i],3.0))

for i in range(0,n2):
  x2v[i] = 0.5*(x2f[i+1] + x2f[i])

for i in range(0,n3):
  x3v[i] = 0.5*(x3f[i+1] + x3f[i])

logr=np.log10(x1v)
#take slice

rhoslice=rho[n3/2,:,:]

rhomax=np.max(rhoslice)
rhomin=np.min(rhoslice)

xmin=np.min(logr)
xmax=np.max(logr)

ymin=np.min(x2v)
ymax=np.max(x2v)

filename='rhoslice.png'
labelname='$\\rho/\\rho_0$'

logrticks=[logr[n1/6],logr[n1/3],logr[n1/2],logr[n1*2/3],logr[n1*5/6]]
rticks=["%2.1f"%x1v[n1/6],"%2.1f"%x1v[n1/3],"%2.1f"%x1v[n1/2],"%2.1f"%x1v[n1*2/3],"%2.1f"%x1v[n1*5/6]]

for i in range(len(rticks)):
   rticks[i]="$"+rticks[i]+"$"

MakeImage(rhoslice, rhomin, rhomax, xmin, xmax, ymin, ymax,logrticks,rticks,labelname,filename)



