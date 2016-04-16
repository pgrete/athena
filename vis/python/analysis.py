
from athena_yt import *

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
import yt

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
    plt.subplots_adjust(left=0.14,right=0.95,top=0.9,bottom=0.15)
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
def MakeImage(data, minval, maxval, xmin, xmax, ymin, ymax,labelname,filename):
    plots, axes = plt.subplots(figsize=(12,9),dpi=300)
    im = axes.imshow(data,cmap='RdGy_r', norm=LogNorm(vmin = minval, vmax=maxval), origin='lower', extent=[xmin,xmax,ymin,ymax])
    cbaxes = plots.add_axes([0.85,0.1,0.03,0.8])
    cbar=plots.colorbar(im,cax=cbaxes)
    cbar.set_label(labelname, size=20)
    axes.set_aspect('auto')
    plt.savefig(filename)
    plt.close(plots)


def loadhdf(filename):
    f=h5py.File(filename,'r')
    # to see all the available keys f.attrs.keys()
    nblock=f.attrs['TotalMeshBlock']
    block_size=f.attrs['MeshBlockSize']
    root_grid_size=f.attrs['RootGridSize']
    maxlevel=f.attrs['MaxLevel']
    cycle=f.attrs['NCycle']
    time=f.attrs['Time']
    nvariable=f.attrs['NVariables']
    #for quantities in each block
    Er=np.array(f[u'MeshBlock0'][u'Er'])

vol = lambda x1m,x1p,x2m,x2p,x3m,x3p: 1.0/3.0 * (x1p**3-x1m**3) * abs(np.cos(x2m)-np.cos(x2p)) * (x3p-x3m)



filename='PolarCap.out1.00157.athdf'

data=yt_loadathdf(filename)

# for volume Rendering
mi, ma=data.all_data().quantities.extrema('rho')

mi=log10(mi.value)
ma=log10(ma.value)

# create a transfer function

tf=yt.ColorTransferFunction((mi,ma))
tf.add_layers(6,w=0.01)

# define the properties and size of the Camera viewpoint

L=[0.3, 0.3, 0.2]  # the view point
c=[0., 0.0, 5.0]  # the focual point
w=1.5*data.domain_width[0] # width
Npixels=512

cam=data.camera(c,L,w,Npixels,tf,fields=['rho'], north_vector=[0,0,1],\
              steady_north=True, sub_samples=5, log_fields=[False])

cam.transfer_function.map_to_colormap(mi,ma,scale=1.0,colormap='algae')
image=cam.snapshot(fn='test.png')
cam.save_image(image,transparent=True)
