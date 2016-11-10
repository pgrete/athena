
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


# prepare function for polar map

def setup_axes(fig, rect, theta, radius):
 
    # PolarAxes.PolarTransform takes radian. However, we want our coordinate
    # system in degree
    tr = Affine2D() + PolarAxes.PolarTransform()
 
    # Find grid values appropriate for the coordinate (degree).
    # The argument is an approximate number of grids.
    grid_locator1 = angle_helper.LocatorD(2)
 
    # And also use an appropriate formatter:
    tick_formatter1 = angle_helper.FormatterDMS()
 
    # set up number of ticks for the r-axis
    grid_locator2 = MaxNLocator(4)
 
    # the extremes are passed to the function
    grid_helper = floating_axes.GridHelperCurveLinear(tr,
                                extremes=(theta[0], theta[1], radius[0], radius[1]),
                                grid_locator1=grid_locator1,
                                grid_locator2=grid_locator2,
                                tick_formatter1=tick_formatter1,
                                tick_formatter2=None,
                                )
 
    ax1 = floating_axes.FloatingSubplot(fig, rect, grid_helper=grid_helper)
    fig.add_subplot(ax1)
 
    # adjust axis
    # the axis artist lets you call axis with
    # "bottom", "top", "left", "right"
    ax1.axis["left"].set_axis_direction("bottom")
    ax1.axis["right"].set_axis_direction("top")
 
    ax1.axis["bottom"].set_visible(False)
    ax1.axis["top"].set_axis_direction("bottom")
    ax1.axis["top"].toggle(ticklabels=True, label=True)
    ax1.axis["top"].major_ticklabels.set_axis_direction("top")
    ax1.axis["top"].label.set_axis_direction("top")
 
    ax1.axis["left"].label.set_text("R")
    ax1.axis["top"].label.set_text(ur"$\theta$")
 
    # create a parasite axes
    aux_ax = ax1.get_aux_axes(tr)
 
    aux_ax.patch = ax1.patch # for aux_ax to have a clip path as in ax
    ax1.patch.zorder=0.9 # but this has a side effect that the patch is
                         # drawn twice, and possibly over some other
                         # artists. So, we decrease the zorder a bit to
                         # prevent this.
 
    return ax1, aux_ax


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



filename='disk.out1.00160.athdf'

data=yt.load(filename)


filename='rhoslice.png'
labelname='$\\rho/\\rho_0$'

logrticks=[logr[n1/6],logr[n1/3],logr[n1/2],logr[n1*2/3],logr[n1*5/6]]
rticks=["%2.1f"%x1v[n1/6],"%2.1f"%x1v[n1/3],"%2.1f"%x1v[n1/2],"%2.1f"%x1v[n1*2/3],"%2.1f"%x1v[n1*5/6]]

for i in range(len(rticks)):
   rticks[i]="$"+rticks[i]+"$"

MakeImage(rhoslice, rhomin, rhomax, xmin, xmax, ymin, ymax,logrticks,rticks,labelname,filename)



