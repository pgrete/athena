from athena_read import *

import sys
sys.settrace
import matplotlib
#matplotlib.use('Agg')
import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from pylab import *
import struct
import array
import os

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

filename='Radshock.out4.00007.athdf'

data=athdf(filename)

Er =data[u'Er']
Fr1=data[u'Fr1']
Fr2=data[u'Fr2']
Fr3=data[u'Fr3']

Er0=data[u'Er0']
Fr01=data[u'Fr01']
Fr02=data[u'Fr02']
Fr03=data[u'Fr03']

Pr11=data[u'Pr11']
Pr12=data[u'Pr12']
Pr13=data[u'Pr13']
Pr21=data[u'Pr21']
Pr22=data[u'Pr22']
Pr23=data[u'Pr23']
Pr31=data[u'Pr31']
Pr32=data[u'Pr32']
Pr33=data[u'Pr33']

sigma_a=data[u'Sigma_a']

sigma_s=data[u'Sigma_s']

rho=data[u'rho']

press=data[u'press']

v1=data[u'vel1']

v2=data[u'vel2']

v3=data[u'vel3']


dim=rho.shape

n3=dim[0]
n2=dim[1]
n1=dim[2]

x1f=data[u'x1f']

if n2 > 1:
  x2f=data[u'x2f']

if n3 > 1:
  x3f=data[u'x3f']

x1v=np.zeros(n1)
for i in range(0,n1):
  x1v[i]=0.5*(x1f[i]+x1f[i+1])


#################################################


crat=1.732e3

Frini=np.loadtxt('rF_4096.dat')
posxini=np.loadtxt('x_4096.dat')

Frini=Frini/crat

Frend=Fr1[0][3][:]

outname='radFr.pdf'

PloTwoLines(x1v,posxini,Frend,Frini, min(x1v),max(x1v),'${\\bf x}$',-0.2, 0.0015, '${\\bf F_r}$', outname, 0)
