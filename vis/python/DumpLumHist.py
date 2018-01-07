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


vol_func = lambda rm,rp,thetam,thetap: \
           (1.0/3.0)*(rp**3-rm**3) * abs(np.cos(thetam)-np.cos(thetap)) * 2.0*np.pi

rarea_func = lambda r,thetam,thetap: \
           r**2.0 * abs(np.cos(thetam)-np.cos(thetap)) * 2.0*np.pi

#calculate Mdot at 5, 10, 15, 20 r_s

files=sorted(glob.glob('LumHistory/Lum*txt'))

num_file=len(files)

#plot the luminosity histoyr



#time, (Mdot, Mdotin, Mdotout)*3, meanrho, meanEr,
ncol=12
Lumhist=np.zeros((num_file,ncol))
count=0


for filename in files:
  print filename
  data=np.loadtxt(filename)
  dim=np.shape(data)


  for n in range(ncol):
    #sum the flux from top and bottom
    Lumhist[count,n]=data[4197,n]+data[3994,n]
  count=count+1

np.savetxt('Lum_history_40r_g.txt',Lumhist,fmt='%4.5e')
