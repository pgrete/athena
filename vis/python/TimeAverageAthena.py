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

from ReduceAthenaData import *

from mpi4py import MPI


# setup latex fonts
#rcParams['text.usetex']=True
#rcParams['text.latex.unicode']=True


#matplotlib.rc('font', family='serif', serif='cm10')
#matplotlib.rc('text', usetex=True)
#matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']


crat=5694.76
vol_func = lambda rm,rp,thetam,thetap,phim,phip: \
           (1.0/3.0)*(rp**3-rm**3) * abs(np.cos(thetam)-np.cos(thetap)) * (phip-phim)

ni=4865
no=4866

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

jobload=(no-ni+1)/nprocs
residual=(no-ni+1)%nprocs

if rank == 0:
  fi=ni
  fo=fi+jobload+residual
else:
  fi=ni+jobload*rank+residual
  fo=fi+jobload

for i in range(fi,fo):
  print i
  filename='disk.out1.'+'{:05d}'.format(i)+'.athdf'
  data=ReduceData(filename)
  quantities=data.keys()
  dim=data[quantities[0]].shape
  if i==fi:
     newdata=data
  else:
     for q in quantities:
       if q != 'x1f' and q != 'x2f' and q != 'x3f' and q != 'x1v' and q != 'x2v' and q != 'x3v':
          newdata[q] = newdata[q] + data[q]


# average over all MPI process

ntot=no-ni+1

recvbuf=np.zeros(dim)

meandata=data

for q in quantities:
    if q != 'x1f' and q != 'x2f' and q != 'x3f' and q != 'x1v' and q != 'x2v' and q != 'x3v':
        comm.Reduce(newdata[q], recvbuf, op=MPI.SUM, root=0)
        if rank==0:
          meandata[q] = recvbuf/ntot

if rank ==0:
  outputname='time_average_'+'{:05d}'.format(ni)+'_'+'{:05d}'.format(no)+'.athdf'

  outf=h5py.File(outputname, 'w')

  for q in quantities:
    outf.create_dataset(q,data=meandata[q],dtype='>f8',shape=meandata[q].shape)

  outf.close()
