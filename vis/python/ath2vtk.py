"""
Read .athdf data files and write new ones as single block at constant refinement level.

Note: Requires h5py.

Note: Only works for 3D data.
"""

# Python modules
import h5py
import numpy as np

#from mpi4py import MPI
import struct


# Athena++ modules
import athena_read
from sys import byteorder

# read in 2D average data
file='disk.out1.04740'

filename=file+'.athdf'

l0=1
d0=1
t0=1
v0=1
pt=1
arad = 7.5659141e-15
c = 2.9979246e+10

leveln=None

quantities=['rho','vel1','vel2','vel3','Er','pgas']

f=h5py.File(filename, 'r')

attributes = f.attrs.items()
attrs = dict(attributes)
level = f.attrs['MaxLevel']

time = f.attrs['Time']
subsample = False
if leveln is not None:
  if level > leveln:
    subsample = True
  level = leveln



print "Real Athena++ athdf file from level {:d}".format(level)

data = athena_read.athdf(filename, quantities=quantities, level=level, subsample=subsample)

# Determine new grid size
nx1 = attrs['RootGridSize'][0] * 2**level
nx2 = attrs['RootGridSize'][1] * 2**level
nx3 = attrs['RootGridSize'][2] * 2**level

f.close()

#####################################
#write the new VTK files


print "Start to write VTK file converted from athdf\n"

outfile=open(file+'.vtk', 'w')
outfile.write("# vtk DataFile Version 2.0\n")
outstr = "# PRIMITIVE vars at time= {:e}, level= {:d}, domain= {:d}\n".format(time,level,0)
outfile.write(outstr)
outfile.write("BINARY\n")
outfile.write("DATASET RECTILINEAR_GRID\n")
outfile.write("DIMENSIONS {:d} {:d} {:d}\n".format(nx1+1,nx2+1,nx3+1))

outfile.write("X_COORDINATES {:d} float\n".format(nx1+1))
array=data['x1f']
#array.tofile(outfile)
#np.savetxt(outfile,array)
myfmt='>'+'f'*len(array)

bin=struct.pack(myfmt,*array)
outfile.write(bin)



outfile.write("\nY_COORDINATES {:d} float\n".format(nx2+1))
array=data['x2f']
#array.tofile(outfile)
#np.savetxt(outfile,array)
myfmt='>'+'f'*len(array)
bin=struct.pack(myfmt,*array)
outfile.write(bin)


#outfile.write(bytearray(lst))
outfile.write("\nZ_COORDINATES {:d} float\n".format(nx3+1))
array=data['x3f']
#array.tofile(outfile)
#np.savetxt(outfile,array)
myfmt='>'+'f'*len(array)
bin=struct.pack(myfmt,*array)
outfile.write(bin)

outfile.write("\nCELL_DATA {:d}".format(nx3*nx2*nx1))

for key in quantities:
  outfile.write("\nSCALARS {:s} float\n".format(key))
  outfile.write("LOOKUP_TABLE default\n")
  array=data[key].ravel()
  myfmt='>'+'f'*len(array)
  bin=struct.pack(myfmt,*array)
  outfile.write(bin)


outfile.close()
