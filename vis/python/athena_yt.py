#!/usr/bin/python

import numpy as np

import h5py
import yt



def yt_loadathdf(filename, quantities=None, Prat=1, Crat=1,geometry='cartesian',flag='yt'):
    f=h5py.File(filename,'r')
    # to see all the available keys f.attrs.keys()
    nblock=f.attrs['TotalMeshBlock']
    block_size=f.attrs['MeshBlockSize']
    root_grid_size=f.attrs['RootGridSize']
    maxlevel=f.attrs['MaxLevel']
    cycle=f.attrs['NCycle']

    time=f.attrs['Time']

#    nvariable=f.attrs['NVariables']

    nx=block_size[0]
    ny=block_size[1]
    nz=block_size[2]

    x1f=np.array(f[u'MeshBlock0'][u'x1f'])
    x2f=np.array(f[u'MeshBlock0'][u'x2f'])
    x3f=np.array(f[u'MeshBlock0'][u'x3f'])

#    location=f[u'MeshBlock0'].attrs[u'LogicalLocation']
#    id=f[u'MeshBlock0'].attrs[u'GlobalID']

    if quantities is None:
       quantities=f[u'MeshBlock0'].keys()
    quantities=[str(q) for q in quantities \
                if q != 'x1f' and q != 'x2f' and q != 'x3f']

    grid_data = [dict() for x in range(nblock)]

    blockid=0
    root_lx1=x1f[0]
    root_rx1=x1f[nx]
    root_lx2=x2f[0]
    root_rx2=x2f[ny]
    root_lx3=x3f[0]
    root_rx3=x3f[nz]

    block_size=block_size

    for block in f.itervalues():
        level=block.attrs['Level'][0]
        x1f=np.array(block[u'x1f'])
        x2f=np.array(block[u'x2f'])
        x3f=np.array(block[u'x3f'])
        root_lx1=min(x1f[0],root_lx1)
        root_rx1=max(x1f[nx],root_rx1)
        root_lx2=min(x2f[0],root_lx2)
        root_rx2=max(x2f[ny],root_rx2)
        root_lx3=min(x3f[0],root_lx3)
        root_rx3=max(x3f[nz],root_rx3)
        left=[x1f[0],x2f[0],x3f[0]]
        right=[x1f[nx],x2f[ny],x3f[nz]]
        grid_data[blockid]['left_edge']=left
        grid_data[blockid]['right_edge']=right
        grid_data[blockid]['level']=level
        # the block size is (nx, ny, nz)
        # we need (nz, ny , nx) corresponding to 3D array

        grid_data[blockid]['dimensions']=block_size
        for q in quantities:
            grid_data[blockid][q]=np.reshape(np.ravel(np.array(block[q]),order='C'),(nx,ny,nz),order='F')
            if q == 'Er' or q=='Er0' or q=='Pr11' or q=='Pr12' \
               or q=='Pr13' or q=='Pr21' or q=='Pr22' or \
               q=='Pr23' or q=='Pr31' or q=='Pr32' or q=='Pr33':
               grid_data[blockid][q] = grid_data[blockid][q] *Prat
            elif q=='Fr01' or q=='Fr02' or q=='Fr03' or q=='Fr1' \
               or q=='Fr2' or q=='Fr3':
               grid_data[blockid][q] = grid_data[blockid][q] *Prat/Crat
        print blockid,nblock
        blockid=blockid+1

    # close the file
    f.close()


    # the unit
    field_units=dict()
    for q in quantities:
      if q=='Er' or q=='Er0' or q=='Pr11' or q=='Pr12' \
         or q=='Pr13' or q=='Pr21' or q=='Pr22' or \
            q=='Pr23' or q=='Pr31' or q=='Pr32' or q=='Pr33' \
            or q=='press':
        field_units[q]='code_mass/(code_time**2*code_length)'
      elif q=='Fr01' or q=='Fr02' or q=='Fr03' or q=='Fr1' \
          or q=='Fr2' or q=='Fr3':
        field_units[q]='code_mass*code_length/code_time'
      elif q=='vel1' or q=='vel2' or q=='vel3':
        field_units[q]='code_length/code_time'
      elif q=='rho':
        field_units[q]='code_mass/code_length**3'

      bbox=np.array([[root_lx1, root_rx1],[root_lx2, root_rx2],[root_lx3, root_rx3]])
      if flag=='yt':
          ds=yt.load_amr_grids(grid_data,block_size,field_units=field_units,sim_time=time,bbox=bbox,geometry=geometry)
          return ds
      else:
          return grid_data


