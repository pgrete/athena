import h5py
import numpy as np


from mpi4py import MPI



def DumpOneBlock(i,bc,x1f,x2f,x3f,data,levels,outf,xdmf,outputname):
    if bc[i].size >1:
       ob=bc[i][0]
       xb=bc[i][1]
       yb=bc[i][2]
       zb=bc[i][4]
       x = np.concatenate([x1f[ob,:-1], x1f[xb,:]])
       y = np.concatenate([x2f[ob,:-1], x2f[yb,:]])
       z = np.concatenate([x3f[ob,:-1], x3f[zb,:]])
    else:
       ob = int(bc[i])
       x = x1f[bc[i],:]
       y = x2f[bc[i],:]
       z = x3f[bc[i],:]
    name='MeshBlock'+str(i)
    outf.create_group(name)
    newnx=len(x)-1
    newny=len(y)-1
    newnz=len(z)-1
    outf[name].attrs.create('Level',levels[ob],dtype='>i4')
    outf[name].attrs.create('BlockSize',[newnx,newny,newnz],dtype='>i4')
    outf[name].create_dataset('x1f',data=x,dtype='>f4',shape=x.shape)
    outf[name].create_dataset('x2f',data=y,dtype='>f4',shape=y.shape)
    outf[name].create_dataset('x3f',data=z,dtype='>f4',shape=z.shape)
      
    xdmf.write('  <Grid Name="{0}" GridType="Uniform">\n'.format(name))
    xdmf.write('    <Time Type="Single" Value="{0}"/>\n'.format(data.attrs['Time']))
    xdmf.write('    <Topology TopologyType="3DRectMesh"' \
                      + ' NumberOfElements="{0} {1} {2}"/>\n'.format(newnz+1,newny+1,newnx+1))
    xdmf.write( '      <Geometry GeometryType="VXVYVZ">\n')
    xdmf.write(('      <DataItem Dimensions="{0}" NumberType="Float" Precision="4"' + \
                      ' Format="HDF">{1}:/{2}/x1f</DataItem>\n').format(\
                        newnx+1,outputname,name))
    xdmf.write(('      <DataItem Dimensions="{0}" NumberType="Float" Precision="4"'+ \
                      ' Format="HDF">{1}:/{2}/x2f</DataItem>\n').format(\
                        newny+1,outputname,name))
    xdmf.write(('      <DataItem Dimensions="{0}" NumberType="Float" Precision="4"'+ \
                      ' Format="HDF">{1}:/{2}/x3f</DataItem>\n').format(\
                        newnz+1,outputname,name))
    xdmf.write('    </Geometry>\n')
    var_offset=0
    for dataset_name,num_vars in zip(data.attrs['DatasetNames'],data.attrs['NumVariables']):
       for var_num in range(num_vars):
          var_name=data.attrs['VariableNames'][var_offset+var_num]
          if var_name in quantities:
    #data array shape follows c convection
             if bc[i].size >1:
                dataxx_y1=np.dstack((data[dataset_name][var_num][bc[i][0]],\
                           data[dataset_name][var_num][bc[i][1]]))
                dataxx_y2=np.dstack((data[dataset_name][var_num][bc[i][2]],\
                           data[dataset_name][var_num][bc[i][3]]))
                data_z1=np.hstack((dataxx_y1,dataxx_y2))
        
                dataxx_y1=np.dstack((data[dataset_name][var_num][bc[i][4]],\
                           data[dataset_name][var_num][bc[i][5]]))
                dataxx_y2=np.dstack((data[dataset_name][var_num][bc[i][6]],\
                           data[dataset_name][var_num][bc[i][7]]))
                data_z2=np.hstack((dataxx_y1,dataxx_y2))
        
                newblockdata=np.vstack((data_z1,data_z2))
             else:
                newblockdata=data[dataset_name][var_num][bc[i]]
        
             outf[name].create_dataset(var_name,data=newblockdata,dtype='>f4',shape=(newnz,newny,newnx))
    # xdmf file information


             xdmf.write(('    <Attribute Name="{0}" AttributeType="Scalar" Center="Cell">\n').format(var_name))
             xdmf.write(('      <DataItem Dimensions="{0} {1} {2}" NumberType="Float" Precision="4"' + \
                      ' Format="HDF">{3}:/{4}/{5}</DataItem>\n').format(\
                             newnz,newny,newnx,outputname,name,var_name))
             xdmf.write('    </Attribute>\n')

       var_offset += num_vars

    xdmf.write('  </Grid>\n')
        



def ConvertVisData(file_num,quantities):

  inputname='disk.out1.'+'{:05d}'.format(file_num)+'.athdf'
  outputname='vis.'+'{:05d}'.format(file_num)+'.athdf'
  data=h5py.File(inputname,'r')
  attributes = data.attrs.items()
  attrs = dict(attributes)
  outf=h5py.File(outputname, 'w')
# now dump the xdmf file

  xdmf=open(outputname+'.xdmf', 'w')
  xdmf.write('<?xml version="1.0" ?>\n')
  xdmf.write('<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n')
  xdmf.write('<Xdmf Version="2.0">\n')
  xdmf.write('<Domain>\n')
  xdmf.write('<Grid Name="Mesh" GridType="Collection">\n')


  num_blocks=data.attrs['NumMeshBlocks']
  block_size=data.attrs['MeshBlockSize']
#3D location of each block
  log_loc=data['LogicalLocations']
# level for all blocks
  levels=data['Levels']
  x1f=data['x1f']
  x2f=data['x2f']
  x3f=data['x3f']
  nx=block_size[0]
  ny=block_size[1]
  nz=block_size[2]
  nbx,nby,nbz=tuple(np.max(log_loc,axis=0)+1)
  nlevel=data.attrs['MaxLevel']+1
  block_grid=-np.ones((nbx,nby,nbz,nlevel),dtype=np.int)
  block_grid[log_loc[:,0],log_loc[:,1],log_loc[:,2],levels[:]] = np.arange(num_blocks)


  bc=[]
  block_list=np.arange(num_blocks,dtype='int64')
  for i in range(num_blocks):
    if block_list[i] > -1:
       ii, jj, kk=log_loc[i]
       neigh = block_grid[ii:ii+2,jj:jj+2,kk:kk+2,levels[i]]
       if np.all(neigh > -1):
          loc_ids = neigh.transpose().flatten()
          bc.append(loc_ids)
          block_list[loc_ids] = -1
       else:
          bc.append(np.array(i))
          block_list[i] = -1
  num_meshes = len(bc)


  NVariables=sum(data.attrs['NumVariables'])
  for key,val in attributes:
      if key == 'NumMeshBlocks':
            value = num_meshes
      else:
            value = val
      outf.attrs.create(key, value, dtype=val.dtype)
  outf.attrs.create('NVariables',NVariables,dtype='>i4')

  for n in range(num_meshes):
    DumpOneBlock(n,bc,x1f,x2f,x3f,data,levels,outf,xdmf,outputname)


  outf.close

  xdmf.write('</Grid>\n')
  xdmf.write('</Domain>\n')
  xdmf.write('</Xdmf>\n')
  xdmf.close


quantities=['rho','vel1','vel2','B1','B2','B3','Er','Sigma_s','Sigma_a']

ni=4576
no=4866

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

for i in range(ni,no+1,nprocs*4):
  file_num=i+rank*4
  print file_num, rank
  ConvertVisData(file_num,quantities)


