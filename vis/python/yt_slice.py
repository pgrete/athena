import yt
from mpi4py import MPI

#yt.enable_parallelism()

crat=5694.76

rhomin=1.e-4
rhomax=50.0

Ermin=0.1
Ermax=2.e4

def makeslice(file_num,variable,label,lim1,lim2,width=60):
    filename='disk.out1.'+'{:05d}'.format(file_num)+'.athdf'
    data=yt.load(filename)
    xcenter=data.domain_center[0].d
    left=data.domain_left_edge[0].d
    shift=xcenter-width*0.5
    slc = yt.SlicePlot(data, 'phi', variable)
    slc.pan([-shift,0])
    slc.set_width(width,180)
    slc.set_figure_size(4.5)
    slc.set_cmap(field=variable, cmap='RdGy_r')

    slc.set_colorbar_label(variable,label)
    slc.set_zlim(field=variable, zmin=lim1, zmax=lim2)
    slc.set_xlabel('$r/r_s$')
    slc.set_ylabel('$z/r_s$')
    slc.annotate_text((0.35, 0.96),'time='+"%4.2f"%(data.current_time*crat)+' $r_s/c$', coord_system='figure',text_args={'color':'black'})
    slc.set_font({'family':'sans-serif', 'style':'italic',
                  'weight':'bold', 'size':24, 'color':'black'})
#slc.annotate_streamlines('velocity_logspherical_logr','velocity_logspherical_theta')
#slc.annotate_quiver('velocity_logspherical_logr', 'velocity_logspherical_theta',40)
    slc.set_font_size(25)
    outputname='disk.'+'{:05d}'.format(file_num)
    slc.save(outputname)

ni=1442
no=4855

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nprocs = comm.Get_size()

for i in range(ni,no+1,nprocs):
  file_num=i+rank
  print file_num, rank
  makeslice(file_num,'density','$\\rho/\\rho_0$',rhomin,rhomax,width=140)
  makeslice(file_num,'Er','$E_r/a_rT_0^4$',Ermin,Ermax,width=140)

