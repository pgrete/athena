import yt
from joblib import Parallel, delayed

#yt.enable_parallelism()

crat=5694.76
rhomin=0.1
rhomax=1.e3

def makeslice(ni,no,width=60):
  for i in range(ni,no+1):
    print i
    filename='disk.out1.'+'{:05d}'.format(i)+'.athdf'
    data=yt.load(filename)
    xcenter=data.domain_center[0].d
    left=data.domain_left_edge[0].d
    shift=xcenter-width*0.5
    slc = yt.SlicePlot(data, 'phi', 'Er')
    slc.pan([-shift,0])
    slc.set_width(width,100)
    slc.set_figure_size(4)
    slc.set_cmap(field="Er", cmap='RdGy_r')

    slc.set_colorbar_label("Er","$E_r/a_rT_0^4$")
    slc.set_zlim(field='Er', zmin=rhomin, zmax=rhomax)
    slc.set_xlabel('$r/r_s$')
    slc.set_ylabel('$z/r_s$')
    slc.annotate_text((0.35, 0.96),'time='+"%4.2f"%(data.current_time*crat)+' $r_s/c$', coord_system='figure',text_args={'color':'black'})
    slc.set_font({'family':'sans-serif', 'style':'italic',
                  'weight':'bold', 'size':24, 'color':'black'})
#slc.annotate_streamlines('velocity_logspherical_logr','velocity_logspherical_theta')
#slc.annotate_quiver('velocity_logspherical_logr', 'velocity_logspherical_theta',40)
    slc.set_font_size(25)
    outputname='disk.'+'{:05d}'.format(i)
    slc.save(outputname)


Parallel(n_jobs=12)(delayed(makeslice)(i,i) for i in range(0,240))
