import pyathena_simple as pa
import numpy as np

file_original = "/Users/munangong/chemistry_Athena/athena/ChangGoo_data/MHD_32pc/id0/MHD_32pc.0000.vtk"
file_new = "/Users/munangong/athena/outputs/cgk_32pc_f8_mpi2.out1.00000.vtk"

tol = 1e-20

original = pa.AthenaDataSet(file_original, noread_mpi=False)
new = pa.AthenaDataSet(file_new)
od = original.read_all_data('density')
nd = new.read_all_data('rho')
ov = original.read_all_data('velocity')
nv = new.read_all_data('vel')
op = original.read_all_data('pressure')
np = new.read_all_data('press')


if od.shape != nd.shape:
    raise RuntimeError( "ERROR: original shape = {}, new shape = {}".format(od.shape,
            nd.shape) )
    
md = ( abs(od-nd) ).max()
mv = ( abs(ov-nv) ).max()
mp =  ( abs(op-np) ).min()

print "Athena4.2 file: {}".format(file_original)
print "Athena++ file: {}".format(file_new)
print "Max difference in density = {:10e}".format( md )
print "Max difference in velocity = {:10e}".format( mv )
print "Max difference in pressure = {:10e}".format( mp )

if max(md, mv, mp) < tol:
    print "Pass!"
else:
    print "Fail! Max error = {:.2e}".format( max(md, mv, mp) )
