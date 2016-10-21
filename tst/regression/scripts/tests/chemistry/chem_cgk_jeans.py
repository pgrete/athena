#regression test for equilibrium chemistry and temperature in cgk simulation
#with jeans shielding 
#compare to know solution.

# Modules
import numpy as np                             # standard Python module for numerics
import sys                                     # standard Python module to change path
import scripts.utils.athena as athena          # utilities for running Athena++
import scripts.utils.comparison as comparison  # more utilities explicitly for testing
sys.path.insert(0, '../../vis/python')         # insert path to Python read scripts
import athena_read                             # utilities for reading Athena++ data

def prepare():
  athena.configure('radiation', 'pp', 
      prob='read_vtk',
      chemistry='gow16', 
      cvode_path='/usr/local',
      cxx="g++"
      #cvode_path='/home/munan/install',
      #cxx="icc"
      )
  athena.make()

def run():
  arguments = [ 
          'problem/G0=1.0',
          'time/tlim=3e14',
          'time/nlim=3',
          'radiation/integrator=jeans',
          'chemistry/isH2RVcooling=0',
          'problem/vtkfile=../data/chem_cgk_input.vtk',
          'mesh/nx1=16',
          'mesh/nx2=16',
          'mesh/nx3=32'
          ]
  athena.run('chemistry/athinput.cgk_32pc_chem', arguments)

def analyze():
  err_control = 1e-6
  _,_,_,data_ref = athena_read.vtk('data/chem_cgk_output.vtk')
  _,_,_,data_new = athena_read.vtk('bin/cgk_32pc_chem.block0.out1.00001.vtk')
  species = ["He+", "OHx", "CHx", "CO", "C+", "HCO+", "H2", "H+", "H3+", "H2+", 
             "O+", "Si+", "E", 'rho', 'press']
  ns = len(species)
  err_all = np.zeros(ns)
  for i in xrange(ns):
    s = species[i]
    xs_ref = data_ref[s]
    xs_new = data_new[s]
    err_all[i] = (abs(xs_ref - xs_new) / (abs(xs_ref) + 1e-20) ).max()
  err_max = err_all.max()
  if err_max < err_control:
    return True
  else:
    print "err_max=", err_max
    return False
