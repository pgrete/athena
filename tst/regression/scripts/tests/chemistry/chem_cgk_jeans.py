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
  athena.configure('pp', 
      prob='read_vtk',
      chemistry='gow16', 
      radiation='loc_jeans',
      #cvode_path='/usr/local',
      #cxx="g++"
      cvode_path='/home/munan/install',
      cxx="icc"
      )
  athena.make()

def run():
  arguments = [ 
          'problem/G0=1.0',
          'time/tlim=3e14',
          'time/nlim=3',
          'chemistry/isH2RVcooling=0',
          'chemistry/reltol=1.0e-5',
          'problem/vtkfile=../data/chem_cgk_input.vtk',
          'mesh/nx1=16',
          'mesh/nx2=16',
          'mesh/nx3=32'
          ]
  athena.run('chemistry/athinput.cgk_32pc_chem', arguments)

def analyze():
  err_control = 2e-2
  _,_,_,data_ref = athena_read.vtk('data/chem_cgk_output.vtk')
  _,_,_,data_new = athena_read.vtk('bin/cgk_32pc_chem.block0.out1.00001.vtk')
  species = ["He+", "OHx", "CHx", "CO", "C+", "HCO+", "H2", "H+", "H3+", "H2+", 
             "O+", "Si+", "E"]
  ns = len(species)
  err_all = np.zeros(ns)
  err_dense = np.zeros(ns)
  indx_rho = data_ref["rho"] > 1e-1
  def get_rel_diff(s):
    xs_ref = data_ref[s]
    xs_new = data_new[s]
    err_s = abs(xs_ref - xs_new) / (abs(xs_ref) + 1e-10)
    return err_s
  #density and pressure should exactly agree
  err_rho = get_rel_diff("rho").max()
  err_press = get_rel_diff("press").max()
  print "relative error: rho: {0:.2e}, press: {1:.2e}".format(err_rho, err_press)
  if err_rho > 1e-10 or err_press > 1e-10:
    return False
  #species
  for i in xrange(ns):
    s = species[i]
    err_s = get_rel_diff(s)
    err_all[i] = err_s.max()
    err_dense[i] = err_s[indx_rho].max()
    print "{0}: all: {1:.2e}, dense: {2:.2e}".format(s, err_all[i], err_dense[i])
  err_max = err_dense.max()
  if err_max < err_control:
    return True
  else:
    print "err_max=", err_max
    return False
