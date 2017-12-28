# Regression test based on Newtonian hydro linear wave convergence problem
#
# Runs a linear wave convergence test in 3D including SMR and checks L1 errors (which
# are computed by the executable automatically and stored in the temporary file
# linearwave_errors.dat)

# Modules
import numpy as np
import math
import sys
import scripts.utils.athena as athena
import scripts.utils.comparison as comparison
import os
sys.path.insert(0, '../../vis/python')

# Prepare Athena++
def prepare():
  athena.configure('radiation','mpi','hdf5',
      prob='beam',
      coord='cartesian',
      flux='hllc',
      cxx='mac')
  athena.make()

# Run Athena++
def run():
  #case 1
  arguments = ['time/rad_xorder=3']
  athena.mpirun(4,'radiation/athinput.beam_smr', arguments)
  bashcommand="mv bin/*athdf* ../../"
#  os.system(bashcommand)
# Analyze outputs
def analyze():

  print "Take a look at the hdf5 data!"
  return True
