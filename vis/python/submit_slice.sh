#!/bin/sh
NODES=`cat $COBALT_NODEFILE | wc -l`
PROCS=$((NODES * 12))

cd /projects/MMRadiation/yanfei/AGNGlobal/Data

mpirun -f $COBALT_NODEFILE -n $PROCS python ./yt_slice.py 
