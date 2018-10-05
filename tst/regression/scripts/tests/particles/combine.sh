#!/bin/bash

# Define the base name.
BASE=UniStream
OUTPUT=$BASE.par.dat

# Find the number of blocks.
read nblocks < <(ls $BASE.block*.out1.00000.par.tab | wc -l)

# Find the number of snapshots.
read ntimes < <(ls $BASE.block0.out1.*.par.tab | wc -l)
(( --ntimes ))

for i in $(seq -f "%05g" 0 $((ntimes - 1))); do
	if (( nblocks > 1 )); then
		cat $BASE.block*.out1.$i.par.tab | sort -n | \
		sed "2,${nblocks}d" > $BASE.$i.par.tab
	else
		mv $BASE.block0.out1.$i.par.tab $BASE.$i.par.tab
	fi
	cat $BASE.block*.out1.$i.tab > $BASE.$i.tab
done

find . -name $BASE'.block*.out1.*.tab' -delete
