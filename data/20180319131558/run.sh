#!/bin/sh
\rm -rf craft.difx
\rm -f *.output
mpirun -machinefile machines -np 5 /Users/phi196/code/difx/bin/mpifxcorr craft.input
