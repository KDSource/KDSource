#!/bin/bash

##### Input #####
INPUT='kds_instrument_example'    # McStas istrument name, without .instr
OutDir='example' # Directory where execute simulation and store output files
N=1E6          # Number of particles to 
KDSourceLibPath=$1
### End Input ###

# Execute McStas
rm -rf $OutDir
export LD_LIBRARY_PATH=$KDSourceLibPath
cp ../../mcstas/contrib/KDSource.comp .
mcrun -c $INPUT --dir=$OutDir -n $N
rm KDSource.comp kds_instrument_example.c kds_instrument_example.out
echo Simulation ended successfully