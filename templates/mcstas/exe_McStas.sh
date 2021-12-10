#!/bin/bash

KDS="$(dirname $(dirname $(which kdtool)))" # KDSource path

##### Input #####

INPUT=input    # McStas istrument name, without .instr
DIR_OUT=outdir # Directory where execute simulation and store output files
N=1E6          # Number of particles to simulate

### End Input ###


# Output directory
rm -r $DIR_OUT
mkdir $DIR_OUT
cp $INPUT.instr $DIR_OUT
cd $DIR_OUT

# Execute McStas
export LD_LIBRARY_PATH=$KDS/lib
mcstas $INPUT.instr -I$KDS/mcstas/contrib -I$KDS/include
gcc $INPUT.c -o $INPUT.out -lkdsource -lmcpl -lm -I$KDS/include -L$KDS/lib
./$INPUT.out -n $N | tee bash.out

cd ..