#!/bin/bash

KS="$(dirname $(dirname $(which kstool)))" # KSource path

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
export LD_LIBRARY_PATH=$KS/lib
mcstas $INPUT.instr -I$KS/mcstas/contrib -I$KS/include
gcc $INPUT.c -o $INPUT.out -lksource -lmcpl -lm -I$KS/include -L$KS/lib
./$INPUT.out -n $N | tee bash.out

cd ..