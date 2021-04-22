#!/bin/bash

##### Input #####

INPUT=1_guia_n # 0_Benchmark_n # 
DIR_OUT=$INPUT # _bw0
N=1E6 # 2991283 # 

### End Input ###

# Directorio de outputs
rm -r $DIR_OUT
mkdir $DIR_OUT
cp $INPUT.instr $DIR_OUT
cd $DIR_OUT

# Ejecutar McStas
KS="$(dirname $(dirname $(which ksource)))" # Path de KSource
export LD_LIBRARY_PATH=$KS/lib:$KS/mcpl/lib
mcstas $INPUT.instr -I$KS/mcstas/contrib -I$KS/include -I$KS/mcpl/include
gcc $INPUT.c -o $INPUT.out -lksource -lmcpl -lm -I$KS/include -L$KS/lib -I$KS/mcpl/include -L$KS/mcpl/lib
./$INPUT.out -n $N | tee bash.out

cd ..