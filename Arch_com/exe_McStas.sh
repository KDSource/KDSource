#!/bin/bash

##### Input #####

INPUT=1_guia_n # 0_Benchmark_n # 
DIR_OUT=${INPUT} # _bw0
N=1E6 # 2991283 # 

### End Input ###

# Directorio de outputs
rm -r $DIR_OUT
mkdir $DIR_OUT
cp $INPUT.instr $DIR_OUT
cd $DIR_OUT

# Ejecutar McStas
KSOURCE=~/Documents/Maestria/KSource
mcstas $INPUT.instr -I $KSOURCE/Arch_com -I $KSOURCE/ksource_c
gcc $INPUT.c -o $INPUT.out -I $KSOURCE/ksource_c -lm
./$INPUT.out -n $N > bash.out

cd ..