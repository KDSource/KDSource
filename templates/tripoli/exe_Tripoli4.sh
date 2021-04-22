#!/bin/bash

##### Input #####

TRIPOLI4=/opt/TRIPOLI-4

INPUT=2_bunker_n # 
DIR_OUT=$INPUT # 2_bunker_act # 
EXT_SOURCE=KSource.c

### End Input ###

KS="$(dirname $(dirname $(which ksource)))"

# Directorio de outputs
rm -r $DIR_OUT
mkdir $DIR_OUT
cp $INPUT.t4 $DIR_OUT
cd $DIR_OUT

# Compilar fuente externa:
# gcc source.c lib.c -o source.so -lm -Ilib_folder -shared -fPIC
KS="$(dirname $(dirname $(which ksource)))" # Path de KSource
export LD_LIBRARY_PATH=$(pwd):$KS/lib:$KS/mcpl/lib
gcc ../$EXT_SOURCE -o source.so -lksource -lmcpl -lm -I$KS/include -L$KS/lib -I$KS/mcpl/include -L$KS/mcpl/lib -shared -fPIC

# Ejecutar tripoli4
# /path/tripoli.exe -d data_file -s NJOY -c PATHFILE -o output_file -l english -u -p parallelism_file -t bds
T4EXE=$TRIPOLI4/CODE/bin/linux-intel-64/static_tripoli4
PATHFILE=$TRIPOLI4/Env/t4path.ceav5
$T4EXE -d $INPUT.t4 -s NJOY -c $PATHFILE -o $INPUT.out -l english -u | tee bash.out

cd ..
