#!/bin/bash

TRIPOLI4=/opt/TRIPOLI-4                       # TRIPOLI-4 path
KDS="$(dirname $(dirname $(which kdsource)))" # KDSource path

##### Input #####

INPUT=input            # TRIPOLI input file name, without suffix
DIR_OUT=outdir         # Directory where execute simulation and store output files
EXT_SOURCE=extsource.c # External source C file

### End Input ###


# Output directory
rm -r $DIR_OUT
mkdir $DIR_OUT
cp $INPUT.t4 $EXT_SOURCE $DIR_OUT
cd $DIR_OUT

# Compile external source:
export LD_LIBRARY_PATH=$(pwd):$KDS/lib
gcc $EXT_SOURCE -o source.so -lkdsource -lmcpl -lm -I$KDS/include -L$KDS/lib -shared -fPIC

# Execute tripoli4
T4EXE=$TRIPOLI4/CODE/bin/linux-intel-64/static_tripoli4
PATHFILE=$TRIPOLI4/Env/t4path.ceav5
$T4EXE -d $INPUT.t4 -s NJOY -c $PATHFILE -o $INPUT.out -l english -u | tee bash.out

cd ..
