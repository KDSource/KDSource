#!/bin/bash

##### Input #####

INPUT=2_bunker_n # 4_hall_n # 5_E_n # 
DIR_OUT=$INPUT # 3_bunker_act #  4_hall_act # 5_E_act # 
EXT_SOURCE=KSource.c # Distrib2Tripoli_hall.c # Tally2Tripoli_activ.c # 

### End Input ###

# Directorio de outputs
rm -r $DIR_OUT
mkdir $DIR_OUT
cp $INPUT.t4 $DIR_OUT
cd $DIR_OUT

# Compilar fuente externa:
# gcc lib.c source.c -shared -o source.so -I ~lib_folder -lm -fPIC
KSOURCE=~/Documents/Maestria/KSource/ksource_c
COM=~/Documents/Maestria/KSource/Arch_com
gcc $KSOURCE/ksource.c $KSOURCE/metrics.c $KSOURCE/plists.c $KSOURCE/aux.c $COM/$EXT_SOURCE -o source.so -I $KSOURCE -shared -lm -fPIC
export LD_LIBRARY_PATH=$(pwd)

# Ejecutar tripoli4
# /path/tripoli.exe -d data_file -s NJOY -c PATHFILE -o output_file -l english -u -p parallelism_file -t bds
TRIPOLI4=/opt/TRIPOLI-4/CODE/bin/linux-intel-64/static_tripoli4
PATHFILE=/opt/TRIPOLI-4/Env/t4path.ceav5
$TRIPOLI4 -d $INPUT.t4 -s NJOY -c $PATHFILE -o $INPUT.out -l english -u

cd ..