#!/bin/bash

echo -e "Este es el script de instalacion del paquete KSource.\n"

if [[ $# -eq 0 ]]; then
	echo -e "Uso: install-ksource.sh /path/to/ksourceinstall\n"
	exit 0
fi

KS=$1
read -r -p "Se instalara en $KS, ¿continuar? s/n " response
if [[ ! "$response" =~ ^[sS]$ ]]; then exit 0; fi
echo ""
mkdir $KS

# Examples
cp -r examples $KS/examples

# Templates
cp -r templates $KS/templates

# McStas components
cp -r mcstas $KS/mcstas

# MCPL
cd mcpl
cmake . -DCMAKE_INSTALL_PREFIX=$KS/mcpl -DBUILD_WITHG4=OFF
make install
cd ..
echo -e "\nInstalacion de MCPL completa"

# Python
mkdir $KS/python
cp -r ksource_py $KS/python

# C
mkdir $KS/include
cp ksource_c/*.h $KS/include
mkdir $KS/lib
export LD_LIBRARY_PATH=$KS/mcpl/lib/:$LD_LIBRARY_PATH
gcc ksource_c/*.c -o $KS/lib/libksource.so -lm -lmcpl -I$KS/include -I$KS/mcpl/include -L$KS/mcpl/lib -shared -fPIC

# Binaries
mkdir $KS/bin
cp scripts/ksource.sh $KS/bin/ksource
cp scripts/templates.sh $KS/bin/ksource-templates
export LD_LIBRARY_PATH=$KS/lib/:$LD_LIBRARY_PATH
gcc scripts/resample.c -o $KS/bin/ksource-resample -lksource -lmcpl -lm -I$KS/include -L$KS/lib -I$KS/mcpl/include -L$KS/mcpl/lib

echo -e "\n¡Instalacion exitosa!"