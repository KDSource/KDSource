#!/bin/bash

KSOURCE=~/Documents/Maestria/KSource

ln -s $KSOURCE/ksource_py $1
cp $KSOURCE/{Arch_com/{exe_McStas.sh,exe_Tripoli4.sh},Auxiliares/*} $1