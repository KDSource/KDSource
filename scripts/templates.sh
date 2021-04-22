#!/bin/bash

KS="$(dirname $(dirname "$0"))"

display_usage(){
	echo -e "Uso: ksource templates dest [opciones]\n"
	echo -e "Copia a dest las plantillas para utilizar KSource desde Python.\n"
	echo -e "Opciones:"
	echo -e "\t--mcstas:  copiar tambien plantillas para utilizar McStas."
	echo -e "\t--tripoli: copiar tambien plantillas para utilizar TRIPOLI-4."
	echo -e "\t--all:     copiar todo."
}

if [ $# -eq 0 ]; then 
	echo "Sin argumentos. Use -h o --help para ayuda."
	exit 1
fi

if [[ "$*" == *-h* || "$*" == *--help* ]]; then 
	display_usage
	exit 1
fi

if [ ! -d "$1" ]; then
	echo "$1 no es un directorio."
	exit 1
fi
DEST="$1"

cp $KS/templates/*.ipynb "$DEST"
echo "Plantillas para preproc/postproc en Python (Jupyter Notebook) copiadas."
if [[ "$*" == *--mcstas* || "$*" == *--all* ]]; then
	cp $KS/templates/mcstas/* "$DEST"
	echo "Plantillas para McStas copiadas."
fi 
if [[ "$*" == *--tripoli* || "$*" == *--all* ]]; then
 	cp $KS/templates/tripoli/* "$DEST"
	echo "Plantillas para TRIPOLI-4 copiadas."
fi
