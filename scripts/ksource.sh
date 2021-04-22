#!/bin/bash

KS="$(dirname $(dirname "$0"))"

display_usage(){
	echo -e "Uso: ksource [opciones]\n"
	echo -e "KSource es un codigo de asistencia al calculo de blindajes por el metodo"
	echo -e "Monte Carlo. Implementa el muestreo de particulas por el metodo Kernel"
	echo -e "Density Estimation.\n"
	echo -e "Opciones:"
	echo -e "\ttemplates:   copiar plantillas para calculo de blindajes."
	echo -e "\tresample:    resamplear muestras en base a definicion de ksource."
	echo -e "\t[Cualquier comando de mcpl]"
}	

if [ $# -eq 0 ]; then 
	echo "Sin argumentos. Use -h o --help para ayuda."
	exit 1
fi

if [[ "$*" == *-h* || "$*" == *--help* ]]; then 
	display_usage
	exit 1
fi

if [ $1 == "templates" ]; then
	$KS/bin/ksource-templates "${@:2}"
	exit 0
fi

if [ $1 == "resample" ]; then
	$KS/bin/ksource-resample "${@:2}"
	exit 0
fi

if [ -f "$KS/mcpl/bin/$1" ]; then
	$KS/mcpl/bin/$1 "${@:2}"
	exit 0
fi

echo "Argumentos invalidos. Use -h o --help para ayuda."
exit 1
