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

case "$1" in
	"templates")
		$KS/bin/ksource-templates "${@:2}"
		exit 0
		;;
	"resample")
		export LD_LIBRARY_PATH=$KS/lib:$KS/mcpl/lib
		$KS/bin/ksource-resample "${@:2}"
		exit 0
		;;
	"-h"|"--help")
		display_usage
		exit 0
		;;
	*)
		if [ -f "$KS/mcpl/bin/$1" ]; then
			$KS/mcpl/bin/$1 "${@:2}"
			exit 0
		fi
		;;
esac

echo "Argumentos invalidos. Use -h o --help para ayuda."
exit 1
