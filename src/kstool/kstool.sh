#!/bin/bash

KSDIR="$(dirname $(dirname "$0"))"

display_usage(){
	echo -e "Uso: kstool [opciones]\n"
	echo -e "KSource es un codigo de asistencia al calculo de blindajes por el metodo"
	echo -e "Monte Carlo. Implementa el muestreo de particulas por el metodo Kernel"
	echo -e "Density Estimation.\n"
	echo -e "Opciones:"
	echo -e "\ttemplates:   copiar plantillas para calculo de blindajes."
	echo -e "\tresample:    resamplear muestras en base a archivo ksource (XML)."
	echo -e "\t[Cualquier comando de mcpl]"
}

if [ $# -eq 0 ]; then 
	echo "Sin argumentos. Use -h o --help para ayuda."
	exit 1
fi

case "$1" in
	"templates")
		$KSDIR/bin/kstool-templates "${@:2}"
		exit 0
		;;
	"resample")
		export LD_LIBRARY_PATH=$KSDIR/lib
		$KSDIR/bin/kstool-resample "${@:2}"
		exit 0
		;;
	"-h"|"--help")
		display_usage
		exit 0
		;;
	*)
		if [ -f "$KSDIR/bin/$1" ]; then
			$KSDIR/bin/$1 "${@:2}"
			exit 0
		fi
		;;
esac

echo "Argumentos invalidos. Use -h o --help para ayuda."
exit 1
