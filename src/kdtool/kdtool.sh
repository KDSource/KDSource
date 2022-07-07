#!/bin/bash

KDSDIR="$(dirname $(dirname "$0"))"

display_usage(){
	echo -e "Usage: kdtool [options]\n"
	echo -e "KDSource is a Monte Carlo calculations assistance tool. It implements particles"
	echo -e "density estimation and sampling by means of Kernel Density Estimation method.\n"
	echo -e "Options:"
	echo -e "\tresample:   Resample particles based on a kdsource XML file."
	echo -e "\ttemplates:  Copy templates for Monte Carlo calculations."
	echo -e "\tbeamtest:   Test source with simple beam calculation."
	echo -e "\t[Any MCPL command]"
	echo -e "\t-h, --help: Display usage instructions."
}

if [ $# -eq 0 ]; then 
	echo "No arguments. Use -h or --help for help."
	exit 1
fi

case "$1" in
	"resample")
		export LD_LIBRARY_PATH=$KDSDIR/lib
		$KDSDIR/bin/kdtool-resample "${@:2}"
		exit 0
		;;
	"templates")
		$KDSDIR/bin/kdtool-templates "${@:2}"
		exit 0
		;;
	"beamtest")
		export LD_LIBRARY_PATH=$KDSDIR/lib
		$KDSDIR/bin/kdtool-beamtest "${@:2}"
		exit 0
		;;
	"-h"|"--help")
		display_usage
		exit 0
		;;
	*)
		if [ -f "$KDSDIR/bin/$1" ]; then
			$KDSDIR/bin/$1 "${@:2}"
			exit 0
		fi
		;;
esac

echo "Invalid arguments. Use -h or --help for help."
exit 1
