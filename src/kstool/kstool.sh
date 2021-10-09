#!/bin/bash

KSDIR="$(dirname $(dirname "$0"))"

display_usage(){
	echo -e "Usage: kstool [options]\n"
	echo -e "KSource is a Monte Carlo calculations assistance tool. It implements particles"
	echo -e "density estimation and sampling by means of Kernel Density Estimation method.\n"
	echo -e "Options:"
	echo -e "\tresample:   Resample particles based on a ksource XML file."
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
		export LD_LIBRARY_PATH=$KSDIR/lib
		$KSDIR/bin/kstool-resample "${@:2}"
		exit 0
		;;
	"templates")
		$KSDIR/bin/kstool-templates "${@:2}"
		exit 0
		;;
	"beamtest")
		export LD_LIBRARY_PATH=$KSDIR/lib
		$KSDIR/bin/kstool-beamtest "${@:2}"
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

echo "Invalid arguments. Use -h or --help for help."
exit 1
