#!/bin/bash

KDSDIR="$(dirname $(dirname "$0"))"

display_usage(){
	echo -e "Usage: kdtool templates dest [options]\n"
	echo -e "Copy to dest templates for KDSource usage in Python, or for interacting with"
	echo -e "Monte Carlo codes.\n"
	echo -e "Options:"
	echo -e "\t--mcstas:   Copy templates for using McStas."
	echo -e "\t--tripoli:  Copy templates for using TRIPOLI-4."
	echo -e "\t--all:      Copy all templates."
	echo -e "\t-h, --help: Display usage instructions."
}

opt_mcstas=0
opt_tripoli=0
DEST=""
while (( "$#" )); do
	case "$1" in
		"-h"|"--help")
			display_usage
			exit 0
			;;
		"--mcstas")
			opt_mcstas=1
			shift 1
			;;
		"--tripoli")
			opt_tripoli=1
			shift 1
			;;
		"--all")
			opt_mcstas=1
			opt_tripoli=1
			shift 1
			;;
		*)
			if [[ "$DEST" == "" ]]; then
				DEST="$1"
				shift 1
			else
				echo "Invalid arguments. Use -h or --help for help."
				exit 1
			fi
			;;
	esac
done

if [[ "$DEST" == "" ]]; then
	echo "No destiny. Use -h or --help for help."
	exit 1
fi
if [ ! -d "$DEST" ]; then
	echo "Created $DEST directory"
	mkdir "$DEST"
else
	echo "Using existing $DEST directory"
fi

cp $KDSDIR/templates/*.ipynb "$DEST"
echo "Copied templates for preproc/postproc in Python (Jupyter Notebook)."
if [[ opt_mcstas -eq 1 ]]; then
	cp $KDSDIR/templates/mcstas/* "$DEST"
	echo "Copied templates for McStas."
fi 
if [[ opt_tripoli -eq 1 ]]; then
	cp $KDSDIR/templates/tripoli/* "$DEST"
	echo "Copied templates for TRIPOLI-4."
fi
