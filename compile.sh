#!/usr/bin/env bash

# Description of what the script is doing
#
# 1. Compiles and installs 
#	Foam-Extend-x
#	ParaFEM
#	OpenFPCI 

# Usage Info
show_help() {
    cat << EOF
    Usage: ${0##*/} [ -d WORKING_DIR ] [ -V ] [ -t ] [ -h ]

	-h display this help and exit
	-d WORKINGDIR, write the result to OUTFILE instead of stanradr output
	-l LOGFILE NAme of logfile
	-V Print version of the script
	-t Run tests

EOF
}

optspec="vVhd:l:t"
while getopts "${optspec}" opt; do
    case ${opt} in
        # for options with required arguments, an additional shift is required
        d )
            WORKING_DIR="${OPTARG}"
            ;;
        l )
            LOGFILE="${OPTARG}"
            ;;
        t )
            RUNTEST=0
            ;;
        v )
            verbose=$((verbose+1))
            ;;
        V )
            echo "${version}"
            exit 1
            ;;
        h ) show_help; exit;;

        *) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
    esac
done

# Enter into Source file directory
cd  src

# Run the Foam-Extend Compile Script
./foam-extend.sh

# Run the paraFEM Compile Script
./parafem.sh

# Run the OpenFPCI Compile Script
./openfpci.sh






