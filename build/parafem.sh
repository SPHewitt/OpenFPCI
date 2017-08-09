#!/usr/bin/env bash

#? paragait 0.1.0
#? Copyright (C) 2017 Nicolas Gruel
#? License MIT
#? This is free software: you are free to change and redistribute it.
#? There is NO WARRANTY, to the extent permitted by law.

version=$(grep "^#?"  "$0" | cut -c 4-)

# Usage info
show_help() {
    cat << EOF
    Usage: ${0##*/} [ -d WORKING_DIR ] [ -l LOGFILE ] [ -V ] [ -h ]

       -h              display this help and exit
       -d WORKING_DIR  write the result to OUTFILE instead of standard output.
       -l LOGFILE      Name of the logfile 
       -V              print version of the script
EOF
}

optspec="vVhd:l:"
while getopts "${optspec}" opt; do
    case ${opt} in
        # for options with required arguments, an additional shift is required
        d )
            WORKING_DIR="${OPTARG}"
            ;;
        l )
            LOGFILE="${OPTARG}"
            ;;
        v )
            verbose=$((verbose+1))
            ;;
        V ) 
            echo "${version}"
            exit 1
            ;;
        h ) show_help
            exit;;
        *)  echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
    esac
done

#Define working directory as the one where the script is executed.

if [ -z "${WORKING_DIR}" ]; then
    WORKING_DIR=`pwd`
    echo "Working directory: " ${WORKING_DIR}
fi

if [ ! -d $WORKING_DIR ]; then
    mkdir $WORKING_DIR
fi

# Set logfile
if [ -z "${LOGFILE}" ]; then
    LOGFILE=`echo $0 | sed 's/.sh/.log/'`
    logfile=$WORKING_DIR/$LOGFILE
else
    logfile=$LOGFILE
fi

if [ ! -f $logfile ]; then
    echo "working directory: " $WORKING_DIR > $logfile
fi

echo "Start Parafem compilation" 
echo "Start Parafem compilation" >> $logfile

####################################################################
# Parafem compilation and installation from source code (repository)
####################################################################

cd $WORKING_DIR

# Parafem is provided through sourceforge and a subversion repository
# I prefer to use git to download the code 

if [ ! -d parafem-code ]; then
    #git svn clone https://svn.code.sf.net/p/parafem/code/trunk parafem-code
    svn co https://svn.code.sf.net/p/parafem/code/trunk parafem-code >> $logfile 2>&1  
    cd parafem-code/parafem
else
    cd parafem-code/parafem
    #git svn fetch
    #git rebase git-svn
    svn update >> $logfile 2>&1  
fi

# Compilation for linuxdesktop
MACHINE=linuxdesktop ./make-parafem >> $logfile 2>&1   

# Testing parafem without and with mpi
if [ ! -d test ]; then
 mkdir test
fi

cp examples/5th_ed/p121/demo/p121_demo.mg test/
cd test

echo "Test parafem without mpi" >> $logfile
../bin/p12meshgen p121_demo >> $logfile 2>&1
../bin/p121 p121_demo >> $logfile 2>&1

# mpi test
echo "Test parafem with mpi" >> $logfile
if [ "$(whoami)" == "root" ]; then
    mpirun --allow-run-as-root ../bin/p121 p121_demo >> $logfile 2>&1
else
    mpirun ../bin/p121 p121_demo >> $logfile 2>&1
fi

cd $WORKING_DIR

echo "End of Parafem compilation"
echo "End of Parafem compilation" >> $logfile
