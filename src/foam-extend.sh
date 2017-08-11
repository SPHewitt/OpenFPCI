#!/usr/bin/env bash

#? paragaitdep 0.1.0
#? Copyright (C) 2017 Nicolas Gruel
#? License MIT
#? This is free software: you are free to change and redistribute it.
#? There is NO WARRANTY, to the extent permitted by law.

RUNTEST=0

version=$(grep "^#?"  "$0" | cut -c 4-)

# Usage info
show_help() {
    cat << EOF
    Usage: ${0##*/} [ -d WORKING_DIR ] [ -V ] [ -t ] [ -h ]

       -h display this help and exit
       -d WORKINGDIR  write the result to OUTFILE instead of standard output.
       -l LOGFILE Name of the logfile
       -V print version of the script
       -t Run Foam-dev tests

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

#Define working directory as the one where the script is executed.

if [ -z "${WORKING_DIR}" ]; then
    WORKING_DIR=`pwd`
    echo "Working directory: " ${WORKING_DIR}
fi

if [ ! -d $WORKING_DIR ]; then
    mkdir $WORKING_DIR
fi

if [ -z "${LOGFILE}" ]; then
    LOGFILE=`echo $0 | sed 's/.sh/.log/'`
    logfile=$WORKING_DIR/$LOGFILE
else
    logfile=$LOGFILE
fi

if [ ! -f $logfile ]; then
    echo "working directory: " $WORKING_DIR > $logfile
fi

##########################################
# FOAM-EXTEND compilation and installation
##########################################

echo "Start FOAM-extend compilation"
echo "Start FOAM-extend compilation" >> $logfile

#READ common variable if present
if [ ! -f version.sh ]; then
    source version.sh
else
    FOAMEXTEND_VERSION=4.0
fi
cd $WORKING_DIR

PREFSFILE=etc/prefs.sh

# Define LIBS_SYSTEM where the dependencies will be install in function if the user is root (Docker) or not.
if [ ! "$(whoami)" == "root" ]; then
    DEPS_PATH=$HOME/.local
else
    DEPS_PATH=/usr/local
fi

# Foam extend
cd $WORKING_DIR
if [ ! -d foam-extend-$FOAMEXTEND_VERSION ]; then
    git clone https://git.code.sf.net/p/foam-extend/foam-extend-$FOAMEXTEND_VERSION foam-extend-$FOAMEXTEND_VERSION >> $logfile 2>&1
    cd foam-extend-$FOAMEXTEND_VERSION
    if [ ! -d $HOME/foam ]; then
        mkdir -p $HOME/foam
        ln -s `pwd` $HOME/foam/
    fi
    cd $HOME/foam/foam-extend-$FOAMEXTEND_VERSION


    # Modifification of preference file to use MPI from system

    cp $PREFSFILE-EXAMPLE $PREFSFILE
    #sed -i 's///' $PREFSFILE

    # Use Compiler from system
    sed -i 's/#compilerInstall=System/compilerInstall=System/' $PREFSFILE

    #Use MPI from system
    sed -i 's/#export WM_MPLIB=SYSTEMOPENMPI/export WM_MPLIB=SYSTEMOPENMPI/' $PREFSFILE
    sed -i 's/#export OPENMPI_DIR=path_to_system_installed_openmpi/export OPENMPI_DIR=\/usr/' $PREFSFILE
    sed -i 's/#export OPENMPI_BIN_DIR=$OPENMPI_DIR\/bin/export OPENMPI_BIN_DIR=$OPENMPI_DIR\/bin/' $PREFSFILE

    # Use CMAKE from system
    sed -i 's/#export CMAKE_SYSTEM=1/export CMAKE_SYSTEM=1/' $PREFSFILE
    sed -i 's/#export CMAKE_DIR=path_to_system_installed_cmake/export CMAKE_DIR=\/usr/' $PREFSFILE
    sed -i 's/#export CMAKE_BIN_DIR=$CMAKE_DIR\/bin/export CMAKE_BIN_DIR=$CMAKE_DIR\/bin/' $PREFSFILE

    # Use Python from system
    #sed -i 's/#export PYTHON_SYSTEM=1/export PYTHON_SYSTEM=1/' $PREFSFILE
    #sed -i 's/#export PYTHON_DIR=path_to_system_installed_python/export PYTHON_DIR=\/usr/' $PREFSFILE
    #sed -i 's/#export PYTHON_BIN_DIR=$PYTHON_DIR\/bin/export PYTHON_BIN_DIR=$PYTHON_DIR\/bin/' $PREFSFILE

    #Modifification of bashrc to use openmpi from system (i.e. not version 1.8.8)
    sed -i 's/export WM_THIRD_PARTY_USE_OPENMPI_188/#export WM_THIRD_PARTY_USE_OPENMPI_188/' etc/bashrc
    source etc/bashrc >> $logfile 2>&1
    ./Allwmake.firstInstall <<< "y" >> $logfile 2>&1
else
    cd foam-extend-4.0
    git pull >> $logfile 2>&1
    source etc/bashrc >> $logfile 2>&1
    cd $HOME/foam/foam-extend-$FOAMEXTEND_VERSION
    ./Allwmake -update <<< "y" >> $logfile 2>&1
fi

# Run foam-extend test and tutorials
if [ $RUNTEST ]; then
    echo "Run foam-extend test and tutorials" >> $logfile 2>&1
    # Test if foam-extend properly installed and working
    USERNAME=`whoami`
    source $HOME/foam/foam-extend-$FOAMEXTEND_VERSION/etc/bashrc

    mkdir -p $FOAM_RUN
    mkdir -p $FOAM_TUTORIALS
    cd $FOAM_TUTORIALS

    # Run all the test
    ./Alltest >> $logfile 2>&1
    ./Allrun >> $logfile 2>&1
fi

cd $WORKING_DIR

echo "End of Foam-extend compilation"
echo "End of Foam-extend compilation" >> $logfile
