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
 
# Define WORK as directory where OpenFPCI is contained

if [ -z "${WORK}" ]; then
    WORK=$WORKING_DIR/../../
    echo "Work: " ${WORK}
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
BASHFILE=etc/bashrc
SETTINGSFILE=etc/settings.sh

CFile=wmake/rules/linux64Gcc/c
CppFILE=wmake/rules/linux64Gcc/c++
CppOptFILE=wmake/rules/linux64Gcc/c++Opt
GENERALFILE=wmake/rules/linux64Gcc/general
MPLIBMPICHFILE=wmake/rules/linux64Gcc/mplibMPICH
MAKEFILE=wmake/src/Makefile


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
#    if [ ! -d $HOME/foam ]; then
#        mkdir -p $HOME/foam
#        ln -s `pwd` $HOME/foam/
#    fi
    cd $WORKING_DIR/foam/foam-extend-$FOAMEXTEND_VERSION


    # Modifification of preference file to use MPI from system

    cp $PREFSFILE-EXAMPLE $PREFSFILE
    #sed -i 's///' $PREFSFILE

    # Use Compiler from system
    sed -i 's/#compilerInstall=System/compilerInstall=System/' $PREFSFILE

    # Use CMAKE from system
    sed -i 's/#export CMAKE_SYSTEM=1/export CMAKE_SYSTEM=1/' $PREFSFILE
    sed -i 's/#export CMAKE_DIR=path_to_system_installed_cmake/export CMAKE_DIR=\/usr/' $PREFSFILE
    sed -i 's/#export CMAKE_BIN_DIR=$CMAKE_DIR\/bin/export CMAKE_BIN_DIR=$CMAKE_DIR\/bin/' $PREFSFILE

    # Edit bashrc File

    # Set MPI to MPICH
    sed -i 's/${WM_MPLIB:=OPENMPI}; export WM_MPLIB/${WM_MPLIB:=MPICH}; export WM_MPLIB/' $BASHFILE
    
    # Set Number of Processors to Install
    sed -i 's/export WM_NCOMPPROCS/export WM_NCOMPPROCS=10/' $BASHFILE

    # Include dynamic linking
    sed -i 's/WM_ARCH=linux64/WM_ARCH=linux64; export XTPE_LINK_TYPE=shared/' $BASHFILE

    # Edit Settings.sh
    sed -i 's/mpi_version=mpich-1.2.4/ /' $SETTINGSFILE
    sed -i 's/export MPI_HOME=$WM_THIRD_PARTY_DIR\/$mpi_version/export MPI_ARCH_PATH=$MPICH_DIR /' $SETTINGSFILE
    sed -i 's/export MPI_ARCH_PATH=$MPI_HOME\/platforms\/$WM_OPTIONS/export MPICH_PATH=$MPI_ARCH_PATH /' $SETTINGSFILE
    sed -i 's/export FOAM_MPI_LIBBIN=$FOAM_LIBBIN\/$mpi_version/export FOAM_MPI_LIBBIN=$FOAM_LIBBIN\/mpich/' $SETTINGSFILE
    sed -i 's/unset mpi_version/export MPICH_ROOT=$MPI_ARCH_PATH/' $SETTINGSFILE

    # Edit Wmake files

    # C file
    sed -i 's/cc          = gcc -m64/cc          = cc/' $CFILE

    # C++ file
    sed -i 's/cc          = gcc -m64/CC          = CC/' $CppFILE
    sed -i 's/LINKEXE     = $(CC) $(c++FLAGS) -Xlinker --add-needed -Xlinker --no-as-needed/LINKEXE     = $(CC) $(c++FLAGS) -dynamic -Xlinker --add-needed -Xlinker --no-as-needed/' $CppFILE

    # C++Opt file
    sed -i 's/c++OPT      = -O3/c++OPT      = -O3 -ftree-vectorize/' $CppOptFILE

    # General file
    sed -i 's/-ldl/-ldl -A64/' $GENERALFILE

    # MPLIBMPICH
    sed -i 's/PFLAGS     =/PFLAGS     = -DMPICH_SKIP_MPICXX/' $MPLIBMPICHFILE
    sed -i 's/PINC       = -I$(MPI_ARCH_PATH)\/include/ /' $MPLIBMPICHFILE
    sed -i 's/PLIBS      = -L$(MPI_ARCH_PATH)\/lib -lmpich/ /' $MPLIBMPICHFILE

    # MAKFILE
    sed -  '/include $(RULES)\/$(WM_LINK_LANGUAGE)/a cc = gcc' input $MAKEFILE
    sed -  '/cc = gcc/a CC = g++' input $MAKEFILE

    source etc/bashrc >> $logfile 2>&1
    ./Allwmake.firstInstall <<< "y" >> $logfile 2>&1
else
    cd foam-extend-4.0
    git pull >> $logfile 2>&1
    source etc/bashrc >> $logfile 2>&1
    cd $HOME/foam/foam-extend-$FOAMEXTEND_VERSION
    ./Allwmake -update <<< "y" >> $logfile 2>&1
fi

# Add alias to bashrc file
echo "Adding alias \"fe40\"to bashrc"
echo "Adding alias \"fe40\"to bashrc">> $logfile
echo "alias fe40='source \$HOME/foam/foam-extend-4.0/etc/bashrc'" >> $HOME/.bashrc

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
