#!/usr/bin/env bash

#? paragaitdep 0.1.0
#? Copyright (C) 2017 Nicolas Gruel
#? License MIT
#? This is free software: you are free to change and redistribute it.
#? There is NO WARRANTY, to the extent permitted by law.

# Description of what the script is doing
#
# 1. Check that parafem is installed 
#
#
#


version=$(grep "^#?"  "$0" | cut -c 4-)

# Usage info
show_help() {
    cat << EOF
    Usage: ${0##*/} [ -d WORKING_DIR ] [ -l LOGFILE ] [ -V ] [ -h ]

       -h display this help and exit
       -d WORKING_DIR Name of the working directory
       -p PARAFEM_DIR Name of the directory where parafem is located
       -l LOGFILE Name of the logfile (optional by default openfpci.log)
       -V print version of the script
EOF
}

optspec="vVhd:p:l:"
while getopts "${optspec}" opt; do
    case ${opt} in
        # for options with required arguments, an additional shift is required
        d )
            WORKING_DIR="${OPTARG}"
            ;;
	p )
	    PARAFEM_DIR="${OPTARG}"
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
        h ) show_help; exit;;

        *) echo "$0: error - unrecognized option $1" 1>&2; exit 1;;
    esac
done

#Define working directory as the one where the script is executed.

if [ -z "${WORKING_DIR}" ]; then
    WORKING_DIR=`pwd`
    echo "Working directory: " ${WORKING_DIR}
fi
echo $WORKING_DIR

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

#READ common variable if present
if [ -f version.sh ]; then
    source version.sh
else
    FOAMEXTEND_VERSION=4.0
    PARAFEM_VERSION=5.0.3
fi

echo "Start OpenFPCI compilation" 
echo "Start OpenFPCI compilation" >> $logfile

if [ -z "${PARAFEM_DIR}" ]; then
    PARAFEM_DIR=$WORKING_DIR/parafem-code/parafem
fi
echo "PARAFEM_DIR="$PARAFEM_DIR >> $logfile

if [ ! -f $PARAFEM_DIR/lib/libParaFEM_mpi.${PARAFEM_VERSION}.a ]; then
    echo "Parafem not present please install it"
    exit 1
else
    echo "Parafem lib used: "  $PARAFEM_DIR/lib/libParaFEM_mpi.${PARAFEM_VERSION}.a >> $logfile
fi

############################################################
# Compilation and Installation of OpenFPCI requirement (FSI)
############################################################

# A Fluid Structure Interaction library which contains a framework 
# for easy implementation of new structural models.
# Installations step from:
# https://openfoamwiki.net/index.php/Extend-bazaar/Toolkits/Fluid-structure_interaction#Install_on_foam-extend-4.0

echo $WM_PROJECT_USER_DIR
echo "FSI compilation"
echo "FSI compilation" >> $logfile

cd $HOME/foam/foam-extend-$FOAMEXTEND_VERSION
source etc/bashrc

echo "Fsi Compilation and Installation" >> $logfile
mkdir -p $WM_PROJECT_USER_DIR
cd $WM_PROJECT_USER_DIR

if [ -f $WORKING_DIR/Fsi_40.tar.gz ]; then
    cp $WORKING_DIR/Fsi_40.tar.gz $WM_PROJECT_USER_DIR
else
    echo "Download Fsi"
    echo "Download Fsi" >> $logfile 2>&1
    wget -c https://openfoamwiki.net/images/d/d6/Fsi_40.tar.gz >> $logfile 2>&1  
fi

if [ ! -d FluidSolidInteraction ]; then
    echo "Uncompress Fsi"
    echo "Uncompress Fsi" >> $logfile 2>&1
    tar -xzf Fsi_40.tar.gz >> $logfile 2>&1  
fi

# build the Toolkit
echo "Fsi compilation"
echo "Fsi compilation" >> $logfile
cd $WM_PROJECT_USER_DIR/FluidSolidInteraction/src
./Allwmake >> $logfile 2>&1

# Corrections to the Fuid Structure Interaction library
# Foam Extend 4.0 Updated the fluxRequired methods

cd $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/fluidSolvers/
sed -i '/setRefCell/a \ \ \ \ mesh.schemesDict().setFluxRequired(p_.name());' icoFluid/icoFluid.C
sed -i '/setRefCell/a \ \ \ \ mesh.schemesDict().setFluxRequired(p().name());' consistentIcoFluid/consistentIcoFluid.C

#######################################
# OpenFPCI compilation and installation
#######################################

echo "OpenFPCI compilation"
echo "OpenFPCI compilation" >> $logfile

cd $WORKING_DIR
#if [ ! -d OpenFPCI ]; then
#    git clone git://github.com/SPHewitt/OpenFPCI.git
#fi

cp -r paraFEM $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/solidSolvers/

cd  $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/solidSolvers

OPENMPI_COMPILE_FLAGS="`mpif90 --showme:compile`"
OPENMPI_LINK_FLAGS="`mpif90 --showme:link`"

echo "solidSolvers/paraFEM/DyParaFEMSolid.C" > paraFEM.files
echo "solidSolvers/paraFEM/DyParaFEMSolidSolve.C" >> paraFEM.files
echo "" >> paraFEM.files
cat paraFEM.files $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/Make/files > tmp.files
mv tmp.files $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/Make/files
    
# Modify the options files
#sed -i 's/EXE_INC = -std=c++11 \\/EXE_INC = -std=c++11 \\\n  '"${OPENMPI_COMPILE_FLAGS//\//\\/}"'\\/' $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/Make/options

#sed -i 's/EXE_LIBS = /EXE_LIBS = \\\n    '"${WM_PROJECT_USER_DIR//\//\\/}"'\/FluidSolidInteraction\/src\/fluidSolidInteraction\/solidSolvers\/paraFEM\/dyparafemsubroutines.o \\\n    -L\/'"${PARAFEM_DIR//\//\\/}"'\/lib -lParaFEM_mpi.5.0.3  -L\/'"${PARAFEM_DIR//\//\\/}"'\/lib -larpack_linuxdesktop -lgfortran \\\n    '"${OPENMPI_LINK_FLAGS//\//\\/}"'/'  $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/Make/options

sed -i 's/EXE_INC = -std=c++11 \\/EXE_INC = -std=c++11 \\\n  '"${OPENMPI_COMPILE_FLAGS//\//\\/}"'\\/' $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/solvers/fsiFoam/Make/options

sed -i 's/EXE_LIBS = /EXE_LIBS = \\\n    '"${WM_PROJECT_USER_DIR//\//\\/}"'\/FluidSolidInteraction\/src\/fluidSolidInteraction\/solidSolvers\/paraFEM\/dyparafemsubroutines.o \\\n    -L\/'"${PARAFEM_DIR//\//\\/}"'\/lib -lParaFEM_mpi.5.0.3  -L\/'"${PARAFEM_DIR//\//\\/}"'\/lib -larpack_linuxdesktop -lgfortran \\\n    '"${OPENMPI_LINK_FLAGS//\//\\/}"'/' $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/solvers/fsiFoam/Make/options

cd $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/solidSolvers/paraFEM
echo "gfortran -fcheck=all -c dyparafemsubroutines.f90 -o dyparafemsubroutines.o -I${PARAFEM_DIR}/include/mpi" >> $logfile
gfortran -fcheck=all -c dyparafemsubroutines.f90 -o dyparafemsubroutines.o -I${PARAFEM_DIR}/include/mpi >> $logfile 2>&1  
cd $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/
./Allwmake >> $logfile 2>&1

# Fixed some dependencies as indicated in the wiki to run the tutorial (not tested!)
# https://openfoamwiki.net/index.php/Extend-bazaar/Toolkits/Fluid-structure_interaction#Install_on_foam-extend-4.0

cd ..
find run -name options | while read item
do
  sed -i -e 's=$(WM_PROJECT_DIR)/applications/solvers/FSI=$(WM_PROJECT_USER_DIR)/FluidSolidInteractio#n/src=' $item
  sed -i -e 's=$(WM_THIRD_PARTY_DIR)/packages/eigen3=$(WM_PROJECT_USER_DIR)/FluidSolidInteraction/src#/ThirdParty/eigen3=' $item
done

cd $WORKING_DIR

echo "End of OpenFPCI compilation"
echo "End of OpenFPCI compilation" >> $logfile

##################
# OpenFPCI Testing 
##################

cd $WORKING_DIR

# Copy Hron Turek test Case to user run folder

cp -r HronTurek $FOAM_RUN/. 

cd $FOAM_RUN/HronTurek/fluid

./Allrun






