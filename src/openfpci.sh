#!/usr/bin/env bash

#? OpenFPCI 1.1.0
#? Copyright (C) 2018 Sam Hewitt
#? License MIT
#? This is free software: you are free to change and redistribute it.
#? There is NO WARRANTY, to the extent permitted by law.

# Description of what the script is doing
#
# 1. Check that parafem is installed 
# 2. Check that foam-extend-x.x is installed
# 3. Install FSI Library
# 4. Install OpenFPCI
# 5. Compile and test
#

version=$(grep "^#?"  "$0" | cut -c 4-)

# Usage info
show_help() {
    cat << EOF
    Usage: ${0##*/} [ -d WORKING_DIR ] [ -l LOGFILE ] [ -p PARAFEM_DIR ] [ -f FOAM_DIR ] [ -V ] [ -h ]

       -h             : Display this help and exit
       -d WORKING_DIR : Name of the working directory
       -p PARAFEM_DIR : Parafem path, default: $HOME/parafem-code/parafem
       -f FOAM_DIR    : Foam-Extend path, default: $HOME/foam/foam-extend-4.0
       -l LOGFILE     : Name of the logfile (optional by default openfpci.log)
       -V             : Print version of the script
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
        f )
            FOAM_DIR="${OPTARG}"
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


# Check for ParaFEM Install
if [ -z ${PARAFEM_DIR} ]; then
  PARAFEM_DIR=$HOME/parafem-code/parafem
fi

echo -e "PARAFEM_DIR="$PARAFEM_DIR
echo "PARAFEM_DIR="$PARAFEM_DIR >> $logfile

if [ ! -f $PARAFEM_DIR/lib/libParaFEM_mpi.${PARAFEM_VERSION}.a ]; then
    echo "Parafem not present please install it"
    exit 1
else
    echo "Parafem lib used: "  $PARAFEM_DIR/lib/libParaFEM_mpi.${PARAFEM_VERSION}.a >> $logfile
fi

# Check for Foam-Extend 4.0 Install
if [ -z ${FOAM_DIR} ]; then
  FOAM_DIR=$HOME/foam/foam-extend-$FOAMEXTEND_VERSION
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

cd $FOAM_DIR
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
echo "Corrections to the FSI Library"
echo "Corrections to the FSI Library" >> $logfile

cd $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/fluidSolvers/
sed -i '/setRefCell/a \ \ \ \ mesh.schemesDict().setFluxRequired(p_.name());' icoFluid/icoFluid.C
sed -i '/setRefCell/a \ \ \ \ mesh.schemesDict().setFluxRequired(p().name());' consistentIcoFluid/consistentIcoFluid.C
sed -i '/setRefCell/a \ \ \ \ mesh.schemesDict().setFluxRequired(p().name());' pisoFluid/pisoFluid.C

#######################################
# OpenFPCI compilation and installation
#######################################

echo "OpenFPCI compilation"
echo "OpenFPCI compilation" >> $logfile

cd $WORKING_DIR

cp -r solidSolvers/* $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/solidSolvers/

cd  $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/solidSolvers

OPENMPI_COMPILE_FLAGS="`mpif90 --showme:compile`"
OPENMPI_LINK_FLAGS="`mpif90 --showme:link`"


echo "Modifying fluidSolidInteraction/Make/options"
echo "Modifying fluidSolidInteraction/Make/options" >> $logfile

echo "solidSolvers/paraFEM/smallStrain/femSmallStrain.C" > paraFEM.files
echo "solidSolvers/paraFEM/smallStrain/femSmallStrainSolve.C" >> paraFEM.files
echo "solidSolvers/paraFEM/largeStrain/femLargeStrain.C" >> paraFEM.files
echo "solidSolvers/paraFEM/largeStrain/femLargeStrainSolve.C" >> paraFEM.files
echo "" >> paraFEM.files

cat paraFEM.files $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/Make/files > tmp.files
mv tmp.files $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/Make/files
    
echo "Modifying fsiFoam/Make/options"
echo "Modifying fsiFoam/Make/options" >> $logfile


# Modify $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/solvers/fsiFoam/Make/options
inctext="\    $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/solidSolvers/paraFEM/fem_routines/objectFiles/parafeml.o \\"
sed -i "/EXE_LIBS = /a $inctext \\" $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/solvers/fsiFoam/Make/options

inctext="\    $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/solidSolvers/paraFEM/fem_routines/objectFiles/parafemnl.o \\"
sed -i "/parafeml.o/a $inctext \\" $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/solvers/fsiFoam/Make/options

inctext="\    $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/solidSolvers/paraFEM/fem_routines/objectFiles/parafemutils.o \\"
sed -i "/parafemnl.o/a $inctext \\" $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/solvers/fsiFoam/Make/options

inctext="\    -L$PARAFEM_DIR/lib -lParaFEM_mpi.5.0.3 \\"
sed -i "/parafemutils.o/a $inctext \\" $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/solvers/fsiFoam/Make/options

inctext="\    -L$PARAFEM_DIR/lib -larpack_linuxdesktop \\"
sed -i "/-lParaFEM_mpi.5.0.3/a $inctext \\" $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/solvers/fsiFoam/Make/options

inctext="\    -lgfortran  \\"
sed -i "/arpack_linuxdesktop/a $inctext \\" $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/solvers/fsiFoam/Make/options

inctext="\    $OPENMPI_LINK_FLAGS \\"
sed -i "/lgfortran/a $inctext \\" $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/solvers/fsiFoam/Make/options


# Compile the Fortran Subroutines
cd $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/solidSolvers/paraFEM/fem_routines
echo "gfortran -fcheck=all -c parafeml.f90 -o parafeml.o -I${PARAFEM_DIR}/include/mpi" >> $logfile
gfortran -fcheck=all -c parafeml.f90 -o parafeml.o -I${PARAFEM_DIR}/include/mpi >> $logfile 2>&1  

echo "gfortran -fcheck=all -c parafemnl.f90 -o parafemnl.o -I${PARAFEM_DIR}/include/mpi" >> $logfile
gfortran -fcheck=all -c parafemnl.f90 -o parafemnl.o -I${PARAFEM_DIR}/include/mpi >> $logfile 2>&1

echo "gfortran -fcheck=all -c parafemutils.f90 -o parafemutils.o -I${PARAFEM_DIR}/include/mpi" >> $logfile
gfortran -fcheck=all -c parafemutils.f90 -o parafemutils.o -I${PARAFEM_DIR}/include/mpi >> $logfile 2>&1

mkdir -p objectFiles
mv *.o objectFiles/.

# Fixed some dependencies as indicated in the wiki to run the tutorial (not tested!)
# https://openfoamwiki.net/index.php/Extend-bazaar/Toolkits/Fluid-structure_interaction#Install_on_foam-extend-4.0

cd $WM_PROJECT_USER_DIR/FluidSolidInteraction
find run -name options | while read item
do
  sed -i -e 's=$(WM_PROJECT_DIR)/applications/solvers/FSI=$(WM_PROJECT_USER_DIR)/FluidSolidInteraction/src=' $item
  sed -i -e 's=$(WM_THIRD_PARTY_DIR)/packages/eigen3=$(WM_PROJECT_USER_DIR)/FluidSolidInteraction/src/ThirdParty/eigen3=' $item
done


# Remake the Library including new solidSolvers ans
cd $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/
./Allwmake >> $logfile 2>&1

cd $WORKING_DIR

echo "End of OpenFPCI compilation"
echo "End of OpenFPCI compilation" >> $logfile

##################
# OpenFPCI Testing 
##################

