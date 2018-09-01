#!/usr/bin/env bash

# OpenFPCI 1.1.0
# Author: Sam Hewitt

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
            printf "${version}"
            exit 1
            ;;
        h ) show_help; exit;;

        *) printf "$0: error - unrecognized option $1" 1>&2; exit 1;;
    esac
done

#Define working directory as the one where the script is executed.

if [ -z "${WORKING_DIR}" ]; then
    WORKING_DIR=`pwd`
    printf "Working directory: " ${WORKING_DIR}
fi
printf $WORKING_DIR

if [ ! -d $WORKING_DIR ]; then
    mkdir $WORKING_DIR
fi

if [ -z "${LOGFILE}" ]; then
    LOGFILE=`printf $0 | sed 's/.sh/.log/'`
    logfile=$WORKING_DIR/$LOGFILE
else
    logfile=$LOGFILE
fi

if [ ! -f $logfile ]; then
    printf "working directory: " $WORKING_DIR > $logfile
fi

#READ common variable if present
if [ -f version.sh ]; then
    source version.sh
else
    FOAMEXTEND_VERSION=4.0
    PARAFEM_VERSION=5.0.3
fi

printf "\n\n--------------------------\n"
printf "Start OpenFPCI compilation" 
printf "\n--------------------------\n\n"

printf "\n\n--------------------------\n">> $logfile
printf "Start OpenFPCI compilation"      >> $logfile
printf "\n--------------------------\n\n">> $logfile

cd $WORKING_DIR/..
export OPENFPCI_DIR=`pwd`
printf "Checking OpenFPCI home directory:\n"
printf "OPENFPCI_DIR="$OPENFPCI_DIR"\n"
printf "Checking OpenFPCI home directory:\n" >> $logfile
printf "OPENFPCI_DIR="$OPENFPCI_DIR"\n"       >> $logfile

printf "\nChecking ParaFEM home directory:\n"
printf "PARAFEM_DIR="$PARAFEM_DIR"\n"
printf "\nChecking ParaFEM home directory:\n" >> $logfile
printf "PARAFEM_DIR="$PARAFEM_DIR"\n"       >> $logfile

if [ ! -f $PARAFEM_DIR/lib/libParaFEM_mpi.${PARAFEM_VERSION}.a ]; then
    printf "Parafem not found please install it"
    exit 1
else
    printf "\nParafem library used: "  $PARAFEM_DIR/lib/libParaFEM_mpi.${PARAFEM_VERSION}.a >> $logfile
fi

printf "\nChecking Foam-Extend home directory:\n"
printf "FOAM_DIR="$FOAM_DIR"\n"
printf "\nChecking Foam-Extend home directory:\n" >> $logfile
printf "FOAM_DIR="$FOAM_DIR"\n"                   >> $logfile

# TO DO:
# Place a check to ensure Foam-Extend has been installed


############################################################
# Compilation and Installation of OpenFPCI requirement (FSI)
############################################################

# A Fluid Structure Interaction library which contains a framework 
# for easy implementation of new structural models.
# Installations step from:
# https://openfoamwiki.net/index.php/Extend-bazaar/Toolkits/Fluid-structure_interaction#Install_on_foam-extend-4.0


printf "\n\n----------------------------------------\n"
printf "FSI Library Compilation and Installation"
printf "\n----------------------------------------\n\n"
printf "\n\n----------------------------------------\n" >> $logfile
printf "FSI Library Compilation and Installation"   >> $logfile
printf "\n----------------------------------------\n\n" >> $logfile


cd $FOAM_DIR
source etc/bashrc

mkdir -p $WM_PROJECT_USER_DIR
cd $WM_PROJECT_USER_DIR

if [ -f $WM_PROJECT_USER_DIR/Fsi_40.tar.gz ]; then
    printf "FSI archive file exists in: "$WM_PROJECT_USER_DIR"\n"
    printf "FSI archive file exists in: "$WM_PROJECT_USER_DIR"\n" >> $logfile
else
    printf "Downloading Fsi\n"
    printf "Downloading Fsi\n" >> $logfile 2>&1
    wget -c https://openfoamwiki.net/images/d/d6/Fsi_40.tar.gz >> $logfile 2>&1  
fi

if [ ! -d FluidSolidInteraction ]; then
    printf "\nUncompress Fsi"
    printf "\nUncompress Fsi" >> $logfile 2>&1
    tar -xzf Fsi_40.tar.gz >> $logfile 2>&1  
fi

# Build the Toolkit
printf "\nFirst Fsi library compilation\n"
printf "\nFirst Fsi library compilation\n" >> $logfile
cd $WM_PROJECT_USER_DIR/FluidSolidInteraction/src
./Allwmake >> $logfile 2>&1

# Corrections to the Fuid Structure Interaction library
# Foam Extend 4.0 Updated the fluxRequired methods which need adding
# to the FSI Library
printf "\nCorrections to the FSI Library\n"
printf "\nCorrections to the FSI Library\n" >> $logfile

cd $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/fluidSolvers/

# If setFluxRequired string exists within file ignore
if grep -Fq "setFluxRequired" icoFluid/icoFluid.C; then
    printf "\nNo corrections required\n"
    printf "\nNo corrections required\n"  >> $logfile
else
    printf "\nAdding mesh.schemesDict().setFluxRequired(p_.name()) to fluid solvers\n"
    sed -i '/setRefCell/a \ \ \ \ mesh.schemesDict().setFluxRequired(p_.name());' icoFluid/icoFluid.C
    sed -i '/setRefCell/a \ \ \ \ mesh.schemesDict().setFluxRequired(p().name());' consistentIcoFluid/consistentIcoFluid.C
    sed -i '/setRefCell/a \ \ \ \ mesh.schemesDict().setFluxRequired(p().name());' pisoFluid/pisoFluid.C
fi

#######################################
# OpenFPCI compilation and installation
#######################################
printf "\n\n--------------------\n"
printf "OpenFPCI compilation"
printf "\n--------------------\n\n"
printf "\n\n--------------------\n" >> $logfile
printf "OpenFPCI compilation"             >> $logfile
printf "\n--------------------\n\n" >> $logfile

printf "Creating a soft link from:\n"
printf "OpenFPCI/src/solidSolvers/paraFEM to FluidSolidInteraction/src/fluidSolidInteraction/solidSolvers\n"

# Create a soft link to the solidSolvers
cd $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/solidSolvers
ln -s $WORKING_DIR/solidSolvers/paraFEM paraFEM


OPENMPI_COMPILE_FLAGS="`mpif90 --showme:compile`"
OPENMPI_LINK_FLAGS="`mpif90 --showme:link`"

if grep -Fq "femLargeStrain.C" $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/Make/files; then
    printf ""
else
    printf "\nModifying fluidSolidInteraction/Make/files"
    printf "\nModifying fluidSolidInteraction/Make/files" >> $logfile
    
    printf "solidSolvers/paraFEM/smallStrain/femSmallStrain.C\n" > paraFEM.files
    printf "solidSolvers/paraFEM/smallStrain/femSmallStrainSolve.C\n" >> paraFEM.files
    printf "solidSolvers/paraFEM/largeStrain/femLargeStrain.C\n" >> paraFEM.files
    printf "solidSolvers/paraFEM/largeStrain/femLargeStrainSolve.C\n" >> paraFEM.files
    printf "" >> paraFEM.files

    cat paraFEM.files $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/Make/files > tmp.files
    mv tmp.files $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/fluidSolidInteraction/Make/files
fi

OPTIONS="$WM_PROJECT_USER_DIR/FluidSolidInteraction/src/solvers/fsiFoam/Make/options"
if grep -Fq "parafemutils.o" $OPTIONS; then
    printf ""
else
    printf "\nModifying fsiFoam/Make/options"
    printf "\nModifying fsiFoam/Make/options" >> $logfile
    
    STRING="$OPENMPI_LINK_FLAGS \\"
    STRING=$(printf "%*s%s" 3 '' "$STRING")
    sed -i "/EXE_LIBS = /a \ ${STRING}\\" $OPTIONS

    STRING="-lgfortran \\"
    STRING=$(printf "%*s%s" 3 '' "$STRING")
    sed -i "/EXE_LIBS = /a \ ${STRING}\\" $OPTIONS
    
    STRING="-L$PARAFEM_DIR/lib -larpack_linuxdesktop \\"
    STRING=$(printf "%*s%s" 3 '' "$STRING")
    sed -i "/EXE_LIBS = /a \ ${STRING}\\" $OPTIONS
    
    STRING="-L$PARAFEM_DIR/lib -lParaFEM_mpi.5.0.3 \\"
    STRING=$(printf "%*s%s" 3 '' "$STRING")
    sed -i "/EXE_LIBS = /a \ ${STRING}\\" $OPTIONS
    
    STRING="${OPENFPCI_DIR}src/solidSolvers/paraFEM/fem_routines/objectFiles/parafemnl.o \\"
    STRING=$(printf "%*s%s" 3 '' "$STRING")
    sed -i "/EXE_LIBS = /a \ ${STRING}\\" $OPTIONS
    
    STRING="${OPENFPCI_DIR}src/solidSolvers/paraFEM/fem_routines/objectFiles/parafeml.o \\"
    STRING=$(printf "%*s%s" 3 '' "$STRING")
    sed -i "/EXE_LIBS = /a \ ${STRING}\\" $OPTIONS
    
    STRING="${OPENFPCI_DIR}src/solidSolvers/paraFEM/fem_routines/objectFiles/parafemutils.o \\"
    STRING=$(printf "%*s%s" 3 '' "$STRING")
    sed -i "/EXE_LIBS = /a \ ${STRING}\\" $OPTIONS
    
fi


printf "\nCompiling OpenFPCI Fortran files"
printf "\nCompiling OpenFPCI Fortran files" >> $logfile

# Compile the Fortran Subroutines
cd $OPENFPCI_DIR/src/solidSolvers/paraFEM/fem_routines
printf "\ngfortran -fcheck=all -c parafeml.f90 -o parafeml.o -I${PARAFEM_DIR}/include/mpi"
printf "\ngfortran -fcheck=all -c parafeml.f90 -o parafeml.o -I${PARAFEM_DIR}/include/mpi" >> $logfile
gfortran -fcheck=all -c parafeml.f90 -o parafeml.o -I${PARAFEM_DIR}/include/mpi >> $logfile 2>&1  

printf "\ngfortran -fcheck=all -c parafemnl.f90 -o parafemnl.o -I${PARAFEM_DIR}/include/mpi"
printf "\ngfortran -fcheck=all -c parafemnl.f90 -o parafemnl.o -I${PARAFEM_DIR}/include/mpi" >> $logfile
gfortran -fcheck=all -c parafemnl.f90 -o parafemnl.o -I${PARAFEM_DIR}/include/mpi >> $logfile 2>&1

printf "\ngfortran -fcheck=all -c parafemutils.f90 -o parafemutils.o -I${PARAFEM_DIR}/include/mpi"
printf "\ngfortran -fcheck=all -c parafemutils.f90 -o parafemutils.o -I${PARAFEM_DIR}/include/mpi" >> $logfile
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


printf "\nRecompiling FSI Library"
printf "\nRecompiling FSI Library\n\n" >> $logfile

# Remake the Library including new solidSolvers ans
cd $WM_PROJECT_USER_DIR/FluidSolidInteraction/src/
./Allwmake >> $logfile 2>&1

cd $WORKING_DIR

printf "\n\n--------------------------\n"
printf "\nEnd of OpenFPCI compilation"
printf "\n--------------------------\n\n"
printf "\n\n--------------------------\n"  >> $logfile
printf "End of OpenFPCI compilation"       >> $logfile
printf "\n--------------------------\n\n"  >> $logfile

##################
# OpenFPCI Testing 
##################

