# OpenFPCI

Open source FOAM to ParaFEM Coupling Interface:

## Motivation

The package allows the coupling of a Finite Element Package ParaFEM to a Computational Fluid Dynamics package Foam-extend. The motivation is to provide a framework for which users can model Fluid Structure Interaction problems on High Performance Computing Systems.

## Description

This package requires Foam-Extend, ParaFEM and an external Fluid Structure Interaction library. These can be installed through the following websites.

ParaFEM: https://sourceforge.net/projects/parafem/

Foam-Extend: https://openfoamwiki.net/index.php/Installation/Linux/foam-extend-4.0

FSI Library: https://openfoamwiki.net/index.php/Extend-bazaar/Toolkits/Fluid-structure\_interaction

The library represents a direct coupling between ParaFEM and Foam-Extend, through C++ wrappers around fortran parafem subroutines.

## Installation

Foam-Extend and ParaFEM must be downloaded and compiled prior to the installation of the OpenFPCI link. It is suggested that users follow the links provided above to install the relevent pacakges suitable for their systems. Once these packages have been compiled and installed the OpenFPCI can be downloaded using the following command: 
```
git clone https://github.com/SPHewitt/OpenFPCI
```

The third-party FSI library and OpenFPCI can be downloaded and compiled using the openfpci.sh script present in the src directory. The script assumes that both Foam-Extend and ParaFEM are installed with the system OpenMPI. The following command can be followed to install the software, where the paths should be replaced:

```
echo "export PARAFEM_DIR=path/to/parafem-code/parafem" >> ~/.bashrc
echo "export FOAM_DIR=path/to/foam/foam-extend.x.x" >> ~/.bashrc
. ~/.bashrc
cd OpenFPCI/src
./openfpci.sh
```

The software has been tested with Foam-Extend-4.0 and ParaFEM.5.0.3 on a linux desktop running Ubuntu 16.04. The code has further been tested to The Computational Shared Facility at Manchester, the SGI system in Leeds (Polaris), The XC30 Cray system in Ediburugh (Archer) and the Tianhe2 Machine in Guangzhou China.

The modules used on each of the systens is summarised below: 

### Linux Desktop

On a linux desktop the code uses openmpi-1.6.5 if another mpi package exists i.e. mpich the script will fail.  

Prerequisites:

* OpenMPI: Version 1.6.5 .
* GCC compilers: gcc 5.x .
* Ubuntu 16.04

### Manchester Computational Facility

The application has been tested with gnu/4.9.1 package and system openmpi version 1.6.5. The default compilers on the system are intel so these need to be swapped for the gnu compilers. 

```
module load mpi/gcc/openmpi/1.8-ib

module load compilers/gcc/4.9.0
```

### SGI - N8 Polaris (Leeds)

The application has been tested with gnu/4.9.1 package and system openmpi version 1.6.5. The default compilers on the system are intel so these need to be swapped for the gnu compilers.

```
module swap intel/"version" gnu/4.9.1

module load openmpi/1.6.5

```

### XC30 -  Archer (Edinburgh)

```
module swap PrgEnv-cray PrgEnv-gnu/5.1.29

module swap gcc/6.3.0 gcc/5.3.0

```

## Contact

For any assistance please contact sam.hewitt@manchester.ac.uk
