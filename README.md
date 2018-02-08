# OpenFPCI

OpenFOAM-ParaFEM Coupling Interface:

## Motivation

The package allows the coupling of a Finite Element Package ParaFEM to a Computational Fluid Dynamics package Foam-extend. The motivation is to provide a framework for which users can model Fluid Structure Interaction problems on High Performance Computing Systems.

## Description

This package requires Foam-Extend, ParaFEM and an external Fluid Structure Interaction library. These can be installed through the following websites.

ParaFEM: https://sourceforge.net/projects/parafem/

Foam-Extend: https://openfoamwiki.net/index.php/Installation/Linux/foam-extend-4.0

FSI Library: https://openfoamwiki.net/index.php/Extend-bazaar/Toolkits/Fluid-structure\_interaction

In the src folder there are a series of scripts to install the pacakges above, however it is strongly advised to install both Foam-Extend and ParaFEM by hand. The openfpci.sh script can then be used to install, compile the FSI library along with the addittions required for the linking.

### Downloading the Code

```
git clone https://github.com/SPHewitt/OpenFPCI
```

This repository contatins the files that couple Foam-extend to ParaFEM. The scripts available in the src directory can be used to install Foam-Extend ParaFEM and OpenFPCI.

The software installation has been tested on a number of systems, a linux desktop (Ubuntu 16.04), The Computational Shared Facility at Manchester and the SGI system in Leeds (Polaris). The instructions for these are as follows. The scripts use the system compilers and mpi packages, so correct modules and packages need to be loaded and installed. 

### Compiling the code

The FSI Library does not currently compile with gcc >= 6 so be careful when loading the gnu modules.


### Linux Desktop

Prerequisites:

* OpenMPI: Version 1.6.5  has been tested.
* GCC compilers: gcc 4.x  have been tested.

```
cd OpenFPCI

./compile.sh
```

### Manchester Computational Facility

Instructuctions to be updated soon.

### SGI - N8 Polaris (Leeds)

The application has been tested with gnu/4.9.1 package and system openmpi version 1.6.5. The default compilers on the system are intel so these need to be swapped for the gnu compilers. The intel compilers can be used however instructions for this are not yet available.

```
module swap intel/"version" gnu/4.9.1

module load openmpi/1.6.5

cd OpenFPCI

./compile.sh
```

### XC30 -  Archer (Edinburgh)

Instructions to be updated soon.
