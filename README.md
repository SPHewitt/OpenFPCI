# OpenFPCI

OpenFOAM-ParaFEM Coupling Interface:

## Motivation

The package allows the coupling of a Finite Element Package ParaFEM to a Computational Fluid Dynamics package OpenFOAM-dev. The motivation is to provide a framework for which users can model Fluid Structure Interaction problems on High Performance Computing Systems.

## Installation

The installation of this package requires Foam-Extend, ParaFEM and an external Fluid Structure Interaction library. These can all be installed using the scripts provided or manually through the following websites.

ParaFEM: https://sourceforge.net/projects/parafem/

Foam-Extend: https://openfoamwiki.net/index.php/Installation/Linux/foam-extend-4.0

FSI Library: https://openfoamwiki.net/index.php/Extend-bazaar/Toolkits/Fluid-structure_interaction

### Downloading the Code

```
git clone https://github.com/SPHewitt/OpenFPCI
```

This repository contatins the files that couple OpenFOAM to ParaFEM. The scripts avialable in the build directory can be used to install Foam-Extend ParaFEM and OpenFPCI.

The software installation has been tested on a number of systems, a linux desktop (Ubuntu 16.04), The Computational Shared Facility at Manchester and the SGI system in Leeds (Polaris). The instructions for these are as follows. The Scripts use the sytem compilers and mpi packages, so correct modules and packages need to be loaded and installed. 

The FSI Library does not currently compile with gcc >= 6 so be careful when loading the gnu modules.


### Linux Desktop

Prerequisites:

* OpenMPI: Version 1.6.5 and 1.8.8 have been tested.
* GCC compilers: gcc 4.x to 5.x have been tested.


### Manchester Computational Facility

Prerequisites:

```
module load compilers/gcc/4.9.0

module load mpi/gcc/openmpi/1.6-ib
```

### SGI - N8 Polaris (Leeds)

Prerequisites:

```
module swap intel/"version" gnu/5.0.3

module load openmpi/1.6.5
```

### XC30 -  Archer (Edinburgh)

Instructions to be updated soon.
