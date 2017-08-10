# OpenFPCI

OpenFOAM-ParaFEM Coupling Interface:

## Motivation

The package allows the coupling of a Finite Element Package ParaFEM to a Computational Fluid Dynamics package OpenFOAM-dev. The motivation is to provide a framework for which users can model Fluid Structure Interaction problems on High Performance Computing Systems.

## Installation

The installation of this package requires Foam-Extend, ParaFEM and an external Fluid Structure Interaction library. These can all be installed using the scripts provided or manually through the following websites.

ParaFEM: https://sourceforge.net/projects/parafem/

Foam-Extend: https://openfoamwiki.net/index.php/Installation/Linux/foam-extend-4.0

FSI Library: https://openfoamwiki.net/index.php/Extend-bazaar/Toolkits/Fluid-structure_interaction

As of the 9th August 2017: The Foam Extend and FSI libraries have only been tested upto GCC 5.0.1

### Downloading the Code

```
git clone https://github.com/SPHewitt/OpenFPCI
```

This repository contatins the files that couple OpenFOAM to ParaFEM. The scripts avialable in the build directory can be used to install Foam-Extend ParaFEM and OpenFPCI.
