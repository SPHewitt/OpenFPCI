
This folder contains the CSM3, CFD3 and FSI3 cases from Turek and Hron [1]. These cases represent a set of benchmarks outlined in the paper, to verify and validate newly developed solid solver for the fluid structure interaction library. These three cases are outlined in the validation section of the OpenFPCI wiki.

The Allrun script runs the FSI3 case only. The CSM3 and CFD3 cases use different executables that are not compiled as part of the openfpci.sh installations script. If you want to run these cases the additions made to solvers/fsiFoam/Make/options have to be added to solvers/fluidFoam/Make/options and solvers/solidFoam/Make/options and recompiled.

[1] S.Turek and J.Hron (2006) Proposal for numerical benchmarking of fluid-structure interaction between an elastic object and laminar incompressible flow.
