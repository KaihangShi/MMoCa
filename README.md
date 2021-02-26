# MMoCa
MMoCa is a Multi-purpose Monte Carlo simulation package for molecular systems. I developed this package for my PhD with Keith Gubbins at North Carolina State University (USA). This MC simulation package can perform multiple types of simulations (NVT, NPT and uVT). This program is capable of calculating many properties of the system, for example, local pressure tensor for different geometries (planar and cylindrical), density profile, isosteric heat, and thermodynamic properties that are measures of fluctuations (e.g., heat capacity). The program also leverages techniques such as cell list and coarse-grained surface potential to facilitate simulations. For more information, check the manual in ```/docs/``` folder. The current version is v2.2.

## Requirement and Installation 
The source codes are in ```/src/```. To compile, simply goto the ```/src/``` directory and hit
```make```
Compilers such as ```gfortran``` or ```ifort``` are required. 

## Getting Started
**/examples/**. This directory contains examples of short MMoCa simulations in
NPT and grand canonical ensembles for a number of systems that I have presented in my paper. 
**/docs/**. This directory contains user manual for setting up the input files and logs of program development. 

