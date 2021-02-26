# MMoCa
MMoCa is a Multi-purpose Monte Carlo simulation package for molecular systems. I developed this package for my PhD with Keith Gubbins at North Carolina State University (USA). This MC simulation package can perform multiple types of simulations (NVT, NPT and uVT). This program is capable of calculating many properties of the system, for example, local pressure tensor for different geometries (planar and cylindrical), density profile, isosteric heat, and thermodynamic properties that are measures of fluctuations (e.g., heat capacity). The program also leverages techniques such as cell list and coarse-grained surface potential to facilitate simulations. For more information, check the manual in ```/docs/``` folder. The current version is v2.2.

## Requirement and Installation 
The source codes are in ```/src/```. To compile, simply goto the ```/src/``` directory and hit

```make```

Compiler such as ```gfortran``` or ```ifort``` is required. 

## Getting Started
**/examples/**. This directory contains examples of short MMoCa simulations in NPT and grand canonical ensembles for a number of systems that I have studied during my PhD. <br/>
**/docs/**. This directory contains user manual for setting up the input files and logs of program development. 

## Reference
[1] K. Shi, E.E. Santiso, K.E. Gubbins, Can we define a unique microscopic pressure in inhomogeneous fluids?, J. Chem. Phys. 154 (2021) 084502. https://doi.org/10.1063/5.0044487.<br/>
[2] K. Shi, E.E. Santiso, K.E. Gubbins, Conformal Sites Theory for Adsorbed Films on Energetically Heterogeneous Surfaces, Langmuir. 36 (2020) 1822–1838. https://doi.org/10.1021/acs.langmuir.9b03633.<br/>
[3] K. Shi, E.E. Santiso, K.E. Gubbins, Bottom-Up Approach to the Coarse-Grained Surface Model: Effective Solid–Fluid Potentials for Adsorption on Heterogeneous Surfaces, Langmuir. 35 (2019) 5975–5986. https://doi.org/10.1021/acs.langmuir.9b00440.<br/>
