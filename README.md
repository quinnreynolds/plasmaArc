# README #

Coupled MHD model for studying plasma arc behaviour. 

### What is this repository for? ###

* Dynamics of incompressible and weakly-compressible (low current) thermal plasma arcs used in welding and arc furnace applications.
* master compiles against OF8.

### How do I get set up? ###

To build plasmaArc:

* Clone
* Init and update emInclude submodule
* wmake

Additional requirements:

* Thermophysical properties - compressibleFluidLut thermo module.
* Radiation - greyPlasmaAbsorptionEmission module (please use P1 or other emission-aware models for correct behaviour).
* Boundary conditions - limitedTemperature, fixedCurrentDensity, variableCurrentCurrentDensity.

### Contribution guidelines ###

* TBD

### Who do I talk to? ###

* quinnr@mintek.co.za
