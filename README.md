# README #

Coupled MHD model for studying plasma arc behaviour. 

### What is this repository for? ###

* Dynamics of incompressible and weakly-compressible (low current) thermal plasma arcs used in welding and arc furnace applications.
* master compiles against OF7.

### How do I get set up? ###

To build plasmaArc:

* Clone
* Init and pull eminclude submodule
* wmake

Additional requirements:

* Thermophysical properties - compressibleFluidLut thermo module (note that incompressibleFluidLut will also work but requires older versions of plasmaArc which uses incompressible formulations of the pressure equation time derivatives).
* Radiation - greyPlasmaAbsorptionEmission module.
* Boundary conditions - limitedTemperature, fixedCurrentDensity, fixedCurrentCurrentDensity.

### Contribution guidelines ###

* TBD

### Who do I talk to? ###

* quinnr@mintek.co.za
