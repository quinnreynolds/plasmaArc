# README #

Coupled MHD model for studying plasma arc behaviour. 

### What is this repository for? ###

* Dynamics of weakly compressible thermal plasma arcs used in welding and arc furnace applications.
* Live branches track ESI OpenFOAM releases unless otherwise indicated.

### How do I get set up? ###

To build plasmaArc:

* Clone
* wmake

Additional requirements:

* Thermophysical properties - compressibleFluidLut thermo module.
* Radiation - greyPlasmaAbsorptionEmission module (please use P1 or other emission-aware models for correct behaviour).
* Boundary conditions - limitedTemperature, fixedCurrentDensity, variableCurrentCurrentDensity, etc.

### Contribution guidelines ###

* TBD

### Who do I talk to? ###

* quinnr@mintek.co.za
