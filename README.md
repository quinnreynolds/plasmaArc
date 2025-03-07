# README #

A coupled magnetohydrodynamics solver for OpenFOAM, for studying the dynamic behaviour of weakly compressible thermal plasma arcs used in welding and metallurgical arc furnace applications. Live branches track ESI OpenFOAM releases unless otherwise indicated.

Copyright (C) Quinn Reynolds and Mintek 2025-present

### Installation ###

To build the solvers and support utilities for `plasmaArc`, you will need a working OpenFOAM installation and activated environment for same on your machine. If this is in place, you should be ok with the following after cloning the repository and changing directory to `plasmaArc`:

* `cd code`
* `./Allwclean`
* `./Allwmake`

In addition to the solver itself this repository also contains the material properties and radiation libraries required for it to work, and a basic selection of solver-specific boundary conditions. Note that these will compile at the same time as the solvers when the `Allwclean`, `Allwmake` commands above are issued:

* Thermophysical properties - `compressibleFluidLut` thermo module.
* Radiation - `greyPlasmaAbsorptionEmission` module (P1 or other emission-aware models are recommended for correct behaviour).
* Boundary conditions - `limitedTemperature`, `fixedCurrentDensity`, `fixedLocationCurrentDensity`, and `fixedLocationAlternatingCurrent`.

### Example cases ###

* TBD

### Who do I talk to? ###

* quinnr@mintek.co.za
