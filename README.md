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

See also READMEs in each case subdirectory for more detailed information.

* `IndustrialAC-1`: A large high-current arc in an industrial silicon smelter, using generic time-dependent boundary conditions.
* `IndustrialAC-2`: A large high-current arc in an industrial silicon smelter, using paired AC boundary conditions.
* `PilotDC-1`: A low-current arc in a pilot-scale furnace, operating in a CO atmosphere.
* `PilotDC-2`: A low-current arc in a pilot-scale furnace, operating in an air atmosphere.

### Who do I talk to? ###

* quinnr@mintek.co.za

### References

* Q.G. Reynolds (2021). Toward computational models of arc dynamics in silicon smelters, in *Proceedings of the 14th International Conference on CFD in the Oil & Gas, Metallurgical and Process Industries (CFD2020)*. Trondheim, Norway: SINTEF Academic Press, 2021, pp. 99-106 [https://www.pyrometallurgy.co.za/Mintek/Files/2021Reynolds-CFD2020.pdf]
* H.V. Haraldsson, Q.G. Reynolds, Y.A. Tesfahunegn, and G.A. Saevarsdottir (2025). Modeling of Industrial Electric Arcs Using Different Plasma Gas Compositions, *Metallurgical and Materials Transactions B*, 2025 [https://doi.org/10/1007/s11663-024-03397-4]
