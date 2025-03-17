# README #

Unified thermophysical model for OpenFOAM simulations, using lookup tables to specify rho(T), Cp(T), h(T), tk(T), and mu(T) for weakly compressible fluids. rho(T) is defined at a reference pressure, typically atmospheric, and ideal gas behaviour is then assumed for the pressure dependence of the properties at a given temperature. 

See example cases for how to use this library in conjuction with plasma arc solvers.
