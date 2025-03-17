### Description ###

3D transient simulation of a low-current DC arc in a pilot-scale direct current plasma arc 
furnace, operating with an air atmosphere. The simulation is initialised at stagnant flow 
and a constant plasma temperature of 3000 K. Inlet-outlet boundaries are specified at the 
top of the region, and solid walls at the perimeter. A ground boundary condition ($V=0$) 
is applied at the anode surface, and a fixed current density is imposed on the cathode spot 
boundary (note that when `fixedCurrentDensity` is used, the user must size the conducting 
boundary appropriately to generate the required total current). 0.025 s of arc development 
and dynamic behaviour is simulated.

### Execution ###

* `./clean`
* `./prepare`
* `./runcase`

### Purpose ###

* Demonstrate `plasmaArc` solver on pilot-scale problems and more complex geometries.
* Demonstrate use of `fixedCurrentDensity` boundary condition.

### Comments ###

Execution of a complete simulation run can take several hours even when running in parallel, 
although output can be monitored continuously as it is produced. It is recommended to run this 
case on a modelling workstation or HPC facility where it will not interfere with other work, 
and adjust the number of processors in `decomposeParDict` as appropriate.

The plasma property dataset for a simplified air mixture was calculated using 
[minplascalc](https://github.com/quinnreynolds/minplascalc) v0.7.0.

