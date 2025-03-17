### Description ###

3D transient simulation of a short, high-current AC arc in the space between the electrode tip 
and the pool in a silicon smelting furnace. The simulation is initialised at stagnant flow and 
a constant plasma temperature of 10,000 K. Inlet-outlet boundaries are specified at the 
perimeter of the region. A fixed-potential boundary is used at the pool (anode) surface, and a 
time-varying current boundary is applied at the electrode (cathode) surface. This produces an 
arc jet which travels in the same direction during each AC half-cycle. Two full AC cycles at 
50 Hz are simulated (0.04 s).

### Execution ###

* `./clean`
* `./prepare`
* `./runcase`

### Purpose ###

* Demonstrate `plasmaArc` solver on industrial-scale problems.
* Demonstrate use of `fixedLocationCurrentDensity` boundary condition.

### Comments ###

Execution of a complete simulation run can take several hours even when running in parallel, 
although output can be monitored continuously as it is produced. It is recommended to run this 
case on a modelling workstation or HPC facility where it will not interfere with other work, 
and adjust the number of processors in `decomposeParDict` as appropriate.

The plasma property dataset for an equimolar mixture of SiO and CO was calculated using 
(minplascalc)[https://github.com/quinnreynolds/minplascalc] v0.7.0.

This example was modified from a case originally developed by Hákon Valur Haraldsson in 2024 
during Ph. D. studies at Reykjavík University.
