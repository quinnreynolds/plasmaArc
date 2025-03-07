## fixedLocationCurrentDensity

This boundary condition provides a fixedLocationCurrentDensity condition for electric potential fields. Total current and current density at the emission faces are specified as functions of time in the form of lookup tables. Faces on the boundary equal in total area to current / currentDensity are selected based on their proximity to a specified point in space.

At the attachment spot faces the BC is calculated as follows, with $\sigma$ the electrical conductivity, $\phi$ the electric potential field, and $j_C$ the emission current density at the conducting surface.

$-\sigma \frac{\partial \phi}{\partial \mathbf{n}} = j_C$

For all other faces, a zero-gradient condition is applied.

#### Usage

- `current` is a lookup table specifying the instantaneous total current (A) as a function of time.
- `currentDensity` is a lookup table specifying $j_C$ (A/m2), the current density in the arc attachment spot.
- `referencePosition` is the point in space used to fix the location of the attachment spot.

Example of the boundary condition specification:

`patchName`  
`{`  
`    type                fixedLocationCurrentDensity;`  
`    current             table`   
`                        (`  
`                            (0 1000)`  
`                            (1 1000)`  
`                        );`  
`    currentDensity      table`   
`                        (`  
`                            (0 2e7)`  
`                            (1 2e7)`  
`                        );`  
`    referencePosition   (0 0 0);`  
`    value               uniform 0;`  
`}`  
