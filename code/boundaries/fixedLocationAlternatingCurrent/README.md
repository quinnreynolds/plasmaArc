## fixedLocationAlternatingCurrent

This boundary condition provides a sinusoidal variable current boundary condition for electric potential 
fields. The BC outputs a variable current for the positive part of the sine curve but then switches to a 
fixed value of zero during the negative part. The curve is specified by a given frequency (in Hz) and 
phase shift (in degrees) in the input directory.

This condition can for example be used to simulate an AC arc that reverses direction as the polarity 
between electrodes reverses periodically. This can be simulated by imposing opposite conditions of the 
`fixedLocationCurrentDensity` type on the two electrodes, shifted by 180°.

Faces on the boundary equal in total area to current / currentDensity are selected and assigned to the 
attachment spot at any given time based on the their proximity to the specified `referencePosition` point. 
This has the effect of fixing the arc root in the same place on the boundary surface during each emission 
period.

At the attachment spot faces the BC is calculated as follows, with $\sigma$ the electrical conductivity, 
$\phi$ the electric potential field, and $j_C$ the emission current density at the conducting surface.

$-\sigma \frac{\partial \phi}{\partial \mathbf{n}} = j_C$

For all other faces, a zero-gradient condition is applied.

This boundary condition was developed together with Hákon Valur Haraldsson in 2024 during Ph. D. 
studies at Reykjavík University.

#### Usage

- `maxCurrent` is the peak current or amplitude of the sine curve.
- `currentDensity` is $j_C$, the current density in the attachment spot of the electrode.
- `frequency` is the frequency of the curve, in Hz.
- `theta` is the phase shift of the curve, in degrees.
- `referencePosition` is the point in space used to fix the location of the attachment spot.

Example of the boundary condition specification:

`patchName`  
`{`  
`    type                fixedLocationAlternatingCurrent;`  
`    maxCurrent          5000;`  
`    currentDensity      2e7;`  
`    frequency           50;`  
`    theta               0;`  
`    referencePosition   (0 0 0);`  
`    value               uniform 0;`  
`}`  
