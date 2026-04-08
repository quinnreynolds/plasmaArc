IOdictionary electromagneticProperties
(
    IOobject
    (
        "electromagneticProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

//magnetic permeability - defaults to free space
dimensionedScalar mu0
(
    "mu0",
    dimensionSet(1, 1, -2, 0, 0, -2, 0),
    electromagneticProperties.lookupOrDefault<scalar>("mu0", 1.25664e-06)
);

// Electron energy transport constant - leave at default unless some weird
// physics is happening
dimensionedScalar constantET
(
    "constantET",
    dimensionSet(1, 2, -3, -1, 0, -1, 0),
    electromagneticProperties.lookupOrDefault<scalar>("constantET", 0.000215433)
);

// MVP decay constant - used in 2D cartesian problems only, should be left at
// default 0 in all other cases
dimensionedScalar constantMVP
(
    "constantMVP",
    dimensionSet(0, -2, 0, 0, 0, 0, 0),
    electromagneticProperties.lookupOrDefault<scalar>("constantMVP", scalar(0.))
);

//electric potential field
Info<< "Reading field ePot" << nl;
volScalarField ePot
(
    IOobject
    (
        "ePot",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//magnetic vector potential field
Info<< "Reading field A" << nl;
volVectorField A
(
    IOobject
    (
        "A",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Setting up field ek" << nl;
volScalarField ek
(
    IOobject
    (
        "ek",
        runTime.name(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//magnetic field
Info<< "Setting up field B" << nl;
volVectorField B
(
    IOobject
    (
        "B",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::curl(A)
);

Info<< "Setting up field j" << nl;
volVectorField j
(
    IOobject
    (
        "j",
        runTime.name(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -ek * fvc::grad(ePot)
);

Info<< "Reading field ddtA" << nl;
volVectorField ddtA
(
    IOobject
    (
        "ddtA",
        runTime.name(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("ddtA", dimensionSet(1,1,-3,0,0,-1,0), vector(0.,0.,0.))
);

scalar dtEM = 0.;
