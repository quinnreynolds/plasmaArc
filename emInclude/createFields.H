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
    dimensionedScalar::lookupOrDefault
    (
        "mu0",
        electromagneticProperties,
        dimensionSet(1, 1, -2, 0, 0, -2, 0),
        scalar(1.25664e-06)
    )
);

// Electron energy transport constant - leave at default unless some weird
// physics is happening
dimensionedScalar constantET
(
    dimensionedScalar::lookupOrDefault
    (
        "constantET",
        electromagneticProperties,
        dimensionSet(1, 2, -3, -1, 0, -1, 0),
        scalar(0.000215433)
    )
);

//MVP decay constant - 2D problems only, should be left at default 0 for 3D
dimensionedScalar constantMVP
(
    dimensionedScalar::lookupOrDefault
    (
        "constantMVP",
        electromagneticProperties,
        dimensionSet(0, -2, 0, 0, 0, 0, 0),
        scalar(0.)
    )
);

//electric potential field
Info<< "Reading field ePot\n" << endl;
volScalarField ePot
(
    IOobject
    (
        "ePot",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//magnetic vector potential field
Info<< "Reading field A\n" << endl;
volVectorField A
(
    IOobject
    (
        "A",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Setting up field ek\n" << endl;
volScalarField ek
(
    IOobject
    (
        "ek",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

//magnetic field
Info<< "Setting up field B\n" << endl;
volVectorField B
(
    IOobject
    (
        "B",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    fvc::curl(A)
);

Info<< "Setting up field j\n" << endl;
volVectorField j
(
    IOobject
    (
        "j",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    -ek * fvc::grad(ePot)
);

Info<< "Reading field ddtA\n" << endl;
volVectorField ddtA
(
    IOobject
    (
        "ddtA",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedVector("ddtA", dimensionSet(1,1,-3,0,0,-1,0), vector(0.,0.,0.))
);

scalar dtEM = 0.;
