// Open dictionary for reading

IOdictionary thermoPropertiesDict
(
    IOobject
    (
        "thermophysicalProperties",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

// LUTs for ek and QRad
fileName ekDataFile( thermoPropertiesDict.subDict("mixture").subDict("transport").lookup("ekLookupTable") );
interpolationTable<scalar> ekLookup(ekDataFile);

fileName QRadDataFile( thermoPropertiesDict.subDict("mixture").subDict("transport").lookup("QRadLookupTable") );
interpolationTable<scalar> QRadLookup(QRadDataFile);

// Read EM constants

// MVP decay constant, must be specified in thermophysicalProperties for 2D problems
dimensionedScalar constantMVP
(
    dimensionedScalar::lookupOrDefault
    (
        "constantMVP",
        thermoPropertiesDict.subDict("EM"),
        dimensionSet(0,-2,0,0,0,0,0),
        scalar(0.)
    )
);

// Permeability of free space
dimensionedScalar mu0
(
    thermoPropertiesDict.subDict("EM").lookup("mu0")
);

// Electron energy transport constant
dimensionedScalar constET
(
    thermoPropertiesDict.subDict("EM").lookup("constET")
);

// Set up EM fields

Info<< "Setting up field QRad\n" << endl;
volScalarField QRad
(
    IOobject
    (
        "QRad",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("QRad", dimensionSet(1,-1,-3,0,0,0,0), 1.0)
);

Info<< "Setting up field QRadSp\n" << endl;
volScalarField QRadSp
(
    IOobject
    (
        "QRadSp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("QRadSp", dimensionSet(1,-1,-3,-1,0,0,0), 1.0)
);

Info<< "Reading field ek\n" << endl;
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

Info<< "Setting up field ekSp\n" << endl;
volScalarField ekSp
(
    IOobject
    (
        "ekSp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("ekSp", dimensionSet(-1,-3,3,-1,0,2,0), 1000.)
);

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

Info<< "Setting up field B\n" << endl;
volVectorField B
(
    IOobject
    (
        "B",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
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
