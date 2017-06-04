Info<< "Reading thermophysical properties\n" << endl;

autoPtr<rhoThermo> rThermo
(
    rhoThermo::New(mesh)
);
rhoThermo& thermo = rThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField& he = thermo.he();
volScalarField& T = thermo.T();
volScalarField& p = thermo.p();
const volScalarField& psi = thermo.psi();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "compressibleCreatePhi.H"

dimensionedScalar rhoMax
(
    dimensionedScalar::lookupOrDefault
    (
        "rhoMax",
        pimple.dict(),
        dimDensity,
        GREAT
    )
);

dimensionedScalar rhoMin
(
    dimensionedScalar::lookupOrDefault
    (
        "rhoMin",
        pimple.dict(),
        dimDensity,
        0
    )
);

Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

Info<< "Creating field dpdt\n" << endl;
volScalarField dpdt
(
    IOobject
    (
        "dpdt",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("dpdt", p.dimensions()/dimTime, 0)
);

Info<< "Creating field kinetic energy K\n" << endl;
volScalarField K("K", 0.5*magSqr(U));


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
fileName ekDataFile
(
    thermoPropertiesDict.subDict("mixture").subDict("transport").lookup("ekLookupTable")
);
interpolationTable<scalar> ekLookup(ekDataFile);

fileName QRadDataFile
(
    thermoPropertiesDict.subDict("mixture").subDict("transport").lookup("QRadLookupTable")
);
interpolationTable<scalar> QRadLookup(QRadDataFile);

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
