Info<< "Reading thermophysical properties\n" << endl;

autoPtr<fluidThermo> pThermo
(
    fluidThermo::New(mesh)
);
fluidThermo& thermo = pThermo();
thermo.validate(args.executable(), "h", "e");

volScalarField& T = thermo.T();
volScalarField& p = thermo.p();
const volScalarField& psi = thermo.psi();

Info<< "Creating field references" << nl;
Info<< "    T field max: " << max(T).value() << nl;
Info<< "    p field max: " << max(p).value() << nl;
Info<< "    psi field max: " << max(psi).value() << nl << nl;

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

pressureControl pressureControl(p, rho, pimple.dict(), false);

mesh.setFluxRequired(p.name());

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

Info<< "Creating field MachNo\n" << endl;
volScalarField MachNo
(
    IOobject
    (
        "MachNo",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
),
    mesh,
    dimensionedScalar("MachNo", dimensionSet(0,0,0,0,0,0,0), 0)
);

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

Info<< "Creating LUT for ek" << nl;

scalar ekStartT(readScalar(thermoPropertiesDict.subDict("mixture")
    .subDict("transport").subDict("ekLookupTable").lookup("startT")));
scalar ekDeltaT(readScalar(thermoPropertiesDict.subDict("mixture")
    .subDict("transport").subDict("ekLookupTable").lookup("deltaT")));
List<scalar> ekData(thermoPropertiesDict.subDict("mixture").subDict("transport")
    .subDict("ekLookupTable").lookup("dataTable"));
List<scalar> ekSpData(ekData);
ekSpData[0] = (ekData[1] - ekData[0]) / ekDeltaT;
ekSpData[ekSpData.size()-1] = (ekData[ekSpData.size()-1]
                             - ekData[ekSpData.size()-2]) / ekDeltaT;
for(label dataI = 1; dataI < ekSpData.size()-1; dataI++)
{
    ekSpData[dataI] = (ekData[dataI+1] - ekData[dataI-1]) / (2. * ekDeltaT);
}
Info<< "    ek temperature range: " << ekStartT << ", "
    << ekStartT + ekDeltaT * (ekData.size() - 1) << nl << nl;

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

dimensionedScalar initialMass = fvc::domainIntegrate(rho);

#include "createDpdt.H"

#include "createK.H"

#include "createMRF.H"

const dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT, pimple.dict());
const dimensionedScalar rhoMin("rhoMin", dimDensity, Zero, pimple.dict());

#include "createFvOptions.H"
