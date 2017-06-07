Info<< "Reading thermophysical properties\n" << endl;

autoPtr<rhoReactionThermo> rrThermo
(
    rhoReactionThermo::New(mesh)
);

rhoReactionThermo& thermo = rrThermo();
thermo.validate(args.executable(), "h", "e");

basicMultiComponentMixture& composition = thermo.composition();
const speciesTable& plasmaSpeciesTable = composition.species();
PtrList<volScalarField>& Y = composition.Y();
volScalarField& he = thermo.he();
volScalarField& T = thermo.T();
volScalarField& p = thermo.p();
const volScalarField& psi = thermo.psi();

forAll(plasmaSpeciesTable, i)
{
    Info<< "Component " << i << ": " << plasmaSpeciesTable[i] << nl;
}

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
// NB assumes all LUTs are same length and keyed the same

scalar ekStartT(readScalar(thermoPropertiesDict.subDict(plasmaSpeciesTable[0]).subDict("transport").subDict("ekLookupTable").lookup("startT")));
scalar ekDeltaT(readScalar(thermoPropertiesDict.subDict(plasmaSpeciesTable[0]).subDict("transport").subDict("ekLookupTable").lookup("deltaT")));
List< List<scalar> > ekDataList;
List< List<scalar> > ekSpDataList;
scalar QRadStartT(readScalar(thermoPropertiesDict.subDict(plasmaSpeciesTable[0]).subDict("transport").subDict("QRadLookupTable").lookup("startT")));
scalar QRadDeltaT(readScalar(thermoPropertiesDict.subDict(plasmaSpeciesTable[0]).subDict("transport").subDict("QRadLookupTable").lookup("deltaT")));
List< List<scalar> > QRadDataList;
List< List<scalar> > QRadSpDataList;

forAll(plasmaSpeciesTable, speciesI)
{
    List<scalar> myQRadData(thermoPropertiesDict.subDict(plasmaSpeciesTable[speciesI]).subDict("transport").subDict("QRadLookupTable").lookup("dataTable"));
    List<scalar> myQRadSpData(myQRadData);
    myQRadSpData[0] = (myQRadData[1] - myQRadData[0]) / QRadDeltaT;
    myQRadSpData[myQRadSpData.size()-1] = (myQRadData[myQRadSpData.size()-1] - myQRadData[myQRadSpData.size()-2]) / QRadDeltaT;
    for(label dataI = 1; dataI < myQRadSpData.size()-1; dataI++)
    {
        myQRadSpData[dataI] = (myQRadData[dataI+1] - myQRadData[dataI-1]) / (2. * QRadDeltaT);
    }
    QRadDataList.append(myQRadData);
    QRadSpDataList.append(myQRadSpData);

    List<scalar> myEkData(thermoPropertiesDict.subDict(plasmaSpeciesTable[speciesI]).subDict("transport").subDict("ekLookupTable").lookup("dataTable"));
    List<scalar> myEkSpData(myEkData);
    myEkSpData[0] = (myEkData[1] - myEkData[0]) / ekDeltaT;
    myEkSpData[myEkSpData.size()-1] = (myEkData[myEkSpData.size()-1] - myEkData[myEkSpData.size()-2]) / ekDeltaT;
    for(label dataI = 1; dataI < myEkSpData.size()-1; dataI++)
    {
        myEkSpData[dataI] = (myEkData[dataI+1] - myEkData[dataI-1]) / (2. * ekDeltaT);
    }
    ekDataList.append(myEkData);
    ekSpDataList.append(myEkSpData);
}


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
