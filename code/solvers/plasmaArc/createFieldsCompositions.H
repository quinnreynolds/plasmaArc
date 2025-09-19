scalar calculateCompositionFieldsYN
(
    thermoPropertiesDict.lookupOrDefault("calculateCompositionFieldsYN", false)
);

List<word> plasmaSpecies;
List<scalar> compStartT;
List<scalar> compDeltaT;
List< List<scalar> > compData;
List< autoPtr<volScalarField> > compFields;

if (calculateCompositionFieldsYN)
{
    plasmaSpecies = thermoPropertiesDict.subDict("compositions").lookup("species");
    compStartT.resize(plasmaSpecies.size());
    compDeltaT.resize(plasmaSpecies.size());
    compData.resize(plasmaSpecies.size());
    compFields.resize(plasmaSpecies.size());

    forAll(plasmaSpecies, i)
    {
        word specI = plasmaSpecies[i];
        
        compStartT[i] = readScalar(thermoPropertiesDict.subDict("compositions").subDict(specI).lookup("startT"));
        compDeltaT[i] = readScalar(thermoPropertiesDict.subDict("compositions").subDict(specI).lookup("deltaT"));

        List<scalar> dataTmp(thermoPropertiesDict.subDict("compositions").subDict(specI).lookup("dataTable"));
        compData[i] = dataTmp;
        
        compFields[i].reset
        (
            new volScalarField
            (
                IOobject
                (
                    "X_species_" + specI,
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(specI, dimensionSet(0,0,0,0,0,0,0), 0.)
            )
        );
    }
}
