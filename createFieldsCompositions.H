IOdictionary compDict
(
    IOobject
    (
        "dataTables/dataCompositions.foam",
        runTime.constant(),
        mesh,
        IOobject::MUST_READ_IF_MODIFIED,
        IOobject::NO_WRITE
    )
);

List<word> plasmaSpecies(compDict.lookup("species"));
List<scalar> compStartT(plasmaSpecies.size());
List<scalar> compDeltaT(plasmaSpecies.size());
List< List<scalar> > compData(plasmaSpecies.size());
List< autoPtr<volScalarField> > compFields(plasmaSpecies.size());

forAll(plasmaSpecies, i)
{
    word specI = plasmaSpecies[i];
    
    compStartT[i] = readScalar(compDict.subDict(specI).lookup("startT"));
    compDeltaT[i] = readScalar(compDict.subDict(specI).lookup("deltaT"));

    List<scalar> dataTmp(compDict.subDict(specI).lookup("dataTable"));
    compData[i] = dataTmp;
    
    compFields[i].reset
    (
        new volScalarField
        (
            IOobject
            (
                specI,
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
