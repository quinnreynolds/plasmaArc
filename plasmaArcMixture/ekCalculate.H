//internal field update

forAll(ek(),cellI)
{
    scalar ekTemp = 0.;
    scalar ekSpTemp = 0.;
    forAll(plasmaSpeciesTable, speciesI)
    {
        ekTemp += Y[speciesI].ref()[cellI]
            * scalarLookup(T.ref()[cellI], ekStartT, ekDeltaT,
                           ekDataList[speciesI]);
        ekSpTemp += Y[speciesI].ref()[cellI]
            * scalarLookup(T.ref()[cellI], ekStartT, ekDeltaT,
                           ekSpDataList[speciesI]);
    }
    ek.ref()[cellI] = ekTemp;
    ekSp.ref()[cellI] = ekSpTemp;
}

ek.correctBoundaryConditions();

//boundary field update
forAll(ek.boundaryField(),patchI) {
    forAll(ek.boundaryFieldRef()[patchI],cellI)
    {
        scalar ekSpTemp = 0.;
        forAll(plasmaSpeciesTable, speciesI)
        {
            ekSpTemp += Y[speciesI].boundaryFieldRef()[patchI][cellI]
                * scalarLookup(T.boundaryFieldRef()[patchI][cellI], ekStartT,
                               ekDeltaT, ekSpDataList[speciesI]);
        }
        ekSp.boundaryFieldRef()[patchI][cellI] = ekSpTemp;
    }
}
