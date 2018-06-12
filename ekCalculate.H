//internal field update
forAll(ek(),cellI)
{
    ek.ref()[cellI] = scalarLookup(T.ref()[cellI], ekStartT, ekDeltaT, ekData);
    ekSp.ref()[cellI] =
        scalarLookup(T.ref()[cellI], ekStartT, ekDeltaT, ekSpData);
}

ek.correctBoundaryConditions();

//boundary field update
forAll(ek.boundaryField(),patchI) {
    forAll(ek.boundaryFieldRef()[patchI],cellI)
    {
        ekSp.boundaryFieldRef()[patchI][cellI] =
            scalarLookup(T.boundaryFieldRef()[patchI][cellI], ekStartT, ekDeltaT, ekSpData);
    }
}
