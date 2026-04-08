if (runTime.timeIndex() == 1 || runTime.timeIndex() % emSolverFreq == 0)
{
    // Internal field update for ek and ekSp
    forAll(ek(),cellI)
    {
        ek.primitiveFieldRef()[cellI] = scalarLookup
        (
            T.primitiveField()[cellI], ekStartT, ekDeltaT, ekData
        );
        ekSp.primitiveFieldRef()[cellI] =
            scalarLookup(T.primitiveField()[cellI], ekStartT, ekDeltaT, ekSpData);
    }

    // Boundary field for ek is updated in emEqns.H

    // Boundary field update for ekSp
    forAll(ekSp.boundaryField(),patchI) {
        forAll(ekSp.boundaryFieldRef()[patchI],cellI)
        {
            ekSp.boundaryFieldRef()[patchI][cellI] = scalarLookup
            (
                T.boundaryField()[patchI][cellI],
                ekStartT, ekDeltaT, ekSpData
            );
        }
    }
}
