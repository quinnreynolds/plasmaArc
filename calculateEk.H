if (runTime.timeIndex() == 1 || runTime.timeIndex() % emSolverFreq == 0)
{
    // Internal field update for ek and ekSp
    forAll(ek(),cellI)
    {
        ek.ref()[cellI] = scalarLookup
        (
            T.ref()[cellI], ekStartT, ekDeltaT, ekData
        );
        ekSp.ref()[cellI] =
            scalarLookup(T.ref()[cellI], ekStartT, ekDeltaT, ekSpData);
    }

    // Boundary field for ek is updated in emEqns.H

    // Boundary field update for ekSp 
    forAll(ekSp.boundaryField(),patchI) {
        forAll(ekSp.boundaryFieldRef()[patchI],cellI)
        {
            ekSp.boundaryFieldRef()[patchI][cellI] = scalarLookup
            ( 
                T.boundaryFieldRef()[patchI][cellI],
                ekStartT, ekDeltaT, ekSpData
            );
        }
    }    
}
