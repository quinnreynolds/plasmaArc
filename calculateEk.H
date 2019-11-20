if (runTime.timeIndex() == 1 || runTime.timeIndex() % emSolverFreq == 0)
{
    //internal field update
    forAll(ek(),cellI)
    {
        ek.ref()[cellI] = scalarLookup
        (
            T.ref()[cellI], ekStartT, ekDeltaT, ekData
        );
        ekSp.ref()[cellI] =
            scalarLookup(T.ref()[cellI], ekStartT, ekDeltaT, ekSpData);
    }
}