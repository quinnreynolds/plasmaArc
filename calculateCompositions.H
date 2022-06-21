if (calculateCompositionFieldsYN && runTime.writeTime())
{
    forAll(plasmaSpecies, compI)
    {
        // Internal field update for ek and ekSp
        forAll(T,cellI)
        {
            compFields[compI]->ref()[cellI] = scalarLookup
            (
                T.ref()[cellI], 
                compStartT[compI], 
                compDeltaT[compI], 
                compData[compI]
            );
        }

        forAll(T.boundaryField(),patchI) 
        {
            forAll(T.boundaryField()[patchI],cellI)
            {
                compFields[compI]->boundaryFieldRef()[patchI][cellI] = scalarLookup
                ( 
                    T.boundaryField()[patchI][cellI],
                    compStartT[compI], 
                    compDeltaT[compI], 
                    compData[compI]
                );
            }
        }
    }
}
