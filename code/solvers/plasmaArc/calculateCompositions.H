if (calculateCompositionFieldsYN)
{
    // Update electron species at every time step - needed for CASM boundary conditions.
    // Info<< plasmaSpecies << nl;

    forAll(plasmaSpecies, compI)
    {
        if (plasmaSpecies[compI] == "e")
        {
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

    // Update all species only at write timesteps.
    if (runTime.writeTime())
    {
        forAll(plasmaSpecies, compI)
        {
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
}
