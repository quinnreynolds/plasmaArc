/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2022 Quinn Reynolds, Mintek
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::fixedCurrentDensityFvPatchScalarField

Group
    grpGenericBoundaryConditions

escription
    This boundary condition provides a fixedCurrentDensity condition for
    electric potential fields, calculated as:
        \f[
            d/dx (ePot) = currentDensity / ek[i]
        \f]
    where
    \vartable
        currentDensity  | scalar field of current density on the patch, A/m2
        ek              | scalar field of electrical conductivity
        ePot            | scalar field of electric potential
    \endtable
Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedCurrentDensity;
        currentDensity  uniform 2e7;
        value           uniform 0.;
    }
    \endverbatim

SourceFiles
    fixedCurrentDensityFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedCurrentDensityFvPatchScalarField_H
#define fixedCurrentDensityFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class fixedCurrentDensityFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedCurrentDensityFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private Data

        //- Field of Types, typically defined across patch faces
        //  e.g. total pressure p0_.  Other options include vectorField
        scalarField currentDensity_;


public:

    //- Runtime type information
    TypeName("fixedCurrentDensity");


    // Constructors

        //- Construct from patch and internal field
        fixedCurrentDensityFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedCurrentDensityFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        fixedCurrentDensityFvPatchScalarField
        (
            const fixedCurrentDensityFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy construct
        fixedCurrentDensityFvPatchScalarField
        (
            const fixedCurrentDensityFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedCurrentDensityFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedCurrentDensityFvPatchScalarField
        (
            const fixedCurrentDensityFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedCurrentDensityFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member Functions

        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchScalarField&,
                const labelList&
            );


        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();


        //- Write
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
