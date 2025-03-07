/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
     \\/     M anipulation  |
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
    Foam::limitedTemperatureFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides an limitedTemperature condition for scalar
    fields. It implements a zero-gadient BC while the internal field temperature
    is below(above) Tbound, and a fixed-value BC while it is above(below) 
    Tbound. The upperBoundYN boolean controls whether Tbound acts as an upper 
    (maximum) or lower (minimum) limit.

Usage
    \table
        Property            | Description           | Req'd? | Default
        Tbound              | single scalar value   | yes    |
        upperBoundYN        | boolean               | no     | true
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type                limitedTemperature;
        Tbound              4100.;
        upperBoundYN        true;
    }
    \endverbatim

SourceFiles
    limitedTemperatureFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef limitedTemperatureFvPatchScalarField_H
#define limitedTemperatureFvPatchScalarField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class limitedTemperatureFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class limitedTemperatureFvPatchScalarField
:
    public fixedValueFvPatchScalarField
{
    // Private data

    scalar Tbound_;

    bool upperBoundYN_;


public:

    //- Runtime type information
    TypeName("limitedTemperature");


    // Constructors

        //- Construct from patch and internal field
        limitedTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        limitedTemperatureFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given fixedValueTypeFvPatchField
        //  onto a new patch
        limitedTemperatureFvPatchScalarField
        (
            const limitedTemperatureFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        limitedTemperatureFvPatchScalarField
        (
            const limitedTemperatureFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new limitedTemperatureFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        limitedTemperatureFvPatchScalarField
        (
            const limitedTemperatureFvPatchScalarField&,
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
                new limitedTemperatureFvPatchScalarField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

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
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
