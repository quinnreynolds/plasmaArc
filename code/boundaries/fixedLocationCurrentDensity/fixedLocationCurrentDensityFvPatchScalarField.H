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
    Foam::fixedLocationCurrentDensityFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a fixedLocationCurrentDensity condition for
    electric potential fields. Faces on the boundary equal in total area to 
    current / currentDensity are selected based on their proximity to a 
    specified point in space. At the attachment spot faces the BC is 
    calculated as:

        \f[
            d/dx (ePot) = jRD / ek[i]
        \f]
    
    For all other faces, a zero-gradient condition is applied.

Usage
    \table
        Property             | Description            | Req'd? | Default
        current              | table of scalar values | yes    |
        currentDensity       | table of scalar values | yes    |
        referencePosition    | vector                 | yes    | 
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            fixedLocationCurrentDensity;
        current         table 
                        (
                            (0 1000)
                            (1 1000)
                        );
        currentDensity  table 
                        (
                            (0 2e7)
                            (1 2e7)
                        );
        referencePosition   (0 0 0);
        value           uniform 0.;
    }
    \endverbatim

SourceFiles
    fixedLocationCurrentDensityFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef fixedLocationCurrentDensityFvPatchScalarField_H
#define fixedLocationCurrentDensityFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class fixedLocationCurrentDensityFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedLocationCurrentDensityFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        autoPtr<Function1<scalar>> current_;

        autoPtr<Function1<scalar>> currentDensity_;

        vector referencePosition_;

    // Private Member Functions

        //- Return current time
        scalar t() const;


public:

    //- Runtime type information
    TypeName("fixedLocationCurrentDensity");


    // Constructors

        //- Construct from patch and internal field
        fixedLocationCurrentDensityFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedLocationCurrentDensityFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mixedTypeFvPatchField
        //  onto a new patch
        fixedLocationCurrentDensityFvPatchScalarField
        (
            const fixedLocationCurrentDensityFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedLocationCurrentDensityFvPatchScalarField
        (
            const fixedLocationCurrentDensityFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedLocationCurrentDensityFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedLocationCurrentDensityFvPatchScalarField
        (
            const fixedLocationCurrentDensityFvPatchScalarField&,
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
                new fixedLocationCurrentDensityFvPatchScalarField
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
