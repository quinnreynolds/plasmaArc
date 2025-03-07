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
    Foam::fixedLocationAlternatingCurrentFvPatchScalarField

Group
    grpGenericBoundaryConditions

Description
    This boundary condition provides a sinusoidal veriable current boundary 
    condition for electric potential fields. The BC outputs a veriable current 
    for the positive part of the sine curve but can then be either set to fixedValue
    or zero gradient. 

    The sine curve can be shifted by specifying theta [deg] in the input directory.

    This conditions should be used to simulate an arc that reverses direction in a
    SAF by having an oposite conditions of the same time on the lower emitter that 
    is shifted by 180°. Can also be used in a different maner if needed. 

    Faces on the boundary equal in total area to current / currentDensity are 
    selected based on the their proximity to the geometry's centerline (x=0, y=0)
    
    At the attachment spot faces the BC is calculated as:

        \f[
            d/dx (ePot) = jRD / ek[i]
        \f]
    
    For all other faces, a zero-gradient condition is applied.

Usage

    maxCurrent -> is the peak current or amplitude of the sine curve.
    currentDensity -> is the current density of the electrode, used to 
                       estimate the spot size  
    theta -> the phase shift of the curve  

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type                fixedLocationAlternatingCurrent;
        maxCurrent          5000;
        currentDensity      2e7;
        frequency           50;
        theta               0;
        referencePosition   (0 0 0);
        value               uniform 0;
    }
    \endverbatim

SourceFiles
    fixedLocationAlternatingCurrentFvPatchScalarField.C

Author
    Quinn Reynolds 2023
    based on mixedVariableCurrentDensity by Hákon Valur Haraldsson
\*---------------------------------------------------------------------------*/

#ifndef fixedLocationAlternatingCurrentFvPatchScalarField_H
#define fixedLocationAlternatingCurrentFvPatchScalarField_H

#include "mixedFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class fixedLocationAlternatingCurrentFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class fixedLocationAlternatingCurrentFvPatchScalarField
:
    public mixedFvPatchScalarField
{
    // Private data

        scalar current_;

        scalar currentDensity_;

        scalar frequency_;

        scalar theta_;

        vector referencePosition_;
        
    // Private Member Functions

        //- Return current time
        scalar t() const;


public:

    //- Runtime type information
    TypeName("fixedLocationAlternatingCurrent");


    // Constructors

        //- Construct from patch and internal field
        fixedLocationAlternatingCurrentFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        fixedLocationAlternatingCurrentFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given mixedTypeFvPatchField
        //  onto a new patch
        fixedLocationAlternatingCurrentFvPatchScalarField
        (
            const fixedLocationAlternatingCurrentFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        fixedLocationAlternatingCurrentFvPatchScalarField
        (
            const fixedLocationAlternatingCurrentFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new fixedLocationAlternatingCurrentFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        fixedLocationAlternatingCurrentFvPatchScalarField
        (
            const fixedLocationAlternatingCurrentFvPatchScalarField&,
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
                new fixedLocationAlternatingCurrentFvPatchScalarField(*this,iF)
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
