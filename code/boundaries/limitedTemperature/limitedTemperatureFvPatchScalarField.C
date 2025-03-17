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

\*---------------------------------------------------------------------------*/

#include "limitedTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::limitedTemperatureFvPatchScalarField::
limitedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    Tbound_( 0. ),
    upperBoundYN_(true)
{}


Foam::limitedTemperatureFvPatchScalarField::
limitedTemperatureFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF),
    Tbound_(readScalar(dict.lookup("Tbound"))),
    upperBoundYN_(dict.lookupOrDefault("upperBoundYN", true))
{
    fvPatchScalarField::operator=
    (
        scalarField("value", dict, p.size())
    );
}


Foam::limitedTemperatureFvPatchScalarField::
limitedTemperatureFvPatchScalarField
(
    const limitedTemperatureFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    Tbound_(ptf.Tbound_),
    upperBoundYN_(ptf.upperBoundYN_)
{}


Foam::limitedTemperatureFvPatchScalarField::
limitedTemperatureFvPatchScalarField
(
    const limitedTemperatureFvPatchScalarField& ptf
)
:
    fixedValueFvPatchScalarField(ptf),
    Tbound_(ptf.Tbound_),
    upperBoundYN_(ptf.upperBoundYN_)
{}


Foam::limitedTemperatureFvPatchScalarField::
limitedTemperatureFvPatchScalarField
(
    const limitedTemperatureFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    Tbound_(ptf.Tbound_),
    upperBoundYN_(ptf.upperBoundYN_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::limitedTemperatureFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchScalarField::autoMap(m);
}


void Foam::limitedTemperatureFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchScalarField::rmap(ptf, addr);
}


void Foam::limitedTemperatureFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalarField TNew(patch().size(), Zero);
    TNew = patchInternalField();
    if (upperBoundYN_)
    {
        forAll(TNew, faceI)
        {
            if (TNew[faceI] > Tbound_)
            {
                TNew[faceI] = Tbound_;
            }
        }
    }
    else
    {
        forAll(TNew, faceI)
        {
            if (TNew[faceI] < Tbound_)
            {
                TNew[faceI] = Tbound_;
            }
        }        
    }
    
    fixedValueFvPatchScalarField::operator==
    (
        TNew
    );

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::limitedTemperatureFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("Tbound", Tbound_);
    os.writeEntry("upperBoundYN", upperBoundYN_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        limitedTemperatureFvPatchScalarField
    );
}

// ************************************************************************* //
