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

\*---------------------------------------------------------------------------*/

#include "fixedCurrentDensityFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedCurrentDensityFvPatchScalarField::
fixedCurrentDensityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    currentDensity_(p.size(), Zero)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::fixedCurrentDensityFvPatchScalarField::
fixedCurrentDensityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    currentDensity_("currentDensity", dict, p.size())
{
    refGrad() = Zero;
    valueFraction() = Zero;

    fvPatchScalarField::operator=
    (
        scalarField("value", dict, p.size())
    );
    refValue() = *this;
}


Foam::fixedCurrentDensityFvPatchScalarField::
fixedCurrentDensityFvPatchScalarField
(
    const fixedCurrentDensityFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    currentDensity_(mapper(ptf.currentDensity_))
{}


Foam::fixedCurrentDensityFvPatchScalarField::
fixedCurrentDensityFvPatchScalarField
(
    const fixedCurrentDensityFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    currentDensity_(ptf.currentDensity_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedCurrentDensityFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::map(ptf, m);
    const fixedCurrentDensityFvPatchScalarField& tiptf =
        refCast<const fixedCurrentDensityFvPatchScalarField>(ptf);
    m(currentDensity_, tiptf.currentDensity_);
}


void Foam::fixedCurrentDensityFvPatchScalarField::reset
(
    const fvPatchScalarField& ptf
)
{
    mixedFvPatchScalarField::reset(ptf);
    const fixedCurrentDensityFvPatchScalarField& tiptf =
        refCast<const fixedCurrentDensityFvPatchScalarField>(ptf);
    currentDensity_ = tiptf.currentDensity_;
}


void Foam::fixedCurrentDensityFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);

    const fixedCurrentDensityFvPatchScalarField& tiptf =
        refCast<const fixedCurrentDensityFvPatchScalarField>(ptf);

    currentDensity_.rmap(tiptf.currentDensity_, addr);
}


void Foam::fixedCurrentDensityFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalarField& ekPatch =
        patch().lookupPatchField<volScalarField, scalar>
        (
            "ek"
        );
    
    mixedFvPatchScalarField::refGrad() =
    (
        currentDensity_ / ekPatch
    );

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::fixedCurrentDensityFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "currentDensity", currentDensity_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedCurrentDensityFvPatchScalarField
    );
}

// ************************************************************************* //
