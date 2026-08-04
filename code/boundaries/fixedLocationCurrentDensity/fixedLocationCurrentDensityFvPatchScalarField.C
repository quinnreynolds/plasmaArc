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

#include "fixedLocationCurrentDensityFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "SortableList.H"

#include <cmath>

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::fixedLocationCurrentDensityFvPatchScalarField::t() const
{
    return db().time().timeOutputValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedLocationCurrentDensityFvPatchScalarField::
fixedLocationCurrentDensityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(p, iF),
    current_(),
    currentDensity_(),
    referencePosition_(0,0,0)
{
    refValue() = Zero;
    refGrad() = Zero;
    valueFraction() = Zero;
}


Foam::fixedLocationCurrentDensityFvPatchScalarField::
fixedLocationCurrentDensityFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    current_(Function1<scalar>::New("current", dict)),
    currentDensity_(Function1<scalar>::New("currentDensity", dict)),
    referencePosition_(vector(dict.lookup("referencePosition")))
{
    refGrad() = Zero;
    valueFraction() = Zero;

    fvPatchScalarField::operator = 
    ( 
        scalarField("value", dict, p.size()) 
    );
}


Foam::fixedLocationCurrentDensityFvPatchScalarField::
fixedLocationCurrentDensityFvPatchScalarField
(
    const fixedLocationCurrentDensityFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf, p, iF, mapper),
    current_(ptf.current_.clone()),
    currentDensity_(ptf.currentDensity_.clone()),
    referencePosition_(ptf.referencePosition_)
{}


Foam::fixedLocationCurrentDensityFvPatchScalarField::
fixedLocationCurrentDensityFvPatchScalarField
(
    const fixedLocationCurrentDensityFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    current_(ptf.current_.clone()),
    currentDensity_(ptf.currentDensity_.clone()),
    referencePosition_(ptf.referencePosition_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedLocationCurrentDensityFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}


void Foam::fixedLocationCurrentDensityFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);
}


void Foam::fixedLocationCurrentDensityFvPatchScalarField::updateCoeffs()
{
    // TEMPORARY DIAGNOSTIC - investigating an AMR+load-balancing SIGFPE.
    // Checks whether every processor calls this function (and its
    // collective gatherList/broadcastList section) the same number of
    // times, especially in the timestep right after a redistribution
    // event, and whether the gathered/selected data looks sane. Safe to
    // remove once the issue is root-caused.
    static label callCount_diag = 0;
    ++callCount_diag;
    Pout<< "[fixedLocationCurrentDensity::updateCoeffs] ENTRY "
        << "proc=" << Pstream::myProcNo()
        << " patch=" << patch().name()
        << " call=" << callCount_diag
        << " patchSize=" << patch().size()
        << " updatedAtEntry=" << updated()
        << " time=" << db().time().timeName()
        << endl;

    if (updated())
    {
        return;
    }

    scalar current = current_->value(t());
    scalar currentDensity = currentDensity_->value(t());
    scalar currentSign = current < 0 ? -1 : 1;

    label nProcs = Pstream::nProcs();
    label myProcNo = Pstream::myProcNo();

    const scalarField& ekPatch = patch().lookupPatchField<volScalarField, scalar>("ek"); 
    const vectorField& posPatch = patch().Cf(); 
    const scalarField distPatch(mag(posPatch - referencePosition_));

    scalarListList procDist(nProcs), procArea(nProcs);
    labelListList procLocalIndex(nProcs), procProc(nProcs);
    procDist[myProcNo].resize(patch().size());
    procArea[myProcNo].resize(patch().size());
    procLocalIndex[myProcNo].resize(patch().size());
    procProc[myProcNo].resize(patch().size());
    forAll (distPatch, faceI)
    {
        procDist[myProcNo][faceI] = distPatch[faceI];
        procArea[myProcNo][faceI] = patch().magSf()[faceI];
        procLocalIndex[myProcNo][faceI] = faceI;
        procProc[myProcNo][faceI] = myProcNo;
    }
    
    #if (OPENFOAM >= 2312)
        Pstream::gatherList(procDist);
        Pstream::broadcastList(procDist);
        Pstream::gatherList(procArea);
        Pstream::broadcastList(procArea);
        Pstream::gatherList(procLocalIndex);
        Pstream::broadcastList(procLocalIndex);
        Pstream::gatherList(procProc);
        Pstream::broadcastList(procProc);
    #else 
        Pstream::gatherList(procDist);
        Pstream::scatterList(procDist);
        Pstream::gatherList(procArea);
        Pstream::scatterList(procArea);
        Pstream::gatherList(procLocalIndex);
        Pstream::scatterList(procLocalIndex);
        Pstream::gatherList(procProc);
        Pstream::scatterList(procProc);
    #endif

    // TEMPORARY DIAGNOSTIC - confirms every processor reached this point
    // (i.e. all four gatherList/broadcastList round-trips completed
    // without hanging) and reports the total gathered face count, which
    // should be identical - the sum of every processor's local patch
    // size - on every rank if the collective calls stayed in sync.
    {
        label totalGathered = 0;
        forAll(procDist, procI)
        {
            totalGathered += procDist[procI].size();
        }
        Pout<< "[fixedLocationCurrentDensity::updateCoeffs] "
            << "POST-GATHER proc=" << Pstream::myProcNo()
            << " call=" << callCount_diag
            << " nProcs=" << nProcs
            << " totalGatheredFaces=" << totalGathered
            << endl;
    }

    DynamicList<scalar> faceDist, faceArea;
    DynamicList<label> faceLocalIndex, faceProc;
    forAll (procDist, procI)
    {
        faceDist.append(procDist[procI]);
        faceArea.append(procArea[procI]);
        faceLocalIndex.append(procLocalIndex[procI]);
        faceProc.append(procProc[procI]);
    }
    SortableList<scalar> sortedDist(faceDist);

    scalar totalArea = mag(current) / currentDensity;

    scalar runningArea = 0;
    
    scalarField elPotenGradient(patch().size(), 0); 

    forAll(sortedDist, n)
    {
        label ni = sortedDist.indices()[n];
        runningArea += faceArea[ni];
        if (faceProc[ni] == myProcNo)
        {
            scalar faceI = faceLocalIndex[ni];

            // TEMPORARY DIAGNOSTIC - guard the two ways this indexing/
            // division could go wrong: faceI landing outside ekPatch's
            // current size (a stale/mismatched "ek" boundary field after
            // redistribution), or dividing by a near-zero conductivity.
            if (faceI < 0 || faceI >= ekPatch.size())
            {
                Pout<< "[fixedLocationCurrentDensity::updateCoeffs] "
                    << "OUT-OF-BOUNDS proc=" << Pstream::myProcNo()
                    << " call=" << callCount_diag
                    << " faceI=" << faceI
                    << " ekPatchSize=" << ekPatch.size()
                    << " patchSize=" << patch().size()
                    << endl;
            }
            else if (mag(ekPatch[faceI]) < SMALL)
            {
                Pout<< "[fixedLocationCurrentDensity::updateCoeffs] "
                    << "NEAR-ZERO-EK proc=" << Pstream::myProcNo()
                    << " call=" << callCount_diag
                    << " faceI=" << faceI
                    << " ekPatchValue=" << ekPatch[faceI]
                    << endl;
            }
            else
            {
                elPotenGradient[faceI] = currentSign * currentDensity / ekPatch[faceI];
            }
        }
        if (runningArea > totalArea)
        {
            break;
        }
    }

    // TEMPORARY DIAGNOSTIC - confirm the final refGrad() payload is finite
    // before it goes into the ePot matrix assembly.
    {
        bool allFinite = true;
        forAll(elPotenGradient, i)
        {
            if (!std::isfinite(elPotenGradient[i]))
            {
                allFinite = false;
                break;
            }
        }
        Pout<< "[fixedLocationCurrentDensity::updateCoeffs] EXIT "
            << "proc=" << Pstream::myProcNo()
            << " call=" << callCount_diag
            << " elPotenGradientFinite=" << allFinite
            << endl;
    }

    this->refGrad() = elPotenGradient;
    this->valueFraction() = 0.0;

    mixedFvPatchScalarField::updateCoeffs();
}


void Foam::fixedLocationCurrentDensityFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    current_->writeData(os);
    currentDensity_->writeData(os);
    os.writeEntry("referencePosition", referencePosition_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedLocationCurrentDensityFvPatchScalarField
    );
}

// ************************************************************************* //
