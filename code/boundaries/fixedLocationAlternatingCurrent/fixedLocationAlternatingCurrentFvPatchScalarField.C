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
    ANelPoten WARRANTelPoten; without even the implied warranty of MERCHANTABILITelPoten or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    elPotenou should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fixedLocationAlternatingCurrentFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "SortableList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::fixedLocationAlternatingCurrentFvPatchScalarField::t() const
{
    return db().time().timeOutputValue();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fixedLocationAlternatingCurrentFvPatchScalarField::
fixedLocationAlternatingCurrentFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchField(p, iF),
    current_(0),
    currentDensity_(0),
    frequency_(0),
    theta_(0),
    referencePosition_(0,0,0)
{
    this->refValue() = Zero;
    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


Foam::fixedLocationAlternatingCurrentFvPatchScalarField::
fixedLocationAlternatingCurrentFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    mixedFvPatchScalarField(p, iF),
    current_(readScalar(dict.lookup("maxCurrent"))),
    currentDensity_(readScalar(dict.lookup("currentDensity"))),
    frequency_(readScalar(dict.lookup("frequency"))),
    theta_(readScalar(dict.lookup("theta"))),
    referencePosition_(vector(dict.lookup("referencePosition")))
{
    if (dict.found("value"))
    {
        fvPatchScalarField::operator=
        (
            scalarField("value", dict, p.size())
        );
    }
    else
    {
        fvPatchScalarField::operator=(this->refValue());
    }

    this->refGrad() = Zero;
    this->valueFraction() = 0.0;
}


Foam::fixedLocationAlternatingCurrentFvPatchScalarField::
fixedLocationAlternatingCurrentFvPatchScalarField
(
    const fixedLocationAlternatingCurrentFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    mixedFvPatchScalarField(ptf),
    current_(ptf.current_),
    currentDensity_(ptf.currentDensity_),
    frequency_(ptf.frequency_),
    theta_(ptf.theta_),
    referencePosition_(ptf.referencePosition_)
{}


Foam::fixedLocationAlternatingCurrentFvPatchScalarField::
fixedLocationAlternatingCurrentFvPatchScalarField
(
    const fixedLocationAlternatingCurrentFvPatchScalarField& ptf
)
:
    mixedFvPatchScalarField(ptf),
    current_(ptf.current_),
    currentDensity_(ptf.currentDensity_),
    frequency_(ptf.frequency_),
    theta_(ptf.theta_),
    referencePosition_(ptf.referencePosition_)
{}


Foam::fixedLocationAlternatingCurrentFvPatchScalarField::
fixedLocationAlternatingCurrentFvPatchScalarField
(
    const fixedLocationAlternatingCurrentFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    mixedFvPatchScalarField(ptf, iF),
    current_(ptf.current_),
    currentDensity_(ptf.currentDensity_),
    frequency_(ptf.frequency_),
    theta_(ptf.theta_),
    referencePosition_(ptf.referencePosition_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fixedLocationAlternatingCurrentFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    mixedFvPatchScalarField::autoMap(m);
}


void Foam::fixedLocationAlternatingCurrentFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    mixedFvPatchScalarField::rmap(ptf, addr);
}

// Simplified and commetns added for clarity by HVH.
void Foam::fixedLocationAlternatingCurrentFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }
    
    const scalar pi = constant::mathematical::pi;
   
    scalar currentDensity = currentDensity_; 

    scalar current = current_ * sin(2 * pi * frequency_ * t() +  theta_ * (2 * pi) / 360);

    // checks if it is negative or possitive    
    scalar currentSign = current <= 0 ? -1 : 1;  
    
    label nProcs = Pstream::nProcs();
    label myProcNo = Pstream::myProcNo();

    if(current > 0)
    {
        const scalarField& ekPatch = patch().lookupPatchField<volScalarField, scalar>("ek"); 
        // creates a sclarfield pointer from the electric field patch field ek
        const vectorField& posPatch = patch().Cf(); 
        const scalarField distPatch(mag(posPatch - referencePosition_));

        scalarListList procDist(nProcs), procArea(nProcs);
        labelListList procLocalIndex(nProcs), procProc(nProcs);
        // List structures to hold the data from each processor

        procDist[myProcNo].resize(patch().size());
        procArea[myProcNo].resize(patch().size());
        procLocalIndex[myProcNo].resize(patch().size());
        procProc[myProcNo].resize(patch().size());
        // Static resize to avoid problems with DynamicList and length-zero lists
        
        forAll (distPatch, faceI)
        {
            procDist[myProcNo][faceI] = distPatch[faceI];
            procArea[myProcNo][faceI] = patch().magSf()[faceI];
            procLocalIndex[myProcNo][faceI] = faceI;
            procProc[myProcNo][faceI] = myProcNo;
        }
        // Load the data into the lists, processor-local
        
        #if OpenFOAM_VERSION >= 2312
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
        // Send and receive the list structures between all processors

        DynamicList<scalar> faceDist, faceArea;
        DynamicList<label> faceLocalIndex, faceProc;
        forAll (procDist, procI)
        {
            faceDist.append(procDist[procI]);
            faceArea.append(procArea[procI]);
            faceLocalIndex.append(procLocalIndex[procI]);
            faceProc.append(procProc[procI]);
        }
        // Recombine the list structures into flat lists

        SortableList<scalar> sortedDist(faceDist);  
        // the list of distance data is sorted

        scalar totalArea = mag(current) / currentDensity;
        // the total needed area is calculated

        scalar runningArea = 0;
        // the scalar that holds the are already looked at in the following "for loop"
        
        scalarField elPotenGradient(patch().size(), 0); 

        forAll(sortedDist, n)
        {
            label ni = sortedDist.indices()[n];
            runningArea += faceArea[ni];
            if (faceProc[ni] == myProcNo) // check to see if face ni is on *this* processor
            {
                scalar faceI = faceLocalIndex[ni];
                elPotenGradient[faceI] = currentSign * currentDensity / ekPatch[faceI];
            }
            if (runningArea > totalArea)
            {
                break;
            }
        }

        this->refGrad() = elPotenGradient;
        this->valueFraction() = 0.0;
    }

    else
    {
        this->refGrad() = 0;
        this->refValue() = 0;
        this->valueFraction() = 1;
    }

    mixedFvPatchScalarField::updateCoeffs();
    
}

void Foam::fixedLocationAlternatingCurrentFvPatchScalarField::write
(
    Ostream& os
) const
{
    fvPatchScalarField::write(os);
    os.writeEntry("maxCurrent", current_);
    os.writeEntry("currentDensity", currentDensity_);
    os.writeEntry("frequency", frequency_);
    os.writeEntry("theta", theta_);
    os.writeEntry("referencePosition", referencePosition_);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * Build Macro Function  * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        fixedLocationAlternatingCurrentFvPatchScalarField
    );
}

// ************************************************************************* //
