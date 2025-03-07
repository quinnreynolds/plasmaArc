/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "greyPlasmaAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace radiation
    {
        defineTypeNameAndDebug(greyPlasmaAbsorptionEmission, 0);

        addToRunTimeSelectionTable
        (
            absorptionEmissionModel,
            greyPlasmaAbsorptionEmission,
            dictionary
        );
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::greyPlasmaAbsorptionEmission::greyPlasmaAbsorptionEmission
(
    const dictionary& dict,
    const fvMesh& mesh
)
:
    absorptionEmissionModel(dict, mesh),
    coeffsDict_(dict.optionalSubDict("plasmaRadiationData")),
    aStartT_(readScalar(coeffsDict_.subDict("absorptionCoeffs").lookup("startT"))),
    aDeltaT_(readScalar(coeffsDict_.subDict("absorptionCoeffs").lookup("deltaT"))),
    aData_(coeffsDict_.subDict("absorptionCoeffs").lookup("dataTable")),
    eStartT_(readScalar(coeffsDict_.subDict("emissionCoeffs").lookup("startT"))),
    eDeltaT_(readScalar(coeffsDict_.subDict("emissionCoeffs").lookup("deltaT"))),
    eData_(coeffsDict_.subDict("emissionCoeffs").lookup("dataTable"))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::greyPlasmaAbsorptionEmission::~greyPlasmaAbsorptionEmission()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::greyPlasmaAbsorptionEmission::aCont(const label bandI) const
{
    tmp<volScalarField> ta
    (
        new volScalarField
        (
            IOobject
            (
                "aCont" + name(bandI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("a", dimless/dimLength, 0.0)
        )
    );

    scalar aInvDT = 1. / aDeltaT_;
    scalarField& a = ta.ref().primitiveFieldRef();
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    //scalarField& T = tT.ref().primitiveFieldRef();

    forAll(a, cellI)
    {
        scalar Tloc = T[cellI];
        label iT = label((Tloc - aStartT_) * aInvDT);
        if (iT > (aData_.size()-2))
        {
            a[cellI] = aData_[aData_.size()-1];
        }
        else if (iT < 0)
        {
            a[cellI] = aData_[0];
        }
        else
        {
            scalar deltaFrac = (Tloc - aStartT_) * aInvDT - scalar(iT);
            a[cellI] = aData_[iT] + deltaFrac * (aData_[iT+1] - aData_[iT]);
        }
    }

    ta.ref().correctBoundaryConditions();
    
    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::greyPlasmaAbsorptionEmission::eCont(const label bandI) const
{
    tmp<volScalarField> te
    (
        new volScalarField
        (
            IOobject
            (
                "eCont" + name(bandI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("e", dimless/dimLength, 0.0)
        )
    );

    scalar eInvDT = 1. / eDeltaT_;
    scalarField& e = te.ref().primitiveFieldRef();
    const volScalarField& T = mesh_.lookupObject<volScalarField>("T");
    //scalarField& T = tT.ref().primitiveFieldRef();

    forAll(e, cellI)
    {
        scalar Tloc = T[cellI];
        label iT = label((Tloc - eStartT_) * eInvDT);
        if (iT > (eData_.size()-2))
        {
            e[cellI] = eData_[eData_.size()-1];
        }
        else if (iT < 0)
        {
            e[cellI] = eData_[0];
        }
        else
        {
            scalar deltaFrac = (Tloc - eStartT_) * eInvDT - scalar(iT);
            e[cellI] = eData_[iT] + deltaFrac * (eData_[iT+1] - eData_[iT]);
        }
    }

    te.ref().correctBoundaryConditions();

    return te;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::greyPlasmaAbsorptionEmission::ECont(const label bandI) const
{
    tmp<volScalarField> tE
    (
      new volScalarField
      (
          IOobject
          (
              "ECont" + name(bandI),
              mesh_.time().timeName(),
              mesh_,
              IOobject::NO_READ,
              IOobject::NO_WRITE
          ),
          mesh_,
          dimensionedScalar("E", dimMass/dimLength/pow3(dimTime), 0.0)
      )
    );

    return tE;
}


// ************************************************************************* //
