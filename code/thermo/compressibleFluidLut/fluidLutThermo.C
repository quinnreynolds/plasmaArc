/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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

Description
    Thermodynamics (density-based) specified by lookup tables, for weakly 
    compressible fluids.

Q Reynolds 2015-2022

\*---------------------------------------------------------------------------*/

#include "fluidLutThermo.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class equationOfState>
Foam::fluidLutThermo<equationOfState>::fluidLutThermo(Istream& is)
:
    equationOfState(is),
    CpStartT_(readScalar(is)),
    CpDeltaT_(readScalar(is)),
    CpData_(is)
{
    is.check("fluidLutThermo::fluidLutThermo(Istream& is)");
    hData_ = CpData_;
    hData_[0] = CpStartT_ * CpData_[0];
    for (label i=1; i<CpData_.size(); i++)
    {
        hData_[i] = hData_[i-1] + 0.5*(CpData_[i]+CpData_[i-1])*CpDeltaT_;
    }
}


template<class equationOfState>
Foam::fluidLutThermo<equationOfState>::fluidLutThermo(const dictionary& dict)
:
    equationOfState(dict),
    CpStartT_(readScalar(dict.subDict("thermodynamics").subDict("CpLookupTable").lookup("startT"))),
    CpDeltaT_(readScalar(dict.subDict("thermodynamics").subDict("CpLookupTable").lookup("deltaT"))),
    CpData_(dict.subDict("thermodynamics").subDict("CpLookupTable").lookup("dataTable"))
{
    hData_ = CpData_;
    hData_[0] = CpStartT_ * CpData_[0];
    for (label i=1; i<CpData_.size(); i++)
    {
        hData_[i] = hData_[i-1] + 0.5*(CpData_[i]+CpData_[i-1])*CpDeltaT_;
    }    
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class equationOfState>
void Foam::fluidLutThermo<equationOfState>::write(Ostream& os) const
{
    equationOfState::write(os);

    dictionary dict("thermodynamics");
    dict.add("CpLookupTable", dictionary("CpLookupTable"));
    dict.subDict("CpLookupTable").add("startT", CpStartT_);
    dict.subDict("CpLookupTable").add("deltaT", CpDeltaT_);
    dict.subDict("CpLookupTable").add("dataTable", CpData_);
    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class equationOfState>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const fluidLutThermo<equationOfState>& ct
)
{
    os  << static_cast<const equationOfState&>(ct)
        << token::SPACE << ct.CpStartT_
        << token::SPACE << ct.CpDeltaT_
        << token::SPACE << ct.CpData_;

    os.check
    (
        "Ostream& operator<<(Ostream& os, const fluidLutThermo& ct)"
    );

    return os;
}


// ************************************************************************* //
