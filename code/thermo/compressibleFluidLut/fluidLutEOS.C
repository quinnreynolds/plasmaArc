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
    Equation of state (density-based) specified by lookup tables, for weakly 
    compressible fluids.

    Q Reynolds 2015-2022

\*---------------------------------------------------------------------------*/

#include "fluidLutEOS.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
Foam::fluidLutEOS<Specie>::fluidLutEOS(Istream& is)
:
    Specie(is),
    pRef_(readScalar(is)),
    rhoStartT_(readScalar(is)),
    rhoDeltaT_(readScalar(is)),
    rhoData_(is)
{
    is.check
    (
        "fluidLutEOS<Specie>::"
        "fluidLutEOS(Istream& is)"
    );
}


template<class Specie>
Foam::fluidLutEOS<Specie>::fluidLutEOS(const dictionary& dict)
:
    Specie(dict),
    pRef_(readScalar(dict.subDict("equationOfState").lookup("pRef"))),
    rhoStartT_(readScalar(dict.subDict("equationOfState").subDict("rhoLookupTable").lookup("startT"))),
    rhoDeltaT_(readScalar(dict.subDict("equationOfState").subDict("rhoLookupTable").lookup("deltaT"))),
    rhoData_(dict.subDict("equationOfState").subDict("rhoLookupTable").lookup("dataTable"))
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
void Foam::fluidLutEOS<Specie>::write(Ostream& os) const
{
    Specie::write(os);
    dictionary dict("equationOfState");
    dict.add("pRef", rhoStartT_);
    dict.add("rhoLookupTable", dictionary("rhoLookupTable"));
    dict.subDict("rhoLookupTable").add("startT", rhoStartT_);
    dict.subDict("rhoLookupTable").add("deltaT", rhoDeltaT_);
    dict.subDict("rhoLookupTable").add("dataTable", rhoData_);

    os  << indent << dict.dictName() << dict;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class Specie>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const fluidLutEOS<Specie>& pg
)
{
    os  << static_cast<const Specie&>(pg)
        << token::SPACE << pg.pRef_
        << token::SPACE << pg.rhoStartT_
        << token::SPACE << pg.rhoDeltaT_
        << token::SPACE << pg.rhoData_;

    os.check
    (
        "Ostream& operator<<"
        "(Ostream& os, const fluidLutEOS<Specie>& st)"
    );
    return os;
}


// ************************************************************************* //
