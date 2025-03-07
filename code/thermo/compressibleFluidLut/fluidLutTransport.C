/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

#include "fluidLutTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
Foam::fluidLutTransport<Thermo>::fluidLutTransport(Istream& is)
:
    Thermo(is),
    muStartT_(readScalar(is)),
    muDeltaT_(readScalar(is)),
    muData_(is),
    tkStartT_(readScalar(is)),
    tkDeltaT_(readScalar(is)),
    tkData_(is)
{
    is.check("fluidLutTransport::fluidLutTransport(Istream& is)");
}


template<class Thermo>
Foam::fluidLutTransport<Thermo>::fluidLutTransport(const dictionary& dict)
:
    Thermo(dict),
    muStartT_(readScalar(dict.subDict("transport").subDict("muLookupTable").lookup("startT"))),
    muDeltaT_(readScalar(dict.subDict("transport").subDict("muLookupTable").lookup("deltaT"))),
    muData_(dict.subDict("transport").subDict("muLookupTable").lookup("dataTable")),
    tkStartT_(readScalar(dict.subDict("transport").subDict("tkLookupTable").lookup("startT"))),
    tkDeltaT_(readScalar(dict.subDict("transport").subDict("tkLookupTable").lookup("deltaT"))),
    tkData_(dict.subDict("transport").subDict("tkLookupTable").lookup("dataTable"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
void Foam::fluidLutTransport<Thermo>::fluidLutTransport::write(Ostream& os) const
{
    os  << this->name() << endl;
    os  << token::BEGIN_BLOCK  << incrIndent << nl;

    Thermo::write(os);

    dictionary dict("transport");
    dict.add("muLookupTable", dictionary("muLookupTable"));
    dict.subDict("muLookupTable").add("startT", muStartT_);
    dict.subDict("muLookupTable").add("deltaT", muDeltaT_);
    dict.subDict("muLookupTable").add("dataTable", muData_);
    dict.add("tkLookupTable", dictionary("tkLookupTable"));
    dict.subDict("tkLookupTable").add("startT", tkStartT_);
    dict.subDict("tkLookupTable").add("deltaT", tkDeltaT_);
    dict.subDict("tkLookupTable").add("dataTable", tkData_);
    os  << indent << dict.dictName() << dict;

    os  << decrIndent << token::END_BLOCK << nl;
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Thermo>
Foam::Ostream& Foam::operator<<(Ostream& os, const fluidLutTransport<Thermo>& ct)
{
    operator<<(os, static_cast<const Thermo&>(ct));
    os  << token::SPACE
        << token::SPACE << ct.muStartT_
        << token::SPACE << ct.muDeltaT_
        << token::SPACE << ct.muData_
        << token::SPACE << ct.tkStartT_
        << token::SPACE << ct.tkDeltaT_
        << token::SPACE << ct.tkData_;


    os.check
    (
        "Ostream& operator<<(Ostream&, const fluidLutTransport&)"
    );

    return os;
}


// ************************************************************************* //
