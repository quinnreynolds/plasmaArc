/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::fluidLutTransport

Description
    Transport properties based on lookup tables for weakly compressible fluids.

SourceFiles
    fluidLutTransportI.H
    fluidLutTransport.C

Q Reynolds 2015-2022

\*---------------------------------------------------------------------------*/

#ifndef fluidLutTransport_H
#define fluidLutTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class fluidLutTransport;

template<class Thermo>
inline fluidLutTransport<Thermo> operator+
(
    const fluidLutTransport<Thermo>&,
    const fluidLutTransport<Thermo>&
);

template<class Thermo>
inline fluidLutTransport<Thermo> operator-
(
    const fluidLutTransport<Thermo>&,
    const fluidLutTransport<Thermo>&
);

template<class Thermo>
inline fluidLutTransport<Thermo> operator*
(
    const scalar,
    const fluidLutTransport<Thermo>&
);

template<class Thermo>
inline fluidLutTransport<Thermo> operator==
(
    const fluidLutTransport<Thermo>&,
    const fluidLutTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const fluidLutTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                           Class fluidLutTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class fluidLutTransport
:
    public Thermo
{
    // Private data

        scalar muStartT_;
        scalar muDeltaT_;
        List<scalar> muData_;
        scalar tkStartT_;
        scalar tkDeltaT_;
        List<scalar> tkData_;


    // Private Member Functions

        //- Construct from components
        inline fluidLutTransport
        (
            const Thermo& t,
            const scalar muStartT,
            const scalar muDeltaT,
            const List<scalar>& muData,
            const scalar tkStartT,
            const scalar tkDeltaT,
            const List<scalar>& tkData
        );


public:

    // Constructors

        //- Construct as named copy
        inline fluidLutTransport(const word&, const fluidLutTransport&);

        //- Construct from Istream
        fluidLutTransport(Istream&);

        //- Construct from dictionary
        fluidLutTransport(const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<fluidLutTransport> clone() const;

        // Selector from Istream
        inline static autoPtr<fluidLutTransport> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<fluidLutTransport> New(const dictionary& dict);


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "fluidLutTransport<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/ms]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/mK]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Thermal diffusivity of enthalpy [kg/ms]
        inline scalar alphah(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member operators

        inline void operator=(const fluidLutTransport&);
        inline void operator+=(const fluidLutTransport&);
        inline void operator-=(const fluidLutTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend fluidLutTransport operator+ <Thermo>
        (
            const fluidLutTransport&,
            const fluidLutTransport&
        );

        friend fluidLutTransport operator- <Thermo>
        (
            const fluidLutTransport&,
            const fluidLutTransport&
        );

        friend fluidLutTransport operator* <Thermo>
        (
            const scalar,
            const fluidLutTransport&
        );

        friend fluidLutTransport operator== <Thermo>
        (
            const fluidLutTransport&,
            const fluidLutTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const fluidLutTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluidLutTransportI.H"

#ifdef NoRepository
#   include "fluidLutTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
