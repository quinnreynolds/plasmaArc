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

Class
    Foam::fluidLutEOS

Description
    Equation of state (density-based) specified by lookup tables, for weakly 
    compressible fluids.

SourceFiles
    fluidLutEOSI.H
    fluidLutEOS.C

Q Reynolds 2015-2022

\*---------------------------------------------------------------------------*/

#ifndef fluidLutEOS_H
#define fluidLutEOS_H

#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Specie> class fluidLutEOS;

template<class Specie>
inline fluidLutEOS<Specie> operator+
(
    const fluidLutEOS<Specie>&,
    const fluidLutEOS<Specie>&
);

template<class Specie>
inline fluidLutEOS<Specie> operator-
(
    const fluidLutEOS<Specie>&,
    const fluidLutEOS<Specie>&
);

template<class Specie>
inline fluidLutEOS<Specie> operator*
(
    const scalar,
    const fluidLutEOS<Specie>&
);

template<class Specie>
inline fluidLutEOS<Specie> operator==
(
    const fluidLutEOS<Specie>&,
    const fluidLutEOS<Specie>&
);

template<class Specie>
Ostream& operator<<
(
    Ostream&,
    const fluidLutEOS<Specie>&
);


/*---------------------------------------------------------------------------*\
                    Class fluidLutEOS Declaration
\*---------------------------------------------------------------------------*/

template<class Specie>
class fluidLutEOS
:
    public Specie
{
    // Private data

        //density-temperature lookup table
        scalar pRef_;
        scalar rhoStartT_;
        scalar rhoDeltaT_;
        List<scalar> rhoData_;


public:

    // Constructors

        //- Construct from components
        inline fluidLutEOS
        (
            const Specie& sp,
            const scalar pRef,
            const scalar rhoStartT,
            const scalar rhoDeltaT,
            const List<scalar>& rhoData
        );

        //- Construct from fluidLutEOS
        inline fluidLutEOS(const fluidLutEOS& sp);

        //- Construct from Istream
        fluidLutEOS(Istream&);

        //- Construct from dictionary
        fluidLutEOS(const dictionary& dict);

        //- Construct as named copy
        inline fluidLutEOS
        (
            const word& name,
            const fluidLutEOS&
        );

        //- Construct and return a clone
        inline autoPtr<fluidLutEOS> clone() const;

        // Selector from Istream
        inline static autoPtr<fluidLutEOS> New(Istream& is);

        // Selector from dictionary
        inline static autoPtr<fluidLutEOS> New
        (
            const dictionary& dict
        );


    // Member functions

        //- Return the instantiated type name
        static word typeName()
        {
            return
                "fluidLutEOS<"
              + word(Specie::typeName_()) + '>';
        }


        // Fundamental properties

            //- Is the equation of state is incompressible i.e. rho != f(p)
            static const bool incompressible = true;

            //- Is the equation of state is isochoric i.e. rho = const
            static const bool isochoric = false;

            //- Return density [kg/m^3]
            inline scalar rho(scalar p, scalar T) const;

            //- Return enthalpy departure [J/kg]
            inline scalar H(const scalar p, const scalar T) const;

            //- Return Cp departure [J/(kg K]
            inline scalar Cp(scalar p, scalar T) const;

            //- Return internal energy departure [J/kg]
            inline scalar E(const scalar p, const scalar T) const;

            //- Return Cv departure [J/(kg K]
            inline scalar Cv(scalar p, scalar T) const;

            //- Return entropy [J/kg/K]
            inline scalar S(const scalar p, const scalar T) const;

            //- Return compressibility rho/p [s^2/m^2]
            inline scalar psi(scalar p, scalar T) const;

            //- Return compression factor []
            inline scalar Z(scalar p, scalar T) const;

            //- Return (Cp - Cv) [J/(kmol K]
            inline scalar CpMCv(scalar p, scalar T) const;


        // IO

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator=(const fluidLutEOS&);
        inline void operator+=(const fluidLutEOS&);
        inline void operator-=(const fluidLutEOS&);
        inline void operator*=(const scalar);


    // Friend operators

        friend fluidLutEOS operator+ <Specie>
        (
            const fluidLutEOS&,
            const fluidLutEOS&
        );

        friend fluidLutEOS operator- <Specie>
        (
            const fluidLutEOS&,
            const fluidLutEOS&
        );

        friend fluidLutEOS operator* <Specie>
        (
            const scalar s,
            const fluidLutEOS&
        );

        friend fluidLutEOS operator== <Specie>
        (
            const fluidLutEOS&,
            const fluidLutEOS&
        );


    // Ostream Operator

        friend Ostream& operator<< <Specie>
        (
            Ostream&,
            const fluidLutEOS&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluidLutEOSI.H"

#ifdef NoRepository
#   include "fluidLutEOS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
