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
    Foam::fluidLutThermo

Description
    Thermodynamics (density-based) specified by lookup tables, for weakly 
    compressible fluids.

SourceFiles
    fluidLutThermoI.H
    fluidLutThermo.C

Q Reynolds 2015-2022
\*---------------------------------------------------------------------------*/

#ifndef fluidLutThermo_H
#define fluidLutThermo_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class EquationOfState> class fluidLutThermo;

template<class EquationOfState>
inline fluidLutThermo<EquationOfState> operator+
(
    const fluidLutThermo<EquationOfState>&,
    const fluidLutThermo<EquationOfState>&
);

template<class EquationOfState>
inline fluidLutThermo<EquationOfState> operator-
(
    const fluidLutThermo<EquationOfState>&,
    const fluidLutThermo<EquationOfState>&
);

template<class EquationOfState>
inline fluidLutThermo<EquationOfState> operator*
(
    const scalar,
    const fluidLutThermo<EquationOfState>&
);

template<class EquationOfState>
inline fluidLutThermo<EquationOfState> operator==
(
    const fluidLutThermo<EquationOfState>&,
    const fluidLutThermo<EquationOfState>&
);

template<class EquationOfState>
Ostream& operator<<
(
    Ostream&,
    const fluidLutThermo<EquationOfState>&
);


/*---------------------------------------------------------------------------*\
                           Class fluidLutThermo Declaration
\*---------------------------------------------------------------------------*/

template<class EquationOfState>
class fluidLutThermo
:
    public EquationOfState
{
    // Private data

        scalar CpStartT_;
        scalar CpDeltaT_;
        List<scalar> CpData_;
        List<scalar> hData_;


    // Private Member Functions

    //- Construct from components
        inline fluidLutThermo
        (
            const EquationOfState& st,
            const scalar CpStartT,
            const scalar CpDeltaT,
            const List<scalar>& CpData
        );


public:

    // Constructors

        //- Construct from Istream
        fluidLutThermo(Istream&);

        //- Construct from dictionary
        fluidLutThermo(const dictionary& dict);

        //- Construct as named copy
        inline fluidLutThermo(const word&, const fluidLutThermo&);

        //- Construct and return a clone
        inline autoPtr<fluidLutThermo> clone() const;

        //- Selector from Istream
        inline static autoPtr<fluidLutThermo> New(Istream& is);

        //- Selector from dictionary
        inline static autoPtr<fluidLutThermo> New(const dictionary& dict);


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "fluidLutThermo<" + EquationOfState::typeName() + '>';
        }

        //- Limit the temperature to be in the range Tlow_ to Thigh_
        inline scalar limit(const scalar T) const;


        // Fundamental properties

            //- Heat capacity at constant pressure [J/(kg K)]
            inline scalar Cp(const scalar p, const scalar T) const;

            //- Absolute Enthalpy [J/kg]
            inline scalar Ha(const scalar p, const scalar T) const;

            //- Sensible enthalpy [J/kg]
            inline scalar Hs(const scalar p, const scalar T) const;

            //- Chemical enthalpy [J/kg]
            inline scalar Hc() const;

            //- Entropy [J/(kg K)]
            inline scalar S(const scalar p, const scalar T) const;

            #include "HtoEthermo.H"

        // I-O

            //- Write to Ostream
            void write(Ostream& os) const;


    // Member operators

        inline void operator+=(const fluidLutThermo&);
        

    // Friend operators

        friend fluidLutThermo operator+ <EquationOfState>
        (
            const fluidLutThermo&,
            const fluidLutThermo&
        );

        friend fluidLutThermo operator- <EquationOfState>
        (
            const fluidLutThermo&,
            const fluidLutThermo&
        );

        friend fluidLutThermo operator* <EquationOfState>
        (
            const scalar,
            const fluidLutThermo&
        );

        friend fluidLutThermo operator== <EquationOfState>
        (
            const fluidLutThermo&,
            const fluidLutThermo&
        );


    // IOstream Operators

        friend Ostream& operator<< <EquationOfState>
        (
            Ostream&,
            const fluidLutThermo&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fluidLutThermoI.H"

#ifdef NoRepository
#   include "fluidLutThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
