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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Specie>
inline Foam::fluidLutEOS<Specie>::fluidLutEOS
(
    const Specie& sp,
    const scalar pRef,
    const scalar rhoStartT,
    const scalar rhoDeltaT,
    const List<scalar>& rhoData
)
:
    Specie(sp),
    pRef_(pRef),
    rhoStartT_(rhoStartT),
    rhoDeltaT_(rhoDeltaT),
    rhoData_(rhoData)
{}


template<class Specie>
inline Foam::fluidLutEOS<Specie>::fluidLutEOS
(
    const fluidLutEOS& ipg
)
:
    Specie(ipg),
    pRef_(ipg.pRef_),
    rhoStartT_(ipg.rhoStartT_),
    rhoDeltaT_(ipg.rhoDeltaT_),
    rhoData_(ipg.rhoData_)
{}


template<class Specie>
inline Foam::fluidLutEOS<Specie>::fluidLutEOS
(
    const word& name,
    const fluidLutEOS<Specie>& ipg
)
:
    Specie(name, ipg),
    pRef_(ipg.pRef_),
    rhoStartT_(ipg.rhoStartT_),
    rhoDeltaT_(ipg.rhoDeltaT_),
    rhoData_(ipg.rhoData_)
{}


template<class Specie>
inline Foam::autoPtr<Foam::fluidLutEOS<Specie> >
Foam::fluidLutEOS<Specie>::clone() const
{
    return autoPtr<fluidLutEOS<Specie> >
    (
        new fluidLutEOS<Specie>(*this)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::fluidLutEOS<Specie> >
Foam::fluidLutEOS<Specie>::New
(
    Istream& is
)
{
    return autoPtr<fluidLutEOS<Specie> >
    (
        new fluidLutEOS<Specie>(is)
    );
}


template<class Specie>
inline Foam::autoPtr<Foam::fluidLutEOS<Specie> >
Foam::fluidLutEOS<Specie>::New
(
    const dictionary& dict
)
{
    return autoPtr<fluidLutEOS<Specie> >
    (
        new fluidLutEOS<Specie>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::scalar Foam::fluidLutEOS<Specie>::rho
(
    scalar p, 
    scalar T
) const
{
    label iT = label((T - rhoStartT_) / rhoDeltaT_);
    if (iT > (rhoData_.size()-2))
    {
        return rhoData_[rhoData_.size()-1];
    }
    if (iT < 0)
    {
        return rhoData_[0];
    }
    scalar deltaFrac = (T - rhoStartT_) / rhoDeltaT_ - scalar(iT);
    scalar rho0 = rhoData_[iT] + deltaFrac * (rhoData_[iT+1] - rhoData_[iT]);

    return rho0 * p / pRef_;
}


template<class Specie>
inline Foam::scalar Foam::fluidLutEOS<Specie>::psi
(
    scalar p, 
    scalar T
) const
{
    label iT = label((T - rhoStartT_) / rhoDeltaT_);
    if (iT > (rhoData_.size()-2))
    {
        return rhoData_[rhoData_.size()-1];
    }
    if (iT < 0)
    {
        return rhoData_[0];
    }
    scalar deltaFrac = (T - rhoStartT_) / rhoDeltaT_ - scalar(iT);
    scalar rho0 = rhoData_[iT] + deltaFrac * (rhoData_[iT+1] - rhoData_[iT]);

    return rho0 / pRef_;
}


template<class Specie>
inline Foam::scalar Foam::fluidLutEOS<Specie>::H
(
    scalar p, 
    scalar T
) const
{
    notImplemented("fluidLutEOS::H");

    return 0;
}


template<class Specie>
inline Foam::scalar Foam::fluidLutEOS<Specie>::Cp
(
    scalar p, 
    scalar T
) const
{
    notImplemented("fluidLutEOS::Cp");

    return 0;
}


template<class Specie>
inline Foam::scalar Foam::fluidLutEOS<Specie>::E
(
    scalar p, 
    scalar T
) const
{
    notImplemented("fluidLutEOS::E");
    
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::fluidLutEOS<Specie>::Cv
(
    scalar p, 
    scalar T
) const
{
    notImplemented("fluidLutEOS::Cv");
    
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::fluidLutEOS<Specie>::S
(
    scalar p, 
    scalar T
) const
{
    notImplemented("fluidLutEOS::S");
    
    return 0;
}


template<class Specie>
inline Foam::scalar Foam::fluidLutEOS<Specie>::Z
(
    scalar p, 
    scalar T
) const
{
    notImplemented("fluidLutEOS::Z");

    return 0;
}


template<class Specie>
inline Foam::scalar Foam::fluidLutEOS<Specie>::CpMCv
(
    scalar p, 
    scalar T
) const
{
    notImplemented("fluidLutEOS::CpMCv");

    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// TODO: No checking is done to see whether the LUTs from both operands are 
// identically keyed, i.e. the same keyset is used across all species.

// TODO: Mixing rules are simplistically based on mass fraction blending - may 
// not be appropriate for some properties.

template<class Specie>
inline void Foam::fluidLutEOS<Specie>::operator=
(
    const fluidLutEOS<Specie>& ifl
)
{
    Specie::operator=(ifl);

    pRef_ = ifl.pRef_;
    rhoStartT_ = ifl.rhoStartT_;
    rhoDeltaT_ = ifl.rhoDeltaT_;
    rhoData_ = ifl.rhoData_;
}

template<class Specie>
inline void Foam::fluidLutEOS<Specie>::operator+=
(
    const fluidLutEOS<Specie>& ifl
)
{
    scalar Y1 = this->Y();
    Specie::operator+=(ifl);
    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = ifl.Y()/this->Y();

        pRef_ = Y1*pRef_ + Y2*ifl.pRef_;
        // rhoStartT_ = rhoStartT_;
        // rhoDeltaT_ = rhoDeltaT_;

        forAll(rhoData_, dataI)
        {
            rhoData_[dataI] = Y1*rhoData_[dataI] + Y2*ifl.rhoData_[dataI];
        }
    }
}


template<class Specie>
inline void Foam::fluidLutEOS<Specie>::operator-=
(
    const fluidLutEOS<Specie>& ifl
)
{
    // NOT IMPLEMENTED YET
}


template<class Specie>
inline void Foam::fluidLutEOS<Specie>::operator*=
(
    const scalar s
)
{
    Specie::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Specie>
inline Foam::fluidLutEOS<Specie> Foam::operator+
(
    const fluidLutEOS<Specie>& ifl1,
    const fluidLutEOS<Specie>& ifl2
)
{
    Specie sp
    (
        static_cast<const Specie&>(ifl1) + static_cast<const Specie&>(ifl2)
    );

    if (mag(sp.Y()) < SMALL)
    {
        return fluidLutEOS<Specie>
        (
            sp,
            ifl1.pRef_,
            ifl1.rhoStartT_,
            ifl1.rhoDeltaT_,
            ifl1.rhoData_
        );
    }
    else
    {
        const scalar Y1 = ifl1.Y()/sp.Y();
        const scalar Y2 = ifl2.Y()/sp.Y();

        List<scalar> fracLut(ifl1.rhoData_.size());
        forAll(ifl2.rhoData_, dataI)
        {
            fracLut[dataI] = Y1*ifl1.rhoData_[dataI] + Y2*ifl2.rhoData_[dataI];
        }

        return fluidLutEOS<Specie>
        (
            sp,
            Y1*ifl1.pRef_ + Y2*ifl2.pRef_,
            ifl1.rhoStartT_,
            ifl1.rhoDeltaT_,
            fracLut
        );
    }
}


template<class Specie>
inline Foam::fluidLutEOS<Specie> Foam::operator-
(
    const fluidLutEOS<Specie>& ipg1,
    const fluidLutEOS<Specie>& ipg2
)
{
    // NOT IMPLEMENTED YET
}


template<class Specie>
inline Foam::fluidLutEOS<Specie> Foam::operator*
(
    const scalar s,
    const fluidLutEOS<Specie>& ifl
)
{
    return fluidLutEOS<Specie>
    (
        s*static_cast<const Specie&>(ifl),
        ifl.pRef_,
        ifl.rhoStartT_,
        ifl.rhoDeltaT_,
        ifl.rhoData_
    );
}


template<class Specie>
inline Foam::fluidLutEOS<Specie> Foam::operator==
(
    const fluidLutEOS<Specie>& pg1,
    const fluidLutEOS<Specie>& pg2
)
{
    // NOT IMPLEMENTED YET
}


// ************************************************************************* //
