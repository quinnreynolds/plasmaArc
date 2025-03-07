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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::fluidLutTransport<Thermo>::fluidLutTransport
(
    const Thermo& t,
    const scalar muStartT,
    const scalar muDeltaT,
    const List<scalar>& muData,
    const scalar tkStartT,
    const scalar tkDeltaT,
    const List<scalar>& tkData
)
:
    Thermo(t),
    muStartT_(muStartT),
    muDeltaT_(muDeltaT),
    muData_(muData),
    tkStartT_(tkStartT),
    tkDeltaT_(tkDeltaT),
    tkData_(tkData)
{}


template<class Thermo>
inline Foam::fluidLutTransport<Thermo>::fluidLutTransport
(
    const word& name,
    const fluidLutTransport& ct
)
:
    Thermo(name, ct),
    muStartT_(ct.muStartT_),
    muDeltaT_(ct.muDeltaT_),
    muData_(ct.muData_),
    tkStartT_(ct.tkStartT_),
    tkDeltaT_(ct.tkDeltaT_),
    tkData_(ct.tkData_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::fluidLutTransport<Thermo> >
Foam::fluidLutTransport<Thermo>::clone() const
{
    return autoPtr<fluidLutTransport<Thermo> >
    (
        new fluidLutTransport<Thermo>(*this)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::fluidLutTransport<Thermo> >
Foam::fluidLutTransport<Thermo>::New
(
    Istream& is
)
{
    return autoPtr<fluidLutTransport<Thermo> >
    (
        new fluidLutTransport<Thermo>(is)
    );
}


template<class Thermo>
inline Foam::autoPtr<Foam::fluidLutTransport<Thermo> >
Foam::fluidLutTransport<Thermo>::New
(
    const dictionary& dict
)
{
    return autoPtr<fluidLutTransport<Thermo> >
    (
        new fluidLutTransport<Thermo>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::fluidLutTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    label iT = label((T - muStartT_) / muDeltaT_);
    if (iT > (muData_.size()-2))
    {
        return muData_[muData_.size()-1];
    }
    if (iT < 0)
    {
        return muData_[0];
    }
    scalar deltaFrac = (T - muStartT_) / muDeltaT_ - scalar(iT);

    return muData_[iT] + deltaFrac * (muData_[iT+1] - muData_[iT]);
}


template<class Thermo>
inline Foam::scalar Foam::fluidLutTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    label iT = label((T - tkStartT_) / tkDeltaT_);
    if (iT > (tkData_.size()-2))
    {
        return tkData_[tkData_.size()-1];
    }
    if (iT < 0)
    {
        return tkData_[0];
    }
    scalar deltaFrac = (T - tkStartT_) / tkDeltaT_ - scalar(iT);

    return tkData_[iT] + deltaFrac * (tkData_[iT+1] - tkData_[iT]);
}


template<class Thermo>
inline Foam::scalar Foam::fluidLutTransport<Thermo>::alphah
(
    const scalar p,
    const scalar T
) const
{
    return kappa(p, T) / this->Cp(p, T);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// TODO: No checking is done to see whether the LUTs from both operands are 
// identically keyed, i.e. the same keyset is used across all species.

// TODO: Mixing rules are simplistically based on mass fraction blending - may 
// not be appropriate for some properties.

template<class Thermo>
inline void Foam::fluidLutTransport<Thermo>::operator=
(
    const fluidLutTransport<Thermo>& ifl
)
{
    Thermo::operator=(ifl);

    muStartT_ = ifl.muStartT_;
    muDeltaT_ = ifl.muDeltaT_;
    muData_ = ifl.muData_;
    tkStartT_ = ifl.tkStartT_;
    tkDeltaT_ = ifl.tkDeltaT_;
    tkData_ = ifl.tkData_;
}


template<class Thermo>
inline void Foam::fluidLutTransport<Thermo>::operator+=
(
    const fluidLutTransport<Thermo>& ifl
)
{
    scalar Y1 = this->Y();
    Thermo::operator+=(ifl);

    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        scalar Y2 = ifl.Y()/this->Y();

        forAll(muData_, dataI)
        {
            muData_[dataI] = Y1*muData_[dataI] + Y2*ifl.muData_[dataI];
        }
        forAll(tkData_, dataI)
        {
            tkData_[dataI] = Y1*tkData_[dataI] + Y2*ifl.tkData_[dataI];
        }
    }
}


template<class Thermo>
inline void Foam::fluidLutTransport<Thermo>::operator-=
(
    const fluidLutTransport<Thermo>& ifl
)
{
    // NOT IMPLEMENTED YET
}


template<class Thermo>
inline void Foam::fluidLutTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::fluidLutTransport<Thermo> Foam::operator+
(
    const fluidLutTransport<Thermo>& ifl1,
    const fluidLutTransport<Thermo>& ifl2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(ifl1) + static_cast<const Thermo&>(ifl2)
    );

    if (mag(t.Y()) < SMALL)
    {
        return fluidLutTransport<Thermo>
        (
            t,
            ifl1.muStartT_,
            ifl1.muDeltaT_,
            ifl1.muData_,
            ifl1.tkStartT_,
            ifl1.tkDeltaT_,
            ifl1.tkData_
        );
    }
    else
    {
        scalar Y1 = ifl1.Y()/t.Y();
        scalar Y2 = ifl2.Y()/t.Y();

        List<scalar> fracLutMu(ifl1.muData_.size());
        forAll(ifl1.muData_, dataI)
        {
            fracLutMu[dataI] = Y1*ifl1.muData_[dataI] + Y2*ifl2.muData_[dataI];
        }
        List<scalar> fracLutTk(ifl1.tkData_.size());
        forAll(ifl1.tkData_, dataI)
        {
            fracLutTk[dataI] = Y1*ifl1.tkData_[dataI] + Y2*ifl2.tkData_[dataI];
        }

        return fluidLutTransport<Thermo>
        (
            t,
            ifl1.muStartT_,
            ifl1.muDeltaT_,
            fracLutMu,
            ifl1.tkStartT_,
            ifl1.tkDeltaT_,
            fracLutTk
        );
    }
}


template<class Thermo>
inline Foam::fluidLutTransport<Thermo> Foam::operator-
(
    const fluidLutTransport<Thermo>& ifl1,
    const fluidLutTransport<Thermo>& ifl2
)
{
    // NOT IMPLEMENTED YET
}


template<class Thermo>
inline Foam::fluidLutTransport<Thermo> Foam::operator*
(
    const scalar s,
    const fluidLutTransport<Thermo>& ifl
)
{
    return fluidLutTransport<Thermo>
    (
        s*static_cast<const Thermo&>(ifl),
        ifl.muStartT_,
        ifl.muDeltaT_,
        ifl.muData_,
        ifl.tkStartT_,
        ifl.tkDeltaT_,
        ifl.tkData_
    );
}


template<class Thermo>
inline Foam::fluidLutTransport<Thermo> Foam::operator==
(
    const fluidLutTransport<Thermo>& ct1,
    const fluidLutTransport<Thermo>& ct2
)
{
    return ct2 - ct1;
}


// ************************************************************************* //
