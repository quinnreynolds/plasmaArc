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

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::fluidLutThermo<EquationOfState>::fluidLutThermo
(
    const EquationOfState& st,
    const scalar CpStartT,
    const scalar CpDeltaT,
    const List<scalar>& CpData
)
:
    EquationOfState(st),
    CpStartT_(CpStartT),
    CpDeltaT_(CpDeltaT),
    CpData_(CpData)
{
    hData_ = CpData_;
    hData_[0] = CpStartT_ * CpData_[0];
    for (label i=1; i<CpData_.size(); i++)
    {
        hData_[i] = hData_[i-1] + 0.5*(CpData_[i]+CpData_[i-1])*CpDeltaT_;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::fluidLutThermo<EquationOfState>::fluidLutThermo
(
    const word& name,
    const fluidLutThermo& ct
)
:
    EquationOfState(name, ct),
    CpStartT_(ct.CpStartT_),
    CpDeltaT_(ct.CpDeltaT_),
    CpData_(ct.CpData_)
{
    hData_ = CpData_;
    hData_[0] = CpStartT_ * CpData_[0];
    for (label i=1; i<CpData_.size(); i++)
    {
        hData_[i] = hData_[i-1] + 0.5*(CpData_[i]+CpData_[i-1])*CpDeltaT_;
    }
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::fluidLutThermo<EquationOfState> >
Foam::fluidLutThermo<EquationOfState>::clone() const
{
    return autoPtr<fluidLutThermo<EquationOfState> >
    (
        new fluidLutThermo<EquationOfState>(*this)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::fluidLutThermo<EquationOfState> >
Foam::fluidLutThermo<EquationOfState>::New(Istream& is)
{
    return autoPtr<fluidLutThermo<EquationOfState> >
    (
        new fluidLutThermo<EquationOfState>(is)
    );
}


template<class EquationOfState>
inline Foam::autoPtr<Foam::fluidLutThermo<EquationOfState> >
Foam::fluidLutThermo<EquationOfState>::New(const dictionary& dict)
{
    return autoPtr<fluidLutThermo<EquationOfState> >
    (
        new fluidLutThermo<EquationOfState>(dict)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::scalar Foam::fluidLutThermo<EquationOfState>::limit
(
    const scalar T
) const
{
    return T;
}


template<class EquationOfState>
inline Foam::scalar Foam::fluidLutThermo<EquationOfState>::Cp
(
    const scalar p,
    const scalar T
) const
{
    label iT = label((T - CpStartT_) / CpDeltaT_);
    if (iT > (CpData_.size()-2))
    {
        return CpData_[CpData_.size()-1];
    }
    if (iT < 0)
    {
        return CpData_[0];
    }
    scalar deltaFrac = (T - CpStartT_) / CpDeltaT_ - scalar(iT);

    return (CpData_[iT] + deltaFrac * (CpData_[iT+1] - CpData_[iT]));
}


template<class EquationOfState>
inline Foam::scalar Foam::fluidLutThermo<EquationOfState>::Hs
(
    const scalar p, 
    const scalar T
) const
{
    label iT = label((T - CpStartT_) / CpDeltaT_);
    if (iT > (hData_.size()-2))
    {
        label iN = hData_.size()-1;
        scalar TN = CpStartT_ + scalar(iN)*CpDeltaT_;
        return hData_[iN] + CpData_[iN]*(T-TN);
    }
    if (iT < 0)
    {
        return hData_[0] + CpData_[0]*(T-CpStartT_);
    }
    scalar Ti = CpStartT_ + scalar(iT)*CpDeltaT_;
    scalar quadCoeff = 0.5*(CpData_[iT+1] - CpData_[iT]) / CpDeltaT_;
    return (hData_[iT] + CpData_[iT]*(T-Ti) + quadCoeff*(T-Ti)*(T-Ti));
}


template<class EquationOfState>
inline Foam::scalar Foam::fluidLutThermo<EquationOfState>::Ha
(
    const scalar p, 
    const scalar T
) const
{
    notImplemented("fluidLutThermo::Ha");

    return 0;
}


template<class EquationOfState>
inline Foam::scalar Foam::fluidLutThermo<EquationOfState>::Hc() const
{
    notImplemented("fluidLutThermo::Hc");

    return 0.;
}


template<class EquationOfState>
inline Foam::scalar Foam::fluidLutThermo<EquationOfState>::S
(
    const scalar p, 
    const scalar T
) const
{
    notImplemented("fluidLutThermo::S");

    return 0;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// TODO: No checking is done to see whether the LUTs from both operands are 
// identically keyed, i.e. the same keyset is used across all species.

// TODO: Mixing rules are simplistically based on mass fraction blending - may 
// not be appropriate for some properties.

template<class EquationOfState>
inline void Foam::fluidLutThermo<EquationOfState>::operator+=
(
    const fluidLutThermo<EquationOfState>& ifl
)
{
    scalar Y1 = this->Y();
    EquationOfState::operator+=(ifl);
    if (mag(this->Y()) > SMALL)
    {
        Y1 /= this->Y();
        const scalar Y2 = ifl.Y()/this->Y();

        forAll(CpData_, dataI)
        {
            CpData_[dataI] = Y1*CpData_[dataI] + Y2*ifl.CpData_[dataI];
        }
        forAll(hData_, dataI)
        {
            hData_[dataI] = Y1*hData_[dataI] + Y2*ifl.hData_[dataI];
        }
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class EquationOfState>
inline Foam::fluidLutThermo<EquationOfState> Foam::operator+
(
    const fluidLutThermo<EquationOfState>& ifl1,
    const fluidLutThermo<EquationOfState>& ifl2
)
{
    EquationOfState eofs
    (
        static_cast<const EquationOfState&>(ifl1)
      + static_cast<const EquationOfState&>(ifl2)
    );

    if (mag(eofs.Y()) < SMALL)
    {
        return fluidLutThermo<EquationOfState>
        (
            eofs,
            ifl1.CpStartT_,
            ifl1.CpDeltaT_,
            ifl1.CpData_
        );
    }
    {
        const scalar Y1 = ifl1.Y()/eofs.Y();
        const scalar Y2 = ifl2.Y()/eofs.Y();

        List<scalar> fracLutCp(ifl1.CpData_.size());
        forAll(ifl1.CpData_, dataI)
        {
            fracLutCp[dataI] = Y1*ifl1.CpData_[dataI] + Y2*ifl2.CpData_[dataI];
        }
        List<scalar> fracLutH(ifl1.hData_.size());
        forAll(ifl1.hData_, dataI)
        {
            fracLutH[dataI] = Y1*ifl1.hData_[dataI] + Y2*ifl2.hData_[dataI];
        }

        return fluidLutThermo<EquationOfState>
        (
            eofs,
            ifl1.CpStartT_,
            ifl1.CpDeltaT_,
            fracLutCp
        );
    }
}


template<class EquationOfState>
inline Foam::fluidLutThermo<EquationOfState> Foam::operator-
(
    const fluidLutThermo<EquationOfState>& ifl1,
    const fluidLutThermo<EquationOfState>& ifl2
)
{
    // NOT IMPLEMENTED YET
}


template<class EquationOfState>
inline Foam::fluidLutThermo<EquationOfState> Foam::operator*
(
    const scalar s,
    const fluidLutThermo<EquationOfState>& ifl
)
{
    return fluidLutThermo<EquationOfState>
    (
        s*static_cast<const EquationOfState&>(ifl),
        ifl.CpStartT_,
        ifl.CpDeltaT_,
        ifl.CpData_
    );
}


template<class EquationOfState>
inline Foam::fluidLutThermo<EquationOfState> Foam::operator==
(
    const fluidLutThermo<EquationOfState>& ifl1,
    const fluidLutThermo<EquationOfState>& ifl2
)
{
    // NOT IMPLEMENTED YET
}


// ************************************************************************* //
