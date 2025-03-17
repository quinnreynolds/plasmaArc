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

\*---------------------------------------------------------------------------*/

#include "rhoThermo.H"
#include "makeThermo.H"

#include "makeReactionThermo.H"
#include "rhoReactionThermo.H"

#include "specie.H"
#include "fluidLutThermo.H"
#include "sensibleEnthalpy.H"
#include "sensibleInternalEnergy.H"
#include "thermo.H"

#include "fluidLutTransport.H"

#include "fluidLutEOS.H"

#include "heRhoThermo.H"
#include "pureMixture.H"
#include "multiComponentMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    fluidLutTransport,
    sensibleEnthalpy,
    fluidLutThermo,
    fluidLutEOS,
    specie
);

makeThermos
(
    rhoThermo,
    heRhoThermo,
    pureMixture,
    fluidLutTransport,
    sensibleInternalEnergy,
    fluidLutThermo,
    fluidLutEOS,
    specie
);

makeReactionThermos
(
    rhoThermo,
    rhoReactionThermo,
    heRhoThermo,
    multiComponentMixture,
    fluidLutTransport,
    sensibleEnthalpy,
    fluidLutThermo,
    fluidLutEOS,
    specie
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
