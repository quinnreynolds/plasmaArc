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

Application
    plasmaArcFoam

Description
    Transient MHD solver for laminar or turbulent flow of compressible or
    incompressible, non-isochoric plasma fluids.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient simulations.

    Based on rhoPimpleFoam and myRhoPimpleFoam.


Q Reynolds Feb 2015

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "rhoThermo.H"
#include "turbulentFluidThermoModel.H"
#include "bound.H"
#include "pimpleControl.H"
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"

#include "directionMixedFvPatchFields.H"
#include "zeroGradientFvPatchField.H"
#include "interpolationTable.H"

#include "scalarLookup.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"

    pimpleControl pimple(mesh);

    #include "createTimeControls.H"
    #include "createRDeltaT.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "eminclude/createFields.H"
    #include "createMRF.H"
    #include "eminclude/readSolverControls.H"
    #include "createFvOptions.H"

    Info<< "\nInitialising surface normals and fraction tensor BCs for A...\n"
        << endl;

    #include "eminclude/setDirectionMixedBC.H"

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        //update EM transport fields and solve
        #include "ekQradCalculate.H"
        #include "eminclude/emEqns.H"

        if (pimple.nCorrPIMPLE() <= 1)
        {
            #include "rhoEqn.H"
        }

        //Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"
            #include "TEqn.H"

            //Pressure corrector loop
            while (pimple.correct())
            {
                if (pimple.consistent())
                {
                    #include "pcEqn.H"
                }
                else
                {
                    #include "pEqn.H"
                }
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        Info<< "Pressure min/max = " << min(p).value() << " / "
            << max(p).value() << " Pa" << nl;

        Info<< "Voltage = " << max(ePot).value() << " V, "
            << "Temperature = " << max(T).value() << " K, "
            << "|U| = " << max(mag(U)).value() << " m/s"
            << nl;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
