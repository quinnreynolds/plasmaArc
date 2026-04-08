/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    plasmaArc

Description
    Transient MHD solver for laminar or turbulent flow of compressible or
    incompressible plasma fluids.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient simulations.

    Based on rhoPimpleFoam.

\*---------------------------------------------------------------------------*/

// Select PIMPLE control for postProcess.H -> createControl.H dispatch
#define PIMPLE_CONTROL

#include "Time.H"
#include "fvMesh.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "argList.H"
#include "timeSelector.H"
#include "fluidThermo.H"
#include "compressibleMomentumTransportModels.H"
#include "fluidThermoThermophysicalTransportModel.H"
#include "bound.H"
#include "pimpleControl.H"
#include "constrainHbyA.H"
#include "constrainPressure.H"
#include "CorrectPhi.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "directionMixedFvPatchFields.H"
#include "zeroGradientFvPatchField.H"

// Individual fvm operator headers (replaces fvCFD.H)
#include "fvmDdt.H"
#include "fvmDiv.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"

using namespace Foam;

#include "scalarLookup.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createPimpleControl.H"
    const bool LTS = false;
    const bool correctPhi = false;
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "emInclude/createFields.H"
    #include "emInclude/readSolverControls.H"
    autoPtr<surfaceVectorField> rhoUf;
    #include "calculateCompositions.H"

    Info<< "\nInitialising surface normals and fraction tensor BCs for A...\n"
        << endl;

    #include "emInclude/setDirectionMixedBC.H"

    turbulence->validate();

    #include "createTimeControls.H"
    #include "compressibleCourantNo.H"
    #include "setInitialDeltaT.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        // Store divrhoU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        autoPtr<volScalarField> divrhoU;
        if (correctPhi)
        {
            divrhoU.reset
            (
                new volScalarField
                (
                    "divrhoU",
                    fvc::div(fvc::absolute(phi, rho, U))
                )
            );
        }

        #include "readTimeControls.H"
        #include "compressibleCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.name() << nl << endl;

        //update EM transport fields and solve
        #include "calculateEk.H"
        #include "emInclude/emEqns.H"

        //Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstPimpleIter())
            {
                #include "rhoEqn.H"
            }

            #include "UEqn.H"
            #include "EEqn.H"

            //Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.finalPimpleIter())
            {
                turbulence->correct();
                thermophysicalTransport->correct();
            }
        }

        rho = thermo.rho();

        #include "calculateMachNo.H"

        #include "calculateCompositions.H"

        runTime.write();

        Info<< "Voltage = " << gMin(ePot) << "/" << gMax(ePot) << " V, "
            << "Temperature = " << gMax(T) << " K, "
            << "|U| = " << gMax(magU) << " m/s, "
            << "Ma = " << gMax(MachNo) << nl;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
