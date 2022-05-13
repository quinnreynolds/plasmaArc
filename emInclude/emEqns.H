// Solve for EM fields (dynamics)

// TODO: Decoupling approach used here is questionable at best. Revise?

if (runTime.timeIndex() == 1 || runTime.timeIndex() % emSolverFreq == 0)
{
    // Time derivative correction factor
    dtEM += runTime.deltaT().value();
    scalar dtCorrect = runTime.deltaT().value() / dtEM;

    // Explicit source and coupling terms (for stability)
    volVectorField UcrossB = U ^ B;
    volVectorField ekGrad = fvc::grad(ek);
    volScalarField dotTerm = ekGrad & ddtA;

    for (int nEM = 0; nEM < nMaxEMLoops; nEM++) {
        ek.correctBoundaryConditions();

        // Set up solution matrix for ePot scalar field
        fvScalarMatrix ePotEqn
        (
            fvm::laplacian(ek, ePot) - fvc::div(ek * UcrossB) + dotTerm
        );
        ePotEqn.relax();

        // Matrix solve and convergence calc
        scalar residualEPot = ePotEqn.solve().initialResidual();

        // Electric potential current density
        volVectorField jPot = fvc::reconstruct
        (
            -fvc::interpolate(ek) * fvc::snGrad(ePot) * mesh.magSf()
        );

        // Set up solution matrices for A vector field
        fvVectorMatrix AEqn
        (
            mu0 * (ek * (dtCorrect * fvm::ddt(A) - UcrossB) - jPot)
          - fvm::laplacian(A)
          + fvm::Sp(constantMVP, A)
        );
        AEqn.relax();

        // Matrix solve and convergence calc
        Vector<scalar> residualA(AEqn.solve().initialResidual());

        // Magnetic field
        B = fvc::curl(A);

        Vector<scalar> bMag(gMax(B) - gMin(B));
        scalar bMagNorm = max(max(bMag.x(), bMag.y()), bMag.z());
        scalar toleranceAx = toleranceA * bMagNorm / bMag.x();
        scalar toleranceAy = toleranceA * bMagNorm / bMag.y();
        scalar toleranceAz = toleranceA * bMagNorm / bMag.z();

        ddtA = dtCorrect * fvc::ddt(A);

        // Total current density
        j = jPot + ek * (-ddtA + UcrossB);

        // Convergence tests
        if
        (
            residualEPot <= toleranceEPot &&
            residualA.x() <= toleranceAx &&
            residualA.y() <= toleranceAy &&
            residualA.z() <= toleranceAz
        )
        {
            break;
        }
        if (nEM == nMaxEMLoops-1)
        {
            Info << "Warning: Maximum EM loops reached - solution may " <<
                    "not be converged." << endl;
        }
    }
    dtEM = 0.;
}
else
{
    dtEM += runTime.deltaT().value();
}
