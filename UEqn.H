volVectorField USource = j ^ B;

tmp<fvVectorMatrix> UEqn
(
    fvm::ddt(rho, U)
        + fvm::div(phi, U)
        + turbulence->divDevRhoReff(U)
        - USource
    ==
        fvOptions(rho, U)

);

fvVectorMatrix& tUEqn = UEqn.ref();

tUEqn.relax();

if (pimple.momentumPredictor())
{
    solve(UEqn() == -fvc::grad(p));

    fvOptions.correct(U);

    K = 0.5*magSqr(U);
}
