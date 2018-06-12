{
    volScalarField& he = thermo.he();
    volScalarField& Cp = thermo.Cp();

    volScalarField qSourceSp = -ekSp * (j & j) / (ek * ek * Cp);
    volScalarField qSourceSc = ((j & j) / ek - qSourceSp * T) / Cp;

    surfaceScalarField electronTransportPhi
        = (constantET / Cp) * linearInterpolate(j) & mesh.Sf();

    Info<< "Getting here 1" << nl;

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he)
      + fvm::div(phi, he)
      + fvm::div(electronTransportPhi, he)
      + fvc::ddt(rho, K)
      + fvc::div(phi, K)
      + ( he.name() == "e"
            ? fvc::div(fvc::absolute(phi / fvc::interpolate(rho), U), p)
            : -dpdt
        )
      ==
        fvm::laplacian(turbulence->alphaEff(), he)
      + fvm::Sp(qSourceSp, he)
      + qSourceSc
      + fvOptions(rho, he)
    );

    Info<< "Getting here 2" << nl;

    EEqn.relax();

    Info<< "Getting here 3" << nl;

    fvOptions.constrain(EEqn);

    Info<< "Getting here 4" << nl;

    EEqn.solve();

    Info<< "Getting here 5" << nl;

    fvOptions.correct(he);

    Info<< "Getting here 6" << nl;

    //he = thermo.he(p, T);

    Info<< "Getting here 7" << nl;

    thermo.correct();
}
