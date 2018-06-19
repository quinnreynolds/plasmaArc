{
    volScalarField& he = thermo.he();

    volScalarField qSourceSp = -ekSp * (j & j) / (ek * ek);
    volScalarField qSourceSc = (j & j) / ek - qSourceSp * T;

    surfaceScalarField electronTransportPhi
        = constantET * linearInterpolate(j) & mesh.Sf();

    Info<< "Getting here 1" << nl;

    fvScalarMatrix TEqn
    (
        thermo.Cp() * (
            fvm::ddt(rho, T)
          + fvm::div(phi, T)
        )
      + fvm::div(electronTransportPhi, T)
      + fvc::ddt(rho, K)
      + fvc::div(phi, K)
      + ( he.name() == "e"
            ? fvc::div(fvc::absolute(phi / fvc::interpolate(rho), U), p)
            : -dpdt
        )
      ==
        fvm::laplacian(turbulence->alphaEff()*thermo.Cp(), T)
      + fvm::Sp(qSourceSp, T)
      + qSourceSc
      + fvOptions(rho, he)
    );

    Info<< "Getting here 2" << nl;

    TEqn.relax();

    Info<< "Getting here 3" << nl;

    fvOptions.constrain(TEqn);

    Info<< "Getting here 4" << nl;

    TEqn.solve();

    Info<< "Getting here 5" << nl;

    fvOptions.correct(T);

    Info<< "Getting here 6" << nl;

    he = thermo.he(p, T);

    Info<< "Getting here 7" << nl;

    thermo.correct();
}
