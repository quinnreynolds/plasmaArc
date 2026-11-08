{
    volScalarField& he = thermo.he();

    volScalarField qSourceSp(-ekSp * (j & j) / (ek * ek * thermo.Cp()));
    volScalarField qSourceSc((j & j) / ek - qSourceSp * he);

    // Electrons carry enthalpy against the conventional current j, so the
    // transporting flux is -j/Cp. The sign is kept in the flux (not on the
    // matrix) so the convection scheme upwinds along the true transport
    // direction.
    surfaceScalarField electronTransportPhi
    (
        "electronTransportPhi",
      - constantET * linearInterpolate(j / thermo.Cp()) & mesh.Sf()
    );

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, he)
      + fvm::div(phi, he)
      + fvm::div(electronTransportPhi, he)
      + fvc::ddt(rho, K)
      + fvc::div(phi, K)
      + (
            he.name() == "e"
          ? fvc::div
            (
                fvc::absolute(phi/fvc::interpolate(rho), U),
                p,
                "div(phiv,p)"
            )
          : -dpdt
        )
      - fvm::laplacian(turbulence->alphaEff(), he)
      ==
        fvm::Sp(qSourceSp, he)
      + qSourceSc
      + fvOptions(rho, he)
    );

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    // tk = thermo.kappa();
    tk = turbulence->kappaEff();
}
