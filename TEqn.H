{
    volScalarField& he = thermo.he();

    volScalarField qSourceSp = -ekSp * (j & j) / (ek * ek);
    volScalarField qSourceSc = (j & j) / ek - qSourceSp * T;

    surfaceScalarField electronTransportPhi
        = constantET * linearInterpolate(j) & mesh.Sf();

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
        Qdot
      + fvm::laplacian(turbulence->alphaEff()*thermo.Cp(), T)
      + fvm::Sp(qSourceSp, T)
      + qSourceSc
      + fvOptions(rho*thermo.Cp(), T)
    );

    TEqn.relax();

    fvOptions.constrain(TEqn);

    TEqn.solve();

    fvOptions.correct(T);

    he = thermo.he(p, T);

    thermo.correct();
}
