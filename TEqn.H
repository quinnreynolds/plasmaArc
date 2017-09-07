{
    volScalarField& he = thermo.he();

    volScalarField qSourceSp = -ekSp * (j & j) / (ek * ek) - QRadSp;
    volScalarField qSourceSc = (j & j) / ek - QRad - qSourceSp * T;

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
      //- fvm::laplacian(thermo.alphaEff(turbulence->mut())*thermo.Cp(), T)
      - fvm::laplacian(turbulence->alphaEff()*thermo.Cp(), T)
      - fvm::Sp(qSourceSp, T)
      - qSourceSc
    );

    Info<< max(thermo.Cp()) << nl;

    TEqn.relax();

    TEqn.solve();

    he = thermo.he(p, T);

    thermo.correct();
}
