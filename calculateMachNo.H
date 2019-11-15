volScalarField magU = mag(U);
scalar adiabaticIndex = 1;

MachNo = magU / (Foam::sqrt(adiabaticIndex * p / rho));
