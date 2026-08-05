// TEMPORARY DIAGNOSTIC - not a permanent feature.
// Checks ek, ePot, UcrossB and dotTerm - the four fields that go into
// building ePotEqn (fvm::laplacian(ek, ePot) - fvc::div(ek*UcrossB) +
// dotTerm) - for non-finite (NaN/Inf) values, both internal and boundary,
// immediately before each ePotEqn construction. Aborts with the exact
// field/cell/patch/face pinpointed on the first bad value found, rather
// than letting it propagate to the later, harder-to-interpret SIGFPE in
// the linear solver. Safe to remove once the AMR+load-balancing crash is
// root-caused.
{
    auto checkFieldFinite = [&](const auto& fld)
    {
        forAll(fld.primitiveField(), celli)
        {
            if (!std::isfinite(mag(fld.primitiveField()[celli])))
            {
                FatalErrorInFunction
                    << "Non-finite value in " << fld.name()
                    << " internal field at cell " << celli
                    << ", value=" << fld.primitiveField()[celli]
                    << ", nEM=" << nEM
                    << ", time=" << runTime.timeName()
                    << abort(FatalError);
            }
        }

        const auto& bf = fld.boundaryField();
        forAll(bf, patchi)
        {
            forAll(bf[patchi], facei)
            {
                if (!std::isfinite(mag(bf[patchi][facei])))
                {
                    FatalErrorInFunction
                        << "Non-finite value in " << fld.name()
                        << " boundary field, patch "
                        << mesh.boundary()[patchi].name()
                        << ", face " << facei
                        << ", value=" << bf[patchi][facei]
                        << ", nEM=" << nEM
                        << ", time=" << runTime.timeName()
                        << abort(FatalError);
                }
            }
        }
    };

    checkFieldFinite(ek);
    checkFieldFinite(ePot);
    checkFieldFinite(UcrossB);
    checkFieldFinite(dotTerm);
}
