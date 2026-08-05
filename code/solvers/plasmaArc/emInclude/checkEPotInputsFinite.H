// TEMPORARY DIAGNOSTIC - not a permanent feature.
// Checks ek, ePot, UcrossB and dotTerm - the four fields that go into
// building ePotEqn (fvm::laplacian(ek, ePot) - fvc::div(ek*UcrossB) +
// dotTerm) - for non-finite (NaN/Inf) values, both internal and boundary,
// immediately before each ePotEqn construction. Collects every anomaly
// found (rather than aborting on the first) so we can tell whether this
// is a single isolated defect or something more widespread - printing
// detail is capped per field at checkFiniteMaxPrint, but the true total
// count is always reported. Aborts with a summary once the full scan
// across all four fields is complete. Safe to remove once the
// AMR+load-balancing crash is root-caused.
{
    const label checkFiniteMaxPrint = 20;
    label checkFiniteTotalBad = 0;

    auto checkFieldFinite = [&](const auto& fld)
    {
        label nBad = 0;

        forAll(fld.primitiveField(), celli)
        {
            if (!std::isfinite(mag(fld.primitiveField()[celli])))
            {
                ++nBad;
                if (nBad <= checkFiniteMaxPrint)
                {
                    Pout<< "[checkEPotInputsFinite] Non-finite " << fld.name()
                        << " internal field at cell " << celli
                        << ", value=" << fld.primitiveField()[celli]
                        << ", nEM=" << nEM
                        << ", time=" << runTime.timeName()
                        << endl;
                }
            }
        }

        const auto& bf = fld.boundaryField();
        forAll(bf, patchi)
        {
            forAll(bf[patchi], facei)
            {
                if (!std::isfinite(mag(bf[patchi][facei])))
                {
                    ++nBad;
                    if (nBad <= checkFiniteMaxPrint)
                    {
                        Pout<< "[checkEPotInputsFinite] Non-finite "
                            << fld.name()
                            << " boundary field, patch "
                            << mesh.boundary()[patchi].name()
                            << ", face " << facei
                            << ", value=" << bf[patchi][facei]
                            << ", nEM=" << nEM
                            << ", time=" << runTime.timeName()
                            << endl;
                    }
                }
            }
        }

        if (nBad > checkFiniteMaxPrint)
        {
            Pout<< "[checkEPotInputsFinite] ... and "
                << (nBad - checkFiniteMaxPrint)
                << " more non-finite " << fld.name()
                << " values (proc=" << Pstream::myProcNo() << ")"
                << endl;
        }

        checkFiniteTotalBad += nBad;
    };

    checkFieldFinite(ek);
    checkFieldFinite(ePot);
    checkFieldFinite(UcrossB);
    checkFieldFinite(dotTerm);

    if (checkFiniteTotalBad > 0)
    {
        FatalErrorInFunction
            << "Found " << checkFiniteTotalBad
            << " non-finite value(s) across ek/ePot/UcrossB/dotTerm"
            << " (proc=" << Pstream::myProcNo() << ")"
            << ", nEM=" << nEM
            << ", time=" << runTime.timeName()
            << abort(FatalError);
    }
}
