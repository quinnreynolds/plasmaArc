// TEMPORARY DIAGNOSTIC - not a permanent feature.
// Checks ePotEqn's diag()/upper()/lower()/source()/internalCoeffs()/
// boundaryCoeffs() for non-finite (NaN/Inf) values. Included twice -
// once right after construction, once right after relax() - with
// checkEPotEqnStage set beforehand to say which. Collects every anomaly
// found (rather than aborting on the first) so we can tell whether this
// is a single isolated defect or something more widespread - printing
// detail is capped per component at checkMaxPrint, but the true total
// count is always reported. Also warns (without counting as an error)
// on a near-zero diagonal entry, since that alone is enough to
// destabilise PCG even without any value actually being non-finite. For
// source() specifically, dumps full forensic detail (see
// dumpEPotCellForensics.H) for the first few bad cells. Aborts with a
// summary once the full scan across all six components is complete.
// Safe to remove once the AMR+load-balancing crash is root-caused.
{
    const label checkMaxPrint = 20;
    const label checkForensicDumpMax = 3;
    label checkTotalBad = 0;
    label checkForensicDumpCount = 0;

    if (ePotEqn.hasDiag())
    {
        const scalarField& diag = ePotEqn.diag();
        label nBad = 0;
        forAll(diag, celli)
        {
            if (!std::isfinite(diag[celli]))
            {
                ++nBad;
                if (nBad <= checkMaxPrint)
                {
                    Pout<< "[checkEPotEqnFinite] Non-finite diag() ["
                        << checkEPotEqnStage << "] at cell " << celli
                        << ", value=" << diag[celli]
                        << ", nEM=" << nEM
                        << ", time=" << runTime.timeName()
                        << endl;
                }
            }
            else if (mag(diag[celli]) < SMALL)
            {
                Pout<< "[checkEPotEqnFinite] WARNING near-zero diag ["
                    << checkEPotEqnStage << "] "
                    << "proc=" << Pstream::myProcNo()
                    << " cell=" << celli
                    << " value=" << diag[celli]
                    << " nEM=" << nEM
                    << " time=" << runTime.timeName()
                    << endl;
            }
        }
        if (nBad > checkMaxPrint)
        {
            Pout<< "[checkEPotEqnFinite] ... and " << (nBad - checkMaxPrint)
                << " more non-finite diag() values" << endl;
        }
        checkTotalBad += nBad;
    }

    if (ePotEqn.hasUpper())
    {
        const scalarField& upper = ePotEqn.upper();
        label nBad = 0;
        forAll(upper, i)
        {
            if (!std::isfinite(upper[i]))
            {
                ++nBad;
                if (nBad <= checkMaxPrint)
                {
                    Pout<< "[checkEPotEqnFinite] Non-finite upper() ["
                        << checkEPotEqnStage << "] at index " << i
                        << ", value=" << upper[i]
                        << ", nEM=" << nEM
                        << ", time=" << runTime.timeName()
                        << endl;
                }
            }
        }
        if (nBad > checkMaxPrint)
        {
            Pout<< "[checkEPotEqnFinite] ... and " << (nBad - checkMaxPrint)
                << " more non-finite upper() values" << endl;
        }
        checkTotalBad += nBad;
    }

    if (ePotEqn.hasLower())
    {
        const scalarField& lower = ePotEqn.lower();
        label nBad = 0;
        forAll(lower, i)
        {
            if (!std::isfinite(lower[i]))
            {
                ++nBad;
                if (nBad <= checkMaxPrint)
                {
                    Pout<< "[checkEPotEqnFinite] Non-finite lower() ["
                        << checkEPotEqnStage << "] at index " << i
                        << ", value=" << lower[i]
                        << ", nEM=" << nEM
                        << ", time=" << runTime.timeName()
                        << endl;
                }
            }
        }
        if (nBad > checkMaxPrint)
        {
            Pout<< "[checkEPotEqnFinite] ... and " << (nBad - checkMaxPrint)
                << " more non-finite lower() values" << endl;
        }
        checkTotalBad += nBad;
    }

    {
        const scalarField& src = ePotEqn.source();
        label nBad = 0;
        forAll(src, celli)
        {
            if (!std::isfinite(src[celli]))
            {
                ++nBad;
                if (nBad <= checkMaxPrint)
                {
                    Pout<< "[checkEPotEqnFinite] Non-finite source() ["
                        << checkEPotEqnStage << "] at cell " << celli
                        << ", value=" << src[celli]
                        << ", nEM=" << nEM
                        << ", time=" << runTime.timeName()
                        << endl;
                }
                if (checkForensicDumpCount < checkForensicDumpMax)
                {
                    ++checkForensicDumpCount;
                    #include "dumpEPotCellForensics.H"
                }
            }
        }
        if (nBad > checkMaxPrint)
        {
            Pout<< "[checkEPotEqnFinite] ... and " << (nBad - checkMaxPrint)
                << " more non-finite source() values" << endl;
        }
        checkTotalBad += nBad;
    }

    forAll(ePotEqn.internalCoeffs(), patchi)
    {
        const scalarField& intCoeffs = ePotEqn.internalCoeffs()[patchi];
        label nBadInt = 0;
        forAll(intCoeffs, facei)
        {
            if (!std::isfinite(intCoeffs[facei]))
            {
                ++nBadInt;
                if (nBadInt <= checkMaxPrint)
                {
                    Pout<< "[checkEPotEqnFinite] Non-finite internalCoeffs() ["
                        << checkEPotEqnStage << "] patch "
                        << mesh.boundary()[patchi].name()
                        << " face " << facei
                        << ", value=" << intCoeffs[facei]
                        << ", nEM=" << nEM
                        << ", time=" << runTime.timeName()
                        << endl;
                }
            }
        }
        if (nBadInt > checkMaxPrint)
        {
            Pout<< "[checkEPotEqnFinite] ... and " << (nBadInt - checkMaxPrint)
                << " more non-finite internalCoeffs() values on patch "
                << mesh.boundary()[patchi].name() << endl;
        }
        checkTotalBad += nBadInt;

        const scalarField& bCoeffs = ePotEqn.boundaryCoeffs()[patchi];
        label nBadB = 0;
        forAll(bCoeffs, facei)
        {
            if (!std::isfinite(bCoeffs[facei]))
            {
                ++nBadB;
                if (nBadB <= checkMaxPrint)
                {
                    Pout<< "[checkEPotEqnFinite] Non-finite boundaryCoeffs() ["
                        << checkEPotEqnStage << "] patch "
                        << mesh.boundary()[patchi].name()
                        << " face " << facei
                        << ", value=" << bCoeffs[facei]
                        << ", nEM=" << nEM
                        << ", time=" << runTime.timeName()
                        << endl;
                }
            }
        }
        if (nBadB > checkMaxPrint)
        {
            Pout<< "[checkEPotEqnFinite] ... and " << (nBadB - checkMaxPrint)
                << " more non-finite boundaryCoeffs() values on patch "
                << mesh.boundary()[patchi].name() << endl;
        }
        checkTotalBad += nBadB;
    }

    if (checkTotalBad > 0)
    {
        FatalErrorInFunction
            << "Found " << checkTotalBad
            << " non-finite ePotEqn coefficient(s) ["
            << checkEPotEqnStage << "]"
            << " (proc=" << Pstream::myProcNo() << ")"
            << ", nEM=" << nEM
            << ", time=" << runTime.timeName()
            << abort(FatalError);
    }
}
