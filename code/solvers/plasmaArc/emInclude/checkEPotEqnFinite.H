// TEMPORARY DIAGNOSTIC - not a permanent feature.
// Checks ePotEqn's diag()/upper()/lower()/source()/internalCoeffs()/
// boundaryCoeffs() for non-finite (NaN/Inf) values. Included twice -
// once right after construction, once right after relax() - with
// checkEPotEqnStage set beforehand to say which. Aborts with the exact
// location pinpointed on the first bad value found. Also warns (without
// aborting) on a near-zero diagonal entry, since that alone is enough to
// destabilise PCG even without any value actually being non-finite.
// Safe to remove once the AMR+load-balancing crash is root-caused.
{
    if (ePotEqn.hasDiag())
    {
        const scalarField& diag = ePotEqn.diag();
        forAll(diag, celli)
        {
            if (!std::isfinite(diag[celli]))
            {
                FatalErrorInFunction
                    << "Non-finite ePotEqn.diag() [" << checkEPotEqnStage
                    << "] at cell " << celli
                    << ", value=" << diag[celli]
                    << ", nEM=" << nEM
                    << ", time=" << runTime.timeName()
                    << abort(FatalError);
            }
            if (mag(diag[celli]) < SMALL)
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
    }

    if (ePotEqn.hasUpper())
    {
        const scalarField& upper = ePotEqn.upper();
        forAll(upper, i)
        {
            if (!std::isfinite(upper[i]))
            {
                FatalErrorInFunction
                    << "Non-finite ePotEqn.upper() [" << checkEPotEqnStage
                    << "] at index " << i
                    << ", value=" << upper[i]
                    << ", nEM=" << nEM
                    << ", time=" << runTime.timeName()
                    << abort(FatalError);
            }
        }
    }

    if (ePotEqn.hasLower())
    {
        const scalarField& lower = ePotEqn.lower();
        forAll(lower, i)
        {
            if (!std::isfinite(lower[i]))
            {
                FatalErrorInFunction
                    << "Non-finite ePotEqn.lower() [" << checkEPotEqnStage
                    << "] at index " << i
                    << ", value=" << lower[i]
                    << ", nEM=" << nEM
                    << ", time=" << runTime.timeName()
                    << abort(FatalError);
            }
        }
    }

    {
        const scalarField& src = ePotEqn.source();
        forAll(src, celli)
        {
            if (!std::isfinite(src[celli]))
            {
                FatalErrorInFunction
                    << "Non-finite ePotEqn.source() [" << checkEPotEqnStage
                    << "] at cell " << celli
                    << ", value=" << src[celli]
                    << ", nEM=" << nEM
                    << ", time=" << runTime.timeName()
                    << abort(FatalError);
            }
        }
    }

    forAll(ePotEqn.internalCoeffs(), patchi)
    {
        const scalarField& intCoeffs = ePotEqn.internalCoeffs()[patchi];
        forAll(intCoeffs, facei)
        {
            if (!std::isfinite(intCoeffs[facei]))
            {
                FatalErrorInFunction
                    << "Non-finite ePotEqn.internalCoeffs() ["
                    << checkEPotEqnStage << "] patch "
                    << mesh.boundary()[patchi].name()
                    << " face " << facei
                    << ", value=" << intCoeffs[facei]
                    << ", nEM=" << nEM
                    << ", time=" << runTime.timeName()
                    << abort(FatalError);
            }
        }

        const scalarField& bCoeffs = ePotEqn.boundaryCoeffs()[patchi];
        forAll(bCoeffs, facei)
        {
            if (!std::isfinite(bCoeffs[facei]))
            {
                FatalErrorInFunction
                    << "Non-finite ePotEqn.boundaryCoeffs() ["
                    << checkEPotEqnStage << "] patch "
                    << mesh.boundary()[patchi].name()
                    << " face " << facei
                    << ", value=" << bCoeffs[facei]
                    << ", nEM=" << nEM
                    << ", time=" << runTime.timeName()
                    << abort(FatalError);
            }
        }
    }
}
