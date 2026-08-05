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
                // TEMPORARY DIAGNOSTIC - forensic dump of this cell and
                // its faces before aborting, since the four scalar/vector
                // inputs (ek, ePot, UcrossB, dotTerm) were already
                // confirmed finite here by checkEPotInputsFinite.H: if
                // the assembled source() is bad anyway, the problem is in
                // the discretisation geometry (volume, face areas,
                // owner/neighbour addressing, interpolation weights),
                // none of which are registered fields the earlier
                // (now-removed) blastAMR-side diagnostics could see.
                Pout<< "[checkEPotEqnFinite] FORENSIC DUMP cell=" << celli
                    << " proc=" << Pstream::myProcNo()
                    << " nEM=" << nEM
                    << " time=" << runTime.timeName()
                    << nl
                    << "  V=" << mesh.V()[celli]
                    << " C=" << mesh.C()[celli]
                    << " ek=" << ek[celli]
                    << " ePot=" << ePot[celli]
                    << " UcrossB=" << UcrossB[celli]
                    << " dotTerm=" << dotTerm[celli]
                    << endl;

                const cell& c = mesh.cells()[celli];
                const label nInternalFaces = mesh.nInternalFaces();
                forAll(c, fi)
                {
                    const label facei = c[fi];
                    const bool isInternal = facei < nInternalFaces;
                    const label own = mesh.faceOwner()[facei];

                    Pout<< "  face=" << facei
                        << " internal=" << isInternal
                        << " owner=" << own
                        << " magSf=" << mag(mesh.faceAreas()[facei])
                        << " Cf=" << mesh.faceCentres()[facei];

                    if (isInternal)
                    {
                        const label nei = mesh.faceNeighbour()[facei];

                        // Angle between the face normal and the owner-to-
                        // neighbour vector: nonOrthDeltaCoeffs() ~ 1/cosAngle,
                        // so this blows up as cosAngle -> 0 (face normal
                        // nearly perpendicular to the cell-to-cell line)
                        // even when magSf/deltaCoeffs/volumes all look fine.
                        const vector d(mesh.C()[nei] - mesh.C()[own]);
                        const vector n(mesh.faceAreas()[facei]/mag(mesh.faceAreas()[facei]));
                        const scalar cosAngle((d & n)/(mag(d) + SMALL));

                        Pout<< " neighbour=" << nei
                            << " deltaCoeffs=" << mesh.deltaCoeffs()[facei]
                            << " nonOrthDeltaCoeffs=" << mesh.nonOrthDeltaCoeffs()[facei]
                            << " nonOrthCorrVec=" << mesh.nonOrthCorrectionVectors()[facei]
                            << " cosAngle=" << cosAngle
                            << " weights=" << mesh.weights()[facei]
                            << " V_owner=" << mesh.V()[own]
                            << " V_neighbour=" << mesh.V()[nei]
                            << " ek_owner=" << ek[own]
                            << " ek_neighbour=" << ek[nei]
                            << endl;
                    }
                    else
                    {
                        const label patchi =
                            mesh.boundaryMesh().whichPatch(facei);
                        const label patchFacei =
                            facei - mesh.boundary()[patchi].start();
                        const polyPatch& pp = mesh.boundaryMesh()[patchi];

                        Pout<< " patch=" << mesh.boundary()[patchi].name()
                            << " patchType=" << pp.type()
                            << " patchFacei=" << patchFacei
                            << " ek_boundary="
                            << ek.boundaryField()[patchi][patchFacei]
                            // deltaCoeffs/weights/non-orthogonal correction
                            // are also defined on boundary (incl. processor)
                            // faces - fvm::laplacian uses them there too,
                            // and they were missing from the first version
                            // of this dump.
                            << " deltaCoeffs="
                            << mesh.deltaCoeffs().boundaryField()[patchi][patchFacei]
                            << " nonOrthDeltaCoeffs="
                            << mesh.nonOrthDeltaCoeffs().boundaryField()[patchi][patchFacei]
                            << " nonOrthCorrVec="
                            << mesh.nonOrthCorrectionVectors().boundaryField()[patchi][patchFacei]
                            << " weights="
                            << mesh.weights().boundaryField()[patchi][patchFacei];

                        if (isA<processorPolyPatch>(pp))
                        {
                            const auto& procPatch =
                                refCast<const processorPolyPatch>(pp);
                            Pout<< " myProcNo=" << procPatch.myProcNo()
                                << " neighbProcNo=" << procPatch.neighbProcNo();
                        }

                        Pout<< endl;
                    }
                }

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
