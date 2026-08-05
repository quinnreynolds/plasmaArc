// TEMPORARY DIAGNOSTIC - not a permanent feature.
// Dumps forensic detail for cell "celli" (must be set before inclusion):
// its own volume/centre/field values, full 1-hop face-by-face detail
// (geometry, non-orthogonal correction, and now ek AND ePot on both
// sides), plus a condensed 2-hop neighbour summary. The 2-hop pass
// exists because fvm::laplacian's "corrected" non-orthogonal scheme
// applies its deferred correction using fvc::grad(ePot) interpolated to
// each face, and a Gauss-linear gradient at a cell depends on that
// cell's own neighbours - so the correction contributing to celli's
// source() can depend on data two hops away, which the original (1-hop
// only) dump never inspected. Does not abort - just prints; the caller
// decides when to stop.
{
    Pout<< "[dumpEPotCellForensics] cell=" << celli
        << " proc=" << Pstream::myProcNo()
        << " nEM=" << nEM
        << " time=" << runTime.timeName()
        << " stage=" << checkEPotEqnStage
        << nl
        << "  V=" << mesh.V()[celli]
        << " C=" << mesh.C()[celli]
        << " ek=" << ek[celli]
        << " ePot=" << ePot[celli]
        << " UcrossB=" << UcrossB[celli]
        << " dotTerm=" << dotTerm[celli]
        << endl;

    const label nInternalFaces = mesh.nInternalFaces();
    labelHashSet twoHopCells;

    const cell& c = mesh.cells()[celli];
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

            // Angle between the face normal and the owner-to-neighbour
            // vector: nonOrthDeltaCoeffs() ~ 1/cosAngle, so this blows up
            // as cosAngle -> 0 even when magSf/deltaCoeffs/volumes all
            // look fine.
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
                << " ePot_owner=" << ePot[own]
                << " ePot_neighbour=" << ePot[nei]
                << endl;

            const label other = (own == celli) ? nei : own;
            twoHopCells.insert(other);
        }
        else
        {
            const label patchi = mesh.boundaryMesh().whichPatch(facei);
            const label patchFacei = facei - mesh.boundary()[patchi].start();
            const polyPatch& pp = mesh.boundaryMesh()[patchi];

            Pout<< " patch=" << mesh.boundary()[patchi].name()
                << " patchType=" << pp.type()
                << " patchFacei=" << patchFacei
                << " ek_boundary=" << ek.boundaryField()[patchi][patchFacei]
                << " ePot_boundary=" << ePot.boundaryField()[patchi][patchFacei]
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
                const auto& procPatch = refCast<const processorPolyPatch>(pp);
                Pout<< " myProcNo=" << procPatch.myProcNo()
                    << " neighbProcNo=" << procPatch.neighbProcNo();
            }

            Pout<< endl;
        }
    }

    // Condensed 2-hop summary: for each direct (1-hop) neighbour, list
    // its own direct neighbours with just enough detail to spot an
    // extreme-but-technically-finite value that checkEPotInputsFinite.H's
    // pure isfinite() check wouldn't have flagged.
    labelHashSet printedTwoHop;
    forAllConstIter(labelHashSet, twoHopCells, hopIter)
    {
        const label oneHopCelli = hopIter.key();
        const cell& c2 = mesh.cells()[oneHopCelli];
        forAll(c2, fi2)
        {
            const label facei2 = c2[fi2];
            if (facei2 >= nInternalFaces)
            {
                continue;
            }

            const label own2 = mesh.faceOwner()[facei2];
            const label nei2 = mesh.faceNeighbour()[facei2];
            const label twoHopCelli = (own2 == oneHopCelli) ? nei2 : own2;

            if
            (
                twoHopCelli == celli
             || twoHopCells.found(twoHopCelli)
             || !printedTwoHop.insert(twoHopCelli)
            )
            {
                continue;
            }

            Pout<< "  2hop via=" << oneHopCelli
                << " cell=" << twoHopCelli
                << " V=" << mesh.V()[twoHopCelli]
                << " ek=" << ek[twoHopCelli]
                << " ePot=" << ePot[twoHopCelli]
                << endl;
        }
    }
}
