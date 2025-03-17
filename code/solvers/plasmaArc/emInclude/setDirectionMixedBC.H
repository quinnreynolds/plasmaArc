// Set mixed boundary conditions for A
// dA.n/dn = 0, A.t=0
{
    vectorField patchNormals;
    forAll(A.boundaryField(), patchI) {
        if (A.boundaryField()[patchI].type() == "directionMixed") {
            patchNormals = mesh.Sf().boundaryField()[patchI] / 
                mag(mesh.Sf().boundaryField()[patchI]);
            directionMixedFvPatchVectorField& dMPatch = 
                refCast<directionMixedFvPatchVectorField>
                (
                    A.boundaryFieldRef()[patchI]
                );
            forAll(dMPatch.valueFraction(), faceI) {
                dMPatch.valueFraction()[faceI] = 
                    I - symm(patchNormals[faceI]*patchNormals[faceI]);
            }
        }
    }
}
