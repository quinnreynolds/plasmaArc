//internal field update
forAll(QRad(),cellI)
{
    QRad.ref()[cellI] = scalarLookup(T.ref()[cellI], QRadStartT, QRadDeltaT, QRadData);
    ek.ref()[cellI] = scalarLookup(T.ref()[cellI], ekStartT, ekDeltaT, ekData);
    QRadSp.ref()[cellI] = scalarLookup(T.ref()[cellI], QRadStartT, QRadDeltaT, QRadSpData);
    ekSp.ref()[cellI] = scalarLookup(T.ref()[cellI], ekStartT, ekDeltaT, ekSpData);
}

ek.correctBoundaryConditions();

//boundary field update
forAll(QRad.boundaryField(),patchI) {
    forAll(QRad.boundaryFieldRef()[patchI],cellI)
    {
        QRad.boundaryFieldRef()[patchI][cellI] =
            scalarLookup(T.boundaryFieldRef()[patchI][cellI], QRadStartT, QRadDeltaT, QRadData);
        QRadSp.boundaryFieldRef()[patchI][cellI] =
            scalarLookup(T.boundaryFieldRef()[patchI][cellI], QRadStartT, QRadDeltaT, QRadSpData);
        ekSp.boundaryFieldRef()[patchI][cellI] =
            scalarLookup(T.boundaryFieldRef()[patchI][cellI], ekStartT, ekDeltaT, ekSpData);
    }
}
