//TODO mixing rules for multicomponent plasma mixtures and switch
// calculation accordingly (or separate solver...)

//internal field update
forAll(QRad(),cellI)
{
    QRad.ref()[cellI] = QRadLookup.interpolate(T.ref()[cellI]);
    ek.ref()[cellI] = ekLookup.interpolate(T.ref()[cellI]);
    QRadSp.ref()[cellI] = QRadSpLookup.interpolate(T.ref()[cellI]);
    ekSp.ref()[cellI] = ekSpLookup.interpolate(T.ref()[cellI]);
}

ek.correctBoundaryConditions();

//boundary field update
forAll(QRad.boundaryField(),patchI) {
    forAll(QRad.boundaryFieldRef()[patchI],cellI)
    {
        QRad.boundaryFieldRef()[patchI][cellI] =
            QRadLookup.interpolate(T.boundaryFieldRef()[patchI][cellI]);
        QRadSp.boundaryFieldRef()[patchI][cellI] =
            QRadSpLookup.interpolate(T.boundaryFieldRef()[patchI][cellI]);
        ekSp.boundaryFieldRef()[patchI][cellI] =
            ekSpLookup.interpolate(T.boundaryFieldRef()[patchI][cellI]);
    }
}
