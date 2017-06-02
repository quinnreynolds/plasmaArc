//internal field update
forAll(QRad(),cellI)
{
    QRad.ref()[cellI] = QRadLookup(T.ref()[cellI]);
    ek.ref()[cellI] = ekLookup(T.ref()[cellI]);
    QRadSp.ref()[cellI] = QRadLookup.rateOfChange(T.ref()[cellI]);
    ekSp.ref()[cellI] = ekLookup.rateOfChange(T.ref()[cellI]);
}

ek.correctBoundaryConditions();

//boundary field update
forAll(QRad.boundaryField(),patchI) {
    forAll(QRad.boundaryFieldRef()[patchI],cellI)
    {
        QRad.boundaryFieldRef()[patchI][cellI] =
            QRadLookup(T.boundaryFieldRef()[patchI][cellI]);
        QRadSp.boundaryFieldRef()[patchI][cellI] =
            QRadLookup.rateOfChange(T.boundaryFieldRef()[patchI][cellI]);
        ekSp.boundaryFieldRef()[patchI][cellI] =
            ekLookup.rateOfChange(T.boundaryFieldRef()[patchI][cellI]);
    }
}
