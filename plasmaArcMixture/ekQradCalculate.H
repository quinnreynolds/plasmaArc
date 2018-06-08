//internal field update

forAll(QRad(),cellI)
{
    scalar QRadTemp = 0.;
    scalar ekTemp = 0.;
    scalar QRadSpTemp = 0.;
    scalar ekSpTemp = 0.;
    forAll(plasmaSpeciesTable, speciesI)
    {
        QRadTemp += Y[speciesI].ref()[cellI]
            * scalarLookup(T.ref()[cellI], QRadStartT, QRadDeltaT, QRadDataList[speciesI]);
        ekTemp += Y[speciesI].ref()[cellI]
            * scalarLookup(T.ref()[cellI], ekStartT, ekDeltaT, ekDataList[speciesI]);
        QRadSpTemp += Y[speciesI].ref()[cellI]
            * scalarLookup(T.ref()[cellI], QRadStartT, QRadDeltaT, QRadSpDataList[speciesI]);
        ekSpTemp += Y[speciesI].ref()[cellI]
            * scalarLookup(T.ref()[cellI], ekStartT, ekDeltaT, ekSpDataList[speciesI]);
    }
    QRad.ref()[cellI] = QRadTemp;
    ek.ref()[cellI] = ekTemp;
    QRadSp.ref()[cellI] = QRadSpTemp;
    ekSp.ref()[cellI] = ekSpTemp;
}

ek.correctBoundaryConditions();

//boundary field update
forAll(QRad.boundaryField(),patchI) {
    forAll(QRad.boundaryFieldRef()[patchI],cellI)
    {
        scalar QRadTemp = 0.;
        scalar QRadSpTemp = 0.;
        scalar ekSpTemp = 0.;
        forAll(plasmaSpeciesTable, speciesI)
        {
            QRadTemp += Y[speciesI].boundaryFieldRef()[patchI][cellI]
                * scalarLookup(T.boundaryFieldRef()[patchI][cellI], QRadStartT, QRadDeltaT, QRadDataList[speciesI]);
            QRadSpTemp += Y[speciesI].boundaryFieldRef()[patchI][cellI]
                * scalarLookup(T.boundaryFieldRef()[patchI][cellI], QRadStartT, QRadDeltaT, QRadSpDataList[speciesI]);
            ekSpTemp += Y[speciesI].boundaryFieldRef()[patchI][cellI]
                * scalarLookup(T.boundaryFieldRef()[patchI][cellI], ekStartT, ekDeltaT, ekSpDataList[speciesI]);
        }
        QRad.boundaryFieldRef()[patchI][cellI] = QRadTemp;
        QRadSp.boundaryFieldRef()[patchI][cellI] = QRadSpTemp;
        ekSp.boundaryFieldRef()[patchI][cellI] = ekSpTemp;
    }
}
