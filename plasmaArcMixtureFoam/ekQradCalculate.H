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
            * QRadLookupList[speciesI](T.ref()[cellI]);
        ekTemp += Y[speciesI].ref()[cellI]
            * ekLookupList[speciesI](T.ref()[cellI]);
        QRadSpTemp += Y[speciesI].ref()[cellI]
            * QRadLookupList[speciesI].rateOfChange(T.ref()[cellI]);
        ekSpTemp += Y[speciesI].ref()[cellI]
            * ekLookupList[speciesI].rateOfChange(T.ref()[cellI]);
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
                * QRadLookupList[speciesI](T.boundaryFieldRef()[patchI][cellI]);
            QRadSpTemp += Y[speciesI].boundaryFieldRef()[patchI][cellI]
                * QRadLookupList[speciesI].rateOfChange(T.boundaryFieldRef()[patchI][cellI]);
            ekSpTemp += Y[speciesI].boundaryFieldRef()[patchI][cellI]
                * ekLookupList[speciesI].rateOfChange(T.boundaryFieldRef()[patchI][cellI]);
        }
        QRad.boundaryFieldRef()[patchI][cellI] = QRadTemp;
        QRadSp.boundaryFieldRef()[patchI][cellI] = QRadSpTemp;
        ekSp.boundaryFieldRef()[patchI][cellI] = ekSpTemp;
    }
}
