
scalar scalarLookup(scalar Temperature, scalar myStartT, scalar myDeltaT, List<scalar>& myData)
{
    label iT = label((Temperature - myStartT) / myDeltaT);

    // Clamp bounds
    if (iT > (myData.size()-2))
    {
        return myData[myData.size()-1];
    }
    if (iT < 0)
    {
        return myData[0];
    }

    scalar deltaFrac = (Temperature - myStartT) / myDeltaT - scalar(iT);

    return myData[iT] + deltaFrac * (myData[iT+1] - myData[iT]);
}
