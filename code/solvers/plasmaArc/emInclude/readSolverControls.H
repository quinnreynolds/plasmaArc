// EM solver loop controls
const dictionary& ePotDict
    = mesh.solutionDict().subDict("solvers").subDict("ePot");

const dictionary& ADict
    = mesh.solutionDict().subDict("solvers").subDict("A");

// TODO: This isn't particularly elegant
scalar toleranceEPot = 10. * ePotDict.lookupOrDefault<scalar>("tolerance",1e-6);
scalar toleranceA = 10. * ADict.lookupOrDefault<scalar>("tolerance",1e-6);

const dictionary& emDict = mesh.solutionDict().subDict("EM");
const int nMaxEMLoops = emDict.lookupOrDefault<int>("nMaxEMLoops", 5);
const int emSolverFreq = emDict.lookupOrDefault<int>("emSolverFreq", 1);
// TODO: Add documentation explaining use of wedgeGeometryYN
const Switch wedgeGeometryYN = emDict.lookupOrDefault<Switch>("wedgeGeometryYN", false);
