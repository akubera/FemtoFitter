
#ifdef __ROOTCLING__

#pragma link off all class;

#pragma link C++ class Data1D;
#pragma link C++ class Data3D;
#pragma link C++ class Value;

#pragma link C++ class Fitter3D<FitterGaussOSL>;
#pragma link C++ class FitterGaussOSL::FitParams;
#pragma link C++ class FitterGaussOSL::FitResult;

#pragma link C++ class Fitter3D<FitterGaussFull>;
#pragma link C++ class FitterGaussFull::FitParams;
#pragma link C++ class FitterGaussFull::FitResult;

#pragma link C++ class Fitter3D<FitterLevy>;
#pragma link C++ class FitterLevy::FitParams;
#pragma link C++ class FitterLevy::FitResult;

#pragma link C++ function apply_momentum_resolution_correction;

#pragma link C++ class FsiGamov;
#pragma link C++ class FsiGamov::Kcalc;
#pragma link C++ class FsiStatic;
#pragma link C++ class FsiStatic::Kcalc;
#pragma link C++ class FsiKFile;
#pragma link C++ class FsiKFile::KCalc;
#pragma link C++ class FsiCalculator;

#pragma link C++ class CoulombHist;

#endif
