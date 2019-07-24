
#ifdef __ROOTCLING__

#pragma link off all class;

#pragma link C++ class Data1D;
#pragma link C++ class Data3D;
#pragma link C++ class Value;

#pragma link C++ class Fitter1D<FitterGauss1D>;
#pragma link C++ class FitterGauss1D::FitParams;
#pragma link C++ class FitterGauss1D::FitResult;

#pragma link C++ class Fitter1D<Fitter1DGaussLin>;
#pragma link C++ class Fitter1DGaussLin::FitParams;
#pragma link C++ class Fitter1DGaussLin::FitResult;

#pragma link C++ class Fitter1D<Fitter1DGaussPolybg>;
#pragma link C++ class Fitter1DGaussPolybg::FitParams;
#pragma link C++ class Fitter1DGaussPolybg::FitResult;

#pragma link C++ class Fitter1D<FitterLevy1D>;
#pragma link C++ class FitterLevy1D::FitParams;
#pragma link C++ class FitterLevy1D::FitResult;

#pragma link C++ class Fitter1D<Fitter1DLevyPolyBg>;
#pragma link C++ class Fitter1DLevyPolyBg::FitParams;
#pragma link C++ class Fitter1DLevyPolyBg::FitResult;

#pragma link C++ class Fitter3D<FitterGaussOSL>;
#pragma link C++ class FitterGaussOSL::FitParams;
#pragma link C++ class FitterGaussOSL::FitResult;

#pragma link C++ class Fitter3D<FitterGaussFull>;
#pragma link C++ class FitterGaussFull::FitParams;
#pragma link C++ class FitterGaussFull::FitResult;

#pragma link C++ class Fitter3D<FitterLevy3D>;
#pragma link C++ class FitterLevy3D::FitParams;
#pragma link C++ class FitterLevy3D::FitResult;

#pragma link C++ class Fitter3D<FitterGauss3DLcmsOS>;
#pragma link C++ class FitterGauss3DLcmsOS::FitParams;
#pragma link C++ class FitterGauss3DLcmsOS::FitResult;

#pragma link C++ class Fitter3D<FitterGauss3DLcmsOL>;
#pragma link C++ class FitterGauss3DLcmsOL::FitParams;
#pragma link C++ class FitterGauss3DLcmsOL::FitResult;

#pragma link C++ class Fitter3D<FitterLevyFull>;
#pragma link C++ class FitterLevyFull::FitParams;
#pragma link C++ class FitterLevyFull::FitResult;

#pragma link C++ function apply_momentum_resolution_correction;

#pragma link C++ class FsiGamov;
#pragma link C++ class FsiGamov::Kcalc;
#pragma link C++ class FsiStatic;
#pragma link C++ class FsiStatic::Kcalc;
#pragma link C++ class FsiKFile;
#pragma link C++ class FsiKFile::KCalc;
#pragma link C++ class FsiCalculator;

#pragma link C++ class MrcRatio1D;
#pragma link C++ class MrcMatrix1D;
#pragma link C++ class MrcMatrix1DJesse;
#pragma link C++ class MrcTransform1D;

#pragma link C++ class MrcRatio3D;
#pragma link C++ class MrcHypercube3D;

#pragma link C++ class ParamHints;

#pragma link C++ class CoulombHist;

#endif
