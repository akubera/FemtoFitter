
#ifdef __ROOTCLING__

#pragma link off all class;

#pragma link C++ class Data1D;
#pragma link C++ class Data3D;
#pragma link C++ class Value;

#pragma link C++ class Fitter1D<Fitter1DGauss>;
#pragma link C++ class Fitter1DGauss::FitParams;
#pragma link C++ class Fitter1DGauss::FitResult;

#pragma link C++ class Fitter1D<Fitter1DGaussLin>;
#pragma link C++ class Fitter1DGaussLin::FitParams;
#pragma link C++ class Fitter1DGaussLin::FitResult;

#pragma link C++ class Fitter1D<Fitter1DGaussPolyBg>;
#pragma link C++ class Fitter1DGaussPolyBg::FitParams;
#pragma link C++ class Fitter1DGaussPolyBg::FitResult;

#pragma link C++ class Fitter1D<Fitter1DLevy>;
#pragma link C++ class Fitter1DLevy::FitParams;
#pragma link C++ class Fitter1DLevy::FitResult;

#pragma link C++ class Fitter1D<Fitter1DLevyPolyBg>;
#pragma link C++ class Fitter1DLevyPolyBg::FitParams;
#pragma link C++ class Fitter1DLevyPolyBg::FitResult;

#pragma link C++ class Fitter3D<Fitter3DGaussLcms>;
#pragma link C++ class Fitter3DGaussLcms::FitParams;
#pragma link C++ class Fitter3DGaussLcms::FitResult;

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

#pragma link C++ class Mrc1DRatio;
#pragma link C++ class Mrc1DRatioMixed;
#pragma link C++ class Mrc1DMatrix;
#pragma link C++ class Mrc1DMatrixJesse;
#pragma link C++ class MrcTransform1D;

#pragma link C++ class Mrc3DRatio;
#pragma link C++ class Mrc3DRatioMixed;
#pragma link C++ class MrcHypercube3D;

#pragma link C++ class ParamHints;

#pragma link C++ class CoulombHist;

#endif
