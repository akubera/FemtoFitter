///
/// \class Fitter1DGauss.cpp
///


#include "Fitter1DGauss.hpp"



double Fitter1DGauss::resid_calc_chi2_mrc(const FitResult &fr)
{
  if (mrc == nullptr) {
    std::cerr << "mrc is null\n";
    return NAN;
  }

  auto params = fr.as_params();
  return Fitter1D::resid_calc_mrc(params, *mrc, CalcChi2::resid_func);
}
