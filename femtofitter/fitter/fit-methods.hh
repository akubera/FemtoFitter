///
/// \file femtofitter/fitter/fit-methods.hh
///


#pragma once

#ifndef FITMETHODS_HH
#define FITMETHODS_HH

#define DECLARE_FIT_METHODS(Fitter) \
  FitResult fit_chi2() { return Fitter::fit_chi2(); } \
  FitResult fit_chi2_mrc() { return Fitter::fit_chi2_mrc(); } \
  FitResult fit_pml() { return Fitter::fit_pml(); } \
  FitResult fit_pml_mrc() { return Fitter::fit_pml_mrc(); }

  // FitResult fit_pml_mrc_quick() { return Fitter::fit_pml_mrc_quick(); }


#define DECLARE_RESID_METHODS(Fitter) \
  double residual_chi2(const FitResult &r) const \
    { return residual_chi2(r.as_params()); } \
  double residual_chi2(const FitParams &p) const \
    { return Fitter::resid_calc(p, CalcChi2::resid_func); } \
  double residual_pml(const FitParams &p) const \
    { return Fitter::resid_calc(p, CalcLoglike::resid_func); } \
  double residual_pml(const FitResult &p) const \
    { return residual_pml(p.as_params()); } \
  double residual_chi2_mrc(const FitResult &r) const \
    { return residual_chi2_mrc(r.as_params()); } \
  double residual_chi2_mrc(const FitParams &p) const \
    { return Fitter::resid_calc_mrc(p, *mrc, CalcChi2::resid_func); } \
  double residual_pml_mrc(const FitParams &p) const \
    { return Fitter::resid_calc_mrc(p, *mrc, CalcLoglike::resid_func); }


#define DECLARE_FILL_METHODS(HistType) \
  void fill(HistType &h, const FitResult &r, UInt_t npoints=1) const \
    { fill(h, r.as_params(), npoints); } \
  void fill(HistType &h, const FitParams &p, UInt_t npoints=1) const \
    { p.fill(h, *fsi, npoints); } \
  void fill_smeared_fit(HistType &h, const FitResult &r) \
    { fill_smeared_fit(h, r.as_params()); } \
  void fill_smeared_fit(HistType &h, const FitParams &p) \
    { Super::fill_smeared_fit(h, p); }



#endif
