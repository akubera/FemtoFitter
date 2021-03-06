///
/// \file femtofitter/fitter/fit-methods.hh
///


#pragma once

#ifndef FITMETHODS_HH
#define FITMETHODS_HH

#define IMPL_FIT_METHOD(__setup_name) \
  if (fsi == nullptr) { throw std::runtime_error("Fitter missing FsiObject"); }\
  TMinuit minuit; minuit.SetPrintLevel(-1); \
  __setup_name(minuit);\
  return do_fit_minuit(minuit);

#define IMPL_BG_FIT_METHOD(__setup_name) \
  if (fsi == nullptr) { throw std::runtime_error("Fitter missing FsiObject"); }\
  TMinuit minuit; minuit.SetPrintLevel(-1); \
  __setup_name(minuit, bglo, bghi);\
  return do_fit_minuit(minuit);


#define DECLARE_FIT_METHODS(Fitter) \
  FitResult fit_chi2() { return Fitter::fit_chi2(); } \
  FitResult fit_chi2_mrc() { return Fitter::fit_chi2_mrc(); } \
  FitResult fit_pml() { return Fitter::fit_pml(); } \
  FitResult fit_pml_mrc() { return Fitter::fit_pml_mrc(); }

  // FitResult fit_pml_mrc_quick() { return Fitter::fit_pml_mrc_quick(); }



#define DECLARE_BG_MINUIT_SETUP_METHODS() \
  void setup_chi2_fitter(TMinuit &minuit, double bglo, double bghi)     \
    { setup_minuit(minuit, bglo, bghi); set_chi2_func(minuit); }        \
  void setup_chi2_mrc_fitter(TMinuit &minuit, double bglo, double bghi) \
    { setup_minuit(minuit, bglo, bghi); set_chi2_mrc_func(minuit); }    \
  void setup_pml_fitter(TMinuit &minuit, double bglo, double bghi)      \
    { setup_minuit(minuit, bglo, bghi); set_pml_func(minuit); }         \
  void setup_pml_mrc_fitter(TMinuit &minuit, double bglo, double bghi)  \
    { setup_minuit(minuit, bglo, bghi); set_pml_mrc_func(minuit); }


#define DECLARE_BG_FIT_METHODS() \
  FitResult fit_chi2(double bglo, double bghi)     \
    { IMPL_BG_FIT_METHOD(setup_chi2_fitter) }      \
  FitResult fit_chi2_mrc(double bglo, double bghi) \
    { IMPL_BG_FIT_METHOD(setup_chi2_mrc_fitter) }  \
  FitResult fit_pml(double bglo, double bghi)      \
    { IMPL_BG_FIT_METHOD(setup_pml_fitter) }       \
  FitResult fit_pml_mrc(double bglo, double bghi)  \
    { IMPL_BG_FIT_METHOD(setup_pml_mrc_fitter) }


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


#define DECLARE_METHODS_GET_CF(HistType) \
  std::unique_ptr<HistType> get_cf(const FitResult &r) const \
    { return get_cf(r.as_params()); } \
  std::unique_ptr<HistType> get_cf(const FitParams &p) const \
    { return mrc ? get_smeared_cf(p) : get_unsmeared_cf(p); } \
  \
  std::unique_ptr<HistType> get_smeared_cf(const FitResult &r) const \
    { return get_smeared_cf(r.as_params()); } \
  std::unique_ptr<HistType> get_smeared_cf(const FitParams &p) const \
    { std::unique_ptr<HistType> cf(static_cast<HistType*>(data.src->num->Clone())); \
      cf->Reset(); cf->SetDirectory(nullptr); cf->SetStats(0); \
      mrc->FillSmearedFit(*cf, p, *fsi, 1); \
      return cf; } \
  \
  std::unique_ptr<HistType> get_unsmeared_cf(const FitResult &r) const \
    { return get_unsmeared_cf(r.as_params()); } \
  std::unique_ptr<HistType> get_unsmeared_cf(const FitParams &p) const \
    { std::unique_ptr<HistType> cf(static_cast<HistType*>(data.src->num->Clone())); \
      cf->Reset(); cf->SetDirectory(nullptr); cf->SetStats(0); \
      fill(*cf, p); \
      return cf; } \


#define DECLARE_FILL_METHODS(HistType) \
  void fill(HistType &h, const FitResult &r, UInt_t npoints=1) const \
    { fill(h, r.as_params(), npoints); } \
  void fill(HistType &h, const FitParams &p, UInt_t npoints=1) const \
    { p.fill(h, *fsi, npoints); } \
  void fill_smeared_fit(HistType &h, const FitResult &r) \
    { fill_smeared_fit(h, r.as_params()); } \
  void fill_smeared_fit(HistType &h, const FitParams &p) \
    { Super::fill_smeared_fit(h, p); }


#define DECLARE_FILL_METHODS_3D() \
  void fill(TH3 &h, const FitParams &p) const \
    { _fill(h, p); } \
  void fill(TH3 &h, const FitResult &r) const \
    { _fill(h, r.as_params()); } \
  void fill_smeared_fit(TH3 &h, const FitParams &p) const \
    { _fill_smeared_fit(h, p); } \
  void fill_smeared_fit(TH3 &h, const FitResult &r) \
    { fill_smeared_fit(h, r.as_params()); }

#endif
