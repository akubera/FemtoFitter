///
/// \file fitter/Fitter1DGaussLin.hpp
///

#pragma once

#ifndef FITTER_FITTER1DGAUSSLIN_HPP
#define FITTER_FITTER1DGAUSSLIN_HPP

#include "Fitter1D.hpp"

#include <TF1.h>
#include <TFitResult.h>


/// \class Fitter1DGaussLin
/// \brief Gaussian 1D fit with linear background
///
struct Fitter1DGaussLin : public Fitter1D<Fitter1DGaussLin> {

  using Super = Fitter1D<Fitter1DGaussLin>;

  struct FitParams;

  static std::string GetName()
    { return "LinGauss1D"; }

  static unsigned char CountParams()
    { return 4; }

  static double
  gauss(double qinv, double RinvSq, double lam, double slope, double K=1.0, double norm=1.0)
    {
      const double
        E = qinv * qinv * RinvSq,
        gauss = 1.0 + std::exp(-E/HBAR_C_SQ),
        result = (1.0 - lam) + lam * K * gauss;

      return (norm + qinv * slope) * result;
    }

  enum {
    DATA_PARAM_IDX = 0,

    LAM_PARAM_IDX = 1,
    R_PARAM_IDX = 2,
    NORM_PARAM_IDX = 3,
    SLOPE_PARAM_IDX = 4,
  };

  /// \class FitResult
  /// \brief Result of the fit
  ///
  struct FitResult : FitResult1D<FitResult, Fitter1DGaussLin> {

    Value norm,
          lam,
          radius,
          slope;

    FitResult()
      : norm(1.0, 0.0)
      , lam(1.0, 0.0)
      , radius(1.0, 0.0)
      , slope(0.0, 0.0)
      { }

    FitResult(const FitResult &orig) = default;

    FitResult(TMinuit &minuit)
      : norm(minuit, NORM_PARAM_IDX)
      , lam(minuit, LAM_PARAM_IDX)
      , radius(minuit, R_PARAM_IDX)
      , slope(minuit, SLOPE_PARAM_IDX)
      { }

    virtual ~FitResult()
      { }

    std::map<std::string, double>
    as_map() const
      {
        #define OUT(__name) {#__name, __name.value}, { # __name "_err", __name.error}

        return {
          OUT(radius),
          OUT(lam),
          OUT(norm),
          OUT(slope),
        };

        #undef OUT
      }

    std::string
    print() const
      {
        std::string result;
        result += Form("Radius: %f ± %f\n", radius.value, radius.error);
        result += Form("Lambda: %f ± %f\n", lam.value, lam.error);
        result += Form("Norm: %f ± %f\n", norm.value, norm.error);
        result += Form("Slope: %f ± %f\n", slope.value, slope.error);
        return result;
      }

    std::string
    __repr__() const
      {
        return Form("<Fitter1DGaussLin::FitResult radius=%f lambda=%f norm=%f slope=%g>",
                    radius.value, lam.value, norm.value, slope.value);
      }

    PyObject*
    as_dict() const
      {
        auto *dict = PyDict_New();
        PyDict_SetItemString(dict, "radius", PyFloat_FromDouble(radius.value));
        PyDict_SetItemString(dict, "radius_err", PyFloat_FromDouble(radius.error));
        PyDict_SetItemString(dict, "lam", PyFloat_FromDouble(lam.value));
        PyDict_SetItemString(dict, "lam_err", PyFloat_FromDouble(lam.error));
        PyDict_SetItemString(dict, "norm", PyFloat_FromDouble(norm.value));
        PyDict_SetItemString(dict, "norm_err", PyFloat_FromDouble(norm.error));
        PyDict_SetItemString(dict, "slope", PyFloat_FromDouble(slope.value));
        PyDict_SetItemString(dict, "slope_err", PyFloat_FromDouble(slope.error));
        return dict;
      }

    double evaluate(const double q, const double K) const
      { return Fitter1DGaussLin::gauss(q, radius * radius, lam, slope, K, norm); }

    FitParams as_params() const;

    void FillMinuit(TMinuit &minuit) const override
      {
        int errflag = 0;
        // minuit.mnparm(NORM_PARAM_IDX, "Norm", norm.value, norm.error * norm.error * 4, 0.0, 0.0, errflag);
        minuit.mnparm(NORM_PARAM_IDX, "Norm", norm.value, 0.005, 0.0, 0.0, errflag);
        minuit.mnparm(LAM_PARAM_IDX, "Lam", lam.value, .01, 0.0, 0.0, errflag);
        minuit.mnparm(R_PARAM_IDX, "Radius", radius.value, 0.2, 0.0, 0.0, errflag);
        minuit.mnparm(SLOPE_PARAM_IDX, "Slope", slope.value, 0.2, 0.0, 0.0, errflag);
      }

    void Normalize(TH1 &h) const override
      {
        h.Scale(1.0 / norm.value);
      }
  };

  struct FitParams : FitParam1D<FitParams> {
    double norm,
           lam,
           radius,
           slope;

    FitParams(const double *par)
      : norm(par[NORM_PARAM_IDX])
      , lam(par[LAM_PARAM_IDX])
      , radius(par[R_PARAM_IDX])
      , slope(par[SLOPE_PARAM_IDX])
      { }

    FitParams(const FitParams &) = default;

    FitParams(const FitResult &res)
      : norm(res.norm)
      , lam(res.lam)
      , radius(res.radius)
      , slope(res.slope)
      { }

    double Rinv() const
      {
        return radius;
      }

    bool is_invalid() const
      {
        #define INVALID(__X) (__X <= 0) || std::isnan(__X)
        return INVALID(norm)
            || INVALID(lam)
            || INVALID(radius)
            || std::isnan(slope);
        #undef INVALID
      }

    std::string
    print() const
      {
        std::string result;
        result += Form("Radius: %f \n", radius);
        result += Form("Lambda: %f\n", lam);
        result += Form("Norm: %f\n", norm);
        result += Form("Slope: %f\n", slope);
        return result;
      }

    std::string
    __repr__() const
      {
        return Form("<Fitter1DGaussLin::FitParam radius=%g lambda=%g norm=%g slope=%g>",
                    radius, lam, norm, slope);
      }

    PyObject*
    as_dict() const
      {
        auto *dict = PyDict_New();
        PyDict_SetItemString(dict, "radius", PyFloat_FromDouble(radius));
        PyDict_SetItemString(dict, "lam", PyFloat_FromDouble(lam));
        PyDict_SetItemString(dict, "norm", PyFloat_FromDouble(norm));
        PyDict_SetItemString(dict, "slope", PyFloat_FromDouble(slope));
        return dict;
      }

    double evaluate(const double q, const double K) const
      {
        return Fitter1DGaussLin::gauss(q, radius * radius, lam, slope, K, norm);
      }

    void Normalize(TH1 &h) const
      {
        const double x = h.GetXaxis()->GetBinCenter(h.GetXaxis()->GetLast());
        double background = norm + x * slope;

        h.Scale(1.0 / background);
      }
  };

  int
  setup_minuit(TMinuit &minuit) const override
    {
      return setup_minuit(minuit, 0.18, 1.0);
    }

  TFitResultPtr fit_background(double bglo, double bghi, TH1 *cf_input=nullptr) const
    {
      auto cf = cf_input ? nullptr : data.cf_src();
      TF1 fg("linbg", "[0] + x * [1]");
      auto bg = (cf_input ? cf_input : cf.get())->Fit(&fg, "SQ0", "", bglo, bghi);
      return bg;
    }

  int
  setup_minuit(TMinuit &minuit, double bglo, double bghi) const
    {
      auto bg = fit_background(bglo, bghi);
      const double *bg_params = bg->GetParams();

      int errflag = 0;
      minuit.mnparm(NORM_PARAM_IDX, "Norm", 0.2, 0.02, 0.0, 0.0, errflag);
      minuit.mnparm(LAM_PARAM_IDX, "Lam", 0.5, 0.01, 0.0, 0.0, errflag);
      minuit.mnparm(R_PARAM_IDX, "Radius", 1.0, 0.2, 0.0, 0.0, errflag);
      minuit.mnparm(SLOPE_PARAM_IDX, "Slope", bg_params[1], 0.0, 0.0, 0.0, errflag);
      minuit.FixParameter(SLOPE_PARAM_IDX);

      const double this_dbl = static_cast<double>((intptr_t)this);
      minuit.mnparm(DATA_PARAM_IDX, "DATA_PTR", this_dbl, 0, 0, INTPTR_MAX, errflag);
      minuit.FixParameter(DATA_PARAM_IDX);

      if (errflag != 0) {
        std::cerr << "Error setting paramters: " << errflag << "\n";
        throw std::runtime_error("Could not set Minuit parameters.");
      }

      return errflag;
    }

  void setup_pml_fitter(TMinuit &minuit, double bglo, double bghi)
    {
      setup_minuit(minuit, bglo, bghi);
      set_pml_func(minuit);
    }

  void setup_pml_mrc_fitter(TMinuit &minuit, double bglo, double bghi)
    {
      setup_minuit(minuit, bglo, bghi);
      set_pml_mrc_func(minuit);
    }

  void
  set_and_fix_variable(TMinuit &minuit, std::string name, double val)
    {
      int idx = (name == "R") || (name == "radius") ? R_PARAM_IDX
              : (name == "lam" || name == "Lam") ? LAM_PARAM_IDX
              : -1;

      if (idx < 0) {
        std::cerr << "Unknown parameter '" << name << "'\n";
        return;
      }

      double stepsize = (idx == LAM_PARAM_IDX) ? 0.1 : 1.0;

      int errflag = 0;
      minuit.mnparm(idx, name, val, stepsize, 0.0, 0.0, errflag);
      minuit.FixParameter(idx);
    }

  Fitter1DGaussLin(const TH1 &num, const TH1 &den, double limit)
    : Fitter1D(num, den, limit)
    { }

  Fitter1DGaussLin(TDirectory &tdir, double limit)
    : Fitter1D(tdir, limit)
    { }

  Fitter1DGaussLin(const Data1D &dat)
    : Fitter1D(dat)
    { }

  virtual ~Fitter1DGaussLin() = default;

  FitResult fit_chi2()
    { return Fitter1D::fit_chi2(); }

  FitResult fit_chi2_mrc()
    { return Fitter1D::fit_chi2_mrc(); }

  FitResult fit_pml()
    { return Fitter1D::fit_pml(); }

  /// Fit with log-likelihood method and momentum-correction smearing
  FitResult fit_pml_mrc()
    { return Fitter1D::fit_pml_mrc(); }

  FitResult fit_pml_mrc_quick()
    { return Fitter1D::fit_pml_mrc_quick(); }

  FitResult fit_pml(double bglo, double bghi)
    {
      if (fsi == nullptr) {
        throw std::runtime_error("Fitter missing Fsi object");
      }

      TMinuit minuit;
      minuit.SetPrintLevel(-1);

      setup_pml_fitter(minuit, bglo, bghi);
      return do_fit_minuit(minuit);
    }

  FitResult fit_pml_mrc(double bglo, double bghi)
    {
      if (mrc == nullptr) {
        throw std::runtime_error("Fitter missing Mrc1D object");
      }

      if (fsi == nullptr) {
        throw std::runtime_error("Fitter missing Fsi object");
      }

      TMinuit minuit;
      minuit.SetPrintLevel(-1);

      setup_pml_mrc_fitter(minuit, bglo, bghi);
      return do_fit_minuit(minuit);
    }

  // void fit_with_random_inits(TMinuit &minuit, FitResult &res, int);

  void fill(TH1 &h, const FitParams &p, UInt_t npoints=1) const
    {
      p.fill(h, *fsi, npoints);
    }

  void fill(TH1 &h, const FitResult &p, UInt_t npoints=1) const
    {
      fill(h, p.as_params(), npoints);
    }

  void fill_smeared_fit(TH1 &h, const FitResult &fr)
    {
      fill_smeared_fit(h, fr.as_params());
    }

  void fill_smeared_fit(TH1 &h, const FitParams &p)
    {
      Fitter1D::fill_smeared_fit(h, p);
    }

  double resid_calc_chi2_mrc(const FitParams &params) const
    {
      return Fitter1D::resid_calc_chi2_mrc(params);
    }

  double resid_calc_chi2_mrc(const FitResult &fr) const
    {
      auto params = fr.as_params();
      return Fitter1D::resid_calc_chi2_mrc(params);
    }

  std::unique_ptr<TH1> get_cf(const FitParams &p) const
    {
      std::unique_ptr<TH1> cf(static_cast<TH1*>(data.src->num->Clone()));
      cf->Reset();
      fill(*cf, p);
      return cf;
    }
};

inline auto
Fitter1DGaussLin::FitResult::as_params() const -> FitParams
{
  return FitParams(*this);
}

#endif
