///
/// \file fitter/Fitter1DGauss.hpp
///

#pragma once

#ifndef FITTER_FITTER1DGAUSS_HPP
#define FITTER_FITTER1DGAUSS_HPP

#include "Fitter1D.hpp"

/// \class Fitter1DGauss
/// \brief Gaussian 1D fit
///
struct Fitter1DGauss : public Fitter1D<Fitter1DGauss> {

  using Super = Fitter1D<Fitter1DGauss>;
  struct FitParams;
  struct FitResult;

  static std::string GetName()
    { return "Gauss1D"; }

  static unsigned char CountParams()
    { return 3; }

  static double
  gauss(double qinv, double RinvSq, double lam, double K=1.0, double norm=1.0)
    {
      const double
        E = qinv * qinv * RinvSq,
        gauss = 1.0 + std::exp(-E/HBAR_C_SQ),
        result = (1.0 - lam) + lam * K * gauss;

      return norm * result;
    }

  enum {
    DATA_PARAM_IDX = 0,

    LAM_PARAM_IDX = 1,
    R_PARAM_IDX = 2,
    NORM_PARAM_IDX = 3,
  };

  /// \class FitResult
  /// \brief Result of the fit
  ///
  struct FitResult : FitResult1D<FitResult, Fitter1DGauss> {

    Value norm,
          lam,
          radius;

    FitResult()
      : norm(1.0, 0.0)
      , lam(1.0, 0.0)
      , radius(1.0, 0.0)
      { }

    FitResult(const FitResult &orig) = default;

    FitResult(TMinuit &minuit)
      : norm(minuit, NORM_PARAM_IDX)
      , lam(minuit, LAM_PARAM_IDX)
      , radius(minuit, R_PARAM_IDX)
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
          OUT(norm)
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
        return result;
      }

    std::string
    __repr__() const
      {
        return Form("<Fitter1DGauss::FitResult radius=%f lambda=%f norm=%f>",
                    radius.value, lam.value, norm.value);
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
        return dict;
      }

    double evaluate(const double q, const double K) const
      { return Fitter1DGauss::gauss(q, radius * radius, lam, K, norm); }

    FitParams as_params() const;

    void FillMinuit(TMinuit &minuit) const override
      {
        int errflag = 0;
        minuit.mnparm(LAM_PARAM_IDX, "Lam", lam.value, .01, 0.0, 0.0, errflag);
        minuit.mnparm(R_PARAM_IDX, "Radius", radius.value, 0.2, 0.0, 0.0, errflag);
        minuit.mnparm(NORM_PARAM_IDX, "Norm", norm.value, 0.005, 0.0, 0.0, errflag);
      }

    void Normalize(TH1 &h) const override
      {
        h.Scale(1.0 / norm.value);
      }
  };

  struct FitParams : FitParam1D<FitParams> {
    double norm,
           lam,
           radius;

    FitParams(const double *par)
      : norm(par[NORM_PARAM_IDX])
      , lam(par[LAM_PARAM_IDX])
      , radius(par[R_PARAM_IDX])
      { }

    FitParams(const FitParams &) = default;

    FitParams(const FitResult &res)
      : norm(res.norm)
      , lam(res.lam)
      , radius(res.radius)
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
            || INVALID(radius);
        #undef INVALID
      }

    std::string
    print() const
      {
        std::string result;
        result += Form("Radius: %f \n", radius);
        result += Form("Lambda: %f\n", lam);
        result += Form("Norm: %f\n", norm);
        return result;
      }

    std::string
    __repr__() const
      {
        return Form("<Fitter1DGauss::FitParam radius=%f lambda=%f norm=%f>",
                    radius, lam, norm);
      }

    PyObject*
    as_dict() const
      {
        auto *dict = PyDict_New();
        PyDict_SetItemString(dict, "radius", PyFloat_FromDouble(radius));
        PyDict_SetItemString(dict, "lam", PyFloat_FromDouble(lam));
        PyDict_SetItemString(dict, "norm", PyFloat_FromDouble(norm));
        return dict;
      }

    double evaluate(const double q, const double K) const
      {
        return Fitter1DGauss::gauss(q, radius * radius, lam, K, norm);
      }

  };

  int
  setup_minuit(TMinuit &minuit) const override
    {
      int errflag = 0;
      minuit.mnparm(NORM_PARAM_IDX, "Norm", 0.2, 0.02, 0.0, 0.0, errflag);
      minuit.mnparm(LAM_PARAM_IDX, "Lam", 0.5, 0.01, 0.0, 0.0, errflag);
      minuit.mnparm(R_PARAM_IDX, "Radius", 1.0, 0.2, 0.0, 0.0, errflag);

      const double this_dbl = static_cast<double>((intptr_t)this);
      minuit.mnparm(DATA_PARAM_IDX, "DATA_PTR", this_dbl, 0, 0, INTPTR_MAX, errflag);
      minuit.FixParameter(DATA_PARAM_IDX);

      if (errflag != 0) {
        std::cerr << "Error setting paramters: " << errflag << "\n";
        throw std::runtime_error("Could not set Minuit parameters.");
      }

      return errflag;
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

  Fitter1DGauss(const TH1 &num, const TH1 &den, double limit)
    : Fitter1D(num, den, limit)
    { }

  Fitter1DGauss(TDirectory &tdir, double limit)
    : Fitter1D(tdir, limit)
    { }

  Fitter1DGauss(const Data1D &dat)
    : Fitter1D(dat)
    { }

  virtual ~Fitter1DGauss() = default;

  double resid_calc_chi2(const FitResult &fr)
    {
      auto params = fr.as_params();
      return Fitter1D::resid_calc(params, CalcChi2::resid_func);
    }

  double resid_calc_chi2(const FitParams &params)
    {
      return Fitter1D::resid_calc(params, CalcChi2::resid_func);
    }

  double resid_calc_pml(const FitParams &params)
    {
      return Fitter1D::resid_calc(params, CalcLoglike::resid_func);
    }

  double calc_chi2_residual(const FitResult &fr)
    {
      return resid_calc_chi2(fr);
    }

  double calc_chi2_residual_mrc(const FitResult &fr)
    {
      return resid_calc_chi2_mrc(fr);
    }

  double resid_calc_chi2_mrc(const FitResult &fr);

  DECLARE_FIT_METHODS(Fitter1D)

  // void fit_with_random_inits(TMinuit &minuit, FitResult &res, int);

  DECLARE_FILL_METHODS(TH1)

};

inline auto
Fitter1DGauss::FitResult::as_params() const -> FitParams
{
  return FitParams(*this);
}

#endif
