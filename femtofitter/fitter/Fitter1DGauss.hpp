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

  using Super = Fitter1D<Fitter1DGauss>;
  struct FitParams;
  struct FitResult;

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

    FitResult(PyObject *pyobj);

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

    /// Build python list of tuples from data
    PyObject* __iter__() const override;

    /// Build python dictionary from data
    PyObject* as_dict() const;

    void FillMinuit(TMinuit &minuit) const override
      {
        int errflag = 0;
        minuit.mnparm(LAM_PARAM_IDX, "Lam", lam.value, .01, 0.0, 0.0, errflag);
        minuit.mnparm(R_PARAM_IDX, "Radius", radius.value, 0.2, 0.0, 0.0, errflag);
        minuit.mnparm(NORM_PARAM_IDX, "Norm", norm.value, 0.005, 0.0, 0.0, errflag);
      }
  };

  struct FitParams : FitParam1D<FitParams> {
    double norm,
           lam,
           radius;

    double evaluate(const double q, const double K) const
      {
        return Fitter1DGauss::gauss(q, radius * radius, lam, K, norm);
      }

    void Normalize(TH1 &h) const
      {
        h.Scale(1.0 / norm);
      }

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

    PyObject* as_dict() const;
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
  set_and_fix_variable(TMinuit &minuit, TString name, double val)
    {
      name.ToLower();
      int idx = (name == "r") || (name == "radius") ? R_PARAM_IDX
              : (name == "lam") ? LAM_PARAM_IDX
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

  DECLARE_FIT_METHODS(Fitter1D)
  DECLARE_RESID_METHODS(Fitter1D);
  DECLARE_FILL_METHODS(TH1)

};

#endif
