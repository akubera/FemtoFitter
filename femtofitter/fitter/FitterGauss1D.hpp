///
/// \file fitter/FitterGauss1D.hpp
///

#pragma once

#ifndef FITTER_FITTERGAUSS1D_HPP
#define FITTER_FITTERGAUSS1D_HPP

#include "Fitter1D.hpp"
#include "Value.hpp"
#include "math/constants.hh"


/// \class FitterGauss1D
/// \brief Gaussian 1D fit
///
struct FitterGauss1D : public Fitter1D<FitterGauss1D> {

  using Super = Fitter1D<FitterGauss1D>;

  struct FitParams;
  struct FitInput;

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
  struct FitResult : FitResult1D<FitResult, FitterGauss1D> {

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
        return Form("<FitterGauss1D::FitResult radius=%f lambda=%f norm=%f>",
                    radius.value, lam.value, norm.value);
      }

    double evaluate(const double q, const double K) const
      { return FitterGauss1D::gauss(q, radius * radius, lam, K, norm); }

    FitParams as_params() const;

    void FillMinuit(TMinuit &minuit) const override
      {
        int errflag = 0;
        // minuit.mnparm(NORM_PARAM_IDX, "Norm", norm.value, norm.error * norm.error * 4, 0.0, 0.0, errflag);
        minuit.mnparm(NORM_PARAM_IDX, "Norm", norm.value, 0.005, 0.0, 0.0, errflag);
        minuit.mnparm(LAM_PARAM_IDX, "Lam", lam.value, .01, 0.0, 0.0, errflag);
        minuit.mnparm(R_PARAM_IDX, "Radius", radius.value, 0.2, 0.0, 0.0, errflag);
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
        return Form("<FitterGauss1D::FitParam radius=%f lambda=%f norm=%f>",
                    radius, lam, norm);
      }

    double evaluate(const double q, const double K) const
      {
        return FitterGauss1D::gauss(q, radius * radius, lam, K, norm);
      }

  };

  int
  setup_minuit(TMinuit &minuit) const override
    {
      int errflag = 0;
      minuit.mnparm(NORM_PARAM_IDX, "Norm", 0.2, 0.02, 0.0, 0.0, errflag);
      minuit.mnparm(LAM_PARAM_IDX, "Lam", 0.5, 0.01, 0.0, 0.0, errflag);
      minuit.mnparm(R_PARAM_IDX, "Radius", 5, 0.2, 0.0, 0.0, errflag);

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

  FitterGauss1D(const TH1 &num, const TH1 &den, double limit)
    : Fitter1D(num, den, limit)
    { }

  FitterGauss1D(TDirectory &tdir, double limit)
    : Fitter1D(tdir, limit)
    { }

  FitterGauss1D(const Data1D &data)
    : Fitter1D(data)
    { }

  virtual ~FitterGauss1D() = default;

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

  // void fit_with_random_inits(TMinuit &minuit, FitResult &res, int);

  void fill(const FitParams &p, TH1 &h, UInt_t npoints=1) const
    {
      p.fill(h, *fsi, npoints);
    }

  void fill(const FitResult &p, TH1 &h, UInt_t npoints=1) const
    {
      fill(p.as_params(), h, npoints);
    }

  void fill_smeared_fit(TH1 &h, const FitResult &fr)
    {
      fill_smeared_fit(h, fr.as_params());
    }

  void fill_smeared_fit(TH1 &h, const FitParams &p)
    {
      Fitter1D::fill_smeared_fit(h, p);
    }
};

#endif
