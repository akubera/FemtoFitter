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
  struct FitResult {
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
  };

  struct FitParams {
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

    double gauss(const double q, const double K) const
      { return FitterGauss1D::gauss(q, radius * radius, lam, K, norm); }

    double evaluate(const double q, const double K) const
      { return gauss(q, K); }

    /// Fill histogram with no FSI factor
    void fill(TH1 &h) const
      {
        const TAxis &xaxis = *h.GetXaxis();
        for (int i=1; i <= xaxis.GetLast(); ++i) {
          double q = xaxis.GetBinCenter(i);
          double k = 1.0;
          double cf = gauss(q, k);
          h.SetBinContent(i, cf);
        }
      }

    /// Fill histogram with average of N-points per bin
    ///
    void fill(TH1 &h, FsiCalculator &fsi, UInt_t npoints=1) const
      {
        const TAxis &xaxis = *h.GetXaxis();
        auto Kfsi = fsi.ForRadius(radius);

        for (int i=1; i <= xaxis.GetLast(); ++i) {
          const double
            qlo = xaxis.GetBinLowEdge(i),
            qhi = xaxis.GetBinUpEdge(i),
            qstep = (qhi - qlo) / npoints,
            qstart = qlo + qstep / 2;

          double sum = 0.0;
          for (double q=qstart; q < qhi; q += qstep) {
            double k = Kfsi(q);
            sum += evaluate(q, k);
          }

          const double cf = sum / npoints;
          h.SetBinContent(i, cf);
        }
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

  // FitResult fit_chi2()
  //   { return Fitter1D::fit_chi2(); }

  FitResult fit_pml()
    { return Fitter1D::fit_pml(); }

  // void fit_with_random_inits(TMinuit &minuit, FitResult &res, int);
};

auto
FitterGauss1D::FitResult::as_params() const -> FitParams
{
  return FitParams(*this);
}

#endif
