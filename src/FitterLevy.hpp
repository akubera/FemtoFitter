///
/// \file FitterLevy.hpp
///

#pragma once


#include <TFile.h>
#include <TH3.h>
#include <TMinuit.h>
#include <TGraph.h>

#include <array>
#include <vector>
#include <string>
#include <memory>
#include <valarray>
#include <iostream>

#include "CoulombHist.hpp"
#include "Fitter.hpp"
#include "Data3D.hpp"



/// \class FitterLevy
/// \brief Fit out-side-long with gaussian parameters
///
struct FitterLevy {
  using PMLCALC = ResidCalculatorPML<FitterLevy>;
  using CalcChi2 = ResidCalculatorChi2<FitterLevy>;

  /// constants used to lookup data from pointer
  enum {
    DATA_PARAM_IDX = 0,

    NORM_PARAM_IDX = 1,
    LAM_PARAM_IDX = 2,
    ROUT_PARAM_IDX = 3,
    RSIDE_PARAM_IDX = 4,
    RLONG_PARAM_IDX = 5,
    ALPHA_PARAM_IDX = 6,
  };

  /// The associated fit data
  Data3D data;

  static double
  gauss(std::array<double, 3> q,
        std::array<double, 3> RSq,
        double lam,
        double alpha,
        double K=1.0,
        double norm=1.0)
  {
    const double
      Eo = std::pow(q[0] * q[0] * RSq[0] / HBAR_C_SQ, alpha/2),
      Es = std::pow(q[1] * q[1] * RSq[1] / HBAR_C_SQ, alpha/2),
      El = std::pow(q[2] * q[2] * RSq[2] / HBAR_C_SQ, alpha/2),

      gauss = 1.0 + std::exp(-(Eo + Es + El)),
      result = (1.0 - lam) + lam * K * gauss;

    return norm * result;
  }

  /// \class FitResult
  /// \brief Values and stderr from minuit results
  ///
  struct FitResult {
    Value lam,
          norm,
          alpha,
          Ro,
          Rs,
          Rl;

    FitResult(const TMinuit &minuit)
      : lam(minuit, LAM_PARAM_IDX)
      , norm(minuit, NORM_PARAM_IDX)
      , alpha(minuit, ALPHA_PARAM_IDX)
      , Ro(minuit, ROUT_PARAM_IDX)
      , Rs(minuit, RSIDE_PARAM_IDX)
      , Rl(minuit, RLONG_PARAM_IDX)
    {
    }

    void print() const
    {
      std::cout << __str__();
    }

    std::string
    __str__() const
    {
      std::vector<char> buff(1000);

      snprintf(buff.data(), buff.size(),
             "Fit-Result:\n"
             "  Ro=%0.4f ± %0.4f \n"
             "  Rs=%0.4f ± %0.4f\n"
             "  Rl=%0.4f ± %0.4f\n"
             "  lam=%0.4f ± %0.4f (%g, %g)\n"
             "  alpha=%0.4f ± %0.4f (%g, %g)\n"
             "  norm=%0.4f ± %0.4f\n"
             " -------------\n", Ro.first, Ro.second,
                                 Rs.first, Rs.second,
                                 Rl.first, Rl.second,
                                 lam.first, lam.second, lam.first - lam.second, lam.first + lam.second,
                                 alpha.first, alpha.second, alpha.first - alpha.second, alpha.first + alpha.second,
                                 norm.first, norm.second);
      return buff.data();
    }

    std::map<std::string, double>
    as_map() const
    {
      #define OUT(__name) {#__name, __name.first}, { # __name "_err", __name.second}

      return {
        OUT(Ro), OUT(Rs), OUT(Rl), OUT(lam), OUT(alpha), OUT(norm)
      };

      #undef OUT
    }
  };

  /// \brief 3D Levy fit parameters
  ///
  ///
  struct FitParams {
    double norm, lam, alpha;
    double Ro, Rs, Rl;
    double gamma {1.0};

    FitParams(double *par)
      : norm(par[NORM_PARAM_IDX])
      , lam(par[LAM_PARAM_IDX])
      , alpha(par[ALPHA_PARAM_IDX])
      , Ro(par[ROUT_PARAM_IDX])
      , Rs(par[RSIDE_PARAM_IDX])
      , Rl(par[RLONG_PARAM_IDX])
    {
    }

    FitParams(const FitResult &res)
      : norm(res.norm)
      , lam(res.lam)
      , alpha(res.alpha)
      , Ro(res.Ro)
      , Rs(res.Rs)
      , Rl(res.Rl)
    {
    }

    bool is_invalid() const
    {
      return Ro < 0
          || Rs < 0
          || Rl < 0
          || lam < 0
          || alpha < 0
          || norm < 0
          || std::isnan(Ro)
          || std::isnan(Rs)
          || std::isnan(Rl)
          || std::isnan(lam)
          || std::isnan(alpha)
          || std::isnan(norm);
    }

    /// Return calculated Rinv: $\sqrt{Ro^2 \gamma + Rs^2 + Rl^2}$
    double PseudoRinv() const
      { return std::sqrt(Ro * Ro * gamma + Rs * Rs + Rl * Rl); }

    double gauss(std::array<double, 3> q, double K) const
      {
        std::array<double, 3> Rsq = {Ro*Ro, Rs*Rs, Rl*Rl};
        return FitterLevy::gauss(q, Rsq, lam, alpha, K, norm);
      }
  };

  /// Utility function for building fitter with tdirectory in file
  /// at specified path
  static std::unique_ptr<FitterLevy>
  From(TFile &file, const std::string &path, double limit=0.0)
  {
    auto tdir = static_cast<TDirectory*>(file.Get(path.c_str()));
    if (!tdir) {
      return nullptr;
    }
    return From(*tdir, limit);
  }

  /// Construct from ("num", "den", "qinv") histograms in tdirectory
  /// limit sets the fit-range.
  ///
  static std::unique_ptr<FitterLevy>
  From(TDirectory &tdir, double limit=0.0)
  {
    TH3 *num = static_cast<TH3*>(tdir.Get("num")),
        *den = static_cast<TH3*>(tdir.Get("den")),
        *qinv = static_cast<TH3*>(tdir.Get("qinv"));

    if (!num || !den || !qinv) {
      return nullptr;
      // throw std::runtime_error("Error loading FitterLevy histograms from path '" + path + "'");
    }
#if __cplusplus > 201103L
    return make_unique<FitterLevy>(*num, *den, *qinv, limit);
#else
    return std::unique_ptr<FitterLevy>(new FitterLevy(*num, *den, *qinv, limit));
#endif
  }

  /// Construct fitter from numerator denominator qinv histograms
  /// and a fit-range limit
  ///
  FitterLevy(TH3 &n, TH3 &d, TH3 &q, double limit=0.0)
    : data(n, d, q, limit)
  {
  }

  /// Number of entries in fitter
  std::size_t size() const
    { return data.size(); }

  static double
  gauss(std::array<double, 3> q, const FitParams &p, double K)
    { return p.gauss(q, K); }

  template <typename ResidFunc>
  double resid_calc(const FitParams &p, ResidFunc resid_calc) const
  {
    double retval = 0;

    double phony_r = p.PseudoRinv();
    auto coulomb_factor = CoulombHist::GetHistWithRadius(phony_r);

    auto &qout = data.qspace[0],
         &qside = data.qspace[1],
         &qlong = data.qspace[2];

#if __cplusplus > 201103L
    auto Kfsi = [&c=coulomb_factor] (double q) {
#else
    auto &c = coulomb_factor;
    auto Kfsi = [&c] (double q) {
#endif
      double coulomb = c.Interpolate(q);
      // if ((rand() * 1.0 / RAND_MAX) < 1e-5) {
      //   printf("K(%g) = %g\n", q, coulomb);
      // }
      return coulomb;
    };

    for (size_t i=0; i<size(); ++i) {
      const double
        qo = qout[i],
        qs = qside[i],
        ql = qlong[i],
        n = data.num[i],
        d = data.den[i],
        q = data.qinv[i];

      const double
        // CF = gauss({qo, qs, ql}, {p.Ro, p.Rs, p.Rl}, p.lam, Kfsi(q), p.norm);
        CF = p.gauss({qo, qs, ql}, Kfsi(q));

      retval += resid_calc(n, d, CF);
    }

    return retval;
  }

  double
  resid_chi2(const FitParams &p) const
    { return resid_calc(p, chi2_calc); }

  double
  resid_pml(const FitParams &p) const
    { return resid_calc(p, loglikelihood_calc); }

  static void
  fit_func(Int_t &,
           Double_t *,
           Double_t &retval,
           Double_t *par,
           Int_t)
  {
    // returned when invalid input or output occurs
    static const double BAD_VALUE = 3e99;
    const auto &data = *(const FitterLevy*)(intptr_t)(par[DATA_PARAM_IDX]);

    FitParams params(par);
    if (params.is_invalid()) {
      retval = BAD_VALUE;
      return;
    }

    retval = data.resid_chi2(params);
    std::cout << "< " << retval << " ("
              << params.Ro << " "
              << params.Rs << " "
              << params.Rl << ") "
              << params.lam << " "
              << params.alpha << " "
              << params.norm << "\n";
  }

  int
  setup_minuit(TMinuit &minuit)
  {
    int errflag = 0;
    minuit.mnparm(NORM_PARAM_IDX, "Norm", 0.25, 0.02, 0.0, 0.0, errflag);
    minuit.mnparm(LAM_PARAM_IDX, "Lam", 0.2, 0.1, 0.0, 1.0, errflag);
    minuit.mnparm(ROUT_PARAM_IDX, "Ro", 2.0, 1.0, 0.0, 0.0, errflag);
    minuit.mnparm(RSIDE_PARAM_IDX, "Rs", 2.0, 1.0, 0.0, 0.0, errflag);
    minuit.mnparm(RLONG_PARAM_IDX, "Rl", 2.0, 1.0, 0.0, 0.0, errflag);
    minuit.mnparm(ALPHA_PARAM_IDX, "Alpha", 1.8, 0.1, 0.0, 0.0, errflag);

    const double this_dbl = static_cast<double>((intptr_t)this);
    minuit.mnparm(DATA_PARAM_IDX, "DATA_PTR", this_dbl, 0, 0, INTPTR_MAX, errflag);

    minuit.FixParameter(DATA_PARAM_IDX);
    if (errflag != 0) {
      std::cerr << "Error setting paramters: " << errflag << "\n";
      throw std::runtime_error("Could not set Minuit parameters.");
    }

    return errflag;
  }

  template <typename F>
  FitResult
  fit(F fcn, double fit_method=1.0)
  {
    TMinuit minuit;
    minuit.SetPrintLevel(-1);
    setup_minuit(minuit);
    minuit.SetFCN(fcn);

    return do_fit_minuit(minuit, fit_method);
  }

  template <typename ResidCalculator_t>
  FitResult
  fit(double fit_method)
  {
    TMinuit minuit;
    minuit.SetPrintLevel(-1);
    setup_minuit(minuit);
    minuit.SetFCN(minuit_f<ResidCalculator_t>);
    return do_fit_minuit(minuit, fit_method);
  }

  FitResult
  do_fit_minuit(TMinuit &minuit, double fit_method)
  {
    double strat_args[] = {1.0};
    double migrad_args[] = {2000.0, fit_method};
    double hesse_args[] = {2000.0, 1.0};

    int errflag;
    minuit.mnexcm("SET STRategy", strat_args, 1, errflag);
    minuit.mnexcm("MIGRAD", migrad_args, 2, errflag);

    strat_args[0] = 2.0;
    minuit.mnexcm("SET STRategy", strat_args, 1, errflag);
    minuit.mnexcm("MIGRAD", migrad_args, 2, errflag);
    minuit.mnexcm("HESSE", hesse_args, 1, errflag);

    return FitResult(minuit);
  }

  FitResult
  fit_pml()
    { return fit<PMLCALC>(0.5); }

  FitResult
  fit_chi2()
    { return fit<CalcChi2>(1.0); }

  FitResult
  fit()
    { return fit(fit_func, 1.0); }

  std::vector<double>
  num_as_vec()
    { return std::vector<double>(std::begin(data.num), std::end(data.num)); }

};
