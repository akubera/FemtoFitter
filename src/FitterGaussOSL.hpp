///
/// \file fitter-gauss-osl.hpp
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



/// \class FitterGaussOSL
/// \brief Fit out-side-long with gaussian parameters
///
struct FitterGaussOSL {
  using PMLCALC = ResidCalculatorPML<FitterGaussOSL>;
  using CalcChi2 = ResidCalculatorChi2<FitterGaussOSL>;

  /// constants used to lookup data from pointer
  enum {
    DATA_PARAM_IDX = 0,

    NORM_PARAM_IDX = 1,
    LAM_PARAM_IDX = 2,
    ROUT_PARAM_IDX = 3,
    RSIDE_PARAM_IDX = 4,
    RLONG_PARAM_IDX = 5,
  };

  static double
  gauss(std::array<double, 3> q,
        std::array<double, 3> RSq,
        double lam,
        double K=1.0,
        double norm=1.0)
  {
    const double
      Eo = q[0] * q[0] * RSq[0],
      Es = q[1] * q[1] * RSq[1],
      El = q[2] * q[2] * RSq[2],
      gauss = 1.0 + std::exp(-(Eo + Es + El) / HBAR_C_SQ),
      result = (1.0 - lam) + lam * K * gauss;

    return norm * result;
  }

  struct Data {
    std::array<std::valarray<double>, 3> qspace;

    std::valarray<double> num,
                          den,
                          qinv;

    size_t size() const
      { return num.size(); }
  };

  /// \brief 3D Gaussian fit parameters
  ///
  ///
  struct FitParams {
    double norm, lam;
    double Ro, Rs, Rl;
    double gamma {1.0};

    FitParams(double *par)
      : norm(par[NORM_PARAM_IDX])
      , lam(par[LAM_PARAM_IDX])
      , Ro(par[ROUT_PARAM_IDX])
      , Rs(par[RSIDE_PARAM_IDX])
      , Rl(par[RLONG_PARAM_IDX])
    {
    }

    bool is_invalid() const
    {
      return Ro < 0
          || Rs < 0
          || Rl < 0
          || lam < 0
          || norm < 0
          || std::isnan(Ro)
          || std::isnan(Rs)
          || std::isnan(Rl)
          || std::isnan(lam)
          || std::isnan(norm);
    }

    /// Return calculated Rinv: $\sqrt{Ro^2 \gamma + Rs^2 + Rl^2}$
    double PseudoRinv() const
      { return std::sqrt(Ro * Ro * gamma + Rs * Rs + Rl * Rl); }

    double gauss(std::array<double, 3> q, double K) const
      {
        std::array<double, 3> Rsq = {Ro*Ro, Rs*Rs, Rl*Rl};
        return FitterGaussOSL::gauss(q, Rsq, lam, K, norm);
      }
  };

  /// \class FitResult
  /// \brief Values and stderr from minuit results
  ///
  struct FitResult {
    Value lam,
          norm,
          Ro,
          Rs,
          Rl;

    FitResult(TMinuit &minuit)
    {
      auto get_param = [&minuit=minuit](int idx, Value &v) {
        minuit.GetParameter(idx, v.first, v.second);
      };

      /*
      auto get_sqrt_param = [&minuit=minuit](int idx, Value &v) {
        minuit.GetParameter(idx, v.first, v.second);
        v.first = std::sqrt(v.first);
        v.second = std::sqrt(v.second) / (2 * v.first);
      };
      */

      get_param(NORM_PARAM_IDX, norm);
      get_param(LAM_PARAM_IDX, lam);
      get_param(ROUT_PARAM_IDX, Ro);
      get_param(RSIDE_PARAM_IDX, Rs);
      get_param(RLONG_PARAM_IDX, Rl);
    }

    void print() const
    {
      printf("Fit-Result:\n"
             "  Ro=%0.4f ± %0.4f \n"
             "  Rs=%0.4f ± %0.4f\n"
             "  Rl=%0.4f ± %0.4f\n"
             "  lam=%0.4f ± %0.4f (%g, %g)\n"
             "  norm=%0.4f ± %0.4f\n"
             " -------------\n", Ro.first, Ro.second,
                                 Rs.first, Rs.second,
                                 Rl.first, Rl.second,
                                 lam.first, lam.second, lam.first - lam.second, lam.first + lam.second,
                                 norm.first, norm.second);
    }

    std::map<std::string, double>
    as_map() const
    {
      #define OUT(__name) {#__name, __name.first}, { # __name "_err", __name.second}

      return {
        OUT(Ro), OUT(Rs), OUT(Rl), OUT(lam), OUT(norm)
      };
        // {"Ro", Ro.first}, {"Ro_err", Ro.second}};

      #undef OUT
    }
  };

  /// The associated fit data
  Data data;

  /// Utility function for building fitter with tdirectory in file
  /// at specified path
  static std::unique_ptr<FitterGaussOSL>
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
  static std::unique_ptr<FitterGaussOSL>
  From(TDirectory &tdir, double limit=0.0)
  {
    TH3 *num = static_cast<TH3*>(tdir.Get("num")),
        *den = static_cast<TH3*>(tdir.Get("den")),
        *qinv = static_cast<TH3*>(tdir.Get("qinv"));

    if (!num || !den || !qinv) {
      return nullptr;
      // throw std::runtime_error("Error loading FitterGaussOSL histograms from path '" + path + "'");
    }

    return make_unique<FitterGaussOSL>(*num, *den, *qinv, limit);
  }

  /// Construct fitter from numerator denominator qinv histograms
  /// and a fit-range limit
  ///
  FitterGaussOSL(TH3 &n, TH3 &d, TH3 &q, double limit=0.0)
  {
    size_t nbinsx = n.GetNbinsX(),
           nbinsy = n.GetNbinsY(),
           nbinsz = n.GetNbinsZ(),
           nbins = nbinsx * nbinsy * nbinsz;

    TAxis *xaxis = n.GetXaxis(),
          *yaxis = n.GetYaxis(),
          *zaxis = n.GetZaxis();

    double low_lim = xaxis->GetBinLowEdge(limit <= 0.0 ? 1 : xaxis->FindBin(-limit)),
          high_lim = xaxis->GetBinUpEdge(limit <= 0.0 ? xaxis->GetNbins() : xaxis->FindBin(limit));


    std::vector<double> qout, qside, qlong, nv, dv, qv;
    qout.reserve(nbins);
    qside.reserve(nbins);
    qlong.reserve(nbins);

    qv.reserve(nbins);
    nv.reserve(nbins);
    dv.reserve(nbins);

    for (size_t k=1; k<=nbinsz; ++k) {
      for (size_t j=1; j<=nbinsy; ++j) {
        for (size_t i=1; i<=nbinsx; ++i) {
          double d_value = d.GetBinContent(i, j, k);
          if (d_value == 0.0) {
            continue;
          }

          const double
            qo = xaxis->GetBinCenter(i),
            qs = yaxis->GetBinCenter(j),
            ql = zaxis->GetBinCenter(k);

          if ((qo < low_lim) || (high_lim < qo) ||
              (qs < low_lim) || (high_lim < qs) ||
              (ql < low_lim) || (high_lim < ql)) {
            continue;
          }

          qout.push_back(qo);
          qside.push_back(qs);
          qlong.push_back(ql);

          qv.push_back(q.GetBinContent(i, j, k));
          nv.push_back(n.GetBinContent(i, j, k));
          dv.push_back(d_value);
        }
      }
    }

    auto count = qout.size();
    data.qspace[0] = std::valarray<double>(qout.data(), count);
    data.qspace[1] = std::valarray<double>(qside.data(), count);
    data.qspace[2] = std::valarray<double>(qlong.data(), count);

    data.num = std::valarray<double>(nv.data(), count);
    data.den = std::valarray<double>(dv.data(), count);
    data.qinv = std::valarray<double>(qv.data(), count);
  }

  /// \class RangedFitter
  /// \brief A fitter with limited fit-range
  ///
  /// redundant...
  ///
  struct RangedFitter {
    double fLimit;

    Data data;

    RangedFitter(FitterGaussOSL &fitter, double limit)
      : fLimit(limit)
    {
      // auto mask = (-limit <= fitter.qspace[0] && fitter.qspace[0] <= limit)
      //          && (-limit <= fitter.qspace[1] && fitter.qspace[1] <= limit)
      //          && (-limit <= fitter.qspace[2] && fitter.qspace[2] <= limit);

      std::valarray<bool> mask(fitter.size());

      // build mask elements outside the limit
      auto &qspace = fitter.data.qspace;
      for (size_t i=0; i < fitter.size(); ++i) {
        mask[i] = (-limit <= qspace[0][i] && qspace[0][i] <= limit)
               && (-limit <= qspace[1][i] && qspace[1][i] <= limit)
               && (-limit <= qspace[2][i] && qspace[2][i] <= limit);
      }

      data.qspace[0] = qspace[0][mask];
      data.qspace[1] = qspace[1][mask];
      data.qspace[2] = qspace[2][mask];

      data.num = fitter.data.num[mask];
      data.den = fitter.data.den[mask];
      data.qinv = fitter.data.qinv[mask];
    }

    /// Number of entries in fitter
    std::size_t size() const
      { return data.size(); }

  };

  /// Number of entries in fitter
  std::size_t size() const
    { return data.size(); }

  static double
  gauss(std::array<double, 3> q, const FitParams &p, double K)
  {
    return p.gauss(q, K);
  }

  template <typename ResidFunc>
  double resid_calc(const FitParams &p, ResidFunc resid_calc) const
  {
    double retval = 0;

    double phony_r = p.PseudoRinv();
    auto coulomb_factor = CoulombHist::GetHistWithRadius(phony_r);

    auto &qout = data.qspace[0],
         &qside = data.qspace[1],
         &qlong = data.qspace[2];

    auto Kfsi = [&c=coulomb_factor] (double q) {
      double coulomb =  c.Interpolate(q);
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
    const auto &data = *(const FitterGaussOSL*)(intptr_t)(par[DATA_PARAM_IDX]);

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
              << params.norm << "\n";
  }

  template <typename F>
  FitResult
  fit(F fcn, double fit_method=1.0)
  {
    TMinuit minuit;

    minuit.SetPrintLevel(-1);

    int errflag = 0;
    minuit.mnparm(NORM_PARAM_IDX, "Norm", 0.25, 0.02, 0.0, 0.0, errflag);
    minuit.mnparm(LAM_PARAM_IDX, "Lam", 0.2, 0.1, 0.0, 1.0, errflag);
    minuit.mnparm(ROUT_PARAM_IDX, "Ro", 2.0, 1.0, 0.0, 0.0, errflag);
    minuit.mnparm(RSIDE_PARAM_IDX, "Rs", 2.0, 1.0, 0.0, 0.0, errflag);
    minuit.mnparm(RLONG_PARAM_IDX, "Rl", 2.0, 1.0, 0.0, 0.0, errflag);

    const double this_dbl = static_cast<double>((intptr_t)this);
    minuit.mnparm(DATA_PARAM_IDX, "DATA_PTR", this_dbl, 0, 0, INTPTR_MAX, errflag);

    minuit.FixParameter(DATA_PARAM_IDX);
    if (errflag != 0) {
      std::cerr << "Error setting paramters: " << errflag << "\n";
      throw std::runtime_error("Could not set Minuit parameters.");
    }

    minuit.SetFCN(fcn);

    double strat_args[] = {1.0};
    double migrad_args[] = {2000.0, fit_method};
    double hesse_args[] = {2000.0, 1.0};

    minuit.mnexcm("SET STRategy", strat_args, 1, errflag);
    minuit.mnexcm("MIGRAD", migrad_args, 2, errflag);

    strat_args[0] = 2.0;
    minuit.mnexcm("SET STRategy", strat_args, 1, errflag);
    minuit.mnexcm("MIGRAD", migrad_args, 2, errflag);

    minuit.mnexcm("HESSE", hesse_args, 1, errflag);

    /*
    TGraph *g = nullptr;

    for (int i=3; i>0; --i) {
      minuit.SetErrorDef(i);
      TGraph *g = (TGraph*) minuit.Contour(10, 3, 4);
      g->Draw(i != 3 ? "SAME" : "");
    }
    */

    return FitResult(minuit);
  }

  template <typename F>
  auto get_fit_func()
  {
    return [] (Int_t &,
           Double_t *,
           Double_t &retval,
           Double_t *par,
           Int_t)
      {
        // returned when invalid input or output occurs
        static const double BAD_VALUE = 3e99;
        const auto &data = *(const FitterGaussOSL*)(intptr_t)(par[DATA_PARAM_IDX]);

        FitParams params(par);
        if (params.is_invalid()) {
          retval = BAD_VALUE;
          return;
        }

        // retval = data.resid_pml(params);
        // retval = resid_calc(data, params);
        retval = F::resid(data, params);
        // std::cout << "< " << retval << " ("
        //           << params.Ro << " "
        //           << params.Rs << " "
        //           << params.Rl << ") "
        //           << params.lam << " "
        //           << params.norm << "\n";
      };
  }

  FitResult
  fit_pml()
  {
    return fit(get_fit_func<PMLCALC>(), 0.5);
  }

  FitResult
  fit_chi2()
  {
    return fit(get_fit_func<CalcChi2>(), 1.0);
  }

  FitResult
  fit()
  {
    return fit(fit_func, 1.0);
    /*
    TMinuit minuit;

    int errflag = 0;
    minuit.mnparm(NORM_PARAM_IDX, "Norm", 0.25, 0.02, 0.0, 0.0, errflag);
    minuit.mnparm(LAM_PARAM_IDX, "Lam", 0.2, 0.1, 0.0, 1.0, errflag);
    minuit.mnparm(ROUT_PARAM_IDX, "Ro", 2.0, 1.0, 0.0, 0.0, errflag);
    minuit.mnparm(RSIDE_PARAM_IDX, "Rs", 2.0, 1.0, 0.0, 0.0, errflag);
    minuit.mnparm(RLONG_PARAM_IDX, "Rl", 2.0, 1.0, 0.0, 0.0, errflag);

    const double this_dbl = static_cast<double>((intptr_t)this);
    minuit.mnparm(DATA_PARAM_IDX, "DATA_PTR", this_dbl, 0, 0, INTPTR_MAX, errflag);

    minuit.FixParameter(DATA_PARAM_IDX);
    if (errflag != 0) {
      std::cerr << "Error setting paramters: " << errflag << "\n";
      throw std::runtime_error("Could not set Minuit parameters.");
    }

    minuit.SetFCN(fit_func);

    double args[] = {.50};
    minuit.mnexcm("SET STRategy", args, 1, errflag);

    double margs[] = {2000.0, 1.0};
    minuit.mnexcm("MIGRAD", margs, 2, errflag);

    // double hargs[] = {2000.0, 1.0};
    // minuit.mnexcm("HESSE", hargs, 1, errflag);


    return FitResult(minuit);
    */
  }

  std::vector<double> num_as_vec()
    { return std::vector<double>(std::begin(data.num), std::end(data.num)); }

};
