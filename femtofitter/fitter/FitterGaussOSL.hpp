///
/// \file fitter/FitterGaussOSL.hpp
///

#pragma once


#include "CoulombHist.hpp"
#include "CalculatorResid.hpp"
#include "Value.hpp"
#include "math/constants.hh"
#include "Fitter3D.hpp"

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


/// \class FitterGaussOSL
/// \brief Fit out-side-long with gaussian parameters
///
struct FitterGaussOSL : public Fitter3D<FitterGaussOSL> {

  /// constants used to lookup data from pointer
  enum {
    DATA_PARAM_IDX = 0,

    NORM_PARAM_IDX = 1,
    LAM_PARAM_IDX = 2,
    ROUT_PARAM_IDX = 3,
    RSIDE_PARAM_IDX = 4,
    RLONG_PARAM_IDX = 5,
  };

  struct FitParams;
  struct FitResult;

  static unsigned char CountParams()
    { return 5; }

  static double
  gauss(const std::array<double, 3> &q,
        const std::array<double, 3> &RSq,
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

  int degrees_of_freedom() const
    { return data.size() - 5; }

  /// \class FitResult
  /// \brief Values and stderr from minuit results
  ///
  struct FitResult {
    Value lam,
          norm,
          Ro,
          Rs,
          Rl;

    FitResult()
      : lam({0, 0})
      , norm({0, 0})
      , Ro({0, 0})
      , Rs({0, 0})
      , Rl({0, 0})
      {}

    FitResult(const FitResult &orig) = default;

    FitResult(TMinuit &minuit)
      : lam(minuit, LAM_PARAM_IDX)
      , norm(minuit, NORM_PARAM_IDX)
      , Ro(minuit, ROUT_PARAM_IDX)
      , Rs(minuit, RSIDE_PARAM_IDX)
      , Rl(minuit, RLONG_PARAM_IDX)
    {
    }

    void print() const
    {
      printf("Fit-Result:\n"
             "  Ro=%0.4f ± %0.4f \n"
             "  Rs=%0.4f ± %0.4f\n"
             "  Rl=%0.4f ± %0.4f\n"
             "  lam=%0.4f ± %0.4f (%g, %g)\n"
             "  norm=%0.4f ± %0.4f\n"
             " -------------\n", Ro.value, Ro.error,
                                 Rs.value, Rs.error,
                                 Rl.value, Rl.error,
                                 lam.value, lam.error, lam.value - lam.error, lam.value + lam.error,
                                 norm.value, norm.error);
    }

    std::map<std::string, double>
    as_map() const
    {
      #define OUT(__name) {#__name, __name.value}, { # __name "_err", __name.error}

      return {
        OUT(Ro), OUT(Rs), OUT(Rl), OUT(lam), OUT(norm)
      };

      #undef OUT
    }

    FitParams as_params() const;
  };

  /// \brief 3D Gaussian fit parameters
  ///
  ///
  struct FitParams {
    double norm, lam;
    double Ro, Rs, Rl;

    FitParams(double *par)
      : norm(par[NORM_PARAM_IDX])
      , lam(par[LAM_PARAM_IDX])
      , Ro(par[ROUT_PARAM_IDX])
      , Rs(par[RSIDE_PARAM_IDX])
      , Rl(par[RLONG_PARAM_IDX])
    {
    }

    FitParams(const FitResult &res)
      : norm(res.norm)
      , lam(res.lam)
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
          || norm < 0
          || std::isnan(Ro)
          || std::isnan(Rs)
          || std::isnan(Rl)
          || std::isnan(lam)
          || std::isnan(norm);
    }

    /// Return calculated Rinv: $\sqrt{Ro^2 \gamma + Rs^2 + Rl^2}$
    double PseudoRinv(double gamma) const
      { return std::sqrt((Ro * Ro * gamma + Rs * Rs + Rl * Rl) / 3.0); }

    double gauss(const std::array<double, 3> &q, double K) const
      { return FitterGaussOSL::gauss(q, {Ro*Ro, Rs*Rs, Rl*Rl}, lam, K, norm); }

    void
    apply_to(TH3 &hist, const Data3D &data)
      {
        const int I = hist.GetNbinsX(),
                  J = hist.GetNbinsY(),
                  K = hist.GetNbinsZ();

        const TAxis
          &qout_ax = *hist.GetXaxis(),
          &qside_ax = *hist.GetYaxis(),
          &qlong_ax = *hist.GetZaxis();

        std::unique_ptr<TH3> qinv { (TH3*)hist.Clone("__qinv") };

        for (const auto &datum : data) {
          int i = qout_ax.FindBin(datum.qo),
              j = qside_ax.FindBin(datum.qs),
              k = qlong_ax.FindBin(datum.ql);

          qinv->SetBinContent(i, j, k, datum.qinv);
        }

        for (int k=1; k<=K; ++k)
        for (int j=1; j<=J; ++j)
        for (int i=1; i<=I; ++i) {
          if (qinv->GetBinContent(i, j, k) == 0.0) {
            double q = qinv->GetBinContent(i-1, j, k)
                     + qinv->GetBinContent(i+1, j, k)
                     + qinv->GetBinContent(i, j, k-1)
                     + qinv->GetBinContent(i, j, k+1)
                     + qinv->GetBinContent(i, j-1, k)
                     + qinv->GetBinContent(i, j+1, k);
            qinv->SetBinContent(i, j, k, q / 6.0);
          }
        }

        return apply_to(hist, *qinv, data.gamma);
      }

    void
    apply_to(TH3 &hist, TH3& qinv, double gamma)
    {
      const int I = hist.GetNbinsX(),
                J = hist.GetNbinsY(),
                K = hist.GetNbinsZ();

      const double phony_r = PseudoRinv(gamma);
      auto coulomb_factor = CoulombHist::GetHistWithRadius(phony_r);

      const TAxis &qout = *hist.GetXaxis(),
                  &qside = *hist.GetYaxis(),
                  &qlong = *hist.GetZaxis();

      for (int k=1; k<=K; ++k)
      for (int j=1; j<=J; ++j)
      for (int i=1; i<=I; ++i) {
        const double
          qo = qout.GetBinCenter(i),
          qs = qside.GetBinCenter(j),
          ql = qlong.GetBinCenter(k),
          q = qinv.GetBinContent(i, j, k),
          K = coulomb_factor.Interpolate(q);

        hist.SetBinContent(i,j,k, hist.GetBinContent(i,j,k) * gauss({qo, qs, ql}, K));
      }
    }
  };

  /// Construct fitter from numerator denominator qinv histograms
  /// and a fit-range limit
  ///
  FitterGaussOSL(TH3 &n, TH3 &d, TH3 &q, double limit=0.0)
    : Fitter3D(n, d, q, limit)
  {
  }

  static std::unique_ptr<FitterGaussOSL>
  FromDirectory(TDirectory &dir, double limit=0.0)
    {
      auto data = Data3D::FromDirectory(dir, limit);
      auto fitter = std::make_unique<FitterGaussOSL>(std::move(data));
      fitter->paramhints = std::make_unique<ParamHints>(dir);
      return fitter;
    }

  FitterGaussOSL(const Data3D &data)
    : Fitter3D(data)
    {
    }

  FitterGaussOSL(Data3D &&data)
    : Fitter3D(std::move(data))
    {
    }

  FitterGaussOSL(std::unique_ptr<Data3D> data)
    : Fitter3D(std::move(data))
    {
    }

  static double
  gauss(const std::array<double, 3>& q, const FitParams &p, double K)
    { return p.gauss(q, K); }

  int
  setup_minuit(TMinuit &minuit) const override
  {
    int errflag = 0;
    minuit.mnparm(NORM_PARAM_IDX, "Norm", 0.25, 0.02, 0.0, 0.0, errflag);
    minuit.mnparm(LAM_PARAM_IDX, "Lam", 0.2, 0.05, 0.0, 1.0, errflag);
    minuit.mnparm(ROUT_PARAM_IDX, "Ro", 2.0, 0.5, 0.0, 0.0, errflag);
    minuit.mnparm(RSIDE_PARAM_IDX, "Rs", 2.0, 0.5, 0.0, 0.0, errflag);
    minuit.mnparm(RLONG_PARAM_IDX, "Rl", 2.0, 0.5, 0.0, 0.0, errflag);

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
    int idx = (name == "Ro") ? ROUT_PARAM_IDX
            : (name == "Rs") ? RSIDE_PARAM_IDX
            : (name == "Rl") ? RLONG_PARAM_IDX
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

  FitResult fit_chi2()
    { return Fitter3D::fit_chi2(); }

  FitResult fit_pml()
    { return Fitter3D::fit_pml(); }

  FitResult fit()
    { return Fitter3D::fit(); }

  void fit_with_random_inits(TMinuit &minuit, FitResult &res);

};


// FitterGaussOSL::FitParams
inline auto
FitterGaussOSL::FitResult::as_params() const -> FitParams
{
  return FitParams(*this);
}

// template<>
void FitterGaussOSL::fit_with_random_inits(TMinuit &minuit, FitResult &res)
{
  int errflag = 0;
  std::cout << "paramhints: " << paramhints.get() << "\n";

  minuit.mnparm(LAM_PARAM_IDX, "Lam", paramhints->GenLam(), 0.01, 0.0, 0.0, errflag);
  minuit.mnparm(ROUT_PARAM_IDX, "Ro", paramhints->GenRo(), 0.1, 0.0, 0.0, errflag);
  minuit.mnparm(RSIDE_PARAM_IDX, "Rs", paramhints->GenRs(), 0.1, 0.0, 0.0, errflag);
  minuit.mnparm(RLONG_PARAM_IDX, "Rl", paramhints->GenRl(), 0.1, 0.0, 0.0, errflag);

  res = do_fit_minuit(minuit, 11);  // -> Impl::FitResult

  // res = FitResult(minuit);
}
