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
#include "Fitter3D.hpp"
#include "Data3D.hpp"


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

      #undef OUT
    }
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
    double PseudoRinv() const
      { return std::sqrt(Ro * Ro * gamma + Rs * Rs + Rl * Rl); }

    double gauss(const std::array<double, 3> &q, double K) const
      { return FitterGaussOSL::gauss(q, {Ro*Ro, Rs*Rs, Rl*Rl}, lam, K, norm); }

    void
    apply_to(TH3 &hist, TH3& qinv)
    {
      const int I = hist.GetNbinsX(),
                J = hist.GetNbinsY(),
                K = hist.GetNbinsZ();

      const double phony_r = PseudoRinv();
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
          // q = qinv.Interpolate(qo, qs, ql),
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

  static double
  gauss(const std::array<double, 3>& q, const FitParams &p, double K)
    { return p.gauss(q, K); }

  int
  setup_minuit(TMinuit &minuit)
  {
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
    return errflag;
  }

};
