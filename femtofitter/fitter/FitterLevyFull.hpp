///
/// \file FitterLevyFull.hpp
///

#pragma once

#include "CoulombHist.hpp"
#include "Fitter3D.hpp"
#include "FitterTraits.hpp"


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


/// \class FitterLevy
/// \brief Fit out-side-long with gaussian parameters
///
struct FitterLevyFull : public Fitter3D<FitterLevyFull> {

  /// constants used to lookup data from pointer
  enum {
    DATA_PARAM_IDX = 0,

    NORM_PARAM_IDX = 1,
    LAM_PARAM_IDX = 2,
    ROUT_PARAM_IDX = 3,
    RSIDE_PARAM_IDX = 4,
    RLONG_PARAM_IDX = 5,
    ALPHAOUT_PARAM_IDX = 6,
    ALPHASIDE_PARAM_IDX = 7,
    ALPHALONG_PARAM_IDX = 8,
  };


  static unsigned char CountParams()
    { return 8; }

  static double
  gauss(const std::array<double, 3> &q,
        const std::array<double, 3> &RSq,
        double lam,
        double alpha_o,
        double alpha_s,
        double alpha_l,
        double K=1.0,
        double norm=1.0)
  {
    const double
      Eo = std::pow(q[0] * q[0] * RSq[0] / HBAR_C_SQ, alpha_o/2),
      Es = std::pow(q[1] * q[1] * RSq[1] / HBAR_C_SQ, alpha_s/2),
      El = std::pow(q[2] * q[2] * RSq[2] / HBAR_C_SQ, alpha_l/2),

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
          alphao,
          alphas,
          alphal,
          Ro,
          Rs,
          Rl;

    FitResult(const TMinuit &minuit)
      : lam(minuit, LAM_PARAM_IDX)
      , norm(minuit, NORM_PARAM_IDX)
      , alphao(minuit, ALPHAOUT_PARAM_IDX)
      , alphas(minuit, ALPHASIDE_PARAM_IDX)
      , alphal(minuit, ALPHALONG_PARAM_IDX)
      , Ro(minuit, ROUT_PARAM_IDX)
      , Rs(minuit, RSIDE_PARAM_IDX)
      , Rl(minuit, RLONG_PARAM_IDX)
    {
    }

    void print() const
      { std::cout << __str__(); }

    /// for printing within a python environment
    auto __str__() const -> std::string
    {
      std::vector<char> buff(1000);

      snprintf(buff.data(), buff.size(),
             "Fit-Result:\n"
             "  Ro=%0.4f ± %0.4f \n"
             "  Rs=%0.4f ± %0.4f\n"
             "  Rl=%0.4f ± %0.4f\n"
             "  lam=%0.4f ± %0.4f (%g, %g)\n"
             "  alpha_out=%0.4f ± %0.4f (%g, %g)\n"
             "  alpha_side=%0.4f ± %0.4f (%g, %g)\n"
             "  alpha_long=%0.4f ± %0.4f (%g, %g)\n"
             "  norm=%0.4f ± %0.4f\n"
             " -------------\n", Ro.value, Ro.error,
                                 Rs.value, Rs.error,
                                 Rl.value, Rl.error,
                                 lam.value, lam.error, lam.value - lam.error, lam.value + lam.error,
                                 alphao.value, alphao.error, alphao.value - alphao.error, alphao.value + alphao.error,
                                 alphas.value, alphas.error, alphas.value - alphas.error, alphas.value + alphas.error,
                                 alphal.value, alphal.error, alphal.value - alphal.error, alphal.value + alphal.error,
                                 norm.value, norm.error);
      return buff.data();
    }

    std::map<std::string, double>
    as_map() const
    {
      #define OUT(__name) {#__name, __name.value}, { # __name "_err", __name.error}

      return {
        OUT(Ro), OUT(Rs), OUT(Rl), OUT(lam), OUT(alphao), OUT(alphas), OUT(alphal), OUT(norm)
      };

      #undef OUT
    }
  };

  /// \brief 3D Levy fit parameters
  ///
  ///
  struct FitParams {
    double norm, lam;
    double alphao, alphas, alphal;
    double Ro, Rs, Rl;
    double gamma {1.0};

    FitParams(const double *par)
      : norm(par[NORM_PARAM_IDX])
      , lam(par[LAM_PARAM_IDX])
      , alphao(par[ALPHAOUT_PARAM_IDX])
      , alphas(par[ALPHASIDE_PARAM_IDX])
      , alphal(par[ALPHALONG_PARAM_IDX])
      , Ro(par[ROUT_PARAM_IDX])
      , Rs(par[RSIDE_PARAM_IDX])
      , Rl(par[RLONG_PARAM_IDX])
    {
    }

    FitParams(const FitResult &res)
      : norm(res.norm)
      , lam(res.lam)
      , alphao(res.alphao)
      , alphas(res.alphas)
      , alphal(res.alphal)
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
          || alphao < 0
          || alphas < 0
          || alphal < 0
          || norm < 0
          || std::isnan(Ro)
          || std::isnan(Rs)
          || std::isnan(Rl)
          || std::isnan(lam)
          || std::isnan(alphao)
          || std::isnan(alphas)
          || std::isnan(alphal)
          || std::isnan(norm);
    }

    /// Return calculated Rinv: $\sqrt{Ro^2 \gamma + Rs^2 + Rl^2}$
    double PseudoRinv(double gamma) const
      { return std::sqrt((Ro * Ro * gamma + Rs * Rs + Rl * Rl) / 3.0); }

    double gauss(const std::array<double, 3> &q, double K) const
      {
        std::array<double, 3> Rsq = {Ro*Ro, Rs*Rs, Rl*Rl};
        return FitterLevyFull::gauss(q, Rsq, lam, alphao, alphas, alphal, K, norm);
      }
  };

  /// Construct fitter from numerator denominator qinv histograms
  /// and a fit-range limit
  ///
  FitterLevyFull(TH3 &n, TH3 &d, TH3 &q, double limit=0.0)
    : Fitter3D(n, d, q, limit)
  {
  }

  FitterLevyFull(const Data3D &data)
    : Fitter3D(data)
    {
    }

  FitterLevyFull(Data3D &&data)
    : Fitter3D(std::move(data))
    {
    }

  static double
  gauss(const std::array<double, 3> &q, const FitParams &p, double K)
    { return p.gauss(q, K); }

  int
  setup_minuit(TMinuit &minuit) const override
  {
    int errflag = 0;
    minuit.mnparm(NORM_PARAM_IDX, "Norm", 0.25, 0.02, 0.0, 0.0, errflag);
    minuit.mnparm(LAM_PARAM_IDX, "Lam", 0.2, 0.1, 0.0, 1.0, errflag);
    minuit.mnparm(ROUT_PARAM_IDX, "Ro", 2.0, 1.0, 0.0, 0.0, errflag);
    minuit.mnparm(RSIDE_PARAM_IDX, "Rs", 2.0, 1.0, 0.0, 0.0, errflag);
    minuit.mnparm(RLONG_PARAM_IDX, "Rl", 2.0, 1.0, 0.0, 0.0, errflag);
    minuit.mnparm(ALPHAOUT_PARAM_IDX, "AlphaOut", 1.0, 0.1, 0.0, 0.0, errflag);
    minuit.mnparm(ALPHASIDE_PARAM_IDX, "AlphaSide", 1.0, 0.1, 0.0, 0.0, errflag);
    minuit.mnparm(ALPHALONG_PARAM_IDX, "AlphaLong", 1.0, 0.1, 0.0, 0.0, errflag);

    const double this_dbl = static_cast<double>((intptr_t)this);
    minuit.mnparm(DATA_PARAM_IDX, "DATA_PTR", this_dbl, 0, 0, INTPTR_MAX, errflag);

    minuit.FixParameter(DATA_PARAM_IDX);
    if (errflag != 0) {
      std::cerr << "Error setting paramters: " << errflag << "\n";
      throw std::runtime_error("Could not set Minuit parameters.");
    }

    return errflag;
  }


  FitResult fit_chi2()
    { return Fitter3D::fit_chi2(); }

  FitResult fit_pml()
    { return Fitter3D::fit_pml(); }

  FitResult fit()
    { return Fitter3D::fit(); }

};
