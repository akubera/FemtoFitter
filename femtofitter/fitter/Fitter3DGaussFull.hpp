///
/// \file Fitter3DGaussFull.hpp
///

#pragma once

#ifndef FITTER3DGAUSSFULL_HPP
#define FITTER3DGAUSSFULL_HPP

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
#include "Fitter3D.hpp"

#include "fit-methods.hh"


/// \class Fitter3DGaussFull
/// \brief Fit full $ q_i q_j R_ij^2 $ gaussian parameters
///
struct Fitter3DGaussFull : public Fitter3D<Fitter3DGaussFull> {

  using Super = Fitter3D<Fitter3DGaussFull>;
  struct FitParams;
  struct FitResult;

  /// constants used to lookup data from pointer
  enum {
    DATA_PARAM_IDX = 0,

    NORM_PARAM_IDX = 1,
    LAM_PARAM_IDX = 2,
    ROUT_PARAM_IDX = 3,
    RSIDE_PARAM_IDX = 4,
    RLONG_PARAM_IDX = 5,
    ROS_PARAM_IDX = 6,
    ROL_PARAM_IDX = 7,
    RSL_PARAM_IDX = 8,
  };

  static std::string GetName()
    { return "Fitter3DGaussFull"; }

  static unsigned char CountParams()
    { return 8; }

  static double
  gauss(std::array<double, 3> q,
        // std::array<double, 6> RSq,
        double Ro, double Rs, double Rl,
        double Ros, double Rol, double Rsl,
        double lam,
        double K=1.0,
        double norm=1.0)
  {
    const double
      qo = q[0],
      qs = q[1],
      ql = q[2],

      E = qo * qo * Ro * Ro
        + qs * qs * Rs * Rs
        + ql * ql * Rl * Rl
        + 2 * qo * qs * std::fabs(Ros) * Ros
        + 2 * qo * ql * std::fabs(Ros) * Rol
        + 2 * qs * ql * std::fabs(Ros) * Rsl,

      gauss = 1.0 + std::exp(-E / HBAR_C_SQ),
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
          Rl,
          Ros,
          Rol,
          Rsl;

    FitResult(TMinuit &minuit)
      : lam(minuit, LAM_PARAM_IDX)
      , norm(minuit, NORM_PARAM_IDX)
      , Ro(minuit, ROUT_PARAM_IDX)
      , Rs(minuit, RSIDE_PARAM_IDX)
      , Rl(minuit, RLONG_PARAM_IDX)
      , Ros(minuit, ROS_PARAM_IDX)
      , Rol(minuit, ROL_PARAM_IDX)
      , Rsl(minuit, RSL_PARAM_IDX)
    {
    }

    void print() const
    {
      printf("Fit-Result:\n"
             "  Ro=%0.4f ± %0.4f \n"
             "  Rs=%0.4f ± %0.4f\n"
             "  Rl=%0.4f ± %0.4f\n"
             "  Ros²=%0.4f ± %0.4f\n"
             "  Rol²=%0.4f ± %0.4f\n"
             "  Rsl²=%0.4f ± %0.4f\n"
             "  lam=%0.4f ± %0.4f (%g, %g)\n"
             "  norm=%0.4f ± %0.4f\n"
             " -------------\n", Ro.value, Ro.error,
                                 Rs.value, Rs.error,
                                 Rl.value, Rl.error,
                                 Ros.value, Ros.error,
                                 Rol.value, Rol.error,
                                 Rsl.value, Rsl.error,
                                 lam.value, lam.error, lam.value - lam.error, lam.value + lam.error,
                                 norm.value, norm.error);
    }

    std::map<std::string, double>
    as_map() const
    {
      #define OUT(__name) {#__name, __name.value}, { # __name "_err", __name.error}

      return {
        OUT(Ro), OUT(Rs), OUT(Rl),
        OUT(Ros), OUT(Rol), OUT(Rsl),
        OUT(lam), OUT(norm)
      };

      #undef OUT
    }

    std::string
    __repr__() const
      {
        return Form("<Fitter3DGaussFull::FitResult Ro=%g Rs=%g Rl=%g Ros=%g Rol=%g Rsl=%g lambda=%g norm=%g>",
                    Ro.value, Rs.value, Rl.value, Ros.value, Rol.value, Rsl.value, lam.value, norm.value);
      }

    PyObject*
    as_dict() const
    {
      #define Add(__name) \
        PyDict_SetItemString(dict, #__name, PyFloat_FromDouble(__name.value));\
        PyDict_SetItemString(dict, #__name "_err", PyFloat_FromDouble(__name.error))

      auto *dict = PyDict_New();
      Add(Ro);
      Add(Rs);
      Add(Rl);
      Add(Ros);
      Add(Rol);
      Add(Rsl);
      Add(lam);
      Add(norm);

      return dict;
      #undef Add
    }

    FitParams as_params() const;
  };

  /// \brief 3D Gaussian fit parameters
  ///
  ///
  struct FitParams : FitParam3D<FitParams> {
    double norm, lam;
    double Ro, Rs, Rl;
    double Ros, Rol, Rsl;
    double gamma {1.0};

    FitParams(double *par)
      : norm(par[NORM_PARAM_IDX])
      , lam(par[LAM_PARAM_IDX])
      , Ro(par[ROUT_PARAM_IDX])
      , Rs(par[RSIDE_PARAM_IDX])
      , Rl(par[RLONG_PARAM_IDX])
      , Ros(par[ROS_PARAM_IDX])
      , Rol(par[ROL_PARAM_IDX])
      , Rsl(par[RSL_PARAM_IDX])
    {
    }

    FitParams(const FitResult &res)
      : norm(res.norm)
      , lam(res.lam)
      , Ro(res.Ro)
      , Rs(res.Rs)
      , Rl(res.Rl)
      , Ros(res.Ros)
      , Rol(res.Rol)
      , Rsl(res.Rsl)
    {
    }

    bool is_invalid() const
    {
      auto invalid = [](double x) { return x < 0 || std::isnan(x); };

      return invalid(Ro)
          || invalid(Rs)
          || invalid(Rl)
          || std::isnan(Ros)
          || std::isnan(Rol)
          || std::isnan(Rsl)
          || invalid(lam)
          || invalid(norm);
    }

    /// Return calculated Rinv: $\sqrt{Ro^2 \gamma + Rs^2 + Rl^2}$
    double PseudoRinv(double gama) const
      { return std::sqrt((Ro * Ro * gama + Rs * Rs + Rl * Rl) / 3.0); }

    double Rinv() const
      { return PseudoRinv(gamma); }

    double evaluate(std::array<double, 3> q, double K) const
      { return gauss(q, K); }

    double gauss(std::array<double, 3> q, double K) const
      {
        // std::array<double, 3> Rsq = {Ro*Ro, Rs*Rs, Rl*Rl};
        // return Fitter3DGaussLcms::gauss(q, Rsq, lam, K, norm);
        return Fitter3DGaussFull::gauss(q, Ro, Rs, Rl, Ros, Rol, Rsl, lam, K, norm);
      }

    /// Multiply histogram with values from this correlation function
    void
    apply_to(TH3 &hist, TH3 &qinv, double gama)
      {
        const int I = hist.GetNbinsX(),
                  J = hist.GetNbinsY(),
                  K = hist.GetNbinsZ();

        const double phony_r = PseudoRinv(gama);
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
            Kq = coulomb_factor.Interpolate(q),

            CF = gauss({qo, qs, ql}, Kq),
            factor = hist.GetBinContent(i,j,k);

          hist.SetBinContent(i,j,k, CF * factor );
        }
      }

    std::string
    __repr__() const
      {
        return Form("<Fitter3DGaussFull::FitParam Ro=%g Rs=%g Rl=%g Ros=%g Rol=%g Rsl=%g lambda=%g norm=%g>",
                    Ro, Rs, Rl, Ros, Rol, Rsl, lam, norm);
      }

    PyObject*
    as_dict() const
    {
      #define Add(__name) \
        PyDict_SetItemString(dict, #__name, PyFloat_FromDouble(__name))

      auto *dict = PyDict_New();
      Add(Ro);
      Add(Rs);
      Add(Rl);
      Add(Ros);
      Add(Rol);
      Add(Rsl);
      Add(lam);
      Add(norm);

      return dict;
      #undef Add
    }
  };

  /// Construct fitter from numeratork, denominator, & qinv histograms
  ///
  /// limit specifies fit-range
  ///
  Fitter3DGaussFull(TH3 &n, TH3 &d, TH3 &q, double limit=0.0)
    : Fitter3D(n, d, q, limit)
    {
    }

  Fitter3DGaussFull(const Data3D &dat)
    : Fitter3D(dat)
    {
    }

  Fitter3DGaussFull(Data3D &&dat)
    : Fitter3D(std::move(dat))
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
    minuit.mnparm(LAM_PARAM_IDX, "Lam", 0.2, 0.01, 0.0, 0.0, errflag);
    minuit.mnparm(ROUT_PARAM_IDX, "Ro", 2.0, 1.0, 0.0, 0.0, errflag);
    minuit.mnparm(RSIDE_PARAM_IDX, "Rs", 2.0, 1.0, 0.0, 0.0, errflag);
    minuit.mnparm(RLONG_PARAM_IDX, "Rl", 2.0, 1.0, 0.0, 0.0, errflag);
    minuit.mnparm(ROS_PARAM_IDX, "Ros", 0.0, 0.05, 0.0, 0.0, errflag);
    minuit.mnparm(ROL_PARAM_IDX, "Rol", 0.0, 0.05, 0.0, 0.0, errflag);
    minuit.mnparm(RSL_PARAM_IDX, "Rsl", 0.0, 0.05, 0.0, 0.0, errflag);

    const double this_dbl = static_cast<double>((intptr_t)this);
    minuit.mnparm(DATA_PARAM_IDX, "DATA_PTR", this_dbl, 0, 0, INTPTR_MAX, errflag);

    minuit.FixParameter(DATA_PARAM_IDX);
    if (errflag != 0) {
      std::cerr << "Error setting paramters: " << errflag << "\n";
      throw std::runtime_error("Could not set Minuit parameters.");
    }

    return errflag;
  }

  DECLARE_FIT_METHODS(Fitter3D);
  DECLARE_RESID_METHODS(Fitter3D);

};

auto
Fitter3DGaussFull::FitResult::as_params() const -> FitParams
{
  return FitParams(*this);
}

#endif
