///
/// \file femtofitter/fitter/Fitter3DLevy.hpp
///

#pragma once

#include "Fitter3D.hpp"
#include "Value.hpp"
#include "math/constants.hh"

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


/// \class Fitter3DLevy
/// \brief Fit out-side-long with gaussian parameters
///
struct Fitter3DLevy : public Fitter3D<Fitter3DLevy> {

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
    ALPHA_PARAM_IDX = 6,
  };

  static unsigned char CountParams()
    { return 6; }

  static double
  levy(const std::array<double, 3> &q,
       const std::array<double, 3> &RSq,
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
  struct FitResult : FitResult3D<FitResult, Fitter3DLevy> {
    Value lam,
          norm,
          alpha,
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

    FitResult(const TMinuit &minuit)
      : lam(minuit, LAM_PARAM_IDX)
      , norm(minuit, NORM_PARAM_IDX)
      , alpha(minuit, ALPHA_PARAM_IDX)
      , Ro(minuit, ROUT_PARAM_IDX)
      , Rs(minuit, RSIDE_PARAM_IDX)
      , Rl(minuit, RLONG_PARAM_IDX)
    {
    }

    virtual ~FitResult() = default;

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
             "  alpha=%0.4f ± %0.4f (%g, %g)\n"
             "  norm=%0.4f ± %0.4f\n"
             " -------------\n", Ro.value, Ro.error,
                                 Rs.value, Rs.error,
                                 Rl.value, Rl.error,
                                 lam.value, lam.error, lam.value - lam.error, lam.value + lam.error,
                                 alpha.value, alpha.error, alpha.value - alpha.error, alpha.value + alpha.error,
                                 norm.value, norm.error);
      return buff.data();
    }

    std::string
    __repr__() const
      {
        return Form("<Fitter3DLevy::FitResult Ro=%g Rs=%g Rl=%g lambda=%g norm=%g>",
                    Ro.value, Rs.value, Rl.value, lam.value, norm.value);

      }

    std::map<std::string, double>
    as_map() const
    {
      #define OUT(__name) {#__name, __name.value}, { # __name "_err", __name.error}

      return {
        OUT(Ro), OUT(Rs), OUT(Rl), OUT(lam), OUT(alpha), OUT(norm)
      };

      #undef OUT
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
      Add(lam);
      Add(alpha);
      Add(norm);

      return dict;
      #undef Add
    }

    void FillMinuit(TMinuit &minuit) const
      {
        int errflag = 0;
        minuit.mnparm(NORM_PARAM_IDX, "Norm", norm.value, 0.02, 0.0, 0.0, errflag);
        minuit.mnparm(LAM_PARAM_IDX, "Lam", lam.value, 0.05, 0.0, 1.0, errflag);
        minuit.mnparm(ALPHA_PARAM_IDX, "alpha", alpha.value, 0.15, 0.0, 0.0, errflag);
        minuit.mnparm(ROUT_PARAM_IDX, "Ro", Ro.value, 0.5, 0.0, 0.0, errflag);
        minuit.mnparm(RSIDE_PARAM_IDX, "Rs", Rs.value, 0.5, 0.0, 0.0, errflag);
        minuit.mnparm(RLONG_PARAM_IDX, "Rl", Rl.value, 0.5, 0.0, 0.0, errflag);
      }

    FitParams as_params() const;
  };

  /// \brief 3D Levy fit parameters
  ///
  ///
  struct FitParams : FitParam3D<FitParams> {
    double norm, lam, alpha;
    double Ro, Rs, Rl;

    FitParams(double *par)
      : norm(par[NORM_PARAM_IDX])
      , lam(par[LAM_PARAM_IDX])
      , alpha(par[ALPHA_PARAM_IDX])
      , Ro(par[ROUT_PARAM_IDX])
      , Rs(par[RSIDE_PARAM_IDX])
      , Rl(par[RLONG_PARAM_IDX])
      { }

    FitParams(const FitResult &res)
      : norm(res.norm)
      , lam(res.lam)
      , alpha(res.alpha)
      , Ro(res.Ro)
      , Rs(res.Rs)
      , Rl(res.Rl)
      { }

    FitParams(PyObject *pyobj)
      {
        std::vector<std::string> missing_keys;

        if (!PyMapping_Check(pyobj)) {
          TPython::Exec(Form("raise TypeError('Object not a collection!')"));
          throw std::runtime_error("Object not a python collection");
        }

        ExtractPythonNumber(pyobj, "norm", norm, missing_keys);
        ExtractPythonNumber(pyobj, "lam", lam, missing_keys);
        ExtractPythonNumber(pyobj, "alpha", alpha, missing_keys);
        ExtractPythonNumber(pyobj, "Ro", Ro, missing_keys);
        ExtractPythonNumber(pyobj, "Rs", Rs, missing_keys);
        ExtractPythonNumber(pyobj, "Rl", Rl, missing_keys);

        if (!missing_keys.empty()) {
          std::string msg = "Python object missing required items:";
          for (const auto &key : missing_keys) {
            msg += " ";
            msg += key;
          }
          TPython::Exec(Form("raise ValueError('%s')", msg.c_str()));
          throw std::runtime_error(msg);
        }
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
    double PseudoRinv(double gamma) const
      { return std::sqrt((Ro * Ro * gamma * gamma + Rs * Rs + Rl * Rl) / 3.0); }

    double evaluate(const std::array<double, 3> &q, double K) const
      {
        std::array<double, 3> Rsq = {Ro*Ro, Rs*Rs, Rl*Rl};
        return Fitter3DLevy::levy(q, Rsq, lam, alpha, K, norm);
      }

    double evaluate(double qo, double qs, double ql, double K) const
      { return evaluate({qo, qs, ql}, K); }

    std::string
    __repr__() const
      {
        return Form("<Fitter3DLevy::FitParam Ro=%g Rs=%g Rl=%g lambda=%g alpha=%g norm=%g>",
                    Ro, Rs, Rl, lam, alpha, norm);
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
      Add(lam);
      Add(alpha);
      Add(norm);

      return dict;
      #undef Add
    }
  };

  /// Construct fitter from numerator denominator qinv histograms
  /// and a fit-range limit
  ///
  Fitter3DLevy(TH3 &n, TH3 &d, TH3 &q, double limit=0.0)
    : Fitter3D(n, d, q, limit)
  {
  }

  Fitter3DLevy(const Data3D &dat)
    : Fitter3D(dat)
  {
  }

  Fitter3DLevy(Data3D &&dat)
    : Fitter3D(std::move(dat))
  {
  }

  static double
  gauss(const std::array<double, 3> &q, const FitParams &p, double K)
    { return p.evaluate(q, K); }

  int
  setup_minuit(TMinuit &minuit) const override
  {
    int errflag = 0;
    minuit.mnparm(NORM_PARAM_IDX, "Norm", 0.25, 0.02, 0.0, 0.0, errflag);
    minuit.mnparm(LAM_PARAM_IDX, "Lam", 0.2, 0.1, 0.0, 1.0, errflag);
    minuit.mnparm(ROUT_PARAM_IDX, "Ro", 2.0, 1.0, 0.0, 0.0, errflag);
    minuit.mnparm(RSIDE_PARAM_IDX, "Rs", 2.0, 1.0, 0.0, 0.0, errflag);
    minuit.mnparm(RLONG_PARAM_IDX, "Rl", 2.0, 1.0, 0.0, 0.0, errflag);
    minuit.mnparm(ALPHA_PARAM_IDX, "Alpha", 1.8, 0.1, 1.0, 3.0, errflag);

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


inline auto
Fitter3DLevy::FitResult::as_params() const -> FitParams
{
  return FitParams(*this);
}
