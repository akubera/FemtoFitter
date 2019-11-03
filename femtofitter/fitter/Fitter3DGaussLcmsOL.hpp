///
/// \file fitter/Fitter3DGaussLcmsOL.hpp
///

#pragma once

#ifndef FITTER_FITTER3DGAUSSLCMSOL_HPP
#define FITTER_FITTER3DGAUSSLCMSOL_HPP

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

#include "Fitter3D.hpp"
#include "Data3D.hpp"

#include "fitter/fit-methods.hh"


/// \class Fitter3DGaussLcmsOL
/// \brief Fit out-side-long with gaussian parameters
///
struct Fitter3DGaussLcmsOL : public Fitter3D<Fitter3DGaussLcmsOL> {

  using Super = Fitter3D<Fitter3DGaussLcmsOL>;
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
    ROL_PARAM_IDX = 6,
  };

  static std::string GetName()
    { return "Fitter3DGaussLcmsOL"; }

  static unsigned char CountParams()
    { return 6; }

  static double
  gauss(const std::array<double, 3> &q,
        const std::array<double, 3> &RSq,
        double Rol,
        double lam,
        double K=1.0,
        double norm=1.0)
  {
    const double
      qo = q[0],
      qs = q[1],
      ql = q[2],
      Eo = qo * qo * RSq[0],
      Es = qs * qs * RSq[1],
      El = ql * ql * RSq[2],
      Eol = 2 * qo * ql * std::fabs(Rol) * Rol,
      gauss = 1.0 + std::exp(-(Eo + Es + El + Eol) / HBAR_C_SQ),
      result = (1.0 - lam) + lam * K * gauss;

    return norm * result;
  }

  /// \class FitResult
  /// \brief Values and stderr from minuit results
  ///
  struct FitResult : FitResult3D<FitResult, Fitter3DGaussLcmsOL> {
    Value lam,
          norm,
          Ro,
          Rs,
          Rl,
          Rol;

    FitResult()
      : lam({0, 0})
      , norm({0, 0})
      , Ro({0, 0})
      , Rs({0, 0})
      , Rl({0, 0})
      , Rol({0, 0})
      {}

    FitResult(TMinuit &minuit)
      : lam(minuit, LAM_PARAM_IDX)
      , norm(minuit, NORM_PARAM_IDX)
      , Ro(minuit, ROUT_PARAM_IDX)
      , Rs(minuit, RSIDE_PARAM_IDX)
      , Rl(minuit, RLONG_PARAM_IDX)
      , Rol(minuit, ROL_PARAM_IDX)
    {
    }

    void FillMinuit(TMinuit &minuit) const override
      {
        int errflag = 0;
        minuit.mnparm(NORM_PARAM_IDX, "Norm", norm.value, 0.02, 0.0, 0.0, errflag);
        minuit.mnparm(LAM_PARAM_IDX, "Lam", lam.value, 0.05, 0.0, 1.0, errflag);
        minuit.mnparm(ROUT_PARAM_IDX, "Ro", Ro.value, 0.5, 0.0, 0.0, errflag);
        minuit.mnparm(RSIDE_PARAM_IDX, "Rs", Rs.value, 0.5, 0.0, 0.0, errflag);
        minuit.mnparm(RLONG_PARAM_IDX, "Rl", Rl.value, 0.5, 0.0, 0.0, errflag);
        minuit.mnparm(ROL_PARAM_IDX, "Rol", Rol.value, 0.5, 0.0, 0.0, errflag);
      }

    void print() const
    {
      printf("Fit-Result:\n"
             "  Ro=%0.4f ± %0.4f \n"
             "  Rs=%0.4f ± %0.4f\n"
             "  Rl=%0.4f ± %0.4f\n"
             "  Rol=%0.4f ± %0.4f\n"
             "  lam=%0.4f ± %0.4f (%g, %g)\n"
             "  norm=%0.4f ± %0.4f\n"
             " -------------\n", Ro.value, Ro.error,
                                 Rs.value, Rs.error,
                                 Rl.value, Rl.error,
                                 Rol.value, Rol.error,
                                 lam.value, lam.error, lam.value - lam.error, lam.value + lam.error,
                                 norm.value, norm.error);
    }

    std::string
    __repr__() const
      {
        return Form("<Fitter3DGaussLcmsOL::FitResult Ro=%g Rs=%g Rl=%g Rol=%g lambda=%g norm=%g>",
                    Ro.value, Rs.value, Rl.value, Rol.value, lam.value, norm.value);
      }

    std::map<std::string, double>
    as_map() const
    {
      #define OUT(__name) {#__name, __name.value}, { # __name "_err", __name.error}

      return {
        OUT(Ro), OUT(Rs), OUT(Rl), OUT(Rol), OUT(lam), OUT(norm)
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
      Add(Rol);
      Add(lam);
      Add(norm);

      return dict;
      #undef Add
    }
  };

  /// \brief 3D Gaussian fit parameters
  ///
  ///
  struct FitParams : FitParam3D<FitParams> {
    double norm, lam;
    double Ro, Rs, Rl, Rol;

    double evaluate(const std::array<double, 3> &q, double K) const
      {
        return Fitter3DGaussLcmsOL::gauss(q, {Ro*Ro, Rs*Rs, Rl*Rl}, Rol, lam, K, norm);
      }

    void Normalize(TH3 &cf) const
      {
        cf.Scale(1.0 / norm);
      }

    FitParams(double *par)
      : norm(par[NORM_PARAM_IDX])
      , lam(par[LAM_PARAM_IDX])
      , Ro(par[ROUT_PARAM_IDX])
      , Rs(par[RSIDE_PARAM_IDX])
      , Rl(par[RLONG_PARAM_IDX])
      , Rol(par[ROL_PARAM_IDX])
    {
    }

    FitParams(const FitResult &res)
      : norm(res.norm)
      , lam(res.lam)
      , Ro(res.Ro)
      , Rs(res.Rs)
      , Rl(res.Rl)
      , Rol(res.Rol)
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
          || std::isnan(Rol)
          || std::isnan(lam)
          || std::isnan(norm);
    }

    /// Return calculated Rinv: $\sqrt{Ro^2 \gamma + Rs^2 + Rl^2}$
    double PseudoRinv(double gamma) const
      { return std::sqrt((Ro * Ro * gamma * gamma + Rs * Rs + Rl * Rl) / 3.0); }

    std::string
    __repr__() const
      {
        return Form("<Fitter3DGaussLcmsOL::FitParam Ro=%g Rs=%g Rl=%g Rol=%g lambda=%g norm=%g>",
                    Ro, Rs, Rl, Rol, lam, norm);
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
      Add(Rol);
      Add(lam);
      Add(norm);

      return dict;
      #undef Add
    }
  };

  /// Construct fitter from numerator denominator qinv histograms
  /// and a fit-range limit
  ///
  Fitter3DGaussLcmsOL(TH3 &n, TH3 &d, TH3 &q, double limit=0.0)
    : Fitter3D(n, d, q, limit)
  {
  }

  static std::unique_ptr<Fitter3DGaussLcmsOL>
  FromDirectory(TDirectory &dir, double limit=0.0)
    {
      auto data = Data3D::FromDirectory(dir, limit);
      if (!data) {
        return nullptr;
      }
      return std::make_unique<Fitter3DGaussLcmsOL>(std::move(data));
    }

  Fitter3DGaussLcmsOL(const Data3D &dat)
    : Fitter3D(dat)
  {
  }

  Fitter3DGaussLcmsOL(Data3D &&dat)
    : Fitter3D(std::move(dat))
  {
  }

  Fitter3DGaussLcmsOL(std::unique_ptr<Data3D> dat)
    : Fitter3D(std::move(dat))
  {
  }

  static double
  gauss(const std::array<double, 3>& q, const FitParams &p, double K)
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
    minuit.mnparm(ROL_PARAM_IDX, "Rol", 0.0, 0.8, 0.0, 0.0, errflag);

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
            : (name == "Rol") ? ROL_PARAM_IDX
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

  DECLARE_FIT_METHODS(Fitter3D);
  DECLARE_RESID_METHODS(Fitter3D);

};

#endif
