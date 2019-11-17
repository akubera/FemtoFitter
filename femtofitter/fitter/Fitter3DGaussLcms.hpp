///
/// \file fitter/Fitter3DGaussLcms.hpp
///

#pragma once

#ifndef FITTER3DGAUSSLCMS_HPP
#define FITTER3DGAUSSLCMS_HPP


#include "CalculatorResid.hpp"
#include "Value.hpp"
#include "math/constants.hh"
#include "Fitter3D.hpp"

#include <TFile.h>
#include <TH3.h>
#include <TMinuit.h>
#include <TGraph.h>

#include <Python.h>

#include <array>
#include <vector>
#include <string>
#include <memory>
#include <valarray>
#include <iostream>


/// \class Fitter3DGaussLcms
/// \brief Fit out-side-long with gaussian parameters
///
struct Fitter3DGaussLcms : public Fitter3D<Fitter3DGaussLcms> {

  static std::string GetName()
    { return "Fitter3DGaussLcms"; }

  static constexpr std::uint8_t CountParams()
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

  using Super = Fitter3D<Fitter3DGaussLcms>;

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
  };

  /// \class FitResult
  /// \brief Values and stderr from minuit results
  ///
  struct FitResult : FitResult3D<FitResult, Fitter3DGaussLcms> {
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

    void FillMinuit(TMinuit &minuit) const override
      {
        int errflag = 0;
        minuit.mnparm(NORM_PARAM_IDX, "Norm", norm.value, 0.02, 0.0, 0.0, errflag);
        minuit.mnparm(LAM_PARAM_IDX, "Lam", lam.value, 0.05, 0.0, 1.0, errflag);
        minuit.mnparm(ROUT_PARAM_IDX, "Ro", Ro.value, 0.5, 0.0, 0.0, errflag);
        minuit.mnparm(RSIDE_PARAM_IDX, "Rs", Rs.value, 0.5, 0.0, 0.0, errflag);
        minuit.mnparm(RLONG_PARAM_IDX, "Rl", Rl.value, 0.5, 0.0, 0.0, errflag);
      }

    FitResult(PyObject *pyobj)
      {
        std::vector<std::string> missing_keys;

        if (!PyMapping_Check(pyobj)) {
          TPython::Exec(Form("raise TypeError('Object not a collection!')"));
          throw std::runtime_error("Object not a python collection");
        }

        ExtractPythonNumber(pyobj, "norm", norm.value, missing_keys);
        ExtractPythonNumber(pyobj, "norm_err", norm.error, missing_keys);
        ExtractPythonNumber(pyobj, "lam", lam.value, missing_keys);
        ExtractPythonNumber(pyobj, "lam_err", lam.error, missing_keys);
        ExtractPythonNumber(pyobj, "Ro", Ro.value, missing_keys);
        ExtractPythonNumber(pyobj, "Ro_err", Ro.error, missing_keys);
        ExtractPythonNumber(pyobj, "Rs", Rs.value, missing_keys);
        ExtractPythonNumber(pyobj, "Rs_err", Rs.error, missing_keys);
        ExtractPythonNumber(pyobj, "Rl", Rl.value, missing_keys);
        ExtractPythonNumber(pyobj, "Rl_err", Rl.error, missing_keys);

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

    std::string
    __repr__() const
      {
        return Form("<Fitter3DGaussLcms::FitResult Ro=%g Rs=%g Rl=%g lambda=%g norm=%g>",
                    Ro.value, Rs.value, Rl.value, lam.value, norm.value);
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
    double Ro, Rs, Rl;

    double evaluate(double qo, double qs, double ql, double K) const
      { return evaluate({qo, qs, ql}, K); }

    double evaluate(const std::array<double, 3> &q, double K) const
      { return Fitter3DGaussLcms::gauss(q, {Ro*Ro, Rs*Rs, Rl*Rl}, lam, K, norm); }

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

    FitParams(PyObject *pyobj)
      {
        std::vector<std::string> missing_keys;

        if (!PyMapping_Check(pyobj)) {
          TPython::Exec(Form("raise TypeError('Object not a collection!')"));
          throw std::runtime_error("Object not a python collection");
        }

        ExtractPythonNumber(pyobj, "norm", norm, missing_keys);
        ExtractPythonNumber(pyobj, "lam", lam, missing_keys);
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

    void Normalize(TH3 &hist) const
      {
        hist.Scale(1.0 / norm);
      }

    /// Return calculated Rinv: $\sqrt{Ro^2 \gamma + Rs^2 + Rl^2}$
    double PseudoRinv(double gamma) const
      { return std::sqrt((Ro * Ro * gamma * gamma + Rs * Rs + Rl * Rl) / 3.0); }

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

    std::string
    __repr__() const
      {
        return Form("<Fitter3DGaussLcms::FitParam Ro=%g Rs=%g Rl=%g lambda=%g norm=%g>",
                    Ro, Rs, Rl, lam, norm);
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
      Add(norm);

      return dict;
      #undef Add
    }
  };

  /// Construct fitter from numerator denominator qinv histograms
  /// and a fit-range limit
  ///
  Fitter3DGaussLcms(TH3 &n, TH3 &d, TH3 &q, double limit=0.0)
    : Fitter3D(n, d, q, limit)
  {
  }

  static std::unique_ptr<Fitter3DGaussLcms>
  FromDirectory(TDirectory &dir, double limit=0.0)
    {
      auto data = Data3D::FromDirectory(dir, limit);
      auto fitter = std::make_unique<Fitter3DGaussLcms>(std::move(data));
      fitter->paramhints = std::make_unique<ParamHints>(dir);
      return fitter;
    }

  Fitter3DGaussLcms(const Data3D &dat)
    : Fitter3D(dat)
    {
    }

  Fitter3DGaussLcms(Data3D &&dat)
    : Fitter3D(std::move(dat))
    {
    }

  Fitter3DGaussLcms(std::unique_ptr<Data3D> dat)
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

  FitResult fit_pml_mrc_quick()
    {
     if (mrc == nullptr) {
       throw std::runtime_error("Fitter missing Mrc1D object");
     }
     if (fsi == nullptr) {
       throw std::runtime_error("Fitter missing Fsi object");
     }
     TMinuit minuit;
     minuit.SetPrintLevel(-1);
     // first fit without smearing (faster)
     setup_pml_mrc_minuit(minuit);
     // auto tmp_res = do_fit_minuit(minuit);
     double strat_args[] = {1.0};
     double migrad_args[] = {2000.0, 1.0};
     int errflag;
     minuit.mnexcm("SET STRategy", strat_args, 1, errflag);
     minuit.mnexcm("MIGRAD", migrad_args, 2, errflag);

     FitResult tmp_res(minuit);
     // then fit with mrc-smearing (slower)
     TMinuit mminuit;
     mminuit.SetPrintLevel(-1);
     setup_pml_mrc_minuit(mminuit);
     tmp_res.FillMinuit(mminuit);
     return do_fit_minuit(mminuit);
   }

  void fit_with_random_inits(TMinuit &minuit, FitResult &res, int rec)
    {
      int errflag = 0;

      minuit.mnparm(NORM_PARAM_IDX, "NORM", 0.14, 1e-1, 0.0, 0.0, errflag);
      minuit.mnparm(LAM_PARAM_IDX, "Lam", paramhints->GenLam(), 1e-1, 0.0, 0.0, errflag);
      minuit.mnparm(ROUT_PARAM_IDX, "Ro", paramhints->GenRo(),  1e0, 0.0, 0.0, errflag);
      minuit.mnparm(RSIDE_PARAM_IDX, "Rs", paramhints->GenRs(), 1e0, 0.0, 0.0, errflag);
      minuit.mnparm(RLONG_PARAM_IDX, "Rl", paramhints->GenRl(), 1e0, 0.0, 0.0, errflag);

      res = do_fit_minuit(minuit, 1.0, rec);
    }

  void fill(TH3 &h, const FitParams &p) const
    {
      p.fill(h, *data.src->qinv, *fsi, data.gamma);
    }

  std::unique_ptr<TH3> get_cf(const FitResult &r) const
    {
      return get_cf(r.as_params());
    }

  std::unique_ptr<TH3> get_cf(const FitParams &p) const
    {
      std::unique_ptr<TH3> cf(static_cast<TH3*>(data.src->num->Clone()));
      cf->Reset();
      cf->SetTitle(Form("Correlation Function (R=<%0.2f,%0.2f,%0.2f,> \\lambda=%0.3f); q_{inv}; CF(q_{inv});", p.Ro, p.Rs, p.Rl, p.lam));
      cf->SetStats(false);
      fill(*cf, p);
      return cf;
    }

  std::unique_ptr<TH3> get_cf_mrc(const FitResult &r) const
    {
      return get_cf_mrc(r.as_params());
    }

  std::unique_ptr<TH3> get_cf_mrc(const FitParams &p) const
    {
      std::unique_ptr<TH3> cf(static_cast<TH3*>(data.src->num->Clone()));
      cf->Reset();
      cf->SetTitle(Form("Correlation Function (R=<%0.2f,%0.2f,%0.2f,> \\lambda=%0.3f); q_{inv}; CF(q_{inv});", p.Ro, p.Rs, p.Rl, p.lam));
      cf->SetStats(false);
      mrc->FillSmearedFit(*cf, p, *data.src->qinv, *fsi, 1);
      return cf;
    }

  DECLARE_FIT_METHODS(Fitter3D);
  DECLARE_RESID_METHODS(Fitter3D);
  // DECLARE_FILL_METHODS(TH3);

};

#endif
