///
/// \file fitter/Fitter3DGaussLcms.hpp
///

#pragma once

#ifndef FITTER3DGAUSSLCMS_HPP
#define FITTER3DGAUSSLCMS_HPP


#include "CoulombHist.hpp"
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

    void FillMinuit(TMinuit &minuit)
      {
        int errflag = 0;
        minuit.mnparm(NORM_PARAM_IDX, "Norm", norm.value, 0.02, 0.0, 0.0, errflag);
        minuit.mnparm(LAM_PARAM_IDX, "Lam", lam.value, 0.05, 0.0, 1.0, errflag);
        minuit.mnparm(ROUT_PARAM_IDX, "Ro", Ro.value, 0.5, 0.0, 0.0, errflag);
        minuit.mnparm(RSIDE_PARAM_IDX, "Rs", Rs.value, 0.5, 0.0, 0.0, errflag);
        minuit.mnparm(RLONG_PARAM_IDX, "Rl", Rl.value, 0.5, 0.0, 0.0, errflag);
      }

    FitParams as_params() const;
  };

  /// \brief 3D Gaussian fit parameters
  ///
  ///
  struct FitParams : FitParam3D<FitParams> {
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
      { return std::sqrt((Ro * Ro * gamma * gamma + Rs * Rs + Rl * Rl) / 3.0); }

    double evaluate(double qo, double qs, double ql, double K) const
      { return evaluate({qo, qs, ql}, K); }

    double evaluate(const std::array<double, 3> &q, double K) const
      { return Fitter3DGaussLcms::gauss(q, {Ro*Ro, Rs*Rs, Rl*Rl}, lam, K, norm); }

    double gauss(const std::array<double, 3> &q, double K) const
      { return evaluate(q, K); }

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
    apply_to(TH3 &hist, const TH3 &qinv, double gamma)
    {
      const int I = hist.GetNbinsX(),
                J = hist.GetNbinsY(),
                K = hist.GetNbinsZ();

      // const double phony_r = PseudoRinv(gamma);
      // auto coulomb_factor = CoulombHist::GetHistWithRadius(phony_r);

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
          Kq = 1.0;

        hist.SetBinContent(i,j,k, hist.GetBinContent(i,j,k) * gauss({qo, qs, ql}, Kq));
      }
    }

    void
    apply_to(TH3 &hist, const TH3 &qinv, FsiCalculator &fsi, double gamma)
    {
      const int I = hist.GetNbinsX(),
                J = hist.GetNbinsY(),
                K = hist.GetNbinsZ();

      const double Rinv = PseudoRinv(gamma);
      auto Kfsi = fsi.ForRadius(Rinv);
      // auto coulomb_factor = CoulombHist::GetHistWithRadius(phony_r);

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
          Kq = Kfsi(q);

        hist.SetBinContent(i,j,k, hist.GetBinContent(i,j,k) * gauss({qo, qs, ql}, Kq));
      }
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

  double residual_chi2(const FitResult &r) const
    {
      return residual_chi2(r.as_params());
    }

  double residual_chi2(const FitParams &p) const
    {
      return Fitter3D::resid_calc(p, ResidCalculatorChi2<Fitter3DGaussLcms>::resid_func);
    }

  double residual_pml(const FitParams &p) const
    {
      return Fitter3D::resid_calc(p, ResidCalculatorPML<Fitter3DGaussLcms>::resid_func);
    }

  double residual_pml(const FitResult &p) const
    {
      return residual_pml(p.as_params());
    }

  double residual_chi2_mrc(const FitResult &r) const
    {
      return residual_chi2_mrc(r.as_params());
    }

  double residual_chi2_mrc(const FitParams &p) const
    {
      return Fitter3D::resid_calc_mrc(p, *mrc, ResidCalculatorChi2<Fitter3DGaussLcms>::resid_func);
    }

  double residual_pml_mrc(const FitParams &p) const
    {
      return Fitter3D::resid_calc_mrc(p, *mrc, ResidCalculatorPML<Fitter3DGaussLcms>::resid_func);
      // return ResidCalculatorPML<Fitter3DGaussLcms>::resid_mrc(*this, p);
    }

  FitResult fit_chi2()
    { return Fitter3D::fit_chi2(); }

  FitResult fit_chi2_mrc()
    { return Fitter3D::fit_chi2_mrc(); }

  FitResult fit_pml()
    { return Fitter3D::fit_pml(); }

  FitResult fit_pml_mrc()
    { return Fitter3D::fit_pml_mrc(); }

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
};

inline auto
Fitter3DGaussLcms::FitResult::as_params() const -> FitParams
{
  return FitParams(*this);
}

#endif
