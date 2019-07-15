///
/// \file Fitter3D.hpp
///

#pragma once

#ifndef FITTER3D_HPP
#define FITTER3D_HPP

#include "CalculatorResid.hpp"
#include "CalculatorFsi.hpp"
#include "CoulombHist.hpp"
#include "Data3D.hpp"
#include "math/fit.hh"

#include "ParamHints.hpp"

#include <TMinuit.h>
#include <TFile.h>
#include <TPython.h>
#include <Python.h>

#include <iostream>
#include <typeinfo>
#include <functional>


/// \class Fitter3D
/// \brief Abstract fitter associated the the Data3D class
///
template <typename Impl>
class Fitter3D {
public:
  using CalcLoglike = ResidCalculatorPML<Impl>;
  using CalcChi2 = ResidCalculatorChi2<Impl>;

  struct FitResult;

  /// The Associated fit data
  Data3D data;

  /// The final-state-interaction calculator
  std::shared_ptr<FsiCalculator> fsi = nullptr;

  /// Used to initialize parameters
  std::unique_ptr<ParamHints> paramhints;

  Fitter3D(TH3 &n, TH3 &d, TH3 &q, double limit, std::shared_ptr<FsiCalculator> fsi_calc=nullptr)
    : data(n, d, q, limit)
    , fsi(fsi_calc)
    { }

  Fitter3D(std::unique_ptr<Data3D> data_)
    : data(std::move(data_))
    // , paramhints(std::make_unique<ParamHints>())
    { }

  Fitter3D(std::shared_ptr<const Data3D> data_)
    : data(*data_)
    { }

  Fitter3D(const Data3D &data_)
    : data(data_)
    { }

  Fitter3D(Data3D &&data_)
    : data(std::move(data_))
    { }

  void SetParamHintsFromDir(const TDirectory &tdir)
    {
      paramhints = std::make_unique<ParamHints>(tdir);
    }

  /// Utility function for building fitter with tdirectory in file
  /// at specified path
  static std::unique_ptr<Impl>
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
  static std::unique_ptr<Impl>
  From(TDirectory &tdir, double limit=0.0)
    {
      TH3 *num = static_cast<TH3*>(tdir.Get("num")),
          *den = static_cast<TH3*>(tdir.Get("den")),
          *qinv = static_cast<TH3*>(tdir.Get("qinv"));

      if (!num || !den || !qinv) {
        std::cerr << "Error loading " << typeid(Impl).name() << " histograms from path "
                  << "'" << tdir.GetName() << "'\n";
        return nullptr;
      }

      return std::make_unique<Impl>(*num, *den, *qinv, limit);
    }

  virtual ~Fitter3D() = default;

  size_t degrees_of_freedom() const
    { return data.size() - Impl::CountParams(); }

  /// Add parameters to minuit object
  virtual int setup_minuit(TMinuit &) const = 0;

  template <typename FitParams>
  double
  resid_chi2_calc(const FitParams &p) const
    {
      double result = 0;

      double pseudo_rinv = p.PseudoRinv(data.gamma);
      // auto coulomb_factor = CoulombHist::GetHistWithRadius(phony_r);
      // auto Kfsi = [&coulomb_factor] (double q) {
      //   return coulomb_factor.Interpolate(q);
      // };


      const std::function<double(double)>
         Kfsi = fsi
              ? fsi->ForRadius(pseudo_rinv)
              : [] (double qinv) { return 1.0; };

      // auto Kfsi_ptr = fsi->ForRadius(pseudo_rinv);
      // auto &Kfsi = *Kfsi_ptr;


      // for (size_t i=0; i < data.size(); ++i) {
      //   const auto &datum = data[i];
      for (const auto &datum : data) {
        // double CF = p.gauss(datum.qspace(), *(fsi)(datum.qinv));
        double CF = p.gauss(datum.qspace(), Kfsi(datum.qinv));
        result += datum.calc_chi2(CF);
      }

      return result;
    }

  template <typename ResidFunc, typename FitParams>
  double resid_calc(const FitParams &p, ResidFunc resid_calc) const
  {
    double retval = 0;

    double phony_r = p.PseudoRinv(data.gamma);
    // auto coulomb_factor = CoulombHist::GetHistWithRadius(phony_r);

    auto Kfsi = fsi->ForRadius(phony_r);
    // std::unique_ptr<FsiQinv> kfsi_ptr = fsi->ForRadius(phony_r);
    // auto &Kfsi = *kfsi_ptr;
    // auto Kfsi = [&coulomb_factor] (double q) {
    //   return coulomb_factor.Interpolate(q);
    // };

    for (const auto &datum : data) {
      const double
        qo = datum.qo,
        qs = datum.qs,
        ql = datum.ql,
        n = datum.num,
        d = datum.den,
        q = datum.qinv,

        CF = p.gauss({qo, qs, ql}, Kfsi(q));

      retval += resid_calc(n, d, CF);
    }

    return retval;
  }

  template <typename FitResult>
  double
  resid_chi2(const FitResult &r) const
  {
    const typename Impl::FitParams &params = r;
    return resid_chi2_calc(params);
  }

  template <typename FitResult>
  double
  resid_pml(const FitResult &r) const
  {
    // auto params = static_cast<const typename Impl::FitParams&>(r);
    const typename Impl::FitParams &params = r;
    return resid_calc(params, loglikelihood_calc);
  }

  /// Automatic Fit Function
  ///
  /// Create minuit function and call minuit_func with the
  /// ResidualCalculation template parameter
  ///
  template <typename ResidCalc_t>
  auto
  fit(double sigma=1.0)  // -> Impl::FitResult
  {
    TMinuit minuit;
    minuit.SetPrintLevel(-1);
    static_cast<Impl*>(this)->setup_minuit(minuit);

    minuit.SetFCN(minuit_func<ResidCalc_t>);

    return do_fit_minuit(minuit, sigma);
  }

  void setup_pml_fitter(TMinuit &minuit)
    {
      static_cast<Impl*>(this)->setup_minuit(minuit);
      minuit.SetFCN(minuit_func<typename Impl::CalcLoglike>);
    }

  void setup_chi2_fitter(TMinuit &minuit)
    {
      static_cast<Impl*>(this)->setup_minuit(minuit);
      minuit.SetFCN(minuit_func<typename Impl::CalcChi2>);
    }

  auto fit_pml()
    {
      TMinuit minuit;
      minuit.SetPrintLevel(-1);
      setup_pml_fitter(minuit);
      return do_fit_minuit(minuit, 1.0);
    }

  auto fit_chi2()
    {
      TMinuit minuit;
      minuit.SetPrintLevel(-1);
      setup_chi2_fitter(minuit);
      return do_fit_minuit(minuit);
    }

  auto fit()
    { return fit_chi2(); }

  auto
  do_fit_minuit(TMinuit &minuit, double sigma=1.0, int recursive_count=0)  // -> Impl::FitResult
  {
    double strat_args[] = {1.0};
    double migrad_args[] = {2000.0, sigma};
    double hesse_args[] = {2000.0, 1.0};

    int errflag;
    minuit.mnexcm("SET STRategy", strat_args, 1, errflag);
    minuit.mnexcm("MIGRAD", migrad_args, 2, errflag);

    strat_args[0] = 2.0;
    minuit.mnexcm("SET STRategy", strat_args, 1, errflag);
    minuit.mnexcm("MIGRAD", migrad_args, 2, errflag);

    minuit.mnexcm("HESSE", hesse_args, 1, errflag);

    auto result = typename Impl::FitResult(minuit);

    for (int pidx=1; pidx < Impl::CountParams(); ++pidx) {
      // Impl
      double val, _err;
      if (minuit.GetParameter(pidx, val, _err)) {
        if (val > 14 || val < 0 || (pidx > 2 && val < 1)) {
          // std::cout << "Bad Fit" << "\n";
          // minuit.mnprin(1, 0.0);

        // minuit.SetParameter(2, paramhints->GenLam(), .01);
        // minuit.SetParameter(3, , .1);
        // minuit.SetParameter(4, paramhints->GenRs(), .1);
        // minuit.SetParameter(5, paramhints->GenRl(), .1);
        /*
        minuit.mnparm(Impl::LAM_PARAM_IDX, "Lam", paramhints->GenLam(), 0.01, 0.0, 0.0, errflag);
        minuit.mnparm(Impl::ROUT_PARAM_IDX, "Ro", paramhints->GenRo(), 0.1, 0.0, 0.0, errflag);
        minuit.mnparm(Impl::RSIDE_PARAM_IDX, "Rs", paramhints->GenRs(), 0.1, 0.0, 0.0, errflag);
        minuit.mnparm(Impl::RLONG_PARAM_IDX, "Rl", paramhints->GenRl(), 0.1, 0.0, 0.0, errflag);
        */

        break;
        }
      }
    }

    return result;
  }

  // template <typename FitRes>
  // void fit_with_random_inits(TMinuit &, FitRes &res);

  std::size_t size() const
    { return data.size(); }

  auto num_as_vec() const -> std::vector<double>
    { return numerator_as_vec(*this); }

  void SetFsi(std::shared_ptr<FsiCalculator> ptr)
    { fsi = ptr; }

  /// Extract a number from a python object and store in dest
  ///
  static
  bool
  ExtractPythonNumber(PyObject *pyobj,
                      const char* key,
                      double &dest,
                      std::vector<std::string> &missing)
    {
      if (!PyMapping_HasKeyString(pyobj, key)) {
        missing.emplace_back(key);
        return false;
      }

      auto *item = PyMapping_GetItemString(pyobj, key);
      if (PyFloat_Check(item)) {
        dest = PyFloat_AS_DOUBLE(item);
        return true;
      }

      if (PyLong_Check(item)) {
        dest = PyLong_AsDouble(item);
        return true;
      }

      missing.emplace_back(key);
      return false;
    }
};

/// \class Fit3DParameters
/// \brief Abstract base class for 3D fitter parameters
struct Fit3DParameters {

  virtual ~Fit3DParameters()
    { }

  virtual void fill(TH3 &, FsiCalculator *fsi=nullptr, UInt_t npoints=1) const = 0;

  virtual void multiply(TH1 &, FsiCalculator *fsi=nullptr, UInt_t npoints=1) const = 0;
};


/// \class FitParam3D
/// \brief Template based superclass for fit parameters
///
template <typename CRTP>
struct FitParam3D : Fit3DParameters {

  virtual ~FitParam3D()
    { }

  void fill(TH3 &h, FsiCalculator *fsi=nullptr, UInt_t npoints=1) const override
    {
      auto &self = static_cast<const CRTP&>(*this);

      auto callback = [&](int i, int j, int k, double cf)
        {
          h.SetBinContent(i, j, k, cf);
        };


      auto loop = [&] (auto fsi_func)
        {
          if (npoints == 1) {
            _loop_over_bins(self, h, fsi_func, callback);
          } else {
            _loop_over_bins(self, h, fsi_func, npoints, callback);
          }
        };

      if (fsi == nullptr) {
        auto no_fsi = [] () { return 1.0; };
        loop(no_fsi);
        // if (npoints == 1) {
        //   _loop_over_bins(self, h, no_fsi, callback);
        // } else {
        //   _loop_over_bins(self, h, no_fsi, npoints, callback);
        // }
      } else {
        auto Kfsi = fsi->ForRadius(self.Rinv());
        loop(Kfsi);
        // if (npoints == 1) {
        //   _loop_over_bins(self, h, Kfsi, callback);
        // } else {
        //   _loop_over_bins(self, h, Kfsi, npoints, callback);
        // }
      }


    }

  void multiply(TH1 &, FsiCalculator *fsi=nullptr, UInt_t npoints=1) const override
    {

    }

private:

  template <typename FsiFuncType, typename FuncType>
  void _loop_over_bins(const CRTP &self, TH1 &h, FsiFuncType Kfsi, FuncType func) const
    {
      const TAxis
        &xaxis = *h.GetXaxis(),
        &yaxis = *h.GetYaxis(),
        &zaxis = *h.GetZaxis();

      for (int k=1; k <= zaxis.GetNbins(); ++k) {
        const double qz = zaxis.GetBinCenter(k);
      for (int j=1; j <= yaxis.GetNbins(); ++j) {
        const double qy = yaxis.GetBinCenter(j);
      for (int i=1; i <= xaxis.GetNbins(); ++i) {
        const double qx = xaxis.GetBinCenter(i);

        // const double k = Kfsi(q);
        const double k = 1.0;
        const double cf = self.evaluate(qx, qy, qz, k);
        func(i, j, k, cf);
      } } }
    }

  template <typename FsiFuncType, typename FuncType>
  void _loop_over_bins(const CRTP &self, TH1 &h, FsiFuncType Kfsi, UInt_t npoints, FuncType func) const
    {
      const TAxis
        &xaxis = *h.GetXaxis(),
        &yaxis = *h.GetYaxis(),
        &zaxis = *h.GetZaxis();

      for (int k=1; k <= zaxis.GetNbins(); ++k) {
        const double
          qzlo = zaxis.GetBinLowEdge(k),
          qzhi = zaxis.GetBinUpEdge(k),
          qzstep = (qzhi - qzlo) / npoints,
          qzstart = qzlo + qzstep / 2;

      for (int j=1; j <= yaxis.GetNbins(); ++j) {
        const double
          qylo = yaxis.GetBinLowEdge(j),
          qyhi = yaxis.GetBinUpEdge(j),
          qystep = (qyhi - qylo) / npoints,
          qystart = qylo + qystep / 2;

      for (int i=1; i <= xaxis.GetNbins(); ++i) {
        const double
          qxlo = xaxis.GetBinLowEdge(i),
          qxhi = xaxis.GetBinUpEdge(i),
          qxstep = (qxhi - qxlo) / npoints,
          qxstart = qxlo + qxstep / 2;

        double sum = 0.0;
        for (double qz=qzstart; qz < qzhi; qz += qzstep) {
        for (double qy=qystart; qy < qyhi; qy += qystep) {
        for (double qx=qxstart; qx < qxhi; qx += qxstep) {
          double k = Kfsi(std::sqrt(qx * qy * qz));
          sum += self.evaluate(qx, qy, qz, k);
        } } }

        const double mean_cf = sum / (npoints * npoints * npoints);
        func(i, j, k, mean_cf);
      } } }
    }

};

#endif
