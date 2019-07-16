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
#include "mrc/Mrc.hpp"

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

  /// The Momentum Resolution Correction
  std::shared_ptr<Mrc3D> mrc = nullptr;

  /// temporary histogram used during fitting
  mutable std::unique_ptr<TH3D> _tmp_cf = nullptr;

  /// Used to initialize parameters
  std::unique_ptr<ParamHints> paramhints = nullptr;

  Fitter3D(const TH3 &n,
           const TH3 &d,
           const TH3 &q,
           double limit)
    : Fitter3D(n, d, q, limit, std::shared_ptr<FsiCalculator>())
    { }

  Fitter3D(const TH3 &n,
	         const TH3 &d,
	         const TH3 &q,
	         double limit,
	         std::shared_ptr<FsiCalculator> fsi_calc)
    : Fitter3D(n, d, q, limit, nullptr, fsi_calc)
    { }

  Fitter3D(const TH3 &n,
	         const TH3 &d,
	         const TH3 &q,
	         double limit,
	         std::shared_ptr<Mrc3D> mrc_ptr,
	         std::shared_ptr<FsiCalculator> fsi_calc=nullptr)
    : data(n, d, q, limit)
    , fsi(fsi_calc)
    , mrc(mrc_ptr)
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
      return From(tdir, {"num", "den", "qinv"}, limit);
    }

  static std::unique_ptr<Impl>
  From(TDirectory &tdir, std::array<const TString, 3> names, double limit=0.0)
    {
      TH3 *num = static_cast<TH3*>(tdir.Get(names[0])),
          *den = static_cast<TH3*>(tdir.Get(names[1])),
          *qinv = static_cast<TH3*>(tdir.Get(names[2]));

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

  template <typename ResidFunc, typename FitParams>
  double resid_calc(const FitParams &p, ResidFunc resid_calc) const
  {
    double retval = 0.0;

    double phony_r = p.PseudoRinv(data.gamma);
    auto Kfsi = fsi->ForRadius(phony_r);

    for (const auto &datum : data) {
      const double
        qo = datum.qo,
        qs = datum.qs,
        ql = datum.ql,
        n = datum.num,
        d = datum.den,
        q = datum.qinv,

        k = Kfsi(q),

        CF = p.gauss({qo, qs, ql}, k);

      retval += resid_calc(n, d, CF);
    }

    return retval;
  }

  template <typename ResidFunc, typename FitParams>
  double resid_calc_mrc(const FitParams &p, Mrc3D &mrc, ResidFunc resid_calc) const
    {
      double retval = 0.0;

      if (_tmp_cf == nullptr) {
        _tmp_cf.reset(static_cast<TH3D*>(data.src->num->Clone()));
      }

      auto &cfhist = *_tmp_cf;
      mrc.FillSmearedFit(cfhist, p, *data.src->qinv, *fsi);

      for (const auto &datum : data) {

        const double
          n = datum.num,
          d = datum.den,

          CF = cfhist.GetBinContent(datum.hist_bin);

        retval += resid_calc(n, d, CF);
      }

      return retval;
    }

  void set_use_chi2_func(TMinuit &minuit) const
    {
      minuit.SetFCN(minuit_func<typename Impl::CalcChi2>);
    }

  void set_use_chi2_mrc_func(TMinuit &minuit) const
    {
      minuit.SetFCN(minuit_func_mrc<typename Impl::CalcChi2>);
    }

  void set_use_pml_func(TMinuit &minuit) const
    {
      minuit.SetFCN(minuit_func<typename Impl::CalcLoglike>);
    }

  void set_use_pml_mrc_func(TMinuit &minuit) const
    {
      minuit.SetFCN(minuit_func_mrc<typename Impl::CalcLoglike>);
    }

  void setup_chi2_minuit(TMinuit &minuit)
    {
      static_cast<Impl*>(this)->setup_minuit(minuit);
      set_use_chi2_func(minuit);
    }

  void setup_chi2_mrc_minuit(TMinuit &minuit)
    {
      static_cast<Impl*>(this)->setup_minuit(minuit);
      set_use_pml_func(minuit);
    }

  void setup_pml_minuit(TMinuit &minuit)
    {
      static_cast<Impl*>(this)->setup_minuit(minuit);
      set_use_pml_func(minuit);
    }

  void setup_pml_mrc_minuit(TMinuit &minuit)
    {
      static_cast<Impl*>(this)->setup_minuit(minuit);
      set_use_pml_mrc_func(minuit);
    }

  auto fit_chi2()
    {
      TMinuit minuit;
      minuit.SetPrintLevel(-1);
      setup_chi2_minuit(minuit);
      return do_fit_minuit(minuit);
    }

  auto fit_chi2_mrc()
    {
      TMinuit minuit;
      minuit.SetPrintLevel(-1);
      setup_chi2_mrc_minuit(minuit);
      return do_fit_minuit(minuit);
    }

  auto fit_pml()
    {
      TMinuit minuit;
      minuit.SetPrintLevel(-1);
      setup_pml_minuit(minuit);
      return do_fit_minuit(minuit);
    }

  auto fit_pml_mrc()
    {
      TMinuit minuit;
      minuit.SetPrintLevel(-1);
      setup_pml_mrc_minuit(minuit);
      return do_fit_minuit(minuit);
    }

  auto fit()
    {
      return fit_pml_mrc();
    }

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
    return result;
  }

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

  template <typename FsiFunc>
  virtual void fill(TH3 &h,  FsiFunc &fsi) const = 0;

  virtual void fill(TH3 &h, const TH3 &qinv, FsiCalculator &fsi) const = 0;
  virtual void multiply(TH3 &h, const TH3 &qinv, FsiCalculator &fsi) const = 0;

  virtual void fill(TH3 &h, const TH3 &qinv, FsiCalculator &fsi, UInt_t npoints) const = 0;
  virtual void multiply(TH3 &h, const TH3 &qinv, FsiCalculator &fsi, UInt_t npoints) const = 0;

};


/// \class FitParam3D
/// \brief Template based superclass for fit parameters
///
template <typename CRTP>
struct FitParam3D : Fit3DParameters {

  virtual ~FitParam3D()
    { }

  // void evalute(double q, double, double k) const
  double evaluate(std::array<double, 3> q, double k) const
    {
      return static_cast<const CRTP*>(this)->gauss(q, k);
    }

  double Rinv() const { return 1.0; }

  void fill(TH3 &h, const TH3 &qinv, FsiCalculator &fsi) const override
    {
      _loop_over_bins(h, qinv, fsi, [&](int i, int j, int k, double cf)
        {
          h.SetBinContent(i, j, k, cf);
        });
    }

  void fill(TH3 &h, const TH3 &qinv, FsiCalculator &fsi) const override
    {
      _loop_over_bins(h, qinv, fsi, [&](int i, int j, int k, double cf)
        {
          h.SetBinContent(i, j, k, cf);
        });
    }

  void multiply(TH3 &h, const TH3 &qinv, FsiCalculator &fsi) const override
    {
      _loop_over_bins(h, qinv, fsi, [&](int i, int j, int k, double cf)
        {
          h.SetBinContent(i, j, k, cf * h.GetBinContent(i, j, k));
        });
    }

  void fill(TH3 &h, const TH3 &qinv, FsiCalculator &fsi, UInt_t npoints) const override
    {
      _loop_over_bins(h, qinv, fsi, npoints, [&](int i, int j, int k, double cf)
        {
          h.SetBinContent(i, j, k, cf);
        });
    }

  void multiply(TH3 &h, const TH3 &qinv, FsiCalculator &fsi, UInt_t npoints) const override
    {
      _loop_over_bins(h, qinv, fsi, npoints, [&](int i, int j, int k, double cf)
        {
          h.SetBinContent(i, j, k, cf * h.GetBinContent(i, j, k));
        });
    }

private:

  template <typename FsiFunc, typename FuncType>
  void _loop_over_bins(TH3 &h, FsiFunc &Kfsi, FuncType func) const
    {
      auto &self = static_cast<const CRTP&>(*this);

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

        const double
          k = Kfsi(qx, qy, qz),
          cf = self.evaluate({qx, qy, qz}, k);
        func(i, j, k, cf);
      } } }
    }

  template <typename FuncType>
  void _loop_over_bins(TH3 &h, const TH3 &qinvh, FsiCalculator &fsi, FuncType func) const
    {
      auto &self = static_cast<const CRTP&>(*this);

      auto Kfsi = fsi.ForRadius(self.Rinv());

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

        const double
          qinv = const_cast<TH3&>(qinvh).Interpolate(qx, qy, qz),
          k = Kfsi(qinv),
          cf = self.evaluate({qx, qy, qz}, k);
        func(i, j, k, cf);
      } } }
    }

  template <typename FuncType>
  void _loop_over_bins(TH3 &h, const TH3 &qinvh, FsiCalculator &fsi, UInt_t npoints, FuncType func) const
    {
      auto &self = static_cast<const CRTP&>(*this);

      auto Kfsi = fsi.ForRadius(self.Rinv());

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
          const double
            qinv = const_cast<TH3&>(qinvh).Interpolate(qx, qy, qz),
            k = Kfsi(qinv),
            cf = self.evaluate({qx, qy, qz}, k);

          sum += cf;
        } } }

        const double mean_cf = sum / (npoints * npoints * npoints);
        func(i, j, k, mean_cf);
      } } }
    }

};

#endif
