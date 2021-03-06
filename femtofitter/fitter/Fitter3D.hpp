///
/// \file Fitter3D.hpp
///

#pragma once

#ifndef FITTER3D_HPP
#define FITTER3D_HPP

#include "CalculatorResid.hpp"
#include "CalculatorFsi.hpp"

#include "Data3D.hpp"
#include "math/constants.hh"
#include "mrc/Mrc.hpp"

#include "fit-methods.hh"

#include "ParamHints.hpp"
#include "PythonInterface.hh"

#include <TMinuit.h>
#include <TFile.h>

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

  mutable std::unique_ptr<TH3> _tmp_fsi = nullptr;

  /// Used to initialize parameters
  std::unique_ptr<ParamHints> paramhints = nullptr;

  Fitter3D(std::unique_ptr<TH3> n,
           std::unique_ptr<TH3> d,
           std::unique_ptr<TH3> q,
           double limit)
    : data(std::move(n), std::move(d), std::move(q), limit)
    { }

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

  /// Add parameters to minuit object
  virtual int setup_minuit(TMinuit &) const = 0;

  template <typename ResidFunc, typename FitParams>
  double resid_calc(const FitParams &p, ResidFunc resid_func) const
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

        CF = p.evaluate({qo, qs, ql}, k),
        res = resid_func(n, d, CF);

      retval += res;
    }

    return retval;
  }

  template <typename ResidFunc, typename FitParams>
  double resid_calc_mrc(const FitParams &p, Mrc3D &mrc3d, ResidFunc resid_func) const
    {
      double retval = 0.0;

      if (_tmp_cf == nullptr) {
        _tmp_cf = std::make_unique<TH3D>();
        data.src->num->Copy(*_tmp_cf);
      }
      auto &cf_hist = *_tmp_cf;

      if (_tmp_fsi == nullptr) {
        _tmp_fsi.reset(static_cast<TH3*>(data.src->qinv->Clone()));
      } else {
        data.src->qinv->Copy(*_tmp_fsi);
      }
      auto &fsi_hist = *_tmp_fsi;

      fsi->FillQinvHist(fsi_hist, p.Ro, p.Rs, p.Rl, data.gamma);
      mrc3d.FillSmearedFit(cf_hist, p, fsi_hist);

      for (const auto &datum : data) {

        const double
          n = datum.num,
          d = datum.den,

          CF = cf_hist.GetBinContent(datum.hist_bin),
          res = resid_func(n, d, CF);

        retval += res;
      }

      return retval;
    }

  template <typename FitParams>
  double resid_calc_chi2_mrc(const FitParams &p) const
    {
      return resid_calc_mrc(p, mrc, Impl::CalcChi2::resid_func);
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
      set_use_chi2_mrc_func(minuit);
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

  void SetFsi(std::shared_ptr<FsiCalculator> ptr)
    { fsi = ptr; }

  void SetMrc(std::shared_ptr<Mrc3D> ptr)
    { mrc = ptr; }

  size_t degrees_of_freedom() const
    { return data.size() - Impl::CountParams(); }

  std::size_t size() const
    { return data.size(); }

  auto num_as_vec() const -> std::vector<double>
    { return numerator_as_vec(*this); }

protected:

  Impl& self()
    {
      return static_cast<Impl&>(*this);
    }

  const Impl& self() const
    {
      return static_cast<const Impl&>(*this);
    }

  template <typename ParamType>
  void _fill(TH3 &h, const ParamType &p) const
    {
      TH3F fsi_hist;
      data.src->qinv->Copy(fsi_hist);
      fsi->FillQinvHist(fsi_hist, p.Ro, p.Rs, p.Rl, data.gamma);

      p.fill(h, fsi_hist);
    }

  template <typename ParamType>
  void _fill_smeared_fit(TH3 &h, const ParamType &p) const
    {
      TH3F fsi_hist;
      data.src->qinv->Copy(fsi_hist);
      fsi->FillQinvHist(fsi_hist, p.Ro, p.Rs, p.Rl, data.gamma);

      mrc->FillSmearedFit(h, p, fsi_hist);
    }

  template <typename ParamType>
  std::unique_ptr<TH3F> _get_cf_hist(const ParamType &p) const
    {
      std::unique_ptr<TH3F> cf(static_cast<TH3F*>(data.src->num->Clone()));
      cf->Reset();
      cf->SetStats(0);
      cf->Sumw2(false);

      cf->SetTitle(self().hist_title_from_params(p));

      if (mrc) {
        self().fill_smeared_fit(*cf, p);
      } else {
        self().fill(*cf, p);
      }

      return cf;
    }
};

/// \class Fit3DParameters
/// \brief Abstract base class for 3D fitter parameters
struct Fit3DParameters {

  using FsiFuncType = std::function<double(double,double,double)>;

  virtual ~Fit3DParameters()
    { }

  virtual void fill(TH3 &h, const TH3 &fsi) const = 0;
  virtual void multiply(TH3 &h, const TH3 &fsi) const = 0;

  // template <typename FsiFunc>
  // void fill(TH3 &h,  FsiFunc &fsi) const;
  virtual void fill(TH3&h, const FsiFuncType &fsi) const = 0;
  virtual void multiply(TH3&h, const FsiFuncType &fsi) const = 0;

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

  double Rinv() const { return 1.0; }

  void fill(TH3 &h, const TH3 &fsi) const override
    {
      _loop_over_bins(h, fsi, [&](int i, int j, int k, double cf)
        {
          h.SetBinContent(i, j, k, cf);
        });
    }

  void multiply(TH3 &h, const TH3 &fsi) const override
    {
      _loop_over_bins(h, fsi, [&](int i, int j, int k, double cf)
        {
          h.SetBinContent(i, j, k, cf * h.GetBinContent(i, j, k));
        });
    }

  void fill(TH3 &h, const FsiFuncType &fsi) const override
    {
      _loop_over_bins(h, fsi, [&](int i, int j, int k, double cf)
        {
          h.SetBinContent(i, j, k, cf);
        });
    }

  void multiply(TH3 &h, const FsiFuncType &fsi) const override
    {
      _loop_over_bins(h, fsi, [&](int i, int j, int k, double cf)
        {
          h.SetBinContent(i, j, k, cf * h.GetBinContent(i, j, k));
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

  std::pair<std::unique_ptr<TH3D>, std::unique_ptr<TH3D>>
  GetSmearedCorrFctnPair(const Mrc3D &mrc, const TH3 &qinv, FsiCalculator &fsi) const
    {
      auto &self = static_cast<const CRTP&>(*this);
      return mrc.GetSmearedCorrFctnPair(self, qinv, fsi);
    }

private:

  template <typename FuncType>
  void _loop_over_bins(TH3 &h, const TH3 &Kfsi, const FuncType &func) const
    {
      auto &self = static_cast<const CRTP&>(*this);

      const TAxis
        &xaxis = *h.GetXaxis(),
        &yaxis = *h.GetYaxis(),
        &zaxis = *h.GetZaxis();

      const int
        Nx = xaxis.GetNbins(),
        Ny = yaxis.GetNbins(),
        Nz = zaxis.GetNbins();

      for (int k=1; k <= Nz; ++k) {
        const double qz = zaxis.GetBinCenter(k);
      for (int j=1; j <= Ny; ++j) {
        const double qy = yaxis.GetBinCenter(j);
      for (int i=1; i <= Nx; ++i) {
        const double qx = xaxis.GetBinCenter(i);

        const double
          K = Kfsi.GetBinContent(i, j, k),
          cf = self.evaluate({qx, qy, qz}, K);

        func(i, j, k, cf);
      } } }
    }

  template <typename FsiFunc, typename FuncType>
  void _loop_over_bins(TH3 &h, const FsiFunc &Kfsi, FuncType func) const
    {
      auto &self = static_cast<const CRTP&>(*this);

      const TAxis
        &xaxis = *h.GetXaxis(),
        &yaxis = *h.GetYaxis(),
        &zaxis = *h.GetZaxis();

      const int
        Nx = xaxis.GetNbins(),
        Ny = yaxis.GetNbins(),
        Nz = zaxis.GetNbins();

      for (int k=1; k <= Nz; ++k) {
        const double qz = zaxis.GetBinCenter(k);
      for (int j=1; j <= Ny; ++j) {
        const double qy = yaxis.GetBinCenter(j);
      for (int i=1; i <= Nx; ++i) {
        const double qx = xaxis.GetBinCenter(i);

        const double
          K = Kfsi(qx, qy, qz),
          cf = self.evaluate({qx, qy, qz}, K);
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
          qinv = qinvh.GetBinContent(i, j, k),
          K = Kfsi(qinv),
          cf = self.evaluate({qx, qy, qz}, K);

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
            K = Kfsi(qinv),
            cf = self.evaluate({qx, qy, qz}, K);

          sum += cf;
        } } }

        const double mean_cf = sum / (npoints * npoints * npoints);
        func(i, j, k, mean_cf);
      } } }
    }

};

/// \class Fit3DParameters
/// \brief Abstract base class for 3D fitter parameters
template <typename CRTP, typename FitterType>
struct FitResult3D {
  using Paramters = typename FitterType::FitParams;

  virtual void FillMinuit(TMinuit &) const = 0;

  virtual ~FitResult3D()
    { }

  Paramters as_params() const
    {
      return Paramters(static_cast<const CRTP&>(*this));
    }

  /// fill histogram with values of the correlation function
  /// represented by this fit result
  void fill_cf(TH3 &h, FsiCalculator &fsi, Mrc3D *mrc=nullptr) const
    {
      if (mrc == nullptr) {
        as_params().fill(h, fsi, 1);
      } else {
        mrc->FillSmearedFit(h, as_params(), fsi, 1);
      }
    }

  double calc_chi2(FitterType &fitter) const
    {
      auto params = as_params();
      // ResidCalculatorChi2
      fitter.resid_calc(params, FitterType::CalcChi2::resid_func);
      return 0.0;
    }

  double evaluate(const std::array<double, 3> &q, double K) const
    {
      return as_params().evaluate(q, K);
    }

  void normalize(TH3 &cf) const
    {
      return as_params().Noramlize(cf);
    }
};

#endif
