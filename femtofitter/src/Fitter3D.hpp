///
/// \file Fitter3D.hpp
///

#pragma once

#ifndef FITTER3D_HPP
#define FITTER3D_HPP

#include "CalculatorResid.hpp"
#include "CoulombHist.hpp"
#include "Data3D.hpp"
#include "math/fit.hh"

#include <TMinuit.h>
#include <TFile.h>

#include <iostream>
#include <typeinfo>


template <typename Impl>
class Fitter3D {
public:
  using CalcLoglike = ResidCalculatorPML<Impl>;
  using CalcChi2 = ResidCalculatorChi2<Impl>;


  /// The Associated fit data
  Data3D data;

  Fitter3D(TH3 &n, TH3 &d, TH3 &q, double limit)
    : data(n, d, q, limit)
  { }

  Fitter3D(std::unique_ptr<Data3D> data_)
    : data(std::move(data_))
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

  /// Add parameters to minuit object
  virtual int setup_minuit(TMinuit &) const = 0;

  template <typename FitParams>
  double
  resid_chi2_calc(const FitParams &p) const
  {
    double result = 0;

    double phony_r = p.PseudoRinv(data.gamma);
    auto coulomb_factor = CoulombHist::GetHistWithRadius(phony_r);
    auto Kfsi = [&coulomb_factor] (double q) {
      return coulomb_factor.Interpolate(q);
    };

    for (size_t i=0; i < data.size(); ++i) {
      const auto &datum = data[i];
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
    auto coulomb_factor = CoulombHist::GetHistWithRadius(phony_r);

    auto Kfsi = [&coulomb_factor] (double q) {
      return coulomb_factor.Interpolate(q);
    };

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
  /// Create minuit function and call minuit_f with the
  /// ResidualCalculation template parameter
  ///
  template <typename ResidCalc_t>
  auto
  fit(double sigma=1.0)  // -> Impl::FitResult
  {
    TMinuit minuit;
    minuit.SetPrintLevel(-1);
    static_cast<Impl*>(this)->setup_minuit(minuit);

    minuit.SetFCN(minuit_f<ResidCalc_t>);

    return do_fit_minuit(minuit, sigma);
  }

  void setup_pml_fitter(TMinuit &minuit)
    {
      static_cast<Impl*>(this)->setup_minuit(minuit);
      minuit.SetFCN(minuit_f<typename Impl::CalcLoglike>);
    }

  void setup_chi2_fitter(TMinuit &minuit)
    {
      static_cast<Impl*>(this)->setup_minuit(minuit);
      minuit.SetFCN(minuit_f<typename Impl::CalcChi2>);
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
  do_fit_minuit(TMinuit &minuit, double sigma=1.0)  // -> Impl::FitResult
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

    return typename Impl::FitResult(minuit);
  }

  std::size_t size() const
    { return data.size(); }

  auto num_as_vec() const -> std::vector<double>
    { return numerator_as_vec(*this); }

};

#endif
