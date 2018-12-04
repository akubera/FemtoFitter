///
/// \file Fitter3D.hpp
///

#pragma once

#ifndef FITTER3D_HPP
#define FITTER3D_HPP

#include <typeinfo>

template <typename Impl>
class Fitter3D {
public:
  using CalcLoglike = ResidCalculatorPML<Impl>;
  using CalcChi2 = ResidCalculatorChi2<Impl>;


  /// The Associated fit data
  Data3D data;

  /// Approximate gamma factor of pairs
  double gamma;

  Fitter3D(TH3 &n, TH3 &d, TH3 &q, double limit)
    : data(n, d, q, limit)
    , gamma(1.0)
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

  template <typename ResidFunc, typename FitParams>
  double resid_calc(const FitParams &p, ResidFunc resid_calc) const
  {
    double retval = 0;

    double phony_r = p.PseudoRinv();
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
    auto params = static_cast<const typename Impl::FitParams&>(r);
    return resid_calc(params, chi2_calc);
  }

  template <typename FitResult>
  double
  resid_pml(const FitResult &r) const
  {
    auto params = static_cast<const typename Impl::FitParams&>(r);
    return resid_calc(params, loglikelihood_calc);
  }


  /// Automatic Fit Function
  ///
  /// Create minuit function and call minuit_f with the
  /// ResidualCalculation template parameter
  ///
  template <typename ResidCalc_t>
  auto fit(double fit_factor) // -> Impl::FitResult
  {
    TMinuit minuit;
    minuit.SetPrintLevel(-1);
    static_cast<Impl*>(this)->setup_minuit(minuit);

    minuit.SetFCN(minuit_f<ResidCalc_t>);

    return do_fit_minuit(minuit, fit_factor);
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
      return do_fit_minuit(minuit, 0.5);
    }

  auto fit_chi2()
    {
      TMinuit minuit;
      minuit.SetPrintLevel(-1);
      setup_chi2_fitter(minuit);
      return do_fit_minuit(minuit, 1.0);
    }

  auto fit()
    { return fit_chi2(); }

  auto do_fit_minuit(TMinuit &minuit, double fit_factor) // -> FitResult
  {
    double strat_args[] = {1.0};
    double migrad_args[] = {2000.0, fit_factor};
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

  static
  void
  fit_func(Int_t &,
           Double_t *,
           Double_t &retval,
           Double_t *par,
           Int_t)
  {
    // returned when invalid input or output occurs
    static const double BAD_VALUE = 3e99;
    const auto &data = *(const Impl*)(intptr_t)(par[Impl::DATA_PARAM_IDX]);

    typename Impl::FitParams params(par);
    if (params.is_invalid()) {
      retval = BAD_VALUE;
      return;
    }

    retval = data.resid_chi2(params);
  }

  std::size_t size() const
    { return data.size(); }

  auto num_as_vec() const -> std::vector<double>
    { return numerator_as_vec(*this); }

};

#endif
