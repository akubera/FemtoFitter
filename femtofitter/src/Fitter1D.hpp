///
/// \file Fitter1D.hpp
///

#pragma once

#ifndef FITTER1D_HPP
#define FITTER1D_HPP

#include "CalculatorResid.hpp"
#include "CalculatorFsi.hpp"

#include "./Data1D.hpp"
#include "mrc/Mrc.hpp"

#include <typeinfo>
#include <TMinuit.h>
#include <TFile.h>
#include <iostream>

#include "CoulombHist.hpp"

/// \class Fitter1D
/// \brief Generic 1D
///
template <typename Impl>
class Fitter1D {
public:
  using CalcLoglike = ResidCalculatorPML<Impl>;
  using CalcChi2 = ResidCalculatorChi2<Impl>;

  /// The associated fit data
  Data1D data;

  /// The final-state-interaction calculator
  std::shared_ptr<FsiCalculator> fsi = nullptr;

  std::shared_ptr<Mrc1D> mrc = nullptr;

  Fitter1D(const TH1 &num, const TH1 &den, double limit)
    : data(num, den, limit)
    {
    }

  Fitter1D(const Data1D &dat)
    : data(dat)
    {
    }

  Fitter1D(TDirectory &tdir, double limit)
    : data(tdir, limit)
    {
    }

  virtual ~Fitter1D() = default;

  template <typename ResidFunc, typename FitParams>
  double resid_calc(const FitParams &p, ResidFunc resid_calc) const
    {
      double retval = 0;

      const std::function<double(double)>
        Kfsi = fsi
             ? fsi->ForRadius(p.radius)
             : [] (double qinv) { return 1.0; };

      for (const auto &datum : data) {
        const double
          n = datum.num,
          d = datum.den,
          q = datum.qinv,

          CF = p.evaluate(q, Kfsi(q));

        retval += resid_calc(n, d, CF);
      }

      return retval;
    }

  template <typename ResidFunc, typename FitParams>
  double resid_calc_mrc(const FitParams &p, Mrc1D &mrc, ResidFunc resid_calc, UInt_t npoints=3) const
    {
      double retval = 0;

      auto *cfhist = static_cast<TH1*>(data.src->num->Clone());
      cfhist->Reset();
      p.fill(*cfhist, *fsi, npoints);
      mrc.Smear(*cfhist);

      for (int i=0; i<cfhist->GetNbinsX(); ++i) {
        if (!data.mask->GetBinContent(i)) {
          continue;
        }

        const auto &datum = data[i];

        const double
          n = datum.num,
          d = datum.den,

          CF = cfhist->GetBinContent(i+1);

        retval += resid_calc(n, d, CF);
      }

      return retval;
    }

  template <typename ResidFunc, typename FitParams>
  double resid_calc_mrc(const FitParams &p, ResidFunc resid_calc) const
    {
      return resid_calc_mrc(p, mrc, resid_calc);
    }

  virtual int setup_minuit(TMinuit &minuit) const = 0;

  void setup_chi2_fitter(TMinuit &minuit)
    {
      static_cast<Impl*>(this)->setup_minuit(minuit);
      minuit.SetFCN(minuit_f<typename Impl::CalcChi2>);
    }

  void setup_pml_fitter(TMinuit &minuit)
    {
      static_cast<Impl*>(this)->setup_minuit(minuit);
      minuit.SetFCN(minuit_func_mrc<typename Impl::CalcLoglike>);
    }

  std::size_t size() const
    { return data.size(); }

  size_t degrees_of_freedom() const
    { return data.size() - Impl::CountParams(); }

  auto num_as_vec() const -> std::vector<double>
    { return numerator_as_vec(*this); }

  auto do_fit_minuit(TMinuit &minuit)
    {
      double strat_args[] = {1.0};
      double migrad_args[] = {2000.0, 1.0};
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

  auto fit_chi2()
    {
      TMinuit minuit;
      minuit.SetPrintLevel(-1);
      setup_chi2_fitter(minuit);
      return do_fit_minuit(minuit);
    }

  auto fit_pml()
    {
      TMinuit minuit;
      minuit.SetPrintLevel(-1);
      setup_pml_fitter(minuit);
      return do_fit_minuit(minuit);
    }

  auto fit_pml_mrc()
    {
      if (mrc == nullptr) {
        throw std::runtime_error("Fitter missing Mrc1D object");
      }

      if (fsi == nullptr) {
        throw std::runtime_error("Fitter missing Fsi object");
      }

      TMinuit minuit;
      minuit.SetPrintLevel(-1);
      setup_pml_fitter(minuit);

      // fit without smearing, first
      do_fit_minuit(minuit);

      minuit.SetFCN(minuit_func_mrc<typename Impl::CalcLoglike>);
      return do_fit_minuit(minuit);
    }

  template <typename FitResult>
  double
  resid_chi2(const FitResult &r) const
    {
      const typename Impl::FitParams &params = r;
      return resid_chi2_calc(params);
    }

  template <typename FitParams>
  double
  resid_chi2_calc(const FitParams &p) const
    {
      return resid_calc<ResidCalculatorChi2, FitParams>(p);
    }
};


template <typename CRTP>
struct FitParam1D {

  /// Fill histogram with no FSI factor
  void fill(TH1 &h) const
    {
      const TAxis &xaxis = *h.GetXaxis();
      for (int i=1; i <= xaxis.GetLast(); ++i) {
        double q = xaxis.GetBinCenter(i);
        double k = 1.0;
        double cf = static_cast<const CRTP*>(this)->evaluate(q, k);
        h.SetBinContent(i, cf);
      }
    }

  /// Fill histogram with average of N-points per bin
  ///
  void fill(TH1 &h, FsiCalculator &fsi, UInt_t npoints=1) const
    {
      auto &self = static_cast<const CRTP&>(*this);

      const TAxis &xaxis = *h.GetXaxis();
      auto Kfsi = fsi.ForRadius(self.Rinv());

      for (int i=1; i <= xaxis.GetLast(); ++i) {
        const double
          qlo = xaxis.GetBinLowEdge(i),
          qhi = xaxis.GetBinUpEdge(i),
          qstep = (qhi - qlo) / npoints,
          qstart = qlo + qstep / 2;

        double sum = 0.0;
        for (double q=qstart; q < qhi; q += qstep) {
          double k = Kfsi(q);
          sum += self.evaluate(q, k);
        }

        const double cf = sum / npoints;
        h.SetBinContent(i, cf);
      }
    }
};

#endif
