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

  /// The Momentum Resolution Correction
  std::shared_ptr<Mrc1D> mrc = nullptr;

  /// "cache" histogram used in fit
  mutable std::unique_ptr<TH1D> _tmp_cf = nullptr;

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
  double resid_calc_mrc(const FitParams &p, Mrc1D &mrc, ResidFunc resid_calc, UInt_t npoints=1) const
    {
      double retval = 0;

      if (_tmp_cf == nullptr) {
        _tmp_cf.reset(static_cast<TH1D*>(data.src->num->Clone()));
      }
      auto *cfhist = _tmp_cf.get();

      // mrc.FillSmearedFit(*cfhist, p, *fsi);
      p.FillAndSmearColMethod(*cfhist, *fsi, mrc);

      for (Int_t i=1; i<=cfhist->GetNbinsX(); ++i) {
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

  void setup_chi2_fitter(TMinuit &minuit) const
    {
      static_cast<const Impl*>(this)->setup_minuit(minuit);
      set_chi2_func(minuit);
    }

  void setup_pml_fitter(TMinuit &minuit) const
    {
      static_cast<const Impl*>(this)->setup_minuit(minuit);
      set_pml_func(minuit);
    }

  void setup_pml_mrc_fitter(TMinuit &minuit) const
    {
      static_cast<const Impl*>(this)->setup_minuit(minuit);
      set_pml_mrc_func(minuit);
    }

  void set_chi2_func(TMinuit &minuit) const
    {
      minuit.SetFCN(minuit_func<typename Impl::CalcChi2>);
    }

  void set_pml_func(TMinuit &minuit) const
    {
      minuit.SetFCN(minuit_func<typename Impl::CalcLoglike>);
    }

  void set_pml_mrc_func(TMinuit &minuit) const
    {
      minuit.SetFCN(minuit_func_mrc<typename Impl::CalcLoglike>);
    }

  size_t size() const
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
      if (fsi == nullptr) {
        throw std::runtime_error("Fitter missing Fsi object");
      }

      TMinuit minuit;
      minuit.SetPrintLevel(-1);
      setup_chi2_fitter(minuit);
      return do_fit_minuit(minuit);
    }

  auto fit_pml()
    {
      if (fsi == nullptr) {
        throw std::runtime_error("Fitter missing Fsi object");
      }

      TMinuit minuit;
      minuit.SetPrintLevel(-1);
      setup_pml_fitter(minuit);
      return do_fit_minuit(minuit);
    }

  /// Creates minuit and uses residual cut 'resid_calc_mrc'
  ///
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

      setup_pml_mrc_fitter(minuit);
      return do_fit_minuit(minuit);
    }

  /// First uses non-smeared fit to quickly get to close result,
  /// then calls pml_mrc to get final result.
  ///
  auto fit_pml_mrc_quick()
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
      setup_pml_fitter(minuit);
      // auto tmp_res = do_fit_minuit(minuit);
      double strat_args[] = {1.0};
      double migrad_args[] = {2000.0, 1.0};

      int errflag;
      minuit.mnexcm("SET STRategy", strat_args, 1, errflag);
      minuit.mnexcm("MIGRAD", migrad_args, 2, errflag);
      auto tmp_res = typename Impl::FitResult(minuit);

      // then fit with mrc-smearing (slower)
      TMinuit mminuit;
      mminuit.SetPrintLevel(-1);
      setup_pml_mrc_fitter(mminuit);
      tmp_res.FillMinuit(mminuit);
      return do_fit_minuit(mminuit);
    }
};


template <typename CRTP>
struct FitParam1D {

  void multiply(TH1 &h, std::shared_ptr<FsiCalculator> fsi, UInt_t npoints=1) const
    {
      multiply(h, *fsi, npoints);
    }

  void multiply(TH1 &h, FsiCalculator &fsi, UInt_t npoints=1) const override
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

      auto Kfsi = fsi.ForRadius(self.Rinv());

      _loop_over_bins(self, h, Kfsi, npoints, [&](int i, double cf) {
        h.SetBinContent(i, cf);
      });
    }

  void FillAndSmearRowMethod(TH1 &h, FsiCalculator &fsi, Mrc1D &mrc) const
    {
      fill(h, fsi, 1);
      mrc.SmearRowMethod(h);
    }

  void FillAndSmearColMethod(TH1 &h, FsiCalculator &fsi, Mrc1D &mrc) const
    {
      mrc.FillUnsmearedDen(h);
      multiply(h, fsi, 1);
      mrc.Smear(h);
      auto den = mrc.GetSmearedDenLike(h);
      h.Divide(den.get());
    }

  /// Multiply histogram contents with average of N-points per bin
  ///
  void multiply(TH1 &h, std::shared_ptr<FsiCalculator> fsi, UInt_t npoints=1) const
    {
      multiply(h, *fsi, npoints);
    }

  void multiply(TH1 &h, FsiCalculator &fsi, UInt_t npoints=1) const
    {
      auto &self = static_cast<const CRTP&>(*this);

      auto Kfsi = fsi.ForRadius(self.Rinv());

      _loop_over_bins(self, h, Kfsi, npoints, [&](int i, double cf) {
        h.SetBinContent(i, h.GetBinContent(i) * cf);
        h.SetBinError(i, h.GetBinError(i) * cf);
      });
    }

private:

  template <typename FsiFuncType, typename FuncType>
  void _loop_over_bins(const CRTP &self, TH1 &h, FsiFuncType Kfsi, UInt_t npoints, FuncType func) const
    {
      const TAxis &xaxis = *h.GetXaxis();

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

        const double mean_cf = sum / npoints;
        func(i, mean_cf);
      }

    }

};

template <typename CRTP>
struct FitResult1D {
  virtual void FillMinuit(TMinuit &) const = 0;
};


#endif
