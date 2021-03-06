///
/// \file Fitter1D.hpp
///

#pragma once

#ifndef FITTER1D_HPP
#define FITTER1D_HPP

#include "CalculatorResid.hpp"
#include "CalculatorFsi.hpp"

#include "Data1D.hpp"
#include "mrc/Mrc.hpp"
#include "Value.hpp"
#include "math/constants.hh"

#include "fit-methods.hh"
#include "PythonInterface.hh"

#include <TMinuit.h>
#include <TFile.h>

#include <iostream>
#include <typeinfo>


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
  double resid_calc(const FitParams &p, ResidFunc resid_func) const
    {
      double retval = 0;

      auto Kfsi = fsi->ForRadius(p.radius);

      for (const auto &datum : data) {
        const double
          n = datum.num,
          d = datum.den,
          q = datum.qinv,
          k = Kfsi(q),

          CF = p.evaluate(q, k);

        retval += resid_func(n, d, CF);
      }

      return retval;
    }

  template <typename ResidFunc, typename FitParams>
  double resid_calc_mrc(const FitParams &p,
                        Mrc1D &mrc1d,
                        ResidFunc resid_func,
                        UInt_t npoints=1) const
    {
      double retval = 0;

      if (_tmp_cf == nullptr) {
        _tmp_cf.reset(static_cast<TH1D*>(data.src->num->Clone("cf_buffer")));
        _tmp_cf->SetDirectory(nullptr);
      }

      auto &cfhist = *_tmp_cf;
      mrc1d.FillSmearedFit(cfhist, p, *fsi);

      for (const auto &datum : data) {

        const double
          n = datum.num,
          d = datum.den,

          CF = cfhist.GetBinContent(datum.hist_bin);

        retval += resid_func(n, d, CF);
      }

      return retval;
    }

  template <typename ResidFunc, typename FitParams>
  double resid_calc_mrc(const FitParams &p, ResidFunc resid_func) const
    {
      return resid_calc_mrc(p, mrc, resid_func);
    }

  virtual int setup_minuit(TMinuit &minuit) const = 0;

  void setup_chi2_fitter(TMinuit &minuit) const
    {
      static_cast<const Impl*>(this)->setup_minuit(minuit);
      set_chi2_func(minuit);
    }

  void setup_chi2_mrc_fitter(TMinuit &minuit) const
    {
      static_cast<const Impl*>(this)->setup_minuit(minuit);
      set_chi2_mrc_func(minuit);
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

  void set_chi2_mrc_func(TMinuit &minuit) const
    {
      minuit.SetFCN(minuit_func_mrc<typename Impl::CalcChi2>);
    }

  void set_pml_func(TMinuit &minuit) const
    {
      minuit.SetFCN(minuit_func<typename Impl::CalcLoglike>);
    }

  void set_pml_mrc_func(TMinuit &minuit) const
    {
      minuit.SetFCN(minuit_func_mrc<typename Impl::CalcLoglike>);
    }

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

  auto fit_chi2_mrc()
    {
      if (mrc == nullptr) {
        throw std::runtime_error("Fitter missing Mrc1D object");
      }

      if (fsi == nullptr) {
        throw std::runtime_error("Fitter missing Fsi object");
      }

      TMinuit minuit;
      minuit.SetPrintLevel(-1);

      setup_chi2_mrc_fitter(minuit);
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

  template <typename Params>
  void fill_smeared_fit(TH1 &h, const Params &p)
    {
      mrc->FillSmearedFit(h, p, *fsi, 1);
    }

  template <typename Params>
  double resid_calc_chi2_mrc(const Params &params) const
    {
      if (mrc == nullptr) {
        std::cerr << "MRC is null\n";
        return NAN;
      }

      return resid_calc_mrc(params, *mrc, Impl::CalcChi2::resid_func);
    }

  size_t size() const
    { return data.size(); }

  size_t degrees_of_freedom() const
    { return data.size() - Impl::CountParams(); }

  auto num_as_vec() const -> std::vector<double>
    { return numerator_as_vec(*this); }

  std::array<std::vector<double>, 3> as_array() const
    {
      std::array<std::vector<double>, 3> result;
      result[0].reserve(data.size());
      result[1].reserve(data.size());
      result[2].reserve(data.size());

      for (const auto &datum : data) {
        result[0].push_back(datum.qinv);
        result[1].push_back(datum.num);
        result[2].push_back(datum.den);
      }

      return result;
    }
};


/// \class Fit1DPrameters
/// \brief Abstract Base Class for 1d fit-parameters
///
struct Fit1DParameters {
  virtual ~Fit1DParameters(){}
  virtual void fill(TH1 &h) const = 0;
  virtual void fill(TH1 &, FsiCalculator &, UInt_t npoints=1) const = 0;
  virtual void multiply(TH1 &, FsiCalculator &, UInt_t npoints=1) const = 0;
};


/// \class FitParam1D
/// \brief Template based superclass for fit parameters
///
template <typename CRTP>
struct FitParam1D : Fit1DParameters {

  virtual ~FitParam1D() = default;

  /// Fill histogram with no FSI factor
  void fill(TH1 &h) const override
    {
      const TAxis &xaxis = *h.GetXaxis();
      for (int i=1; i <= xaxis.GetNbins(); ++i) {
        double q = xaxis.GetBinCenter(i);
        double k = 1.0;
        double cf = static_cast<const CRTP*>(this)->evaluate(q, k);
        h.SetBinContent(i, cf);
      }
    }

  /// Fill histogram with average of N-points per bin
  ///
  void fill(TH1 &h, FsiCalculator &fsi, UInt_t npoints=1) const override
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

  void multiply(TH1 &h, FsiCalculator &fsi, UInt_t npoints=1) const override
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

      for (int i=1; i <= xaxis.GetNbins(); ++i) {
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

template <typename CRTP, typename FitterType>
struct FitResult1D {
  using Paramters = typename FitterType::FitParams;

  virtual ~FitResult1D() = default;

  virtual void FillMinuit(TMinuit &) const = 0;

  Paramters as_params() const
    {
      return Paramters(static_cast<const CRTP&>(*this));
    }

  /// fill histogram with values of the correlation function
  /// represented by this fit result
  void fill_cf(TH1 &h, FsiCalculator &fsi, Mrc1D *mrc=nullptr) const
    {
      if (mrc == nullptr) {
        as_params().fill(h, fsi, 1);
      } else {
        mrc->FillSmearedFit(h, as_params(), fsi, 1);
      }
    }

  virtual void Evaluate(const double q, const double K) const
    {
      as_params().evaluate(q, K);
    }

  virtual void Normalize(TH1 &h) const
    {
      as_params().Normalize(h);
    }

  virtual PyObject* __iter__() const
    {
      auto *dict = static_cast<const CRTP*>(this)->as_dict();
      auto *list = PyDict_Items(dict);

      return PyObject_GetIter(list);
    }

};


#endif
