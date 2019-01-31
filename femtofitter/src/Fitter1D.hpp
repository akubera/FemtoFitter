///
/// \file Fitter1D.hpp
///

#pragma once

#ifndef FITTER1D_HPP
#define FITTER1D_HPP

#include "CalculatorResid.hpp"

#include "./Data1D.hpp"

#include <typeinfo>
#include <TMinuit.h>
#include <TFile.h>
#include <iostream>

#include "CoulombHist.hpp"


template <typename Impl>
class Fitter1D {
public:
  using CalcLoglike = ResidCalculatorPML<Impl>;
  using CalcChi2 = ResidCalculatorChi2<Impl>;

  /// The associated fit data
  Data1D data;

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

  template <typename ResidFunc, typename FitParams>
  double resid_calc(const FitParams &p, ResidFunc resid_calc) const
    {
      double retval = 0;

      // FitParams
      auto coulomb_factor = CoulombHist::GetHistWithRadius(p.radius);

      // auto Kfsi = [&coulomb_factor] (double q) {
      //   return coulomb_factor.Interpolate(q);
      // };
      auto Kfsi = coulomb_factor.Interpolate;

      for (const auto &datum : data) {
        const double
          n = datum.num,
          d = datum.den,
          q = datum.qinv,

          CF = p.gauss(q, Kfsi(q));

        retval += resid_calc(n, d, CF);
      }

      return retval;
    }

  void setup_chi2_fitter(TMinuit &minuit)
    {
      static_cast<Impl*>(this)->setup_minuit(minuit);
      minuit.SetFCN(minuit_f<typename Impl::CalcChi2>);
    }

  void setup_pml_fitter(TMinuit &minuit)
    {
      static_cast<Impl*>(this)->setup_minuit(minuit);
      minuit.SetFCN(minuit_f<typename Impl::CalcLoglike>);
    }

  std::size_t size() const
    { return data.size(); }

  auto num_as_vec() const -> std::vector<double>
    { return data.size(); }

};


#endif
