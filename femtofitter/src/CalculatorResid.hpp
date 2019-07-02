///
/// \file femtofitter/src/CalculatorResid.hpp
///

#pragma once

#ifndef CALCULATORRESID_HPP_
#define CALCULATORRESID_HPP_


#include "math/fit.hh"


/// \brief static minuit function, forwards paramters to
///        ResidCalc_t::resid static method
///
/// Expects a DATA_PARAM_IDX member of ResidCalc_t::Fitter pointing
/// to the data
///
template <typename ResidCalculator>
static void minuit_f(Int_t&, Double_t*, Double_t &retval, Double_t *par, Int_t)
{
  using Fitter_t = typename ResidCalculator::Fitter;
  using FitterParams_t = typename ResidCalculator::FitParams;

  static const double BAD_VALUE = INFINITY;
  const auto &fitter = *(const Fitter_t*)(intptr_t)(par[Fitter_t::DATA_PARAM_IDX]);

  FitterParams_t params(par);
  if (params.is_invalid()) {
    retval = BAD_VALUE;
    return;
  }

  retval = ResidCalculator::resid(fitter, params);
}

/// Fit using masked histograms instead of data objects
template <typename ResidCalculator>
static void minuit_func_mrc(Int_t&, Double_t*, Double_t &retval, Double_t *par, Int_t)
{
  using Fitter_t = typename ResidCalculator::Fitter;
  using FitterParams_t = typename ResidCalculator::FitParams;

  static const double BAD_VALUE = INFINITY;
  const auto &fitter = *(const Fitter_t*)(intptr_t)(par[Fitter_t::DATA_PARAM_IDX]);

  FitterParams_t params(par);
  if (params.is_invalid()) {
    retval = BAD_VALUE;
    return;
  }

  retval = ResidCalculator::resid_mrc(fitter, params);
}


/// \brief ResidCalculatorPML
template <typename Fitter_t>
struct ResidCalculatorPML {
  using Fitter = Fitter_t;
  using FitParams = typename Fitter::FitParams;

  static double resid(const Fitter &f, const FitParams &p)
    { return f.resid_calc(p, loglikelihood_calc); }

  static double resid_mrc(const Fitter &f, const FitParams &p)
    { return f.resid_calc_mrc(p, *f.mrc, loglikelihood_calc); }

  static double resid(const Fitter &f, const typename Fitter_t::FitResult &p)
    { return f.resid_calc(p, loglikelihood_calc); }

  static double resid_mrc(const Fitter &f, const typename Fitter_t::FitResult &p)
    { return f.resid_calc_mrc(p, *f.mrc, loglikelihood_calc); }

};

template <typename Fitter_t>
struct ResidCalculatorChi2 {
  using Fitter = Fitter_t;
  using FitParams = typename Fitter::FitParams;

  static double resid(const Fitter &f, const FitParams &p)
    { return f.resid_chi2(p); }
};


#endif
