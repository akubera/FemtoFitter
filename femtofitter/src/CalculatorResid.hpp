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
template <typename ResidCalc_t>
static void minuit_f(Int_t&, Double_t*, Double_t &retval, Double_t *par, Int_t)
{
  using Fitter_t = typename ResidCalc_t::Fitter;

  static const double BAD_VALUE = 3e99;
  const auto &data = *(const Fitter_t*)(intptr_t)(par[Fitter_t::DATA_PARAM_IDX]);

  typename Fitter_t::FitParams params(par);
  if (params.is_invalid()) {
    retval = BAD_VALUE;
    return;
  }

  retval = ResidCalc_t::resid(data, params);
}


/// \brief ResidCalculatorPML
template <typename Fitter_t>
struct ResidCalculatorPML {
  using Fitter = Fitter_t;
  using FitParams = typename Fitter::FitParams;

  static double resid(const Fitter &f, const FitParams &p)
    { return f.resid_calc(p, loglikelihood_calc); }

  static double resid(const Fitter &f, const typename Fitter_t::FitResult &p)
    { return f.resid_calc(p, loglikelihood_calc); }
};

template <typename Fitter_t>
struct ResidCalculatorChi2 {
  using Fitter = Fitter_t;
  using FitParams = typename Fitter::FitParams;

  static double resid(const Fitter &f, const FitParams &p)
    { return f.resid_chi2(p); }
};


#endif
