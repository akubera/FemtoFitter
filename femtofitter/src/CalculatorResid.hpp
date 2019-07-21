///
/// \file femtofitter/src/CalculatorResid.hpp
///

#pragma once

#ifndef CALCULATORRESID_HPP_
#define CALCULATORRESID_HPP_

#include "math/fit.hh"

#include <cstdint>


/// \brief static minuit function, forwards paramters to
///        ResidCalc_t::resid static method
///
/// Expects a DATA_PARAM_IDX member of ResidCalc_t::Fitter pointing
/// to the data
///
template <typename ResidCalculator>
static void minuit_func(int&, double*, double &retval, double *par, int)
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


/// Fit using momentum-resolution smearing of the correlation function
template <typename ResidCalculator>
static void minuit_func_mrc(int&, double*, double &retval, double *par, int)
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


/// \class ResidCalculatorT
/// \brief Abstraction over residual calcuations
///
template <typename Calc, typename Fitter_t>
struct ResidCalculatorT {
  using Fitter = Fitter_t;
  using FitParams = typename Fitter::FitParams;

  static constexpr auto resid_func = Calc::residual;

  static double resid(const Fitter &f, const FitParams &p)
    { return f.resid_calc(p, resid_func); }

  static double resid_mrc(const Fitter &f, const FitParams &p)
    { return f.resid_calc_mrc(p, *f.mrc, resid_func); }

  static double resid(const Fitter &f, const typename Fitter_t::FitResult &p)
    { return f.resid_calc(p, resid_func); }

  static double resid_mrc(const Fitter &f, const typename Fitter_t::FitResult &p)
    { return f.resid_calc_mrc(p, *f.mrc, resid_func); }

};


struct Chi2ResidualFunctor {
  static constexpr auto residual = chi2_calc;
};

struct LogLikelihoodResidualFunctor {
  static constexpr auto residual = loglikelihood_calc;
};

template <typename Fitter_t>
using ResidCalculatorPML = ResidCalculatorT<LogLikelihoodResidualFunctor, Fitter_t>;

template <typename Fitter_t>
using ResidCalculatorChi2 = ResidCalculatorT<Chi2ResidualFunctor, Fitter_t>;


#endif
