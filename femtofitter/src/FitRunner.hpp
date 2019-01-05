///
/// \file FitRunner.hpp
///


#pragma once

#ifndef FITRUNNER_HPP
#define FITRUNNER_HPP


#include <TMinuit.h>

template <typename T>
struct FitRunner {

  using FitResult = typename T::FitResult;

  double fit_factor;

  FitRunner(double ff)
    : fit_factor(ff)
  {}

  auto operator()(TMinuit &minuit) const -> FitResult
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

    return FitResult(minuit);
  }

};


template <typename T>
struct Chi2FitRunner : FitRunner<Chi2FitRunner<T>> {
  using FitResult = typename T::FitResult;

  Chi2FitRunner()
    : FitRunner<Chi2FitRunner<T>>(1.0)
  {}

};

template <typename T>
struct PmlFitRunner : FitRunner<PmlFitRunner<T>> {
  using FitResult = typename T::FitResult;

  PmlFitRunner()
    : FitRunner<PmlFitRunner<T>>(0.5)
  {}

};

#endif
