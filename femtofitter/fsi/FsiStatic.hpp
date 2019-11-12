///
/// \file femtofitter/fsi/FsiStatic.hpp
///

#pragma once

#ifndef FSI_FSISTATIC_HPP
#define FSI_FSISTATIC_HPP

#include "CalculatorFsi.hpp"

#include <vector>


/// Final state interaction calculator that always returns a single
/// value (default 1.0)
///
/// Static Value Calculator - used primarliy to ignore FSI
///
struct FsiStatic : public FsiCalculator {
  double value;

  FsiStatic(double val=1.0)
    : value(val)
    {
    }

  struct Kcalc : FsiQinv {

    double val;

    Kcalc(double v)
      : val(v)
      {}

    double operator()(double qinv) override
      {
        return val;
      }
  };

  std::function<double(double)> ForRadius(double) override
    {
      return [value=value](double) { return value; };
    }

  void Fill(std::vector<double> &dest, double _R)
    {
      std::fill(dest.begin(), dest.end(), value);
    }

  void FillQinvHist(TH3 &hist, double _Ro, double _Rs, double _Rl, double _gamma) const override
    {
      #pragma omp for
      for (int k=1; k <= hist.GetNbinsZ(); ++k)
      for (int j=1; j <= hist.GetNbinsY(); ++j)
      for (int i=1; i <= hist.GetNbinsX(); ++i) {
        hist.SetBinContent(i, j, k, value);
      }
    }

  std::string ClassName() const override
    { return std::string("FsiStatic[") + Form("%g", value) + "]"; }

  static std::shared_ptr<FsiCalculator> new_shared_ptr(double value)
    { return std::make_shared<FsiStatic>(value); }

  static std::shared_ptr<FsiCalculator> From(double value)
    { return std::make_shared<FsiStatic>(value); }
};

#endif
