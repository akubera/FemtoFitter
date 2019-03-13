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
      return [=](double) { return value; };
    }

  void Fill(std::vector<double> &dest, double _R)
    {
      for (auto &n : dest) {
        n = value;
      }
    }

  std::string ClassName() const override
    { return std::string("FsiStatic[") + Form("%g", value) + "]"; }

  static std::shared_ptr<FsiCalculator> new_shared_ptr(double value)
    { return std::make_shared<FsiStatic>(value); }
};

#endif
