///
/// \file femtofitter/fsi/FsiStatic.hpp
///

#pragma once

#ifndef FSI_FSISTATIC_HPP
#define FSI_FSISTATIC_HPP

#include "CalculatorFsi.hpp"

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

  double
  operator()(double R) const override
    {
      return value;
    }

};


#endif
