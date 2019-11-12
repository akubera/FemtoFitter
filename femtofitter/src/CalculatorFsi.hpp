///
/// \file FsiCalculator.hpp
///

#pragma once

#ifndef FSICALCULATOR_HPP_
#define FSICALCULATOR_HPP_

#include <array>
#include <vector>
#include <string>
#include <cmath>
#include <memory>
#include <functional>

class TH3;


/// \class FsiQinv
/// \brief Returns final-state-interaction from qinv
///
struct FsiQinv {

  virtual ~FsiQinv()
    {}

  virtual double operator()(double qinv) = 0;
};


/// \class FsiCalculator
/// \brief Interface for calculating final state interaction component
///         of correlation function
///
/// FSI should be called with
///
struct FsiCalculator {

  virtual ~FsiCalculator()
    {}

  virtual std::function<double(double)> ForRadius(double Rinv) = 0;

  virtual std::function<double(double)> ForRadius(double Ro,
                                                  double Rs,
                                                  double Rl,
                                                  double gamma=1.0)
    {
      double Rinv = std::sqrt(gamma*gamma*Ro*Ro + Rs*Rs + Rl*Rl);
      return ForRadius(Rinv);
    }

  virtual std::string ClassName() const = 0;
  virtual std::string Describe() const
    {
      return ClassName();
    }

  // Fill histogram with qinv datapoints
  virtual void FillQinvHist(TH3 &hist, double Ro, double Rs, double Rl, double gamma) const = 0;
};

#endif
