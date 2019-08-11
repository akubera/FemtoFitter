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

/// \class FsiQinv
/// \brief Returns final-state-interaction from qinv
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

  // virtual std::unique_ptr<FsiQinv> ForRadius(double Rinv) = 0;

  // virtual double operator()(double R) const = 0;
  virtual std::function<double(double)> ForRadius(double Rinv) = 0;

/*
  virtual double operator()(const std::array<double, 3> qinv,
                            const std::array<double, 3> &R,
                            double gamma) const
    {
      const double
        Ro = R[0],
        Rs = R[1],
        Rl = R[2],
        Rinv = std::sqrt((gamma * gamma * Ro*Ro + Rs*Rs + Rl*Rl) / 3);

      return (*this)(qinv, Rinv);
    }

  virtual double operator()(const std::array<double, 3> &R, double gamma) const
    {
      const double
        Ro = R[0],
        Rs = R[1],
        Rl = R[2],
        Rinv = std::sqrt((gamma * gamma * Ro*Ro + Rs*Rs + Rl*Rl) / 3);

      return (*this)(Rinv);
    }
*/
  // virtual double eval(double qinv) const = 0;

  virtual std::string ClassName() const = 0;
  virtual std::string Describe() const
    {
      return ClassName();
    }
};


#if false
template <typename CRTP>
class FsiCalculatorImpl : public FsiCalculator {
public:

/*
  std::vector<double> calculator()(std::vector<double> &qinv, double Rinv) const
    {
      std::vector<double> K;
      K.reserve(qinv.size());

      std::transform(qinv.begin(), qinv.end(), std::back_inserter(K),
                     [&](double q) { return static_cast<const CRTP&>(*this)(q, Rinv); });
      //                [&](double q) { return (*this)(q, Rinv); });

      return K;
    }
*/

  static std::shared_ptr<CRTP> NewPtr()
    {
      return std::make_shared<CRTP>();
    }

};
#endif


#endif
