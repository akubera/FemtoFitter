///
/// \file src/Fitter.hpp
///

#pragma once


#include <valarray>


/// struct
struct Value {
  double first, second;
  operator double() const
    { return first; }
};


template <typename Fitter>
struct ResidCalculatorPML {
  using FitParams = typename Fitter::FitParams;
  static double resid(const Fitter &f, const FitParams &p)
    { return f.resid_pml(p); }
};

template <typename Fitter>
struct ResidCalculatorChi2 {
  using FitParams = typename Fitter::FitParams;
  static double resid(const Fitter &f, const FitParams &p)
    { return f.resid_chi2(p); }
};


const double
  HBAR_C = .297,
  HBAR_C_SQ = HBAR_C * HBAR_C;


double chi2_calc(double N, double D, double C)
{
  const double
    R = N / D,
    // variance = (1.0 + R) * R * R / N;
    // variance = R * (N + D) / (D*D);
    // variance = R * std::sqrt((1.0/N  + 1.0/D));
    variance = R * std::sqrt((1.0 + R) / N);

  return (N == 0.0)
       ? 0.0
       : (R-C) * (R-C) / variance;
}

double loglikelihood_calc(double A, double B, double C)
{
  const double
    ta = (A == 0.0) ? 0.0 : A * std::log((C/A * (A+B) / (C+1.0))),
    tb = B * std::log((A+B) / B / (C+1.0));

  return -2 * (ta + tb);
}

/// \class Fitter
/// \brief Abstract base class for fitting
///
///
class Fitter {
//     std::atomic<std::unique_ptr<std::valarray<double>>> fRatio {nullptr};

  std::valarray<double> num,
                        den;

  std::valarray<double> fRatio;

public:
  Fitter();

  template <typename Params>
  double
  chi2(const Params &p) const
  {
    auto diff = ratio() - evaluate(p);
    auto e = fRatio * (1.0 + fRatio) / den / den;
    return (diff * diff / e).sum();
  }

  const std::valarray<double>&
  ratio() const
  {
    if (fRatio.size() == 0) {
        const_cast<std::valarray<double>&>(fRatio) = num / den;
    }
    return fRatio;
  }

  // template <typename Params>
  // std::valarray<double> evaluate(Params &const) const;
};


// template <typename T>
// std::valarray<double>
// Fitter::evaluate(Params &const) const
// {
// }
