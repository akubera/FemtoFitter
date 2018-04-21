///
/// \file femtomath.h
///


#ifndef FEMTOMATH_H
#define FEMTOMATH_H


#include <cmath>


const double HBAR_C = 0.19732697, // GeV * fm
             HBAR_C_SQ = HBAR_C * HBAR_C;

/// Return the principle least likelihood $\chi^2_{PML}$
inline
double
loglikelihood_calc(double A, double B, double C)
{
  const double
    ta = (A == 0) ? 0.0 : A * std::log((C/A * (A+B) / (C+1.0))),
    tb = (B == 0) ? 0.0 : B * std::log((A+B) / B / (C+1.0));

  return ta + tb;
}

struct LoglikelihoodCalculator {
  static double Calculate(double a, double b, double c) { return loglikelihood_calc(a, b, c); }
};

/// Calculate the $\chi^2$ of a ratio of numbers, assuming errors based on
/// Poissonian statistics ($\sqrt{N}$)
inline
double chi2_ratio_calc(double n, double d, double ratio)
{
  const double observed = n / d,
               variance = (n*n + n*d)/(d*d*d),
               diff = observed - ratio,
               chi2 = diff*diff/variance;
  return chi2;
}

struct Chi2Calculator {
  static double Calculate(double a, double b, double c) { return chi2_ratio_calc(a, b, c); }
};



#endif
