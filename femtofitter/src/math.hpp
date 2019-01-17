///
/// \file math.hpp
///

#include <cmath>

#define IF_NOT_ZERO(_val, expr) __builtin_expect(_val == 0.0, 0) ? 0.0 : expr


static const double
  HBAR_C = 0.19732697,
  HBAR_C_SQ = HBAR_C * HBAR_C;


/// \brief chi2 of ratio - assuming stderr is sqrt of values
///
inline
double chi2_calc(double N, double D, double C)
{
  const double
    R = N / D,
    variance = R * (1.0 + R) / D;

  return IF_NOT_ZERO(variance, (R-C) * (R-C) / variance);
}

inline
double loglikelihood_calc(double A, double B, double C)
{
  const double
    ta = IF_NOT_ZERO(A, A * std::log((A+B) * C / A / (C+1.0))),
    tb = IF_NOT_ZERO(B, B * std::log((A+B) / B / (C+1.0)));

  return -2.0 * (ta + tb);
}
