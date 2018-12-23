///
/// \file femtofitter/math/fit.hh
///

/// \brief chi2 of ratio - assuming stderr is sqrt of values
///
inline
double chi2_calc(double N, double D, double C)
{
  const double
    R = N / D,
    // variance = (N * N * D + D * D * N) / (D * D * D * D)
    // variance = (1.0 + R) * R * R / N;
    variance = R * (N + D) / (D*D);
    // variance = R * std::sqrt((1.0/N  + 1.0/D));
    // variance = R * std::sqrt((1.0 + R) / N);

  // return (N == 0.0)
  return (variance == 0.0)
       ? 0.0
       : (R-C) * (R-C) / variance;
}


/// \brief Calculation of loglikelihood
///
inline
double loglikelihood_calc(double A, double B, double C)
{
  const double
    ta = (A == 0.0) ? 0.0 : A * std::log((A+B) * C / A / (C+1.0)),
    tb = (B == 0.0) ? 0.0 : B * std::log((A+B) / B / (C+1.0));

  return -2.0 * (ta + tb);
}
