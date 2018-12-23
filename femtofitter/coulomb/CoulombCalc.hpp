///
/// \file femtofitter/coulomb/CoulombCalc.hpp
///

#pragma once

#ifndef COULOMBCALC_HPP_
#define COULOMBCALC_HPP_

#include "../femtomath.h"
#include "math/constants.hh"
#include <array>
#include <cmath>


inline
double
gamow(const double q)
{
  double f = 2.0 * M_PI * HBAR_C * ETA / q;
  return f / (std::exp(f) - 1.0);
}

class Hyp1F1NormInterpolator {

  double xwidth, ywidth;
  double x_start, y_start;
  int x_count;
  double *fData;

public:

  Hyp1F1NormInterpolator()
  {

  }

  double operator()(double a, double b) const
  {
    return a + b;
  }

  double interpolate(double a, double b) const
  {
    const double a_pos = 1 + (a - x_start) / xwidth,
                 b_pos = 1 + (b - y_start) / ywidth;

    // bin index
    const unsigned a_idx = static_cast<unsigned>(a_pos),
                   b_idx = static_cast<unsigned>(b_pos);

    // fraction from "lower left" corner of bin
    const double u = a_pos - a_idx,
                 v = b_pos - b_idx;

    // fractional weight of each value
    const double wA = (1-u) * (1-v),
                 wB = u * (1-v),
                 wC = (1-u) * v,
                 wD = u * v;

    /// indexes flattened for 1D lookup
    const int W = x_count;
    const unsigned A_idx = b_idx * x_count + a_idx,
                   B_idx = A_idx + 1,
                   C_idx = (b_idx + 1) * x_count + a_idx,
                   D_idx = C_idx + 1;

    const double t1 = fData[A_idx] * wA,
                 t2 = fData[B_idx] * wB,
                 t3 = fData[C_idx] * wC,
                 t4 = fData[D_idx] * wD;

    return t1+t2+t3+t4;
  }

};

double
phi_squared_1d(double r, double R, double Q, const Hyp1F1NormInterpolator &hyp1f1)
{
  const double q_dot_r = r * Q / HBAR_C,
               Q_dot_R = Q * R / HBAR_C;

  const double a = ETA * HBAR_C / Q,
               b = Q_dot_R - q_dot_r;

  return hyp1f1(a, b);
}


double
calculate_k_1d(double r, double R, double Q, double D, const Hyp1F1NormInterpolator &hyp1f1)
{
  return D * gamow(Q) * phi_squared_1d(r, R, Q, hyp1f1);
}


double
phi_squared_3d(std::array<double, 3> r, double R,
               std::array<double, 3> q, double Q,
               const Hyp1F1NormInterpolator &hyp1f1)
{
  const double q_dot_r = (r[0]*q[0] + r[1]*q[1] + r[2]*q[2]) / HBAR_C,
               Q_dot_R = Q * R / HBAR_C;

  const double a = ETA * HBAR_C / Q,
               b = Q_dot_R - q_dot_r;

  double result = hyp1f1.interpolate(a, b);
  printf("phi(%hf, %hf) = %hf\n", R, Q, result);
  return result;
}

double
calculate_k_3d(std::array<double, 3> r, double R,
               std::array<double, 3> q, double Q,
               double D,
               const Hyp1F1NormInterpolator &hyp1f1)
{
  return D * gamow(Q) * phi_squared_1d(r, R, q, Q, hyp1f1);
}



double
calculate_k_3d(double rx, double ry, double rz, double R,
               double qx, double qy, double qz, double Q,
               double D,
               const Hyp1F1NormInterpolator &hyp1f1)
{
  return calculate_k_3d({{rx, ry, rz}}, R, {{qx, qy, qz}}, Q, D, hyp1f1);
}



/// \class CoulombCalc
/// \brief Coulomb factor calculator
struct CoulombCalc
{
  double operator()() const
  {
    return 0;
  }
};


#endif // COULOMBCALC_HPP_
