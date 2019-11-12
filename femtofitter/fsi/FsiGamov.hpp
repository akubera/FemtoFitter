///
/// \file fsi/FsiGamov.hpp
///

#pragma once

#ifndef FSIGAMOV_HPP
#define FSIGAMOV_HPP


#include "CalculatorFsi.hpp"
#include "math/constants.hh"


/// \class FsiGamov
/// \brief Coulomb calculation assuming point-source
///
//struct FsiGamov : public FsiCalculatorImpl<FsiGamov> {
struct FsiGamov : public FsiCalculator {

  struct Kcalc : FsiQinv {

    double operator()(double qinv) override
      {
        return FsiGamov::GamovFactor(qinv);
      }
  };

  // std::unique_ptr<FsiQinv> ForRadius(double Rinv) override
  //   {
  //     return std::make_unique<Kcalc>();
  //   }

  std::function<double(double)> ForRadius(double Rinv) override
    {
      return FsiGamov::GamovFactor;
    }

  static double GamovFactor(double qinv)
    {
      double x = 2.0 * M_PI * HBAR_C * ETA_PION / qinv;
      return x / (std::exp(x) - 1);
    }

  void FillQinvHist(TH3 &hist, double _Ro, double _Rs, double _Rl, double _gamma) const override
    {
      #pragma omp for
      for (int k=1; k <= hist.GetNbinsZ(); ++k)
      for (int j=1; j <= hist.GetNbinsY(); ++j)
      for (int i=1; i <= hist.GetNbinsX(); ++i) {
        const double q = hist.GetBinContent(i, j, k);
        hist.SetBinContent(i, j, k, FsiGamov::GamovFactor(q));
      }
    }

  std::string ClassName() const override
    { return "FsiGamov"; }

  static std::shared_ptr<FsiCalculator> new_shared_ptr()
    { return std::make_shared<FsiGamov>(); }

  static std::shared_ptr<FsiCalculator> From()
    { return std::make_shared<FsiGamov>(); }
};

#endif
