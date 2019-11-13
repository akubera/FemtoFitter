///
/// \file femtofitter/fsi/Fsi.hpp
///

#pragma once

#ifndef FSI_FSI_HPP
#define FSI_FSI_HPP


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

  // Fill histogram with qinv datapoints
  virtual void FillQinvHist(TH3 &hist, double Ro, double Rs, double Rl, double gamma);
};



// template <typename Subclass, typename Mixin>
// str

#endif
