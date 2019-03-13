///
/// \file femtofitter/fsi/FsiPrecalc.hpp
///


#ifndef FSI_FSIPRECALC_HPP
#define FSI_FSIPRECALC_HPP

#include <array>
#include <vector>
#include "Data3D.hpp"


/// \class FsiPrecalc
/// \brief Precalculated final state interaction
///
struct FsiPrecalc {

  ///
  std::vector<double> qinv;

  FsiPrecalc(const Data3D &data)
    {
      for (auto &data = )
    }

  double operator()(const std::array<double 3> R, double gamma)
    {
      double a = gamma;
    }

  std::vector<double> operator()()
    {

    }
};



#endif
