///
/// \file femtofitter/data/DataTypes.hpp
///

#pragma once


#ifndef DATA_DATATYPES_HPP
#define DATA_DATATYPES_HPP

#include <cmath>


template <typename DType>
struct DataTypes {

  constexpr double gamma_from_kt(double kt)
    {
      const double rel_kt = kt / DType::mass;
      return std::sqrt(rel_kt * rel_kt + 1.0);
    };

};


struct PionData {
  static constexpr double mass = 0.13957;
  static constexpr double eta = 1.0 / 388.0;
};


#endif
