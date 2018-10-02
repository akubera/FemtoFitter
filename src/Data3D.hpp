///
/// \file Data3D.hpp
///


#pragma once

#include <array>
#include <valarray>


template <typename T> struct data_traits;

/// \class Data3D
/// \brief 3D Correlation function
///
/// Linearlized correlation function data
///
struct Data3D {
  using V = std::valarray<double>;

  std::array<V, 3> qspace;

  V num,
    den,
    qinv;

  size_t size() const
    { return num.size(); }

};

template <>
struct data_traits<Data3D> {
  // static const size_t ndim = decltype(Data3D::qspace)::size_type;
  static const size_t ndim = std::rank<decltype(Data3D::qspace)>::value;
};


template <>
struct std::rank<Data3D> {
  static const size_t value = data_traits<Data3D>::ndim;
};
