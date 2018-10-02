///
/// \file Data3D.hpp
///


#pragma once

#include <array>
#include <vector>
#include <valarray>

class TH3;
class TDirectory;


template <typename T> struct data_traits;

/// \class Data3D
/// \brief 3D Correlation function
///
/// Linearized correlation function data
///
struct Data3D {
  using V = std::valarray<double>;
  // using V = std::vector<double>;

  std::array<V, 3> qspace;

  V num,
    den,
    qinv;

  /// Build out of standard tdirectory;
  Data3D(TDirectory &, double limit=0.0);

  /// Construct from histograms
  ///
  /// The q-space is taken from the numerator.
  /// It is assumed the axes of the three histograms are the same.
  ///
  Data3D(const TH3& num, const TH3& den, const TH3& qinv, double limit=0.0);

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
