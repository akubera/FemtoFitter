///
/// \file Data3D.hpp
///


#pragma once

#include <array>
#include <vector>
#include <valarray>

#include <TH3.h>

class TH3;
class TDirectory;


template <typename T> struct data_traits;

/// \class Data3D
/// \brief 3D Correlation function
///
/// Linearized correlation function data
///
struct Data3D {

  /// Unit of data
  struct Datum {
    double qo,
           qs,
           ql,
           num,
           den,
           qinv;
  };

  std::vector<Datum> data;

  double limit,
         true_limit;

  /// Build out of standard tdirectory;
  static std::unique_ptr<Data3D> FromDirectory(TDirectory &, double limit=0.0);

  /// Construct from histograms
  ///
  /// The q-space is taken from the numerator.
  /// It is assumed the axes of the three histograms are the same.
  ///
  Data3D(const TH3& num, const TH3& den, const TH3& qinv, double limit=0.0);

  size_t size() const
    { return data.size(); }

  auto begin()
    { return data.begin(); }

  auto end()
    { return data.end(); }

  auto begin() const
    { return data.cbegin(); }

  auto end() const
    { return data.cend(); }

};
