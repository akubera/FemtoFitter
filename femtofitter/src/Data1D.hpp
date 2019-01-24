///
/// \file Data1D.hpp
///


#pragma once

#ifndef DATA1D_HPP
#define DATA1D_HPP

#include <array>
#include <memory>
#include <vector>
#include <valarray>

class TDirectory;
class TH1;


/// \class Data1D
/// \brief 1D Correlation function
///
struct Data1D {

  /// Unit of data (qinv, numerator, denominator)
  struct Datum {
    double qinv,
           num,
           den;
  };

  std::vector<Datum> data;

  double limit,
         true_limit,
         gamma;

  /// Build out of standard tdirectory;
  static std::unique_ptr<Data1D> FromDirectory(TDirectory &, double limit=0.0);

  /// Construct from histograms
  ///
  /// The q-space is taken from the numerator.
  /// It is assumed the axes of the three histograms are the same.
  ///
  Data1D(const TH1& num, const TH1& den, double limit=0.0);

  /// Build from TDirectory
  Data1D(TDirectory &tdir, double limit);

  size_t size() const
    { return data.size(); }

  const Datum& operator[](size_t idx) const
    { return data[idx]; }

  auto begin()
    { return data.begin(); }

  auto end()
    { return data.end(); }

  auto begin() const
    { return data.begin(); }

  auto end() const
    { return data.end(); }

};


#endif
