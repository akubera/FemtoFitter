///
/// \file Data1D.hpp
///


#pragma once

#ifndef DATA1D_HPP
#define DATA1D_HPP

#include <TH1.h>
#include <TH1C.h>

#include <array>
#include <memory>
#include <vector>
#include <valarray>

class TDirectory;

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

  /// reference back to source histograms
  struct Source {
    std::shared_ptr<const TH1> num;
    std::shared_ptr<const TH1> den;

    Source(const TH1& n, const TH1& d)
      : num(static_cast<TH1*>(n.Clone()))
      , den(static_cast<TH1*>(d.Clone()))
      {}

    Source(std::unique_ptr<const TH1> n, std::unique_ptr<const TH1> d)
      : num(std::move(n))
      , den(std::move(d))
      {}
  };

  std::vector<Datum> data;

  double limit,
         true_limit,
         gamma;

  std::shared_ptr<Source> src;

  std::unique_ptr<TH1C> mask;

  /// Build out of standard tdirectory;
  static std::unique_ptr<Data1D> From(TDirectory &data, double limit=0.0);
  static std::unique_ptr<Data1D> From(TDirectory &data, const TH1 &mrc, double limit=0.0);

  /// Build from TDirectories of data and momentum resolution correction
  static std::unique_ptr<Data1D> From(TDirectory &data, TDirectory &mrc, double limit=0.0);

  /// Construct from histograms
  ///
  /// The q-space is taken from the numerator.
  /// It is assumed the axes of the three histograms are the same.
  ///
  Data1D(const TH1& num, const TH1& den, double limit=0.0);

  /// Build from TDirectory
  Data1D(TDirectory &tdir, double limit);

  /// copy constructor
  Data1D(const Data1D&);

  /// assignment operator
  Data1D& operator=(const Data1D&) = delete;

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
