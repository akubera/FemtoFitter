///
/// \file Data3D.hpp
///


#pragma once

#include <array>
#include <vector>
#include <valarray>
#include <memory>

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

  double gamma;

  /// Build out of standard tdirectory;
  static std::unique_ptr<Data3D> FromDirectory(TDirectory &, double limit=0.0);

  /// Data3D with applied momentum resolution correction
  static std::unique_ptr<Data3D> FromDirectory(TDirectory &tdir, const TH3 *mrc, double limit=0.0)
    { return FromDirectory(tdir, *mrc, limit); }

  /// Data3D with applied momentum resolution correction
  static std::unique_ptr<Data3D> FromDirectory(TDirectory &, const TH3 &mrc, double limit=0.0);

  /// Default empty constructor - use is discouraged
  Data3D();

  /// Construct from histograms
  ///
  /// The q-space is taken from the numerator.
  /// It is assumed the axes of the three histograms are the same.
  ///
  Data3D(const TH3& num, const TH3& den, const TH3& qinv, double limit=0.0);

  /// Data limit
  Data3D(std::vector<Datum> data, double limit, double true_limit);

  /// Copy Constructor
  Data3D(const Data3D &orig)
    : data(orig.data)
    , limit(orig.limit)
    , true_limit(orig.true_limit)
  {}

  /// Move Constructor
  Data3D(Data3D &&ptr)
    : data(std::move(ptr.data))
    , limit(ptr.limit)
    , true_limit(ptr.true_limit)
  {}

  /// Move value from unique ptr
  Data3D(std::unique_ptr<Data3D> ptr)
    : data(std::move(ptr->data))
    , limit(ptr->limit)
    , true_limit(ptr->true_limit)
  {}

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

  /// get out-side datapoints from quadrants I-III
  std::unique_ptr<Data3D> cowboy_subset() const
    {
      std::vector<Datum> subset;
      for (auto &dat : data) {
        bool accept =
          ((dat.qo > 0.0) && (dat.qs > 0.0)) ||
          ((dat.qo < 0.0) && (dat.qs < 0.0));

        if (accept) {
          subset.emplace_back(dat);
        }
      }

      return std::make_unique<Data3D>(subset, limit, true_limit);
    }

  /// get out-side datapoints from quadrants II-IV
  std::unique_ptr<Data3D> sailor_subset() const
    {
      std::vector<Datum> subset;
      for (auto &dat : data) {
        bool accept =
          ((dat.qo < 0.0) && (dat.qs > 0.0)) ||
          ((dat.qo > 0.0) && (dat.qs < 0.0));

        if (accept) {
          subset.emplace_back(dat);
        }
      }

      return std::make_unique<Data3D>(subset, limit, true_limit);
    }
};
