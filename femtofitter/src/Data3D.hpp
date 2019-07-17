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
           num_err,
           den,
           qinv;

    unsigned hist_bin;

    Datum(double q0, double q1,double q2, double n, double d, double qi, unsigned bin)
      : qo(q0)
      , qs(q1)
      , ql(q2)
      , num(n)
      , num_err(std::sqrt(n))
      , den(d)
      , qinv(qi)
      , hist_bin(bin)
      {
      }

    const double
    calc_chi2(const double model) const
    {
      const double
        ratio = num / den,
        diff = ratio - model,
        variance = num_err * num_err / den / den + ratio * ratio / den;

      return variance == 0.0 ? 0.0 : diff * diff / variance;
    }

    std::array<double, 3>
    qspace() const
      { return {qo, qs, ql}; }

  };

  /// saved copy of source histograms
  struct Source {
    std::shared_ptr<const TH3> num;
    std::shared_ptr<const TH3> den;
    std::shared_ptr<const TH3> qinv;

    Source(const TH3& n, const TH3& d, const TH3& qi)
      : num(static_cast<TH3*>(n.Clone()))
      , den(static_cast<TH3*>(d.Clone()))
      , qinv(static_cast<TH3*>(qi.Clone()))
      { }

    Source(std::shared_ptr<const TH3> n,
           std::shared_ptr<const TH3> d,
           std::shared_ptr<const TH3> q)
      : num(n)
      , den(d)
      , qinv(q)
      { }

    Source(std::unique_ptr<const TH3> n, std::unique_ptr<const TH3> d)
      : num(std::move(n))
      , den(std::move(d))
      , qinv(nullptr)
      { }

    Source(const Source &orig) = default;
  };

  std::vector<Datum> data;

  std::shared_ptr<Source> src;

  double limit,
         true_limit;

  double gamma;

  /// Build out of standard tdirectory;
  static std::unique_ptr<Data3D> FromDirectory(TDirectory &, double limit=0.0);

  static std::unique_ptr<Data3D> FromDirectory(TDirectory &, double limit, double minimum);

  /// Data3D with applied momentum resolution correction
  static std::unique_ptr<Data3D> FromDirectory(TDirectory &tdir, const TH3 *mrc, double limit=0.0)
    { return FromDirectory(tdir, *mrc, limit); }

  /// Data3D with applied momentum resolution correction
  static std::unique_ptr<Data3D> FromDirectory(TDirectory &, const TH3 &mrc, double limit=0.0);

  /// Construct with advanced MRC
  static std::unique_ptr<Data3D> FromDirectory(TDirectory &data, TDirectory &mrc, double limit);

  /// Data3D with minimum ratio
  static std::unique_ptr<Data3D> FromDirectory(TDirectory &, const TH3 &mrc, double limit, double minimum);

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
    , src(orig.src)
    , limit(orig.limit)
    , true_limit(orig.true_limit)
    , gamma(orig.gamma)
  {}

  /// Move Constructor
  Data3D(Data3D &&orig)
    : data(std::move(orig.data))
    , src(std::move(orig.src))
    , limit(orig.limit)
    , true_limit(orig.true_limit)
    , gamma(orig.gamma)
  {}

  /// Move value from unique ptr
  Data3D(std::unique_ptr<Data3D> ptr)
    : data(std::move(ptr->data))
    , src(std::move(ptr->src))
    , limit(ptr->limit)
    , true_limit(ptr->true_limit)
    , gamma(ptr->gamma)
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

  const Datum& operator[](size_t idx) const
    { return data[idx]; }

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

  std::unique_ptr<Data3D> cone_subset(bool same_sign=true, double alpha=0.46) const
    {
      const double ratio = std::tan(alpha) * (same_sign ? -1 : 1);

      std::vector<Datum> subset;
      for (auto &dat : data) {

        bool accept = std::fabs(dat.ql) > 0.01
                   || (same_sign ? (dat.qs * dat.qo > 0.0) : (dat.qs * dat.qo < 0.0))
                   || std::abs(dat.qs / dat.qo) > ratio;

        if (accept) {
          subset.emplace_back(dat);
        }
      }

      return std::make_unique<Data3D>(subset, limit, true_limit);
    }

};
