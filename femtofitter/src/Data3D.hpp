///
/// \file Data3D.hpp
///


#pragma once


#include <iostream>
#include <array>
#include <vector>
#include <valarray>
#include <memory>

#include <TH3.h>
#include <TDirectory.h>

class TH3;


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

    unsigned hist_bin;

    Datum(double q0, double q1,double q2, double n, double d, double qi, unsigned bin)
      : qo(q0)
      , qs(q1)
      , ql(q2)
      , num(n)
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
        variance = ratio / den * (1.0 + ratio);

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
  };

  std::vector<Datum> data;

  std::shared_ptr<Source> src;

  double limit,
         true_limit;

  double gamma;

  static std::unique_ptr<Data3D> From(TDirectory &tdir, double limit=0.0);

  /// Build out of standard tdirectory;
  static std::unique_ptr<Data3D> FromDirectory(TDirectory &tdir, double limit=0.0)
    {
      return FromDirectory(tdir, {"num", "den", "qinv"}, limit);
    }

  static std::unique_ptr<Data3D> FromDirectory(TDirectory &tdir, const std::array<TString, 3> &names, double limit=0.0)
    {
      return From(tdir, names[0], names[1], names[2], limit);
    }


  /// Build from directory and histogram names
  static std::unique_ptr<Data3D> From(TDirectory &tdir,
                                      const TString &num_name,
                                      const TString &den_name,
                                      const TString &qinv_name,
                                      double limit=0.0);

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

  /// Construct with unique pointers
  Data3D(std::unique_ptr<TH3> num, std::unique_ptr<TH3> den, std::unique_ptr<TH3> qinv, double limit=0.0);

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

  Data3D(const Data3D &orig, std::vector<Datum> dat)
    : data(std::move(dat))
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
  {
  }

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

  std::unique_ptr<Data3D> subset_pos_side() const
    {
      std::vector<Datum> subset;
      subset.reserve(size() / 2);
      for (auto &dat : data) {
        if (dat.qs > 1e-4) {
          subset.emplace_back(dat);
        }
      }
      subset.shrink_to_fit();

      return std::make_unique<Data3D>(*this, std::move(subset));
    }

  std::unique_ptr<Data3D> subset_neg_side() const
    {
      std::vector<Datum> subset;
      subset.reserve(size() / 2);
      for (auto &dat : data) {
        if (dat.qs < -1e-4) {
          subset.emplace_back(dat);
        }
      }
      subset.shrink_to_fit();

      return std::make_unique<Data3D>(*this, std::move(subset));
    }

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

  /// Estimate gamma from TDirectory name
  ///
  /// TDirectory name should have form <kTlo>_<kThi>.
  /// Mother directory is also checked.
  /// If sibling TH1 "kTDist" is found, it is used to find more accurate
  /// mean kT of pairs.
  ///
  static double calc_gamma_from_tdir(const TDirectory &tdir, double mass=0.13957);
};


inline
double Data3D::calc_gamma_from_tdir(const TDirectory &tdir, double mass)
{
  double gamma = 0.0;

  TString ktname = tdir.GetName();
  Ssiz_t underscore = ktname.First('_');

  // check parent directory for kt name
  if (underscore == TString::kNPOS || ktname.First('.') == TString::kNPOS) {
    if (const auto *mdir = tdir.GetMotherDir()) {
      ktname = mdir->GetName();
      underscore = ktname.First('_');
    }
  }

  if (underscore == TString::kNPOS || ktname.First('.') == TString::kNPOS) {
    std::cerr << "Warning: Could not determine kt-range from directory\n";
    return gamma;
  }

  const double
    mass_sqr = mass * mass,
    kt_lo = TString(ktname(0, underscore)).Atof(),
    kt_hi = TString(ktname(underscore+1, ktname.Length())).Atof(),
    est_mean_kt = (kt_hi + kt_lo) / 2.0,
    est_mean_kt_sqr = est_mean_kt * est_mean_kt;

  gamma = std::sqrt(1.0 + 4 * est_mean_kt_sqr / mass_sqr);

  if (auto *cent_dir = tdir.GetMotherDir()) {
    if (auto obj = std::unique_ptr<TObject>(cent_dir->Get("kTDist"))) {
      if (auto *kthist = dynamic_cast<TH1*>(obj.get())) {
        kthist->GetXaxis()->SetRangeUser(kt_lo, kt_hi);

        const double
          mean_kt = kthist->GetMean(),
          kt_mass_ratio_sqr = mean_kt * mean_kt / mass_sqr;

        gamma = std::sqrt(1.0 + 4 * kt_mass_ratio_sqr);
      }
    }
  }

  return gamma;
}
